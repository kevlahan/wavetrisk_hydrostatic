"""One-time mechanical migration of scalar record accesses, not solver logic.

Run only on the original direct-array access sites after changing allocation
and the two INOUT call sites. Existing scalar storage API calls are retained.
"""
import argparse
from pathlib import Path
import re

REF=r'block_scalar_tendency\s*\([^()]+\)\s*%\s*(?:patch|bdry|ghost)'
REF_RE=re.compile(REF,re.I)


def close_paren(text,start):
    depth=0;quote=None
    for i in range(start,len(text)):
        c=text[i]
        if quote:
            if c==quote:quote=None
            continue
        if c in ('"',"'"):quote=c
        elif c=='(':depth+=1
        elif c==')':
            depth-=1
            if depth==0:return i
    raise ValueError('Unclosed access: '+text[start:])


def bounds(index):
    depth=0
    for i,c in enumerate(index):
        if c in '([':depth+=1
        elif c in ')]':depth-=1
        elif c==':' and depth==0:return index[:i].strip(),index[i+1:].strip()
    return None


def reads(text):
    text=re.sub(r'\bsize\(\s*('+REF+r')\s*,\s*kind\s*=\s*int64\s*\)',r'int(scalar_extent(\1),int64)',text,flags=re.I)
    text=re.sub(r'\bsize\(\s*('+REF+r')\s*\)',r'scalar_extent(\1)',text,flags=re.I)
    text=re.sub(r'\ballocated\(\s*('+REF+r')\s*\)',r'scalar_is_allocated(\1)',text,flags=re.I)
    result=[];pos=0
    for match in REF_RE.finditer(text):
        if match.start()<pos:continue
        at=match.end()
        while at<len(text) and text[at].isspace():at+=1
        if at==len(text) or text[at]!='(':continue
        end=close_paren(text,at)
        index=text[at+1:end].strip();span=bounds(index)
        replacement=(f'scalar_read_range({match[0]},{span[0]},{span[1]})' if span else
                     f'scalar_read({match[0]},{index})')
        result.extend([text[pos:match.start()],replacement]);pos=end+1
    return ''.join(result)+text[pos:]


def transform(text):
    conditional=re.match(r'\s*if\s*\(',text,re.I)
    if conditional:
        end=close_paren(text,conditional.end()-1)
        return reads(text[:end+1])+' '+transform(text[end+1:].strip())
    match=REF_RE.match(text.strip())
    if match:
        text=text.strip();at=match.end();index=None
        if at<len(text) and text[at]=='(':
            end=close_paren(text,at);index=text[at+1:end];at=end+1
        rest=text[at:].lstrip()
        if rest.startswith('=') and not rest.startswith('=>'):
            rhs=reads(rest[1:].strip());span=bounds(index) if index is not None else ('1',f'scalar_extent({match[0]})')
            if span:return f'call scalar_write_range({match[0]},{span[0]},{span[1]},{rhs})'
            return f'call scalar_write({match[0]},{index},{rhs})'
    return reads(text)


def wrap(text,indent):
    result=[];prefix=indent
    while len(prefix+text)>125:
        quote=None;breaks=[]
        for i,c in enumerate(text[:125-len(prefix)]):
            if quote:
                if c==quote:quote=None
            elif c in ('"',"'"):quote=c
            elif c in (' ', ','):breaks.append(i+1)
        if not breaks:raise ValueError('Cannot safely wrap: '+text)
        at=breaks[-1];result.append(prefix+text[:at].rstrip()+' &\n')
        text=text[at:].lstrip();prefix=indent+'     '
    result.append(prefix+text+'\n')
    return ''.join(result)


def migrate(source):
    lines=source.splitlines(keepends=True);out=[];at=0;changes=0
    while at<len(lines):
        first=at;parts=[lines[at]];at+=1
        if parts[0].lstrip().startswith('!'):
            out.extend(parts);continue
        while parts[-1].rstrip().endswith('&'):
            parts.append(lines[at]);at+=1
        text=' '.join(s.strip().rstrip('&').lstrip('&').strip() for s in parts)
        if not REF_RE.search(text):out.extend(parts);continue
        updated=transform(text)
        if updated==text:out.extend(parts);continue
        indent=re.match(r'\s*',lines[first])[0]
        out.append(wrap(updated,indent));changes+=1
    return ''.join(out),changes


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('path',type=Path)
    a=p.parse_args();text,n=migrate(a.path.read_text());a.path.write_text(text)
    print('Mechanically migrated statements:',n)

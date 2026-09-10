#!/usr/bin/env python3
"""Inject allocation-size inquiries into an isolated, verified block archive.

Counts allocated capacity, including nested owned allocatables and descriptor
storage. Pointer targets and module-local temporary allocations are excluded.
This is a diagnostic build, never a timing executable or a numerical change.
"""
import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import re
import shutil
from build_legacy import tracked_entries, verify
from experiment import digest
from prepare_legacy_profile import replace_once

OWNERS=('domain','comm_mpi','parallel_block','parallel_block_mpi')
TYPE_SOURCES=('shared','patch','dyn_array',*OWNERS)
OPAQUE={'mpi_request','mpi_comm','mpi_datatype','mpi_status'}


def split_top(text):
    result=[];level=0;start=0
    for i,ch in enumerate(text):
        if ch=='(':level+=1
        elif ch==')':level-=1
        elif ch==',' and level==0:result.append(text[start:i].strip());start=i+1
    if level:raise ValueError('Unbalanced declaration: '+text)
    return result+[text[start:].strip()]


def statements(text):
    pending=''
    for line in text.splitlines():
        line=line.split('!',1)[0].strip().lower()
        if not line:continue
        pending+=' '+line.lstrip('&').rstrip('&').strip()
        if line.endswith('&'):continue
        yield pending.strip();pending=''
    if pending:raise ValueError('Unfinished continuation')


@dataclass
class Field:
    name:str
    kind:str
    rank:int
    allocated:bool
    pointer:bool


def declaration(line):
    if '::' not in line:return []
    left,right=line.split('::',1)
    if re.search(r'\b(parameter|external|procedure)\b',left):return []
    match=re.match(r'(type\s*\(\s*(\w+)\s*\)|real(?:\([^)]*\))?|integer(?:\([^)]*\))?|logical(?:\([^)]*\))?|character(?:\([^)]*\))?)(?=\s|,|$)',left)
    if not match:
        if 'allocatable' in left:raise ValueError('Unsupported owning declaration: '+line)
        return []
    kind=match[2] or 'intrinsic'
    dims=re.search(r'dimension\s*\(([^)]*)\)',left)
    result=[]
    for var in split_top(right):
        var=var.split('=')[0].strip()
        m=re.fullmatch(r'(\w+)\s*(?:\((.*)\))?',var)
        if not m:raise ValueError('Unsupported declaration: '+line)
        shape=m[2] if m[2] is not None else dims[1] if dims else None
        result.append(Field(m[1],kind,len(split_top(shape)) if shape is not None else 0,
                            'allocatable' in left,'pointer' in left))
    return result


def specification(text):
    types={};roots=[];current=None;interface=False
    for line in statements(text):
        if line=='contains':break
        if line.startswith(('interface','abstract interface')):interface=True;continue
        if line.startswith('end interface'):interface=False;continue
        if interface:continue
        m=re.fullmatch(r'type\s*(?:,\s*(?:public|private)\s*)?(?:::)?\s+(\w+)',line)
        if m:current=m[1];types[current]=[];continue
        if line.startswith('end type'):current=None;continue
        fields=declaration(line)
        if current:types[current].extend(fields)
        else:roots.extend(fields)
    return types,roots


def generate(owner,roots,types):
    rows=[];body=[];omitted=[];depth_max=0
    def dynamic(kind,trail=()):
        if kind=='intrinsic' or kind in OPAQUE:return False
        if kind not in types:raise ValueError('Unresolved owned type: '+kind)
        if kind in trail:raise ValueError('Recursive owned type: '+kind)
        return any(not f.pointer and (f.allocated or dynamic(f.kind,trail+(kind,))) for f in types[kind])
    def walk(field,expr,key,depth):
        nonlocal depth_max
        if field.pointer:omitted.append(key+' (pointer target)');return
        has_children=dynamic(field.kind)
        if not field.allocated and not has_children:return
        if field.allocated:
            body.append(f'    if (allocated({expr})) then')
            rows.append(key);n=len(rows)
            elements=f'size({expr},kind=mc_i8)' if field.rank else '1_mc_i8'
            body.extend([f'    mc_bytes({n})=mc_bytes({n}) + &',
                         f'         {elements} * &',f'         int(storage_size({expr}),mc_i8)/8_mc_i8'])
        if has_children:
            scalar=expr
            if field.rank:
                indexes=[]
                for d in range(1,field.rank+1):
                    idx='mc_j'+str(depth+d);indexes.append(idx)
                    body.extend([f'    do {idx}=lbound({expr},{d}), &',f'         ubound({expr},{d})'])
                depth_max=max(depth_max,depth+field.rank)
                scalar=expr+'('+','.join(indexes)+')'
            for child in types[field.kind]:walk(child,scalar+'%'+child.name,key+'.'+child.name,depth+field.rank)
            body.extend(['    end do']*field.rank)
        if field.allocated:body.append('    end if')
    for root in roots:walk(root,root.name,owner+'.'+root.name,0)
    if not rows:raise ValueError('No allocation roots in '+owner)
    if len(set(rows))!=len(rows):raise ValueError('Duplicate allocation path')
    lines=[f'  subroutine census_{owner}(mc_unit)',
           '    use iso_fortran_env, only : mc_i8 => int64','    integer, intent(in) :: mc_unit',
           f'    integer(mc_i8) :: mc_bytes({len(rows)})']
    for d in range(1,depth_max+1):lines.append('    integer :: mc_j'+str(d))
    lines+=['    mc_bytes=0_mc_i8',*body]
    for n,key in enumerate(rows,1):
        lines.append(f"    write(mc_unit,'(i0,1x,a)') mc_bytes({n}), &")
        pieces=[key[i:i+65] for i in range(0,len(key),65)]
        lines.extend("         '"+part+"'"+(' // &' if i<len(pieces)-1 else '') for i,part in enumerate(pieces))
    lines.append(f'  end subroutine census_{owner}')
    if any(len(line)>132 for line in lines):raise ValueError('Generated Fortran line exceeds 132 columns')
    return '\n'.join(lines)+'\n',{'allocations':rows,'excluded_pointer_targets':omitted}


DRIVER='''
  subroutine allocation_census(label)
    use iso_fortran_env, only : mc_i8 => int64
    use domain_mod, only : census_domain
    use comm_mpi_mod, only : census_comm_mpi
    use parallel_block_mod, only : census_parallel_block
    character(*), intent(in) :: label
    character(16) :: value
    character(80) :: filename
    integer :: status
    integer, save :: unit=0, sequence=0
    logical, save :: initialized=.false., enabled=.false.
    if (.not.initialized) then
       call get_environment_variable('WAVETRISK_ALLOCATION_CENSUS',value,status=status)
       enabled=status==0.and.trim(value)=='1'
       initialized=.true.
       if (enabled) then
          write(filename,'(a,i0,a)') 'allocation-rank-',rank,'.txt'
          open(newunit=unit,file=trim(filename),status='new',action='write')
       end if
    end if
    if (.not.enabled) return
    sequence=sequence+1
    write(unit,'(a,1x,i0,1x,a)') 'sample',sequence,label
    call census_domain(unit)
    call census_comm_mpi(unit)
    call census_parallel_block(unit)
    call census_parallel_block_mpi(unit)
    flush(unit)
  end subroutine allocation_census
'''


def prepare(repo,baseline,out,ref):
    repo=repo.resolve(strict=True);baseline=baseline.resolve(strict=True);out=out.resolve()
    if (baseline/'.git').exists():raise ValueError('Baseline must be an archived build, not a Git checkout')
    if any(out.is_relative_to(p) or p.is_relative_to(out) for p in (repo,baseline)):
        raise ValueError('Use a new external output directory')
    entries=tracked_entries(repo,ref);verify(baseline,entries)
    sources={name:(baseline/'src'/f'{name}.f90').read_text() for name in TYPE_SOURCES}
    types={};roots={}
    for name,text in sources.items():
        local,roots[name]=specification(text)
        if types.keys() & local.keys():raise ValueError('Duplicate derived type')
        types.update(local)
    generated={};coverage={}
    for name in OWNERS:generated[name],coverage[name]=generate(name,roots[name],types)
    shutil.copytree(baseline,out,symlinks=True)
    for name in OWNERS:
        source=sources[name]
        source=re.sub(r'(?im)^\s*contains\s*$',f'  public :: census_{name}\ncontains',source,count=1)
        anchor=re.search(r'(?im)^\s*end module\b',source)
        if not anchor:raise ValueError('Missing module end')
        source=source[:anchor.start()]+'\n'+generated[name]+source[anchor.start():]
        if name=='parallel_block_mpi':
            source=source.replace('  public :: census_parallel_block_mpi','  public :: allocation_census, census_parallel_block_mpi',1)
            source=replace_once(source,'    if (detail_enabled) call sample_scalar_working_storage',
                                "    call allocation_census('scalar-storage-ready')\n"+
                                '    if (detail_enabled) call sample_scalar_working_storage')
            anchor=re.search(r'(?im)^\s*end module\b',source)
            source=source[:anchor.start()]+'\n'+DRIVER+source[anchor.start():]
        (out/'src'/f'{name}.f90').write_text(source)
    main=out/'src/main.f90';source=main.read_text()
    source=replace_once(source,'module main_mod\n\n',
                        'module main_mod\n  use parallel_block_mpi_mod, only : allocation_census\n\n')
    source=replace_once(source,'  end subroutine time_step',"    call allocation_census('step-end')\n  end subroutine time_step")
    main.write_text(source)
    make=out/'Makefile'
    make.write_text(make.read_text()+'\n# Isolated allocation census dependency\n'+
                    '$(BUILD_DIR)/parallel_block_mpi.o: $(BUILD_DIR)/comm_mpi.o\n')
    manifest={'revision':ref,'baseline_binary_sha256':digest(baseline/'bin/climate'),
              'scope':'Owned allocatable capacity in four modules; pointer aliases excluded; no MPI synchronization',
              'excluded':'Other modules, automatic arrays, allocator metadata, MPI internals, OS/runtime and transient allocation peaks',
              'coverage':coverage,'changed_files':{f'src/{name}.f90':digest(out/'src'/f'{name}.f90') for name in (*OWNERS,'main')}}
    manifest['changed_files']['Makefile']=digest(make)
    (out/'allocation-census-identity.json').write_text(json.dumps(manifest,indent=2)+'\n')
    verify(baseline,entries)
    print(json.dumps({k:len(v['allocations']) for k,v in coverage.items()},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('repo','baseline','out'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--ref',default='1b01510ae655eee0e2aab1a8d911e55c147173f2')
    a=p.parse_args();prepare(a.repo,a.baseline,a.out,a.ref)

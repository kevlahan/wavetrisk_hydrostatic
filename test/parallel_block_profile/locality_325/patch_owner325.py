"""Add only the audited, fixed-case L256 selector to the sealed Stage322 seam."""
from pathlib import Path
import hashlib
import json

BASE_SHA = '0ec557fc4220e6c1b7c09fe53c8af25602403bbd78fd4c7280e3fbde540ced36'
CURRENT_SHA = 'd5519da50e82252d5f7dede09bc99633056e29e982770af559f0c836597ed75a'


def maps(root):
    root = Path(root)
    a, b = [], []
    for name, expected, target in [('current256.json', CURRENT_SHA, a),
                                   ('owner_label_permutation256.json', BASE_SHA, b)]:
        raw = (root / name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != expected:
            raise ValueError('Unsealed map: ' + name)
        target.extend(json.loads(raw))
    if len(a) != 640 or len(b) != 640:
        raise ValueError('Map extent')
    p = {}
    for old, new in zip(a, b):
        if old in p and p[old] != new:
            raise ValueError('Changed coownership')
        p[old] = new
    if set(p) != set(range(256)) or set(p.values()) != set(range(256)):
        raise ValueError('Not a bijection')
    if sum(i // 4 == r for i, r in enumerate(b)) != 400:
        raise ValueError('Unexpected home count')
    return a, b, [p[r] for r in range(256)]


def replace_once(text, old, new):
    if text.count(old) != 1:
        raise ValueError('Source anchor changed: ' + old[:100])
    return text.replace(old, new, 1)


def apply(tree, map_dir=None):
    root = Path(map_dir) if map_dir else Path(__file__).resolve().parent
    a, b, permutation = maps(root)
    path = Path(tree) / 'src/parallel_block_mpi.f90'
    text = path.read_text()
    text = replace_once(text, "       case('C256');value=3\n", "       case('C256');value=3\n       case('L256');value=4\n")
    text = replace_once(text, "off, A160, B256, or C256'", "off, A160, B256, C256, or L256'")
    rows = [','.join(map(str, permutation[i:i+16])) for i in range(0,256,16)]
    array = '    ! Exact Stage324 permutation; never load an arbitrary owner file.\n'
    array += '    integer,parameter :: locality325(256)=[ &\n'
    array += ', &\n'.join('       ' + row for row in rows) + ' ]\n'
    text = replace_once(text, '    if(mode==3) then\n       r=native256(i)\n',
                        array + "    if(i<1.or.i>640.or.mode<1.or.mode>4) error stop 'Stage325 invalid owner lookup'\n"
                        '    if(mode==4) then\n       r=locality325(native256(i)+1)\n'
                        '    else if(mode==3) then\n       r=native256(i)\n')
    text = replace_once(text, '((mode==2.or.mode==3).and.n_process/=256)',
                        '((mode==2.or.mode==3.or.mode==4).and.n_process/=256)')
    anchor = '  subroutine check_owner_bridge322_catalog(require_bridge)\n'
    guard = '''  subroutine check_locality325_permutation()
    ! Once per construction: include this control work in the qualification ledger.
    ! Independent of input ordering; every native group has exactly one new owner.
    integer :: i,old,new,homes
    integer :: mapping(0:255),seen(0:255),counts(0:255)
    mapping=-1;seen=0;counts=0;homes=0
    do i=1,640
       old=selected_owner322(i,3);new=selected_owner322(i,4)
       if(old<0.or.old>255.or.new<0.or.new>255) error stop 'Stage325 owner outside communicator'
       if(mapping(old)/=-1.and.mapping(old)/=new) error stop 'Stage325 split numerical group'
       mapping(old)=new;counts(new)=counts(new)+1
       if(new==(i-1)/4) homes=homes+1
    end do
    if(any(mapping<0).or.any(counts<1).or.homes/=400) error stop 'Stage325 invalid locality map'
    do i=0,255
       seen(mapping(i))=seen(mapping(i))+1
    end do
    if(any(seen/=1).or.count(counts==3)/=192.or.count(counts==1)/=64) &
         error stop 'Stage325 owner permutation is not bijective'
  end subroutine

'''
    text = replace_once(text, anchor, guard + anchor)
    text = replace_once(text, '    call check_owner_bridge322_catalog(.false.)\n',
                        '    call check_owner_bridge322_catalog(.false.)\n'
                        '    if(mode==4) call check_locality325_permutation()\n')
    path.write_text(text)
    return dict(selector='L256', mode=4, current_map_sha256=CURRENT_SHA,
                candidate_map_sha256=BASE_SHA, home_blocks=400, numerical_owners=256)


if __name__ == '__main__':
    import sys
    print(json.dumps(apply(sys.argv[1]), indent=2))

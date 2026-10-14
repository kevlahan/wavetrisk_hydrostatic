"""Extract the accepted step1 address traversal verbatim, without its arithmetic."""
import re

def source(ops):
    body=ops.split('    id = offs(0)\n',1)[1].split('  contains',1)[0]
    body='    id = offs(0)\n'+body
    body=re.sub(r'dom%patch%elts\(p\+1\)%neigh\((NORTH|EAST)\) < 0',r'boundary(\1)',body)
    body=re.sub(r'dom%bdry_patch%elts\(-dom%patch%elts\(p\+1\)%neigh\((NORTH|EAST)\)\+1\)%side > 0',r'boundary(\1)',body)
    body=re.sub(r'dom%bdry_patch%elts\(-dom%patch%elts\(p\+1\)%neigh\((NORTH|EAST)\)\+1\)%side <= 0',r'.not.boundary(\1)',body)
    assert 'dom%' not in body
    offsets=ops.split('    offs(1:N_BDRY)  =',1)[1].split('  end subroutine comp_offs3',1)[0]
    offsets='    offs(1:N_BDRY)  ='+offsets
    return """module rhs_schedule_oracle_mod
 use shared_mod, only : N_BDRY,NORTH,EAST,SOUTH,WEST,NORTHEAST,NORTHWEST,SOUTHEAST,SOUTHWEST
 use patch_mod, only : PATCH_SIZE
 implicit none
 contains
 subroutine legacy_trace(absolute,dims,boundary,events,sw_events,nout)
 integer, intent(in) :: absolute(0:N_BDRY),dims(2,N_BDRY)
 logical, intent(in) :: boundary(4)
 integer, intent(out) :: events(:,:),nout
 logical, intent(out) :: sw_events(:)
 integer, parameter :: LAST=PATCH_SIZE-1
 integer :: offs(0:N_BDRY),id,n,e,s,w,ne,sw,j
 logical :: S_bdry,W_bdry
 offs=absolute;nout=0;events=0;sw_events=.false.
 S_bdry=boundary(SOUTH);W_bdry=boundary(WEST)
"""+offsets+body+"""
 contains
 subroutine comput
 nout=nout+1
 if (nout>size(events,2)) error stop 'legacy trace output extent'
 events(:,nout)=[id,id+n,id+e,id+s,id+w,id+ne,id+sw]+1
 end subroutine
 subroutine comp_SW
 sw_events(nout)=.true.
 end subroutine
 end subroutine
end module
"""

"""Compile actual candidate and baseline assignment code at O0/check and O2."""
import hashlib,json,random,re,subprocess,sys,tarfile,tempfile
from pathlib import Path
from common import save,sha,verify
from partition_evidence import balanced,legacy


def check(kit,out):
    kit=Path(kit).resolve(strict=True);out=Path(out).resolve();out.mkdir(exist_ok=False);verify(kit)
    with tarfile.open(kit/'sources.tar.gz') as t:
        texts={kind:t.extractfile(kind+'/src/parallel_block_mpi.f90').read().decode() for kind in ('stage266','candidate')}
    new=re.search(r'  subroutine assign_balanced_block_owners\b.*?end subroutine assign_balanced_block_owners',texts['candidate'],re.S).group()
    old=texts['stage266'];a=old.index('  n_assigned     = 0\n');b=old.index('  if (any(block_owner_new < 0))',a)
    legacy_code='''subroutine old_assign(weights,n_process,block_owner_new)
 integer,intent(in)::weights(:),n_process
 integer,intent(out)::block_owner_new(:)
 type block
  integer::weight
 end type
 type(block)::block_catalog(size(weights))
 integer::n_block_total,n_assigned,r,load_proposed(n_process)
 real(dp)::balanced_weight,imbalance_goal
 real(dp),parameter::init_goal=.05_dp,incr_goal=1.2_dp
 n_block_total=size(weights)
 block_catalog%weight=weights
 balanced_weight=real(sum(weights),dp)/n_process
'''+old[a:b]+'''end subroutine
'''
    program='''module actual
 use iso_fortran_env,only:int64,real64
 implicit none
 integer,parameter::dp=real64
 contains
'''+new+'\n'+legacy_code+'''end module
program driver
 use actual
 implicit none
 integer::cases,c,n,p,mode
 integer,allocatable::weights(:),owners(:)
 read(*,*) cases,mode
 do c=1,cases
  read(*,*) n,p
  allocate(weights(n),owners(n))
  if(n>0) read(*,*) weights
  if(mode==0) then
   call assign_balanced_block_owners(weights,p,owners)
  else
   call old_assign(weights,p,owners)
  end if
  write(*,'(i0,1x,*(i0,1x))') c,owners
  deallocate(weights,owners)
 end do
end program
'''
    source=out/'actual.f90';source.write_text(program)
    cases=[([],4),([0]*9,4),([1]*640,256),([1]*640,640),([1]*640,160),([1]*640,83),([1,3],5),([2147483647]*7,3),([1000000]+[1]*31,8),([1]*31+[1000000],8),([0,0,1,0,0,9,0,1],4),([2,9,1,5,8,2,1,11],3)]
    rng=random.Random(273)
    for i in range(160):
        n=rng.randrange(1,85);cases.append(([rng.randrange(0,1000) for _ in range(n)],rng.randrange(1,100)))
    cases.extend([([1]*n,p) for n in range(1,20) for p in (1,2,3,7)])
    results=[]
    for opt in ('-O0','-O2'):
        binary=out/('partition-'+opt[1:])
        subprocess.run(['gfortran',opt,'-fcheck=all','-ffree-line-length-132',str(source),'-o',str(binary)],cwd=out,check=True,capture_output=True)
        for mode in (0,1):
            batch=cases if mode==0 else [(w,p) for w,p in cases if len(w)>=p and all(x>0 for x in w) and sum(w)<2147483647]
            data=f'{len(batch)} {mode}\n'+''.join(f'{len(w)} {p}\n'+(' '.join(map(str,w))+'\n' if w else '') for w,p in batch)
            proc=subprocess.run([str(binary)],input=data,text=True,capture_output=True,check=True,timeout=30)
            rows=proc.stdout.splitlines();assert len(rows)==len(batch)
            for index,((weights,ranks),line) in enumerate(zip(batch,rows),1):
                got=list(map(int,line.split()));assert got[0]==index
                expected=(balanced if mode==0 else legacy)(weights,ranks)
                assert got[1:]==expected,(opt,mode,weights,ranks,got,expected)
            results.append(dict(optimization=opt,mode='candidate' if mode==0 else 'baseline',cases=len(batch),passed=True))
        for data in ('1 0\n1 0\n1\n','1 0\n2 1\n1 -1\n'):
            proc=subprocess.run([str(binary)],input=data,text=True,capture_output=True,timeout=5)
            assert proc.returncode!=0
    uniform=balanced([1]*640,256);old=legacy([1]*640,256)
    from collections import Counter
    record=dict(status='passed',source_manifest_sha256=sha(kit/'source-files.json'),extracted_source_sha256=sha(source),cases=results,invalid_cases=4,
        equal_weight_candidate=dict(Counter(Counter(uniform).values())),equal_weight_baseline=dict(Counter(Counter(old).values())),
        limitation='Assignment correctness and work distribution only, not a solver speedup or actual J5 catalogue replay.')
    save(out/'partition-tests.json',record);print(json.dumps(record,indent=2));return record

if __name__=='__main__':check(*sys.argv[1:])

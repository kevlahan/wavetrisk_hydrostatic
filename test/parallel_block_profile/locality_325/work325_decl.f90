  ! Stage325 qualification ledger: counters only, never a timing instrument.
  public :: work325_enabled,work325_record
  integer,parameter :: work325_limit=4096,work325_width=24
  integer :: work325_n=0,work325_rank=-1
  logical :: work325_on=.false.
  type :: Work325_Entry
     character(48) :: tag=''
     integer :: phase=0,width=0
     integer(int64) :: calls=0,total(work325_width)=0,peak(work325_width)=0
     integer(int64),allocatable :: sent(:),received(:),send_peak(:),receive_peak(:)
  end type
  type(Work325_Entry),allocatable :: work325(:)

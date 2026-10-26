"""Three native correctness/work launches; performance scheduling is absent."""
from pathlib import Path
import json,os
STAGE=325
REVISION='locality-qualification-v1'
STEP_RESTRICTIONS=('SLURM_CPUS_PER_TASK','SRUN_CPUS_PER_TASK','SLURM_EXACT','SLURM_EXCLUSIVE','SLURM_TRES_PER_TASK','SLURM_HINT','SLURM_THREADS_PER_CORE')
SCHEDULE=(
 dict(label='checked-L256',kind='checked',owner='L256',ranks=256,steps=4),
 dict(label='qual-C0',kind='candidate',owner='off',ranks=256,steps=4),
 dict(label='qual-L256',kind='candidate',owner='L256',ranks=256,steps=4),
)

def environment(kit,kind,owner):
 if kind not in ('checked','candidate') or owner not in ('off','L256'):raise ValueError('Unsupported Stage325 arm')
 if kind=='checked' and owner!='L256':raise ValueError('Only candidate checked arm is sealed')
 env={k:v for k,v in os.environ.items() if not k.startswith(('WAVETRISK_','MV2_','OMPI_','GFORTRAN_')) and not k.endswith('NUM_THREADS')}
 env.update(json.loads((Path(kit)/(kind+'-environment.json')).read_text()))
 env.update(WAVETRISK_CAMPAIGN312_STEPS='4',WAVETRISK_SELECTIVE_BALANCE='1',WAVETRISK_SUBDIVISION_TARGET_TASKS='256',
 WAVETRISK_OWNER_MAP322=owner,WAVETRISK_EXECUTION_AUDIT316='0',WAVETRISK_WORK_LEDGER325='1',
 WAVETRISK_PLACEMENT322='1',WAVETRISK_WINDOW322='0',WAVETRISK_SAMPLE_RANKS201='[]',
 OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
 # The semantic checker owns its explicitly bounded capture selector. It may not
 # turn on the retired unbounded Stage316 value stream as an implicit fallback.
 from verify325 import capture_environment
 knobs=capture_environment(kind,owner)
 if not isinstance(knobs,dict) or any(not k.startswith('WAVETRISK_') for k in knobs):raise ValueError('Invalid semantic capture environment')
 env.update(knobs)
 if env.get('WAVETRISK_WINDOW322')!='0' or env.get('WAVETRISK_EXECUTION_AUDIT316')!='1' or env.get('WAVETRISK_COMPACT_AUDIT325')!='1':raise ValueError('Require compact semantic mode; timing/full real-value capture forbidden')
 return env

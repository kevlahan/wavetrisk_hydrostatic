"""Reviewed lifecycle dependencies anchored to source files and named routines."""
import json,re,subprocess,sys
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'test/parallel_block_profile'))
from capability_223 import digest,save
ROWS=[
 ('RK publication','src/parallel_block_mpi.f90',['prepare_parallel_block_grid_change'],
  'Final owners -> explicit Domain physical/alias image; retained atmospheric wavelets remain pending.',
  'Preserve the distinct saved final-consumer alias epoch; cannot delete publication while physics and adaptation consume it.'),
 ('Simple physics','src/physics/physics_simple.f90',['physics_simple_step','physics_call'],
  'Reads atmospheric sol, means, topography, geometry, mask and soil; writes velocity, temperature and soil in Domain traversal order.',
  'First bounded migration removes Domain hydrostatic scratch reads/writes. Owner dispatch still needs versioned kinetic-energy neighbors and validation of the package state/first-call contract.'),
 ('Physics initialization','src/physics/init_physics_module.f90',['init_physics'],
  'Package parameters and extra levels are initialized on every rank before stepping, including zero-Domain ranks.',
  'Preserve first-call and per-column package state when changing execution ownership.'),
 ('Kinetic energy ordering','src/diagnostics.f90',['kinetic_energy'],
  'Reads current and west/southwest/south velocity stencils; preceding physics columns can already have updated them.',
  'Compile read-after-write dependencies or establish a separately validated simultaneous-update formulation before parallel dispatch.'),
 ('Other physics','src/main.f90',['apply_domain_physics_consumers'],
  'Vertical diffusion and Held-Suarez remain distinct optional Domain consumers; Simple is the qualification target.',
  'Migrate these only with fixtures exercising their actual branches.'),
 ('Soil wavelets','src/main.f90',['update_domain_soil_wavelet_consumer'],
  'Simple physics writes sol temperature zmin:0; legacy WT writes soil wav_coeff; atmospheric pending decision state must remain unchanged.',
  'Move soil transform operands, coarse restriction and signed aliases together. Surface/soil fields are not atmospheric RK scratch.'),
 ('Adaptation and rebuild','src/main.f90',['time_step'],
  'adapt mutates geometry and wavelets, legacy inverse reconstructs physical fields, build_source_blocks imports the new topology.',
  'Separate metadata mutation from numerical prolongation and restrict/import lifetimes before deleting Domain arrays.'),
 ('Restart','src/main.f90',['restart'],
  'Checkpoint reload supplies Domain physical/reference state before block construction and migration.',
  'Introduce bounded I/O-to-owner import including poles, means, masks and soil; preserve exact reload behavior.'),
 ('Output and checkpoint','src/main.f90',['write_output_consumer','prepare_end_step_domain_consumers'],
  'Audited ready-state consumers serialize the compatibility image.',
  'Replace persistent numerical mirrors with bounded serialization adapters after remaining mutators migrate.'),
]
def audit():
 rows=[]
 for family,file,names,state,next_step in ROWS:
  s=(BASE/file).read_text();anchors=[]
  for name in names:
   m=re.search(r'^\s*(?:subroutine|function)\s+'+name+r'\b',s,re.M|re.I)
   assert m,(file,name)
   anchors.append(dict(routine=name,line=s[:m.start()].count('\n')+1))
  rows.append(dict(family=family,file=file,sha256=digest(BASE/file),anchors=anchors,ownership_and_lifetime=state,next_migration=next_step))
 save(BASE/'docs/parallel-blocks-stage237-dependencies.json',dict(stage=237,checkpoint_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=BASE,text=True).strip(),
  classification='Manual dependency and ordering audit, source anchored; not measured performance attribution',dependencies=rows,
  first_boundary='Simple physics column-local hydrostatics',complete_physics_owner_migration=False,performance_claim=None))
 print('PASS lifecycle audit:',len(rows),'families')
if __name__=='__main__':audit()

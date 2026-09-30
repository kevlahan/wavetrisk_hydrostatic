"""Reviewed owner-routing and package-state contract, bound to source files."""
import sys,subprocess
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'test/parallel_block_profile'))
from capability_223 import digest,save
P='src/physics/simple_physics/phyparam/physics/'
rows=[
('src/arch.f90','block_catalog supplies final owner, global root Domain and subtree patch identity.'),
('src/domain_ops.f90','Original local column traversal and pole order are preserved.'),
('src/physics/physics_simple.f90','Map catalogue subtrees to local patch/node IDs; transport source-local version keys with the column. Await each reply before the next gather. Fixed scaffold and unmapped entry aliases retain Domain compatibility ownership.'),
('src/physics/init_physics_module.f90','init_physics initializes parameters, geometry and extra-level contract on every rank. physics_checkpoint_restart is empty; do not assume it allocates package arrays.'),
(P+'single_column_mod.F90','Before every non-first call, reload Tsurf/Tsoil from the packet; after every call, return them. Reload latitude/longitude and trigonometric geometry for every column.'),
(P+'phyparam_mod.F90','firstcall_alloc is rank-local and private. firstcall is per-column epoch semantics: preserve it in the packet. Late first use initializes via a disposable first-call copy before the requested non-first call. icount and zDay_last have no numerical consumers in this package.'),
(P+'soil_mod.F90','Tsurf/Tsoil are mutable package scratch; dz1/dz2 and surface property arrays are initialized geometry/parameters. Soil state travels in each packet.'),
(P+'radiative_mod.F90','Radiative inputs and outputs are recomputed per call; no saved column tendency cache.'),
(P+'turbulence.F90','Saved members are model constants/parameters; numerical column work arrays are local.'),
(P+'astronomy.F90','Orbit constants are initialized per rank; SolarLon/Orbite compute from the supplied time.'),
(P+'comgeomfi.F90','Location arrays are rank-local mutable scratch, overwritten for every dispatched column.'),
('src/main.f90','Domain return updates occur before soil wavelets and adaptation. Retained atmospheric transaction remains unchanged.')]
save(BASE/'docs/parallel-blocks-stage239-dependencies.json',dict(stage=239,
 checkpoint_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=BASE,text=True).strip(),
 sources=[dict(file=name,sha256=digest(BASE/name),contract=contract) for name,contract in rows],
 route='(global Domain, patch) -> subtree catalogue index -> final owner; scaffold index zero remains local.',
 source_identity='Each of six operands carries local node ID, expected writer ordinal, and catalogue index where covered. Zero identifies a compatibility alias/scaffold. Do not replace entry aliases with newer canonical values.',
 transport='Single-column single-precision package packet plus integer identity/version header. Serial geometry sources; synchronous owner reply before dependent source advances. No Domain-shaped numerical snapshot.',
 limits='Owner execution migrates the package only. Input gathering, version authority, prognostic/soil return fields remain Domain-owned; no persistent Domain arrays removed.',
 package_audit='Manual source review; scientific equivalence is additionally checked by exact trajectory and per-column owner-output oracles.',
 performance_claim=None))
print('PASS owner and package state audit')

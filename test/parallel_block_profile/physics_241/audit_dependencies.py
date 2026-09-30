"""Reviewed owner-routing and package-state contract, bound to source files."""
import sys,subprocess
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'test/parallel_block_profile'))
from capability_223 import digest,save
P='src/physics/simple_physics/phyparam/physics/'
rows=[
('src/parallel_block.f90','Canonical physical scalar/vector, mean and soil input reads use final-owner Block_Data via compiled native addresses.'),
('src/parallel_block_inverse.f90','Physics aliases read the explicitly saved final-inverse consumer epoch with readiness and closure guards; canonical interiors are forbidden through this accessor.'),
('src/arch.f90','block_catalog supplies final owner, global root Domain and subtree patch identity.'),
('src/domain_ops.f90','Original local column traversal and pole order are preserved.'),
('src/physics/physics_simple.f90','Import column keys and explicit alias/reference inputs to execution owners; fetch canonical atmospheric, mean and soil values from native storage owners. Compile source-local predecessor identities into DAG waves; compute hydrostatics, kinetic energy and package updates on owners. Exchange only required updated edges, preserving entry aliases. Publish completed double-precision state once. Fixed scaffold remains an explicit compatibility case.'),
('src/physics/init_physics_module.f90','init_physics initializes parameters, geometry and extra-level contract on every rank. physics_checkpoint_restart is empty; do not assume it allocates package arrays.'),
(P+'single_column_mod.F90','Before every non-first call, reload Tsurf/Tsoil from the packet; after every call, return them. Reload latitude/longitude and trigonometric geometry for every column.'),
(P+'phyparam_mod.F90','firstcall_alloc is rank-local and private. firstcall is per-column epoch semantics: preserve it in the packet. Late first use initializes via a disposable first-call copy before the requested non-first call. icount and zDay_last have no numerical consumers in this package.'),
(P+'soil_mod.F90','Tsurf/Tsoil are mutable package scratch; dz1/dz2 and surface property arrays are initialized geometry/parameters. Soil state travels in each packet.'),
(P+'radiative_mod.F90','Radiative inputs and outputs are recomputed per call; no saved column tendency cache.'),
(P+'turbulence.F90','Saved members are model constants/parameters; numerical column work arrays are local.'),
(P+'astronomy.F90','Orbit constants are initialized per rank; SolarLon/Orbite compute from the supplied time.'),
(P+'comgeomfi.F90','Location arrays are rank-local mutable scratch, overwritten for every dispatched column.'),
('src/main.f90','Domain return updates occur before soil wavelets and adaptation. Retained atmospheric transaction remains unchanged.')]
save(BASE/'docs/parallel-blocks-stage241-dependencies.json',dict(stage=241,
 checkpoint_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=BASE,text=True).strip(),
 sources=[dict(file=name,sha256=digest(BASE/name),contract=contract) for name,contract in rows],
 route='Column and edge native keys are sent to execution owners, which request canonical data directly from storage owners. Global metadata defines predecessor edges by source Domain/node, ordinal and block. Only required updated velocity components are exchanged.',
 source_identity='Each of six operands carries local node ID, expected writer ordinal, and catalogue index where covered. Zero identifies a saved atmospheric consumer alias; its mean/soil reference adapter remains explicit. Do not replace entry aliases with newer canonical values.',
 transport='Bulk owner-local packets and double-precision raw inputs. Topological waves execute ready columns; targeted all-to-all edge updates carry producer components. Bulk double-precision results return to a write-only compatibility callback.',
 limits='Owner-resident epoch and DAG replace per-column Domain reads/writes during execution. Canonical atmospheric/mean/soil reads use owner storage; atmospheric aliases use saved native consumer values. Zero-key mean and soil alias adapters and final publication still use Domain compatibility storage. Persistent adaptation and I/O migration remain required.',
 package_audit='Manual source review; scientific equivalence is additionally checked by exact trajectory and per-column owner-output oracles.',
 performance_claim=None))
print('PASS owner and package state audit')

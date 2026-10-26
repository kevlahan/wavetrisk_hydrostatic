"""Closed, source-bound Stage322 level5 coverage.  This is NOT merely pair equality.
Expected schedules are derived from configured source loops; unknown records fail.
"""
"""Source-grounded RHS lifecycle and actual halo-order coverage (qualification only)."""
import math
import struct
import json
from pathlib import Path
PERMUTATION325=json.loads(Path(__file__).with_name("map325.json").read_text())["permutation"]
def owner325(i,arm):
 old=i//3 if i<576 else 192+i-576
 return PERMUTATION325[old] if arm in ("L256","checked-L256","qual-L256") else old
def blocks325(rank,arm):
 return sum(owner325(i,arm)==rank for i in range(640))


FAMILIES = ('gradient', 'pressure.patch', 'diffusion.producer', 'diffusion.restriction',
            'diffusion.divergence', 'rotation.curl', 'rotation.pentagon',
            'dynamics.divergence', 'dynamics.physics', 'dynamics.source.direct',
            'dynamics.source.restriction')
REAL_FAMILIES = set(FAMILIES)-{'pressure.patch','dynamics.source.restriction'}
LOGICAL_FAMILIES = set(FAMILIES)-{'pressure.patch','diffusion.divergence','rotation.pentagon'}
RHS_PROGRAM_FIELDS = {'program.header','program.surface_visits','program.level','pressure.node',
    'pressure.producer_ordinal','dynamics.injection.shape','dynamics.injection',
    'dynamics.source.action','dynamics.source.operand'}
for family in FAMILIES:
    RHS_PROGRAM_FIELDS.update((family+'.shape',family+'.i'))
    if family in REAL_FAMILIES: RHS_PROGRAM_FIELDS.add(family+'.r')
    if family in LOGICAL_FAMILIES: RHS_PROGRAM_FIELDS.add(family+'.l')
RHS_REQUIRED_FIELDS = {
    'rhs.prepare_check': {'prepare_check.generation','prepare_check.rebuild'},
    'rhs.prepare': RHS_PROGRAM_FIELDS | {'backend.generation','backend.geometry.base',
        'placement.key.shape','placement.key','placement.base','placement.owned',
        'transport.epoch','transport.route.bounds','transport.bound','transport.route_phase',
        'transport.route_level','transport.scratch.shape','backend.adapter_domain',
        'backend.adapter_node','backend.rk_domain','backend.rk_node','backend.owned',
        'backend.mask.domains','backend.surface','retention.blocks',
        'executor.scalar_rhs.shape','executor.velocity_rhs.shape'},
    'rhs.evaluate': {'evaluate.epoch','evaluate.first_after_prepare','execution.shape',
        'execution.order','execution.constants','execution.nu.shape','execution.nu',
        'execution.nu_pentagon','execution.state.begin','execution.state.end',
        'execution.halo.shape','execution.halo.order'},
    'rhs.prepared_state': {'state.scalar.shape','state.scalar','state.mean.shape','state.mean',
        'state.velocity.shape','state.velocity','state.surface'},
    'rhs.retain': {'retain.epoch'},
    'rhs.rk_prepare': {'rk_prepare.epoch','rk_prepare.h'},
    'rhs.rk_install': {'rk_install.epoch','rk_install.validate'},
}
for field in ('lower phi pressure midpoint exner thermodynamic flux pv ke bernoulli vorticity '
              'diffusion_flux laplacian divu lap_divu rotation source physics scalar_rhs '
              'velocity_rhs scalar_physics sw_physics').split():
    RHS_REQUIRED_FIELDS['rhs.prepare'].add('workspace.'+field+'.shape')
for field in 'scalar mean velocity adapter_value rk_value rk_scalar rk_mean rk_velocity'.split():
    RHS_REQUIRED_FIELDS['rhs.prepare'].add('backend.'+field+'.shape')
RHS_CAPTURE_FIELDS = {'program.header','backend.generation','prepare_check.generation','prepare_check.rebuild',
    'evaluate.epoch','evaluate.first_after_prepare','execution.shape','execution.order',
    'execution.state.begin','execution.state.end','execution.halo.shape','execution.halo.order',
    'state.scalar.shape','state.mean.shape','state.velocity.shape','retain.epoch',
    'rk_prepare.epoch','rk_install.epoch','rk_install.validate'}

def _ints(scope,name,kind=None):
    e=scope['values'].get(name)
    if e is None: raise ValueError('Missing captured RHS field: '+name)
    if kind is not None and e.kind!=kind: raise ValueError('Wrong RHS field kind: '+name)
    if e.kind not in (1,2,4): raise ValueError('Wrong RHS integer/logical kind: '+name)
    code={1:'i',2:'q',4:'B'}[e.kind]
    return list(struct.unpack('<'+str(e.count)+code,e.data))

def _descriptor(scope,name):
    records=scope.get('record_descriptors',{})
    value=records.get(name)
    if value is not None: return tuple(value)
    e=scope['values'].get(name)
    if e is not None: return e.kind,e.count
    raise ValueError('Missing RHS payload descriptor: '+name)

def expected_halos(lo,hi,nlayer,order):
    """Exactly mirrors the sealed call schedule, not numerical calculations."""
    if not (0<=lo<=hi<=100 and 1<=nlayer<=1000 and len(order)==3 and all(v in (0,1,2) for v in order)):
        raise ValueError('Invalid RHS schedule dimensions/orders')
    result=[]
    def emit(phase,level,columns,k,repeated=0):
        physical_phase=16 if phase==18 else phase
        physical_level=lo if phase in (16,18,15,7) else level
        result.extend((phase,level,columns,physical_phase,physical_level,repeated,k))
    for k in range(1,nlayer+1):
        if order[1]!=0:
            for level in range(hi,lo-1,-1):
                for phase in (5,1,7): emit(phase,level,1,k)
        if order[0]==2:
            for level in range(hi,lo-1,-1):
                for phase in (5,1,2): emit(phase,level,2,k)
        if order[1]==2:
            for level in range(hi,lo-1,-1):
                for phase in (5,1,2): emit(phase,level,1,k)
        for level in range(hi,lo-1,-1):
            for phase,columns in ((11,1),(12,1),(4,1),(16,2),(18,2),(20,2),(13,2)):
                emit(phase,level,columns,k)
            emit(13,level,1,k,1)
            if level>lo: emit(14,level,2,k)
            if order[2]==2:
                for phase in (3,6,4): emit(phase,level,1,k)
            emit(15,level,1,k)
    return result

def validate_rhs322(scopes,steps,coverage='candidate'):
    """scopes: begin-order dicts with direct tags, captured values and descriptors."""
    if coverage not in ('candidate','checked'): raise ValueError('Unknown RHS coverage variant')
    selected={name:[s for s in scopes if s['name']==name] for name in RHS_REQUIRED_FIELDS}
    for name,items in selected.items():
        for scope in items:
            if not RHS_REQUIRED_FIELDS[name] <= scope['tags']:
                raise ValueError('Incomplete RHS scope '+name+': '+','.join(sorted(RHS_REQUIRED_FIELDS[name]-scope['tags'])))
    evaluations=selected['rhs.evaluate'];preparations=selected['rhs.prepare']
    if len(evaluations)!=4*steps: raise ValueError('Require four RHS evaluations per timestep')
    for name in ('rhs.retain','rhs.rk_prepare','rhs.rk_install','rhs.prepare_check'):
        if len(selected[name])!=len(evaluations): raise ValueError('Incomplete RHS lifecycle: '+name)
    if [s['epoch'] for s in preparations]!=list(range(1,len(preparations)+1)):
        raise ValueError('Missing/reordered RHS preparation epoch')
    if [s['epoch'] for s in evaluations]!=list(range(1,len(evaluations)+1)):
        raise ValueError('Missing/reordered RHS evaluation epoch')
    if evaluations and not preparations: raise ValueError('RHS execution without a prepared program')
    if len(selected['rhs.prepared_state'])!=len(evaluations):
        raise ValueError('Require one installed input state at EVERY RHS evaluation')
    positions={id(s):i for i,s in enumerate(scopes)}
    used=set()
    for ordinal,evaluation in enumerate(evaluations,1):
        earlier=[s for s in preparations if positions[id(s)]<positions[id(evaluation)]]
        if not earlier: raise ValueError('RHS evaluation precedes preparation')
        prepare=earlier[-1];epoch=prepare['epoch']
        generation=_ints(prepare,'backend.generation',2)
        if len(generation)!=4 or generation[1]!=epoch:
            raise ValueError('RHS preparation epoch/header differs')
        ev=_ints(evaluation,'evaluate.epoch',2)
        if ev!=[generation[0],epoch,ordinal,ordinal-1,ordinal-1]:
            raise ValueError('RHS evaluation generation/mask/state/RK linkage differs')
        first=epoch not in used
        check=selected['rhs.prepare_check'][ordinal-1]
        flags=_ints(check,'prepare_check.rebuild',4)
        before=_ints(check,'prepare_check.generation',2)
        if check['epoch']!=generation[0] or check.get('step')!=evaluation.get('step'):
            raise ValueError('RHS preparation check belongs to wrong topology/step')
        if len(flags)!=2 or flags[0]>flags[1] or flags[1]!=int(first):
            raise ValueError('RHS collective rebuild decision differs from actual preparation')
        if len(before)!=2 or before[1]!=(epoch-1 if first else epoch):
            raise ValueError('RHS preparation check mask history differs')
        if positions[id(check)]>=positions[id(evaluation)]:
            raise ValueError('RHS preparation check follows evaluation')
        if first and positions[id(check)]>=positions[id(prepare)]:
            raise ValueError('RHS preparation precedes its rebuild check')
        if _ints(evaluation,'evaluate.first_after_prepare',4)!=[int(first)]:
            raise ValueError('RHS first-after-prepare marker differs')
        children=[s for s in selected['rhs.prepared_state'] if s.get('parent') and
                  tuple(s['parent'][-1])==('rhs.evaluate',ordinal)]
        if len(children)!=1: raise ValueError('Prepared RHS state missing/duplicated/misnested')
        program=_ints(prepare,'program.header',1)
        dims=_ints(evaluation,'execution.shape',1)
        if len(program)!=4 or len(dims)!=2 or dims[0]!=program[1] or dims[0]<0:
            raise ValueError('RHS program/execution dimensions differ')
        n,nlayer=dims;lo,hi=program[-2:]
        order=_ints(evaluation,'execution.order',1)
        if (lo,hi,nlayer,order)!=(5,5,30,[2,2,2]):
            raise ValueError('RHS schedule differs from sealed level5 uniform campaign')
        counts=prepare.get('record_counts',{})
        for name in RHS_REQUIRED_FIELDS['rhs.prepare']:
            expected_count=1
            if counts.get(name,0)!=expected_count:
                raise ValueError('Incomplete/duplicate per-level RHS program record: '+name)
        expected=expected_halos(lo,hi,nlayer,order)
        if _ints(evaluation,'execution.halo.shape',1)!=[7,len(expected)//7]:
            raise ValueError('RHS actual halo dimensions differ from sealed schedule')
        if _ints(evaluation,'execution.halo.order',1)!=expected:
            raise ValueError('RHS actual halo order/width/binding/repeat/layer differs from sealed schedule')
        begin=_ints(evaluation,'execution.state.begin',2);end=_ints(evaluation,'execution.state.end',2)
        if len(begin)!=1 or len(end)!=1 or end[0]-begin[0]!=len(expected)//7:
            raise ValueError('RHS halo state delta differs from actual exchange count')
        if children:
            state=children[0]
            if state['epoch']!=epoch or state.get('step')!=evaluation.get('step'):
                raise ValueError('Prepared RHS input state belongs to wrong epoch/step')
            for name,shape in [('scalar',[n,2,nlayer]),('mean',[n,2,nlayer]),('velocity',[3*n,nlayer])]:
                if _ints(state,'state.'+name+'.shape',1)!=shape:
                    raise ValueError('Prepared RHS input shape differs: '+name)
                if _descriptor(state,'state.'+name)!=(3,math.prod(shape)):
                    raise ValueError('Prepared RHS input payload extent differs: '+name)
            if _descriptor(state,'state.surface')!=(3,n): raise ValueError('Prepared surface extent differs')
            used.add(epoch)
        for name,field in [('rhs.retain','retain.epoch'),('rhs.rk_prepare','rk_prepare.epoch'),('rhs.rk_install','rk_install.epoch')]:
            scope=selected[name][ordinal-1];values=_ints(scope,field,2)
            if scope['epoch']!=ordinal or scope.get('step')!=evaluation.get('step'):
                raise ValueError('RHS lifecycle scope linked to wrong evaluation: '+name)
            if len(values)!=4 or values[:3]!=ev[:3]: raise ValueError('RHS lifecycle epoch fields differ: '+name)
            if name=='rhs.rk_prepare' and values[3]!=ordinal-1: raise ValueError('RK preparation counter differs')
            if name=='rhs.rk_install':
                if values[3]!=ordinal: raise ValueError('RK installation counter differs')
                if _ints(scope,'rk_install.validate',4)!=[int(coverage=='checked')]:
                    raise ValueError('RHS oracle validation selector differs from sealed coverage variant')
            if positions[id(scope)]<=positions[id(evaluation)]: raise ValueError('RHS lifecycle precedes evaluation')
        chain=[evaluation,selected['rhs.retain'][ordinal-1],selected['rhs.rk_prepare'][ordinal-1],selected['rhs.rk_install'][ordinal-1]]
        if [positions[id(s)] for s in chain]!=sorted(positions[id(s)] for s in chain):
            raise ValueError('RHS retain/RK lifecycle reordered')
    if used!=set(range(1,len(preparations)+1)):
        raise ValueError('Prepared RHS epoch never supplied installed state to evaluation')

"""Sealed level5 RK4 transform/inverse coverage, independent of pair equality.

Required per-step final inverse counts follow prepare_native_inverse's four
initial gathers/boundaries/scatters, one initial scalar synchronization, and
(level_end-level_start) transitions: four gathers/scatters, five boundaries,
two scalar kernels, one inner kernel and one outer kernel per transition.
The three provisional RK transforms use the separate resident backend.
"""
import struct
from collections import Counter

TRANSFORM_REQUIRED_FIELDS = {
 'transform.plan': {'masks','extent','placement%key.shape','placement%key','placement%owner',
  'placement%state_key','placement%level','placement%mask','placement%owned',
  'placement%scalar_order','placement%inner_order','placement%outer_order',
  'info','interior_address','boundary_address','save_node','save_index','fixed_node','fixed_index',
  'interior_node','publication_domain','publication_node','publication_row','forward_inner','base',
  'fixed_domain','fixed_physical_node','scalar.program.shape','scalar.program','scalar.coefficients',
  'inner.program.shape','inner.program','inner.coefficients','outer.program.shape','outer.program',
  'outer.coefficients','scalar.shape','vector.shape','scalar_wavelet.shape','vector_wavelet.shape',
  'fixed_scalar.shape','fixed_vector.shape','wire.shape','exchange_scratch.shape'},
 'transform.execute': {'masks','commit'},
 'transform.level': {'scalar.boundary.key','scalar.boundary.sign','vector.boundary.key','vector.boundary.sign'},
 'transform.scalar_adapter': {'extent','plan%sc','plan%sd','plan%rc','plan%rd','plan%source',
  'plan%index','plan%target','plan%sign','plan%send_peer','plan%receive_peer','plan%sent.shape',
  'plan%received.shape','plan%request.shape'},
 'inverse.plan': {'extent','options','node_domain','node_id','operation_first','operation_last',
  'scaffold_slot','consumer_alias_slot','node_flags','restore203_index','restore203_slot',
  'outer.program.shape','outer.program','outer.coefficients','scalar.program.shape','scalar.program',
  'inner.program.shape','inner.program','publication.program.shape','publication.program',
  'value.shape','consumer_alias_value.shape','scaffold_value.shape','restore203.shape','saved203.shape',
  'publication_tile.shape','transfer_send.shape','transfer_recv.shape','alias_send.shape','alias_recv.shape',
  'requests.shape'},
 'inverse.prepare': {'first_level'},
 'inverse.begin': {'options'},
 'inverse.end': {'work'},
 'inverse.node_map': {'slot','key.shape','key'},
 'inverse.level': {'component'},
 'inverse.outer_level': {'component'},
 'inverse.local_alias': {'source','dest','range'},
 'inverse.scalar': {'level_phase','selected_and_active','coefficients'},
 'inverse.inner': {'level','selected','coefficients'},
 'inverse.outer': {'level','operation_range'},
 'inverse.gather': {'family','flags'},
 'inverse.scatter': {'family','flags'},
 'inverse.boundary': {'family_and_range'},
}
TRANSFORM_REQUIRED_FIELDS['transform.vector_adapter']=TRANSFORM_REQUIRED_FIELDS['transform.scalar_adapter']
for label in ('interior','boundary','level_interior','level_boundary','outer_install'):
 TRANSFORM_REQUIRED_FIELDS['inverse.'+label] = {'edgewise','plan%sc','plan%sd','plan%rc','plan%rd','plan%key',
  'plan%slot','plan%address','plan%source_info.shape','plan%target_info.shape','plan%self_slot.shape'}
# 'inverse.boundary' is also an invocation scope. Required fields differ there
# and are selected by presence of family_and_range in required_fields().
for label in ('alias','range_alias'):
 TRANSFORM_REQUIRED_FIELDS['inverse.'+label] = {'plan%sc','plan%sd','plan%rc','plan%rd','plan%source','plan%dest',
  'plan%install_order.shape'}
for label in ('scalar_halo','vector_halo','outer_halo','scalar_initial','vector_initial','scalar_wave_halo',
              'vector_wave_halo','scalar_publication','vector_publication','fixed_scalar_route','fixed_vector_route'):
 TRANSFORM_REQUIRED_FIELDS['transform.'+label] = {'route.allocated','route.epoch','route.header','route.flags',
  'route.sc','route.sd','route.rc','route.rd','route.source','route.target','route.address','route.orientation',
  'route.send_peer','route.receive_peer','route.target_peer','route.allocations'}
TRANSFORM_CAPTURE_FIELDS = {'extent','masks','commit','first_level','options','level_phase','level','operation_range',
 'family','flags','family_and_range','component','range','domain.local_glo_id'}

def required_fields(scope):
 if scope['name']=='inverse.boundary' and not any(p[0]=='inverse.plan' for p in scope['parent']):
  return {'family_and_range'}
 return TRANSFORM_REQUIRED_FIELDS.get(scope['name'],set())

def vals(scope,name):
 e=scope['values'].get(name)
 if e is None: raise ValueError('Missing '+scope['name']+':'+name)
 code={1:'i',2:'q',4:'B'}[e.kind]
 return list(struct.unpack('<'+str(e.count)+code,e.data))

def validate_transform322(scopes,steps,coverage='candidate'):
 if coverage not in ('candidate','checked'):raise ValueError('Unknown sealed coverage variant')
 checked=coverage=='checked'
 for s in scopes:
  missing=required_fields(s)-s['tags']
  if missing: raise ValueError(f"{s['name']} epoch {s['epoch']}: missing fields {sorted(missing)}")
 if not steps:return
 plans={}
 transform_plan=None
 for s in scopes:
  name=s['name']
  if name=='transform.plan':
   ext=vals(s,'extent')
   if len(ext)!=8 or ext[5:]!=[5,5,30]:raise ValueError('Transform plan differs from sealed level5/30 layers')
   transform_plan=(s['epoch'],vals(s,'masks'))
  elif name=='transform.execute':
   if transform_plan!=(s['epoch'],vals(s,'masks')):raise ValueError('Transform execution has no matching preceding plan')
   if vals(s,'commit')!=[int(not checked)]:raise ValueError('Provisional transform commit differs from sealed coverage variant')
  elif name=='inverse.plan':
   ext=vals(s,'extent')
   if len(ext)!=6 or ext[3:]!=[30,5,5]:raise ValueError('Inverse plan differs from sealed level5/30 layers')
   if s['epoch'] in plans:raise ValueError('Duplicate inverse plan generation')
   plans[s['epoch']]=s
   if vals(s,'options')!=[1,0]:raise ValueError('Unexpected inverse implementation/options')
  elif name.startswith('inverse.') and name not in ('inverse.plan',):
   if s['epoch'] not in plans and not any(p[0]=='inverse.plan' for p in s['parent']):
    raise ValueError('Inverse invocation has no preceding plan generation')
 for step in range(1,steps+1):
  seq=[s for s in scopes if s['step']==step]
  # Final inverse runs outside inverse.plan. In particular the same label
  # inverse.boundary also names its compiled base route within the plan.
  runtime=[s for s in seq if not any(p[0]=='inverse.plan' for p in s['parent'])]
  c=Counter(s['name'] for s in runtime)
  ninverse=4 if checked else 1
  expected={'inverse.plan':1,'inverse.prepare':ninverse,'inverse.begin':ninverse,'inverse.end':ninverse,
   'inverse.scalar':0,'inverse.inner':0,'inverse.outer':0,
   'inverse.gather':26 if checked else 5,'inverse.scatter':26 if checked else 5,
   'inverse.boundary':26 if checked else 5,'transform.execute':3}
  for name,n in expected.items():
   if c[name]!=n:raise ValueError(f'Step {step}: require {n} {name} records, found {c[name]}')
  if not 1<=c['transform.plan']<=3:raise ValueError(f'Step {step}: require one to three provisional transform plans')
  kernels=[(s['name'],vals(s,'level_phase' if s['name']=='inverse.scalar' else 'level'))
    for s in runtime if s['name'] in ('inverse.scalar','inverse.inner','inverse.outer')]
  want=[]
  for l in ():want.extend([('inverse.scalar',[l,1]),('inverse.outer',[l]),('inverse.scalar',[l,2]),('inverse.inner',[l])])
  want=want*ninverse
  if kernels!=want:raise ValueError(f'Step {step}: inverse kernel level/phase/order differs')
  generation=next(s['epoch'] for s in runtime if s['name']=='inverse.plan')
  for s in runtime:
   if s['name'].startswith('inverse.') and s['epoch']!=generation:
    raise ValueError(f'Step {step}: inverse invocation generation differs')
  for begin in (s for s in runtime if s['name']=='inverse.begin'):
   if vals(begin,'options')!=[1,int(checked),0]:raise ValueError('Inverse options differ from sealed checked/candidate contract')
 # Each compiled plan must include all final bound per-level transfer plans.
 for s in scopes:
  if s['name'] not in ('inverse.plan','transform.plan'):continue
  children=[x for x in scopes if (s['name'],s['epoch']) in [tuple(p) for p in x['parent']]
            and x['record']>s['record'] and x['end_record']<s['end_record']]
  cc=Counter(x['name'] for x in children)
  if s['name']=='inverse.plan':
   needed={'inverse.interior':1,'inverse.boundary':1,'inverse.level':3,'inverse.level_interior':6,
    'inverse.level_boundary':6,'inverse.outer_level':2,'inverse.outer_install':2,
    'inverse.alias':4,'inverse.local_alias':2,'inverse.range_alias':4}
  else:
   needed={'transform.level':3,'transform.scalar_halo':3,'transform.vector_halo':3,'transform.outer_halo':0,
    'transform.scalar_initial':1,'transform.vector_initial':1,'transform.scalar_wave_halo':1,'transform.vector_wave_halo':1,
    'transform.scalar_publication':1,'transform.vector_publication':1,'transform.fixed_scalar_route':1,
    'transform.fixed_vector_route':1,'transform.scalar_adapter':1,'transform.vector_adapter':1}
  for name,n in needed.items():
   if cc[name]!=n:raise ValueError(f"{s['name']} generation {s['epoch']}: require {n} {name} components")

# Additional scope closure, independent of equality across the two executions.
import json
from pathlib import Path

def insist(ok,code,detail=''):
 if not ok:raise ValueError(code+(': '+str(detail) if detail else ''))

def one(s,k):
 insist(s['record_counts'].get(k)==1,'coverage.field_count',s['name']+'/'+k)
 return vals(s,k)

@__import__('functools').lru_cache(maxsize=1)
def source_schema():
 return json.loads((Path(__file__).with_name('schema322.json')).read_text())

def source_fields(s,schema):
 """Writer inventories close optional allocation branches without learning a run."""
 name=s['name'];writers=schema['writer_fields'];mapping={
 'catalog.executed':'audit_catalog316','migration.routes':'audit_migration316',
 'block':'audit_block316','ghost.plan':'audit_ghost316','writeback.plan':'audit_writeback316',
 'backend.input':'audit_rhs_input316','backend.rk_input':'audit_rhs_input316',
 'inverse.node_map':None,'transform.scalar_adapter':'audit_transform_adapter316',
 'transform.vector_adapter':'audit_transform_adapter316'}
 writer=mapping.get(name)
 if not writer:return
 required=set(writers.get(writer,[]))
 insist(bool(required),'coverage.source_writer',writer)
 # A field serialized only in the allocated branch is mandatory iff shape says allocated.
 for field in sorted(required):
  if field+'.shape' in required:
   shapes=[e for e in s['records'] if e.name==field+'.shape']
   insist(len(shapes)==1,'coverage.shape_count',name+'/'+field)
   shape=vals({'name':name,'values':{field+'.shape':shapes[0]}},field+'.shape')
   if shape==[-1]:
    insist(s['record_counts'].get(field,0)==0,'coverage.unallocated_payload',field)
    continue
   insist(all(x>=0 for x in shape),'coverage.source_shape',field)
   payload=s['values'].get(field)
   insist(payload is not None and payload.count==shape_payload_extent322(name,field,shape),'coverage.source_payload_extent',name+'/'+field)
  if name=='catalog.executed' and field in ('cut.domain','cut.root','catalog.domain_extent322'):
   expected=len(vals(s,'domain.local_glo_id'))
  else:expected=1
  insist(s['record_counts'].get(field,0)==expected,'coverage.required_source_field',name+'/'+field)

def validate_coverage322(scopes,rank,n,steps=4,coverage='candidate',arm=None):
 insist(steps==4,'coverage.steps')
 schema=source_schema()
 names=Counter(s['name'] for s in scopes)
 for s in scopes:
  insist(s['name'] in schema['scopes'],'coverage.unknown_scope',s['name'])
  for e in s['records']:
   kinds=schema['field_kinds'].get(e.name)
   insist(kinds is not None and e.kind in kinds,'coverage.unknown_field',s['name']+'/'+e.name)
  validate_allowed_fields322(s,schema)
  source_fields(s,schema)
  validate_shape_payloads322(s)
  validate_fixed_world_extents322(s,n,arm,rank)
 validate_rhs322(scopes,steps,coverage)
 validate_transform322(scopes,steps,coverage)
 validate_static_components322(scopes)
 validate_program_payloads322(scopes)
 insist(names['initialization']==1,'coverage.initialization')
 insist([s['epoch'] for s in scopes if s['name']=='step']==[1,2,3,4],'coverage.step_order')
 catalogs=[s for s in scopes if s['name']=='catalog.executed']
 insist([s['epoch'] for s in catalogs]==[1,2,3,4,5],'coverage.catalog_epochs')
 for name in ('migration.routes','storage.installed','ghost.plan','writeback.plan'):
  insist(names[name]==5,'coverage.catalog_lifecycle',name)
 for s in catalogs:
  insist(one(s,'physical.shape')==[6,640],'coverage.catalog_shape')
  phys=one(s,'physical');insist(len(phys)==3840,'coverage.catalog_extent')
  for i in range(640):
   row=phys[6*i:6*i+6];expected_owner=i//4 if arm in ('A0','B','checked-B') else owner325(i,arm)
   insist(row==[i,i//4,2+i%4,5,expected_owner,16],'coverage.owner_identity',[i,row])
  domains=one(s,'domain.local_glo_id')
  insist(domains==([rank] if rank<160 else []),'coverage.domain_identity',[rank,domains])
  if arm in ('B','checked-B') and rank>=160:
   insist(not domains,'coverage.empty_rank_domain')
 for step in range(1,5):
  current=[s for s in scopes if s['step']==step]
  count=Counter(s['name'] for s in current)
  for name in ('rhs.prepare_check','rhs.evaluate','rhs.prepared_state','rhs.retain','rhs.rk_prepare','rhs.rk_install'):
   insist(count[name]==4,'coverage.step_lifecycle',f'{step}/{name}')
  insist(count['rhs.prepare']==1,'coverage.prepare_per_step',step)
  insist(count['transform.plan']==1,'coverage.transform_plan_per_step',step)
 # This list is generated ONLY from literal additions supplied in the sealed source schema,
 # while required call counts are configured here, independent of observed equality.
 for name,num in schema.get('stage322_required_scopes',{}).items():
  if name=='inverse.values322' and coverage=='checked':num=48
  if name=='adapt322':
   # Frozen resume=3, level5/tol0 restart: one mask adaptation followed by
   # exactly two prepare_tolerance312 iterations, then one per timestep.
   # Keep the original schema sealed and reject a different schema contract.
   insist(num==5,'coverage.adapt_sealed_count')
   num=steps+3
  insist(names[name]==num,'coverage.stage322_count',name)
 validate_adaptation_lifecycle322(scopes,steps)
 for s in scopes:
  req=schema.get('stage322_scope_fields',{}).get(s['name'])
  if req is not None:
   for field in req:
    insist(s['record_counts'].get(field)==1,'coverage.stage322_field',s['name']+'/'+field)
 validate_transport_calls322(scopes)
 validate_added322(scopes,rank,n,arm,coverage)
 if arm in ('B','checked-B') and rank>=160:
  for s in scopes:
   if s['name']=='block':raise ValueError('coverage.empty_rank_block')
   if s['name']=='rhs.prepare':
    insist(one(s,'retention.blocks')==[0],'coverage.empty_rank_retention')
    insist(one(s,'backend.mask.domains')==[0],'coverage.empty_rank_masks')
    insist(one(s,'program.header')[1]==0,'coverage.empty_rank_operators')
    insist(s['values']['placement.owned'].count==0,'coverage.empty_rank_owned')
   if s['name']=='inverse.plan':
    insist(one(s,'extent')[:2]==[0,0],'coverage.empty_rank_inverse')
   if s['name']=='transform.plan':
    insist(one(s,'extent')[:5]==[0,0,0,0,0],'coverage.empty_rank_transform')
 return {'passed':True,'scopes':len(scopes),'rhs_evaluations':names['rhs.evaluate'],
         'prepared_states':names['rhs.prepared_state'],'source_inventory':schema['schema']}


def validate_transport_calls322(scopes):
 events=[s for s in scopes if s['name']=='transport.event322']
 footer=[s for s in scopes if s['name']=='transport.footer322']
 insist(len(footer)==1,'coverage.transport_footer')
 totals=one(footer[0],'transport.calls322');insist(len(totals)==4,'coverage.transport_totals')
 seen=[0,0,0,0]
 for s in events:
  identity=one(s,'transport.identity322');insist(len(identity)==3 and 1<=identity[0]<=4,'coverage.transport_kind')
  kind=identity[0];seen[kind-1]+=1
  insist(s['epoch']==seen[kind-1],'coverage.transport_sequence')
  insist(len(one(s,'transport.epoch322'))==1,'coverage.transport_epoch')
  counts=one(s,'transport.count322');insist(len(counts)==4 and all(x>=0 for x in counts),'coverage.transport_extent')
  if kind==4:insist(len(one(s,'transport.flags322'))==3,'coverage.transport_flags')
 insist(totals==seen,'coverage.transport_complete')

def children_of(scopes,parent,name=None):
 return [s for s in scopes if s['record']>parent['record'] and s['end_record']<parent['end_record'] and
         (name is None or s['name']==name)]

def repeated(s,k):
 return [vals({'name':s['name'],'values':{k:e}},k) for e in s['records'] if e.name==k]

def payload_extent(s,name,shape):
 shape_values=one(s,shape)
 insist(all(x>=0 for x in shape_values),'coverage.numeric_shape',name)
 insist(s['record_counts'].get(name)==1 and s['values'][name].kind==3 and
        s['values'][name].count==math.prod(shape_values),'coverage.numeric_extent',name)

def validate_added322(scopes,rank,n,arm,coverage):
 byname={name:[s for s in scopes if s['name']==name] for name in {s['name'] for s in scopes}}
 nblocks=4 if arm in ('A0','B','checked-B') and rank<160 else (0 if arm in ('B','checked-B') else blocks325(rank,arm))
 insist(len(byname.get('block',[]))==5*nblocks,'coverage.installed_blocks')
 insist(len(byname.get('rhs.retained322',[]))==16*nblocks,'coverage.retained_blocks')
 for p in byname.get('rhs.prepare',[]):
  nd=one(p,'backend.mask.domains')[0]
  for tag in ('backend.mask.domain','backend.mask.node','backend.mask.edge'):
   insist(p['record_counts'].get(tag,0)==nd,'coverage.mask_count',tag)
  nb=one(p,'retention.blocks')[0];insist(nb==nblocks,'coverage.retention_count')
  for tag in ('retention.block','retention.node','retention.length.shape','retention.length'):
   insist(p['record_counts'].get(tag,0)==nb,'coverage.retention_field_count',tag)
  for tag in ('retention.column','retention.ready.shape','retention.exner.shape','retention.theta.shape','retention.ready'):
   insist(p['record_counts'].get(tag,0)==3*nb,'coverage.retention_column_count',tag)
  insist(repeated(p,'retention.column')==[[i] for _ in range(nb) for i in (1,2,3)],'coverage.retention_column_order')
  bound=one(p,'transport.bound');bounds=one(p,'transport.route.bounds')
  insist(bounds==[5,5] and len(bound)==20,'coverage.route_bounds')
  insist(repeated(p,'transport.route.coordinate')==[[i+1,5] for i,v in enumerate(bound) if v],
         'coverage.route_coordinates')
  for family in FAMILIES:
   header=one(p,family+'.shape');insist(len(header)==4,'coverage.operator_shape',family)
   for width,suffix,kind in zip(header[:3],('i','r','l'),(1,3,4)):
    if width:
     e=p['values'].get(family+'.'+suffix)
     insist(e is not None and e.kind==kind and e.count==width*header[3],'coverage.operator_extent',family)
 for s in byname.get('rhs.output322',[]):
  for field in ('scalar','velocity'):payload_extent(s,'output.'+field,'output.'+field+'.shape')
 for s in byname.get('rhs.rk_values322',[]):
  for field in ('adapter','value','scalar','mean','velocity'):payload_extent(s,'rk.'+field,'rk.'+field+'.shape')
 for s in byname.get('rhs.rk_installed322',[]):
  na,nr,nz=one(s,'installed.extent');insist(nz==30,'coverage.installed_layers')
  for name,count,width in [('installed.adapter',na,3),('installed.resident',nr,3)]:
   rows=repeated(s,name);insist(len(rows)==count and [r[0] for r in rows]==list(range(1,count+1)),
                              'coverage.installed_order',name)
  for name,count,width in [('installed.values',na*30,5),('installed.temperature',nr*30,1)]:
   es=[e for e in s['records'] if e.name==name]
   insist(len(es)==count and all(e.kind==3 and e.count==width for e in es),'coverage.installed_values',name)
 for step in range(1,5):
  bounds=[s for s in byname.get('main.boundary322',[]) if s['step']==step]
  insist([one(s,'main.identity')[0] for s in bounds]==list(range(1,7)),'coverage.main_boundary_order')
  for s in bounds:
   row=one(s,'main.identity');insist(len(row)==14 and row[1:5]==[step,rank,n,int(rank<160)] and row[6:11]==[5,5,5,5,30],
                                  'coverage.main_configuration')
   # After initialized storage, active owners and empty extras are explicit.
   insist(row[5]==nblocks,'coverage.main_owned_blocks')
   fam={1:[1],2:[1,3],3:[1],4:[2],5:[1,2,3],6:[1]}[row[0]]
   kids=children_of(scopes,s,'domain.fields322')
   insist([one(q,'fields.header')[0] for q in kids]==fam,'coverage.main_field_families')
  provisional=[s for s in byname.get('rk.provisional322',[]) if s['step']==step]
  insist([one(s,'rk.stage') for s in provisional]==[[1,4],[2,4],[3,4]],'coverage.rk_stage_order')
  for s in provisional:
   insist([one(q,'fields.header')[0] for q in children_of(scopes,s,'domain.fields322') if q['parent'][-1][0]=='rk.provisional322']==[20,22],
          'coverage.rk_boundary_families')
  final=[s for s in byname.get('rk.final322',[]) if s['step']==step]
  insist(len(final)==1 and one(final[0],'rk.stage')==[4,4],'coverage.rk_final')
  insist([one(q,'fields.header')[0] for q in children_of(scopes,final[0],'domain.fields322') if q['parent'][-1][0]=='rk.final322']==[20,21,22,23],
         'coverage.rk_final_families')
 validate_catalog_geometry322(scopes,rank)
 for s in byname.get('domain.fields322',[]):validate_domain_fields322(s,rank)
 _validate_publication322(scopes,coverage)
 _validate_retention322(scopes,nblocks)
 validate_transform_calls322(scopes,coverage)
 for s in byname.get('rhs.rk_prepare',[]):
  kids=[q for q in children_of(scopes,s,'domain.fields322') if q['parent'][-1][0]=='rhs.rk_prepare']
  insist([one(q,'fields.header')[0] for q in kids]==[10],'coverage.rk_base_family')
 for s in byname.get('transform.values322',[]):
  for field in ('scalar','vector','scalar_wavelet','vector_wavelet'):
   payload_extent(s,'transform.'+field,'transform.'+field+'.shape')
 for s in byname.get('transform.execute',[]):
  valsc=children_of(scopes,s,'transform.values322')
  insist([one(t,'transform.phase')[0] for t in valsc]==[1,2,3,4],'coverage.transform_value_phases')
 for s in byname.get('inverse.values322',[]):
  for field in ('value','scaffold','alias'):payload_extent(s,'inverse.'+field,'inverse.'+field+'.shape')
  if one(s,'inverse.restore.shape')!=[-1]:payload_extent(s,'inverse.restore','inverse.restore.shape')
 for step in range(1,5):
  phases=[one(s,'inverse.phase')[0] for s in byname.get('inverse.values322',[]) if s['step']==step]
  insist(phases==[1,2,3]*(4 if coverage=='checked' else 1),'coverage.inverse_value_phases')
 # Known root path: one startup and one end-step CFL, four physics calls.
 for name in ('cfl322','physics322'):
  for step in range(1,5):
   insist(len([s for s in byname.get(name,[]) if s['step']==step])==1,'coverage.consumer_call_count',name)
 _validate_visits322(scopes,rank)
 validate_capture_extents322(scopes)
 validate_numeric_bindings322(scopes)
 if arm in ('B','checked-B') and rank>=160:
  for s in byname.get('transport.event322',[]):
   insist(one(s,'transport.count322')==[0,0,0,0],'coverage.empty_rank_transport')


def validate_adaptation_lifecycle322(scopes,steps=4):
 """Exact frozen restart and timestep schedule; not a count learned from runs."""
 insist(steps==4,'coverage.adapt_configured_steps')
 adaptations=[s for s in scopes if s['name']=='adapt322']
 expected=[(None,0,(('initialization',0),),[0,0]),
           (None,0,(('initialization',0),),[1,0]),
           (None,0,(('initialization',0),),[1,0])]
 expected.extend((step,step,(('step',step),('step322',step)),[1,1])
                 for step in range(1,steps+1))
 insist(len(adaptations)==len(expected),'coverage.adapt_lifecycle_count')
 for s,(step,epoch,parent,options) in zip(adaptations,expected):
  insist((s['step'],s['epoch'],tuple(tuple(p) for p in s['parent']))==(step,epoch,parent),
         'coverage.adapt_lifecycle_context')
  actual=one(s,'adapt.options');entry=s['values']['adapt.options']
  insist(entry.kind==4 and entry.count==2 and actual==options,
         'coverage.adapt_lifecycle_options')
  insist(repeated(s,'adapt.operation')==[[i] for i in range(1,13)],'coverage.adaptation_order')
  insist(one(s,'adapt.refined')==[0],'coverage.fixed_grid_refinement')


def validate_fixed_world_extents322(scope,n,arm,rank):
 """Closed level5 ghost peer matrices and inverse MPI_Request capacity."""
 def ints(name,width):
  values=one(scope,name);entry=scope['values'][name]
  insist(entry.kind==1 and entry.count==width,'coverage.world_extent_descriptor',name)
  return values
 if scope['name']=='inverse.plan':
  insist(ints('requests.shape',1)==[2*n],'coverage.inverse_request_capacity')
 elif scope['name']=='ghost.plan':
  insist(ints('route_level_start',1)==[5] and ints('route_level_count',1)==[1],
         'coverage.ghost_level_range')
  pre='ghost_exchange_plan%'
  for direction in ('send','recv'):
   stem=pre+direction
   counts=ints(stem+'_record_count',n)
   insist(all(x>=0 for x in counts),'coverage.ghost_level_negative_count')
   prefix=[];total=0
   for count in counts:prefix.append(total);total+=count
   insist(ints(stem+'_record_displ',n)==prefix,'coverage.ghost_record_prefix')
   for suffix,expected in (('count',counts),('displ',prefix)):
    field=stem+'_level_'+suffix
    insist(ints(field+'.shape',2)==[n,1],'coverage.ghost_level_shape',field)
    insist(ints(field,n)==expected,'coverage.ghost_level_record_identity',field)
   insist(ints(stem+'_level_base.shape',1)==[1] and ints(stem+'_level_base',1)==[0],
          'coverage.ghost_level_base')
   if arm in ('B','checked-B'):
    insist(not any(counts[160:]),'coverage.ghost_level_extra_peer')
    if rank>=160:insist(not any(counts),'coverage.ghost_level_empty_rank')


def validate_domain_fields322(s,rank):
 header=one(s,'fields.header');insist(len(header)==4,'coverage.fields_header')
 family,nvar,nz,nd=header
 layers={1:41,2:41,3:41,10:30,20:30,21:30,22:30,23:30,30:30,31:30}
 insist(family in layers and nvar==3 and nz==layers[family] and nd==int(rank<160),'coverage.fields_extent')
 rows=repeated(s,'fields.domain');insist(len(rows)==nd,'coverage.fields_domains')
 for field in ('fields.node_level','fields.node_mask','fields.edge_mask','fields.poles'):
  insist(s['record_counts'].get(field,0)==nd,'coverage.fields_geometry_count',field)
 if not nd:
  insist(len(s['records'])==1,'coverage.empty_rank_field_records');return
 d,g,nn,np,nb=rows[0];insist([d,g]==[1,rank] and nn==321 and np>=4 and nb>=0,'coverage.fields_domain_identity')
 insist(s['values']['fields.node_level'].count==nn and s['values']['fields.node_mask'].count==nn and
        s['values']['fields.edge_mask'].count==3*nn,'coverage.fields_mask_extent')
 patches=repeated(s,'fields.patch');insist(len(patches)==np,'coverage.fields_patch_count')
 insist([v[0] for v in patches]==list(range(np)) and all(len(v)==8 for v in patches),'coverage.fields_patch_identity')
 coordinates=repeated(s,'fields.coordinate')
 insist(len(coordinates)==nvar*nz,'coverage.fields_coordinate_count')
 insist([r[:3] for r in coordinates]==[[1,k,v] for k in range(1,nz+1) for v in range(1,nvar+1)],
        'coverage.fields_coordinate_order')
 payloads=[e for e in s['records'] if e.name=='fields.values']
 boundaries=[e for e in s['records'] if e.name=='fields.boundary_current']
 insist(len(payloads)==len(boundaries)==len(coordinates),'coverage.fields_values_count')
 for row,payload,flag in zip(coordinates,payloads,boundaries):
  variable=row[2];expected_length=(3 if variable==1 else 1)*nn;expected_position=2 if variable==1 else 1
  insist(len(row)==6 and row[3]==expected_position and row[5]==expected_length and payload.kind==3 and payload.count==expected_length and flag.kind==4 and flag.count==1,
         'coverage.fields_values_extent')


def _validate_visits322(scopes,rank):
 fields=[s for s in scopes if s['name']=='domain.fields322']
 insist(bool(fields),'coverage.missing_domain_state')
 geometry=fields[0]
 if rank>=160:count=0
 else:
  patches=repeated(geometry,'fields.patch');active=[p for p in patches if p[1]==5 and p[3]==0]
  insist(len(active)==4,'coverage.uniform_patch_count')
  poles=one(geometry,'fields.poles');insist(len(poles)==2,'coverage.pole_extent')
  count=64+sum(poles)
 for s in scopes:
  if s['name']=='cfl322':
   rows=repeated(s,'cfl.visit')
   insist(one(s,'cfl.visits')==[count] and len(rows)==count,'coverage.cfl_visit_count')
   insist(all(len(row)==10 and row[:2]==[1,rank] and row[5]==5 for row in rows),'coverage.cfl_visit_identity')
   insist(len({tuple(row[:6]) for row in rows})==count,'coverage.cfl_duplicate_visit')
  elif s['name']=='physics322':
   columns=children_of(scopes,s,'physics.column322')
   insist(one(s,'physics.visits')==[count] and len(columns)==count,'coverage.physics_visit_count')
   insist([c['epoch'] for c in columns]==list(range(1,count+1)),'coverage.physics_visit_order')
   rows=[one(c,'physics.column') for c in columns]
   insist(all(len(row)==9 and row[:2]==[1,rank] for row in rows),'coverage.physics_visit_identity')
   insist(len({tuple(row[:8]) for row in rows})==count,'coverage.physics_duplicate_visit')


def _validate_publication322(scopes,coverage):
 for s in scopes:
  if s['name'] in ('inverse.publish322','transform.installed322'):
   fam=31 if s['name']=='inverse.publish322' else 30
   kids=children_of(scopes,s,'domain.fields322')
   insist(len(kids)==1 and one(kids[0],'fields.header')[0]==fam,'coverage.publication_child',s['name'])
 for step in range(1,5):
  pubs=[one(s,'publication.mode') for s in scopes if s['step']==step and s['name']=='inverse.publish322']
  expected=([[1,1,0],[2,1,0]]*4+[[1,1,1]]) if coverage=='checked' else [[1,1,0],[1,2,0],[2,1,0],[1,1,1]]
  insist(pubs==expected,'coverage.inverse_publication_schedule',step)

def _validate_retention322(scopes,nblocks):
 for step in range(1,5):
  prepare=next(s for s in scopes if s['step']==step and s['name']=='rhs.prepare')
  block_rows=repeated(prepare,'retention.block')
  insist(len(block_rows)==nblocks and [row[0] for row in block_rows]==list(range(1,nblocks+1)),
         'coverage.retention_local_order')
  calls=[s for s in scopes if s['step']==step and s['name']=='rhs.retain']
  for stage,call in enumerate(calls):
   emitted=children_of(scopes,call,'rhs.retained322')
   insist(len(emitted)==nblocks,'coverage.retained_per_rhs')
   rows=[one(s,'retained.block') for s in emitted]
   insist([row[:2] for row in rows]==[[catalog,local] for local,catalog in block_rows],
          'coverage.retained_identity_order')
   insist(all(s['epoch']==stage for s in emitted),'coverage.retained_stage')


def validate_capture_extents322(scopes):
 for s in scopes:
  if s['name']=='physics.column322':
   for tag,count in [('physics.coordinates',2),('physics.input',256),('physics.output',131)]:
    e=s['values'].get(tag);insist(e is not None and e.kind==3 and e.count==count,'coverage.physics_payload_extent',tag)
  elif s['name']=='rhs.retained322':
   row=one(s,'retained.block');insist(len(row)==4 and row[2:]==[1,16],'coverage.retained_patch_extent')
   for tag,count in [('retained.scalar',2*41*16),('retained.velocity',3*41*16)]:
    e=s['values'].get(tag);insist(e is not None and e.kind==3 and e.count==count,'coverage.retained_payload_extent',tag)


def allowed_fields322(scope,schema):
 name=scope['name'];writers=schema['writer_fields']
 # Both frozen inverse invocation writers emit level only when changed_level
 # is present, recorded by flags[2]. It is not an allocation-shape field.
 if name in ('inverse.gather','inverse.scatter'):return {'family','flags','level'}
 static={'catalog.executed':'audit_catalog316','migration.routes':'audit_migration316',
 'block':'audit_block316','ghost.plan':'audit_ghost316','writeback.plan':'audit_writeback316',
 'backend.input':'audit_rhs_input316','backend.rk_input':'audit_rhs_input316',
 'transform.scalar_adapter':'audit_transform_adapter316','transform.vector_adapter':'audit_transform_adapter316'}
 if name in static:return set(writers[static[name]])
 if name in ('initialization','step','storage.installed'):return set()
 if name in ('backend.adapter_route','input.scalar_route','input.vector_route','input.mean_route') or name in (
 'transform.scalar_halo','transform.vector_halo','transform.outer_halo','transform.scalar_initial','transform.vector_initial',
 'transform.scalar_wave_halo','transform.vector_wave_halo','transform.scalar_publication','transform.vector_publication',
 'transform.fixed_scalar_route','transform.fixed_vector_route'):
  return set(writers['audit_rhs_exchange316'])
 if name=='rhs.prepare':
  ans=set(RHS_REQUIRED_FIELDS[name])
  for writer in ('audit_rhs_distributed316','audit_rhs_program316','audit_rhs_retention316','audit_rhs_exchange316'):
   ans.update(writers[writer])
  ans.update(('backend.mask.domain','backend.mask.node','backend.mask.edge'))
  return ans
 if name in RHS_REQUIRED_FIELDS:
  ans=set(RHS_REQUIRED_FIELDS[name])
  if name=='rhs.evaluate':ans.update(writers['run_rhs_distributed'])
  return ans
 if name in TRANSFORM_REQUIRED_FIELDS:
  ans=set(required_fields(scope));ans.update(f+'.shape' for f in list(ans) if not f.endswith('.shape'))
  if name=='transform.execute':ans.add('transform.call322')
  if name in ('inverse.interior','inverse.level_interior','inverse.level_boundary','inverse.outer_install') or (
      name=='inverse.boundary' and any(p[0]=='inverse.plan' for p in scope['parent'])):
   ans.update(writers['audit_transfer316'])
  if name in ('inverse.alias','inverse.range_alias'):ans.update(writers['audit_alias316'])
  return ans
 ans=set(schema.get('stage322_scope_fields',{}).get(name,[]))
 extras={
 'domain.fields322':set(writers['audit_domain_fields322']),
 'rhs.rk_installed322':{'installed.adapter','installed.values','installed.resident','installed.temperature'},
 'adapt322':{'adapt.operation'},'cfl322':{'cfl.visit'},
 'rk.provisional322':{'transform.call322'},
 'inverse.values322':{'inverse.restore'},'transform.installed322':set(),
 'transport.event322':{'transport.flags322'}}
 ans.update(extras.get(name,set()))
 return ans


def validate_inverse_level322(scope):
 """Source-bound optional changed_level contract for gather/scatter only."""
 if scope['name'] not in ('inverse.gather','inverse.scatter'):return
 for field,kind,width in (('family',1,2),('flags',4,3)):
  entries=[e for e in scope['records'] if e.name==field]
  insist(scope['record_counts'].get(field,0)==1 and len(entries)==1 and
         entries[0].kind==kind and entries[0].count==width,
         'coverage.inverse_invocation_descriptor',scope['name']+'/'+field)
 flags=one(scope,'flags')
 insist(all(x in (0,1) for x in flags),'coverage.inverse_invocation_flags',scope['name'])
 levels=[e for e in scope['records'] if e.name=='level']
 insist(scope['record_counts'].get('level',0)==flags[2] and len(levels)==flags[2],
        'coverage.inverse_optional_level_presence',scope['name'])
 if levels:
  insist(levels[0].kind==1 and levels[0].count==1,
         'coverage.inverse_optional_level_descriptor',scope['name'])


def validate_allowed_fields322(scope,schema=None):
 if schema is None:schema=source_schema()
 allowed=allowed_fields322(scope,schema)
 for e in scope['records']:
  insist(e.name in allowed,'coverage.field_scope',scope['name']+'/'+e.name)
 validate_inverse_level322(scope)


def validate_catalog_geometry322(scopes,rank):
 """Every domain snapshot is tied to all five committed catalog geometries."""
 catalogs=[s for s in scopes if s['name']=='catalog.executed']
 insist(len(catalogs)==5,'coverage.geometry_catalog_count')
 expected=None
 for catalog in catalogs:
  rows=repeated(catalog,'catalog.domain_extent322')
  insist(len(rows)==int(rank<160),'coverage.catalog_domain_extent_count')
  if rank<160:
   insist(len(rows[0])==5 and rows[0][:3]==[1,rank,321] and rows[0][3]>=4 and rows[0][4]>=0,
          'coverage.catalog_domain_extent')
  if expected is None:expected=rows
  insist(rows==expected,'coverage.geometry_changed')
 for scope in scopes:
  if scope['name']=='domain.fields322':
   insist(repeated(scope,'fields.domain')==expected,'coverage.fields_catalog_geometry')


def expected_transform_calls322(plan,route_scopes,coverage):
 """Level5 source schedule, including each actual kernel, reset and publication."""
 def route(name,epoch=None):
  selected=[s for s in route_scopes if s['name']==name and (epoch is None or ('transform.level',epoch) in s['parent'])]
  insist(len(selected)==1,'coverage.transform_route_count',name)
  return one(selected[0],'route.header')[0]
 scalar=one(plan,'scalar.shape');vector=one(plan,'vector.shape')
 sw=one(plan,'scalar_wavelet.shape');vw=one(plan,'vector_wavelet.shape')
 ns=plan['values']['save_node'].count;ni=plan['values']['interior_node'].count
 program=one(plan,'scalar.program.shape');insist(len(program)==2 and program[0]==33,'coverage.transform_program_shape')
 for shape in (scalar,vector,sw,vw):insist(len(shape)==2,'coverage.transform_runtime_shape')
 rows=[[7,1,0,0],[7,2,0,0]]
 def exchange(name,shape,epoch=None):rows.append([1,route(name,epoch),*shape])
 exchange('transform.scalar_initial',scalar);exchange('transform.vector_initial',vector)
 rows.extend([4,i,0,0] for i in range(1,program[1]+1))
 if coverage!='checked':rows.append([3,2,ni,0])
 rows.extend([[2,1,0,ns],[2,2,0,ns]])
 exchange('transform.scalar_halo',scalar,4);exchange('transform.vector_halo',vector,4)
 exchange('transform.scalar_wave_halo',sw);exchange('transform.vector_wave_halo',vw)
 rows.append([2,1,1,ns]);exchange('transform.scalar_halo',scalar,5)
 rows.extend([[2,1,1,ns],[2,2,1,ns]])
 exchange('transform.scalar_halo',scalar,3);exchange('transform.vector_halo',vector,3)
 if coverage!='checked':rows.append([3,1,ni,0])
 exchange('transform.scalar_publication',scalar);exchange('transform.vector_publication',vector)
 return rows


def validate_transform_calls322(scopes,coverage):
 for step in range(1,5):
  plans=[s for s in scopes if s['name']=='transform.plan' and s['step']==step]
  insist(len(plans)==1,'coverage.transform_call_plan')
  plan=plans[0];routes=children_of(scopes,plan)
  expected=expected_transform_calls322(plan,routes,coverage)
  calls=[s for s in scopes if s['name']=='transform.execute' and s['step']==step]
  insist(len(calls)==3,'coverage.transform_call_invocations')
  for scope in calls:
   insist(repeated(scope,'transform.call322')==expected,'coverage.transform_call_schedule')
  for scope in [s for s in scopes if s['name']=='rk.provisional322' and s['step']==step]:
   expected_install=[[3,1,plan['values']['interior_node'].count,0]] if coverage=='checked' else []
   insist(repeated(scope,'transform.call322')==expected_install,'coverage.transform_checked_install')


def validate_shape_payloads322(scope):
 # Ordinary .shape records capture Fortran allocation extents. The two
 # inverse node maps serialize only their initialized extent[0] prefix; keep
 # their allocation capacities in the stream and validate both separately.
 # RHS operator shapes encode record widths plus count and are checked apart.
 if scope['name']=='inverse.plan':
  extent=one(scope,'extent')
  insist(len(extent)==6 and extent[0]>=0,'coverage.inverse_node_prefix_extent')
  capacities=[]
  for field in ('node_domain','node_id'):
   dims=one(scope,field+'.shape')
   insist(len(dims)==1 and dims[0]>=extent[0],
          'coverage.inverse_node_prefix_capacity',field)
   insist(scope['record_counts'].get(field)==1,'coverage.inverse_node_prefix_payload',field)
   payload=scope['values'][field]
   insist(payload.kind==1 and payload.count==extent[0],
          'coverage.shape_payload_extent',scope['name']+'/'+field)
   capacities.append(dims[0])
  insist(capacities[0]==capacities[1],'coverage.inverse_node_prefix_capacity_pair')
 for tag in scope['record_counts']:
  if not tag.endswith('.shape') or tag[:-6] in FAMILIES:continue
  base=tag[:-6]
  if base not in scope['record_counts']:continue
  shapes=[e for e in scope['records'] if e.name==tag]
  payloads=[e for e in scope['records'] if e.name==base]
  insist(len(shapes)==len(payloads),'coverage.shape_payload_count',scope['name']+'/'+base)
  for shape,payload in zip(shapes,payloads):
   dims=vals({'name':scope['name'],'values':{tag:shape}},tag)
   expected=shape_payload_extent322(scope['name'],base,dims)
   if scope['name']=='inverse.plan' and base in ('node_domain','node_id'):
    expected=extent[0]  # Both capacity/prefix contracts were checked above.
   insist(all(x>=0 for x in dims) and payload.count==expected,
          'coverage.shape_payload_extent',scope['name']+'/'+base)


def validate_static_components322(scopes):
 for plan in [s for s in scopes if s['name'] in ('inverse.plan','transform.plan')]:
  children=children_of(scopes,plan)
  if plan['name']=='transform.plan':
   insist([s['epoch'] for s in children if s['name']=='transform.level']==[3,4,5],
          'coverage.transform_level_order')
  else:
   for name,epochs in [('inverse.level',[3,4,5]),('inverse.outer_level',[4,5]),('inverse.local_alias',[1,2])]:
    selected=[s for s in children if s['name']==name]
    insist([s['epoch'] for s in selected]==epochs,'coverage.inverse_component_scope',name)
    for scope in selected:
     tag='range' if name=='inverse.local_alias' else 'component'
     insist(repeated(scope,tag)==[[1],[2]],'coverage.inverse_component_order',name)


def validate_numeric_bindings322(scopes):
 for step in range(1,5):
  current=[s for s in scopes if s['step']==step]
  plans={name:next(s for s in current if s['name']==name) for name in ('rhs.prepare','transform.plan','inverse.plan')}
  for scope in current:
   pairs=[];plan=None
   if scope['name']=='rhs.output322':
    plan=plans['rhs.prepare'];pairs=[('output.'+f+'.shape','executor.'+f+'_rhs.shape') for f in ('scalar','velocity')]
   elif scope['name']=='rhs.rk_values322':
    plan=plans['rhs.prepare'];pairs=[('rk.'+f+'.shape','backend.'+g+'.shape') for f,g in
      [('adapter','adapter_value'),('value','rk_value'),('scalar','rk_scalar'),('mean','rk_mean'),('velocity','rk_velocity')]]
   elif scope['name']=='transform.values322':
    plan=plans['transform.plan'];pairs=[('transform.'+f+'.shape',f+'.shape') for f in ('scalar','vector','scalar_wavelet','vector_wavelet')]
    insist(one(scope,'transform.phase')[1:]==one(plan,'extent')[:3],'coverage.transform_phase_geometry')
   elif scope['name']=='inverse.values322':
    plan=plans['inverse.plan'];pairs=[('inverse.'+f+'.shape',g+'.shape') for f,g in
      [('value','value'),('scaffold','scaffold_value'),('alias','consumer_alias_value'),('restore','restore203')]]
    extent=one(plan,'extent');insist(one(scope,'inverse.phase')[1:]==[extent[0],extent[2],30],'coverage.inverse_phase_geometry')
   for field,static_field in pairs:
    insist(one(scope,field)==one(plan,static_field),'coverage.numeric_plan_shape',scope['name']+'/'+field)


def shape_payload_extent322(scope,field,shape):
 # audit_block316 emits five semantic integers per Block_Stencil object:
 # storage,id,offset,dims(2). All other ordinary array captures are flat values.
 width=5 if (scope,field)==('block','stencil') else 1
 return width*math.prod(shape)


def validate_program_payloads322(scopes):
 for scope in scopes:
  if scope['name'] not in ('inverse.plan','transform.plan'):continue
  insist(all(c==1 for c in scope['record_counts'].values()),'coverage.plan_duplicate_field',scope['name'])
  widths={'outer':9} if scope['name']=='inverse.plan' else {'scalar':12,'inner':32,'outer':9}
  for family,width in widths.items():
   shape=one(scope,family+'.program.shape');payload=scope['values'].get(family+'.coefficients')
   insist(len(shape)==2 and payload is not None and payload.kind==3 and payload.count==width*shape[1],
          'coverage.program_coefficient_extent',scope['name']+'/'+family)

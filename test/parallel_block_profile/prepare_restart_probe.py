#!/usr/bin/env python3
"""Prepare an external source/build snapshot for read-only checkpoint probes."""
import argparse
from pathlib import Path
import shutil
from prepare_legacy_profile import replace_once

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--source', type=Path, required=True)
parser.add_argument('--build', type=Path, required=True)
parser.add_argument('--out', type=Path, required=True)
parser.add_argument('--block', action='store_true')
parser.add_argument('--phase-step', type=int, action='append', default=[],
                    help='Also snapshot this cumulative step; repeat for several steps')
parser.add_argument('--phase-range', type=int, nargs=2, metavar=('FIRST','LAST'))
parser.add_argument('--physics', action='store_true', help='Read-only packed physics input/output probe')
parser.add_argument('--physics-replay-u', type=int, nargs=5, action='append', default=[],
                    metavar=('STEP','DOMAIN','NODE','LEVEL','BITS'),
                    help='CAUSAL DIAGNOSTIC ONLY: replace one packed U input in the isolated copy')
parser.add_argument('--physics-replay-velocity', nargs=6, action='append', default=[],
                    metavar=('FIELD','STEP','DOMAIN','NODE','LEVEL','BITS'),
                    help='CAUSAL DIAGNOSTIC ONLY: replace one packed U/V/W input; repeat for several')
parser.add_argument('--integrator', choices=('RK3','RK4'),
                    help='Override integrator in the isolated source copy only')
args = parser.parse_args()
if (args.physics_replay_u or args.physics_replay_velocity) and not args.physics:
    raise ValueError('--physics-replay-u requires --physics')
source=args.source.resolve(strict=True)
out=args.out.resolve()
if out.exists() or out.is_relative_to(source):
    raise ValueError('Use a new external directory')
out.mkdir()
shutil.copytree(source/'src',out/'src',symlinks=True)
shutil.copytree(source/'test/climate',out/'test/climate',symlinks=True)
shutil.copytree(args.build,out/'build',symlinks=True)
shutil.copy2(source/'Makefile',out/'Makefile')
shutil.copy2(Path(__file__).with_name('restart_probe.f90'),out/'src/restart_probe.f90')
if args.integrator:
    import re
    path=out/'src/shared.f90'
    text,count=re.subn(r'(timeint_type\s*=\s*)"RK[34]"',
                      lambda m:m[1]+'"'+args.integrator+'"',path.read_text())
    if count!=1: raise ValueError('Expected exactly one default integrator')
    path.write_text(text)
path=out/'Makefile'
path.write_text(replace_once(path.read_text(),'SRC += main.f90',
                            'SRC += restart_probe.f90 '+chr(92)+'\n       main.f90'))
path=out/'src/main.f90'
text=replace_once(path.read_text(),'module main_mod\n\n','module main_mod\n  use restart_probe_mod\n\n')
text=replace_once(text,'    call dump_adapt_mpi (cp_idx)',
                 "    call restart_probe('before-dump')\n    call dump_adapt_mpi (cp_idx)\n"+
                 "    call restart_probe('after-dump')")
text=replace_once(text,'  end subroutine restart',
                 "    call restart_probe('restart-ready')\n  end subroutine restart")
if args.phase_step or args.phase_range:
    condition=' .or. &\n         '.join(f'istep_cumul=={step}' for step in args.phase_step)
    if args.phase_range:
        if args.phase_step: raise ValueError('Use either --phase-step or --phase-range')
        first,last=args.phase_range
        if first>last: raise ValueError('Invalid phase range')
        condition=f'istep_cumul>={first}.and.istep_cumul<={last}'
    for marker, tag in [('    !    Physics split step', 'after-dynamics'),
                        ('    !    Grid adaptation', 'after-physics'),
                        ('    !    Vertical remapping (', 'after-adaptation')]:
        text=replace_once(text,marker,
            f"    if ({condition}) call restart_probe('{tag}')\n"+marker)
    text=replace_once(text,'  end subroutine time_step',
        f"    if ({condition}) call restart_probe('step-end')\n"+
        '  end subroutine time_step')
if args.block:
    text=replace_once(text,'    writeback_before = block_domain_production_writeback_count()\n    if (checkpoint_required) then',
        "    call restart_probe('before-consumer-sync')\n"+
        '    writeback_before = block_domain_production_writeback_count()\n    if (checkpoint_required) then')
    text=replace_once(text,'    compatibility_ready = .true.',
                     "    call restart_probe('after-consumer-sync')\n    compatibility_ready = .true.")
path.write_text(text)
if args.physics:
    path=out/'src/physics/physics_simple.f90'
    text=replace_once(path.read_text(),'  use kind_mod,',
        '  use arch_mod, only : rank, glo_id\n  use shared_mod, only : istep_cumul\n  use kind_mod,')
    text=replace_once(text,'    ! Backwards Euler step on current column',
        "    call physics_probe('input')\n    ! Backwards Euler step on current column")
    text=replace_once(text,'    ! Assign solution at t+h',
        "    call physics_probe('output')\n    ! Assign solution at t+h")
    replay_keys=set()
    replays=[('U',*values) for values in args.physics_replay_u]
    replays += [(values[0],*map(int,values[1:])) for values in args.physics_replay_velocity]
    for field,step,domain,node,level,bits in replays:
        if (field not in ('U','V','W') or min(step,domain,node)<0 or level<1 or
                not 0<=bits<2**32 or (bits & 0x7f800000)==0x7f800000):
            raise ValueError('Unsupported replay address or nonfinite float32 bits')
        signed_bits=bits if bits<2**31 else bits-2**32
        key=(field,step,domain,node,level)
        if key in replay_keys: raise ValueError('Duplicate replay address')
        replay_keys.add(key)
        replay=f'''    ! Explicitly requested causal perturbation, never a production correction.
    if (istep_cumul=={step}.and.glo_id(rank+1,d)=={domain}.and.id=={node}) then
       if ({level}>zlevels) error stop 'physics replay level invalid'
       write(6,'(a,3(i0,1x),2(es24.16,1x))') 'PHYSICS REPLAY {field}: ', &
            istep_cumul,glo_id(rank+1,d),id,phys_{field}({level}),transfer({signed_bits},0.0_sp)
       phys_{field}({level})=transfer({signed_bits},0.0_sp)
    end if
'''
        text=replace_once(text,"    call physics_probe('input')",replay+"    call physics_probe('input')")
    probe=Path(__file__).with_name('physics_probe.inc').read_text()
    text=replace_once(text,'  contains\n    subroutine pack_physics_vars',
        '  contains\n'+probe+'\n    subroutine pack_physics_vars')
    path.write_text(text)
import hashlib,json
(out/'probe-identity.json').write_text(json.dumps({
    'source':str(source),'block':args.block,'phase_step':args.phase_step,'phase_range':args.phase_range,
    'physics':args.physics,'physics_replay_u':args.physics_replay_u,
    'physics_replay_velocity':args.physics_replay_velocity,
    'modified_source_sha256':{str(p.relative_to(out)):hashlib.sha256(p.read_bytes()).hexdigest()
        for p in (out/'src/main.f90',out/'src/restart_probe.f90',out/'src/physics/physics_simple.f90')}
},indent=2)+'\n')
print(out)

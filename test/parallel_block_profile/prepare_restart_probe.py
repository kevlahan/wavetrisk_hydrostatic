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
parser.add_argument('--integrator', choices=('RK3','RK4'),
                    help='Override integrator in the isolated source copy only')
args = parser.parse_args()
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
if args.phase_step:
    condition='.or.'.join(f'istep_cumul=={step}' for step in args.phase_step)
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
print(out)

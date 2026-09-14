"""Freeze additive legacy phase instrumentation; used only by kit preparation."""
from pathlib import Path
import re
import shutil
import subprocess
from prepare_legacy_profile import boundary_patch, replace_once, scope, wrap_calls


def instrument(repo, out):
    helper=Path(__file__).resolve().parent
    shutil.copy2(repo/'src/parallel_block_profile.f90',out/'src/parallel_block_profile.f90')
    shutil.copy2(helper/'legacy_profile_io.f90',out/'src/legacy_profile_io.f90')
    make=out/'Makefile'
    sources=('SRC = parallel_block_profile.f90','      legacy_profile_io.f90','      kind.f90')
    make.write_text(replace_once(make.read_text(),'SRC = kind.f90',(' '+chr(92)+'\n').join(sources)))
    subprocess.run(['git','apply','--check','-'],cwd=out,input=boundary_patch(repo),check=True)
    subprocess.run(['git','apply','-'],cwd=out,input=boundary_patch(repo),check=True)
    for filename in ('main', 'multi_level', 'adapt', 'wavelet', 'time_integr'):
        path = out / f'src/{filename}.f90'
        text = path.read_text()
        text = re.sub(r'^(module \w+\s*)$', r'\1\n  use parallel_block_profile_mod', text, count=1, flags=re.M)
        if filename == 'main':
            text = scope(text, 'time_step', '    ! New time step', 'DP_STEP')
            text = wrap_calls(text, 'dt_step', 'DP_DYNAMICS')
            text = wrap_calls(text, 'physics_simple_step', 'DP_PHYSICS')
            text = wrap_calls(text, 'remap_vertical_coordinates', 'DP_REMAP')
            text = wrap_calls(text, 'write_and_export', 'DP_OUTPUT')
            text = wrap_calls(text, 'restart', 'DP_RESTART')
        elif filename == 'multi_level':
            text = scope(text, 'trend_ml', '    call update_bdry (q, NONE, 967)', 'DP_SHARED')
            for name, region in [('basic_operators', 'DP_BASIC'), ('cal_scalar_trend', 'DP_TENDENCY'),
                                 ('velocity_trend_source', 'DP_NATIVE_VELOCITY'),
                                 ('velocity_trend_grad', 'DP_NATIVE_VELOCITY'),
                                 ('cpt_or_restr_flux', 'DP_SCALAR_REPLAY'),
                                 ('cpt_or_restr_scalar', 'DP_SCALAR_REPLAY')]:
                text = wrap_calls(text, name, region)
        elif filename == 'adapt':
            text = scope(text, 'adapt', '    n_patch_old = grid%patch%length', 'DP_ADAPT')
            text = scope(text, 'WT_after_step', '    call zero_float (wavelet)', 'DP_WT_DRIVER')
            text = wrap_calls(text, 'compress_wavelets', 'DP_COMPRESSION')
        elif filename == 'wavelet':
            text = scope(text, 'inverse_wavelet_transform', '    if (present(jmin_in)) then', 'DP_INVERSE')
        elif filename == 'time_integr':
            text = wrap_calls(text, 'RK_sub_step', 'DP_RK_ASSEMBLE')
        path.write_text(text)

    path = out / 'test/climate/climate.f90'
    text = replace_once(path.read_text(), 'program climate\n  !',
                        'program climate\n  use legacy_profile_io_mod\n  !')
    text = replace_once(text, '  do while (time < time_end)',
                        '  call legacy_profile_begin\n  do while (time < time_end)')
    text = replace_once(text, '  call finalize', '  call legacy_profile_end(rank)\n  call finalize')
    path.write_text(text)

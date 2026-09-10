#!/usr/bin/env python3
"""Instrument a NEW experimental copy of a verified, built main archive.

Never edit a Git checkout or the baseline archive. Exact source anchors and
the existing boundary-only instrumentation patch constrain the transformation.
This is NOT an authoritative legacy binary; retain the unmodified comparator.
"""
import argparse
import difflib
import json
from pathlib import Path
import re
import shutil
import subprocess

from build_legacy import tracked_entries, verify
from experiment import LEGACY_REF, digest


def replace_once(text, old, new):
    if text.count(old) != 1:
        raise ValueError(f"Expected one exact anchor: {old!r}")
    return text.replace(old, new, 1)


def wrap_calls(text, name, region):
    """Wrap only unconditional single-line CALLs; preserve call and order exactly."""
    pattern = re.compile(r'^( +)(call ' + re.escape(name) + r'\b[^\n]*)$', re.M | re.I)
    def wrap(match):
        indent, statement = match.groups()
        if '&' in statement or ';' in statement:
            raise ValueError(f"Unsupported continued/compound call: {statement}")
        return (f'{indent}call detail_enter({region})\n{indent}{statement}\n'
                f'{indent}call detail_leave({region})')
    text, count = pattern.subn(wrap, text)
    if not count:
        raise ValueError(f'No calls to {name}')
    return text


def scope(text, name, first, region):
    """Bracket a known routine without early returns/internal subprograms."""
    pattern = re.compile(r'^  subroutine ' + name + r'\b.*?^  end subroutine ' + name + r'\b', re.M | re.S | re.I)
    matches = list(pattern.finditer(text))
    if len(matches) != 1:
        raise ValueError(f'Expected one routine {name}')
    m = matches[0]
    body = m.group()
    if re.search(r'^\s*(?! !).*\breturn\b|^\s*contains\b', '\n'.join(
            line for line in body.splitlines() if not line.lstrip().startswith('!')), re.M | re.I):
        raise ValueError(f'Unsupported control flow in {name}')
    body = replace_once(body, first, f'    call detail_enter({region})\n' + first)
    body = re.sub(r'(?im)^  end subroutine', f'    call detail_leave({region})\n  end subroutine', body)
    return text[:m.start()] + body + text[m.end():]


def prepare(repo, baseline, out):
    repo, baseline, out = repo.resolve(), baseline.resolve(), out.resolve()
    if out.exists() or out.is_relative_to(repo) or out.is_relative_to(baseline):
        raise ValueError('Output must be a NEW directory outside repository and baseline')
    if (baseline / '.git').exists():
        raise ValueError('Baseline must be an archived build, not a checkout')
    entries = tracked_entries(repo, LEGACY_REF)
    verify(baseline, entries)
    shutil.copytree(baseline, out, symlinks=True)
    if (out / 'legacy-build.json').exists():
        (out / 'legacy-build.json').rename(out / 'baseline-build.json')
    helper = Path(__file__).resolve().parent
    shutil.copy2(repo / 'src/parallel_block_profile.f90', out / 'src/parallel_block_profile.f90')
    shutil.copy2(helper / 'legacy_profile_io.f90', out / 'src/legacy_profile_io.f90')
    make = out / 'Makefile'
    sources = ('SRC = parallel_block_profile.f90', '      legacy_profile_io.f90', '      kind.f90')
    make.write_text(replace_once(make.read_text(), 'SRC = kind.f90', (' ' + chr(92) + '\n').join(sources)))

    # Transplant ONLY uncommitted timer/counter additions to the common boundary
    # routines. The rest of legacy comm_mpi (including ownership) stays untouched.
    patch = subprocess.check_output(['git', '-C', str(repo), 'diff', 'HEAD', '--', 'src/comm_mpi.f90'])
    for line in patch.decode().splitlines():
        if line.startswith('-') and not line.startswith('---') and line[1:].strip():
            raise ValueError('Boundary instrumentation patch removes executable text')
    subprocess.run(['patch', '--dry-run', '-p1'], cwd=out, input=patch, check=True)
    subprocess.run(['patch', '-p1'], cwd=out, input=patch, check=True)

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
    changed = {}
    diffs = []
    for _, _, name in entries:
        path, original = out / name, baseline / name
        if path.is_symlink():
            continue
        if path.read_bytes() != original.read_bytes():
            changed[name] = digest(path)
            diffs.extend(difflib.unified_diff(original.read_text().splitlines(True), path.read_text().splitlines(True),
                                            fromfile='a/' + name, tofile='b/' + name))
    (out / 'experimental-instrumentation.patch').write_text(''.join(diffs))
    (out / 'experimental-identity.json').write_text(json.dumps({
        'revision': LEGACY_REF, 'authoritative': False, 'instrumented_files': changed,
        'baseline_binary_sha256': digest(baseline / 'bin/climate')}, indent=2) + '\n')
    verify(baseline, entries)
    print(out)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, required=True)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    prepare(args.repo, args.baseline, args.out)

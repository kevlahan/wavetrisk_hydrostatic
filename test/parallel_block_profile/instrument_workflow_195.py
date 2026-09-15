"""Add identical coarse observers using pinned, checked source anchors."""
from pathlib import Path
import re
from prepare_legacy_profile import replace_once


def routine(text,name,first,tag):
    pat=re.compile(r'^  subroutine '+name+r'\b.*?^  end subroutine '+name+r'\b',re.M|re.S|re.I)
    matches=list(pat.finditer(text))
    if len(matches)!=1:raise ValueError('Missing unique routine '+name)
    m=matches[0];body=m.group()
    if re.search(r'^\s*return\b',body,re.M|re.I):raise ValueError('Early return '+name)
    body=replace_once(body,first,f'    call wf_enter({tag})\n'+first)
    body=re.sub(r'^  end subroutine',f'    call wf_leave({tag})\n  end subroutine',body,flags=re.M)
    return text[:m.start()]+body+text[m.end():]


def instrument(root,block):
    root=Path(root)
    original={n:(root/n).read_text() for n in ('Makefile','test/climate/climate.f90','src/main.f90','src/adapt.f90','src/time_integr.f90' if block else 'src/multi_level.f90')}
    (root/'src/workflow_profile_195.f90').write_bytes(Path(__file__).with_name('workflow_profile_195.f90').read_bytes())
    p=root/'Makefile';s=p.read_text();s=replace_once(s,'SRC = ','SRC = workflow_profile_195.f90 '+chr(92)+'\n      ');p.write_text(s)
    p=root/'test/climate/climate.f90';s=p.read_text()
    s,n=re.subn(r'^program climate\n','program climate\n  use workflow_profile_195_mod\n',s,count=1,flags=re.M)
    if n!=1:raise ValueError('Missing climate program')
    s=replace_once(s,'  do while (time < time_end)','  call wf_init\n  do while (time < time_end)')
    s=replace_once(s,'     call start_timing ; call time_step; call stop_timing',
                   '     call wf_begin\n     call start_timing ; call time_step; call stop_timing\n     call wf_end')
    s=replace_once(s,'  call finalize','  call wf_finish(rank)\n  call finalize');p.write_text(s)
    for name in ['main','adapt','time_integr' if block else 'multi_level']:
        p=root/f'src/{name}.f90';s=p.read_text()
        s,n=re.subn(r'^(module \w+[^\n]*\n)',r'\1  use workflow_profile_195_mod\n',s,count=1,flags=re.M)
        if n!=1:raise ValueError('Missing module '+name)
        if name=='main':
            # Headings delimit the same physical workflow in both versions.
            for phase,heading in [(1,'Dynamics time step'),(2,'Physics split step'),(3,'Grid adaptation'),(4,'Vertical remapping (after grid adaptation')]:
                lines=s.splitlines(True);hits=[i for i,v in enumerate(lines) if heading in v]
                if len(hits)!=1:raise ValueError('Ambiguous phase heading '+heading)
                lines.insert(hits[0],f'    call wf_phase({phase})\n');s=''.join(lines)
            anchor='    output_profile_start = &' if block else '    if (log_total_mass) call cal_total_mass (.false.)'
            # Existing block header also covers post-step diagnostics/clock advance.
            if block:
                # Only the time_step occurrence; other procedures also use output timers.
                start=s.index('  subroutine time_step');end=s.index('  end subroutine time_step',start)
                body=replace_once(s[start:end],anchor,'    call wf_phase(5)\n'+anchor);s=s[:start]+body+s[end:]
            else:s=replace_once(s,anchor,'    call wf_phase(5)\n'+anchor)
            # Flags do not require clock reads; work includes checkpoint restart.
            s=replace_once(s,'       call remap_vertical_coordinates','       call wf_flag(2)\n       call remap_vertical_coordinates')
            # Mark checkpoint writes in their routine, not initial restart.
            start=s.index('  subroutine write_checkpoint');end=s.index('  end subroutine write_checkpoint',start)
            body=s[start:end];anchor='    call restart'
            body=replace_once(body,anchor,'    call wf_flag(1)\n'+anchor);s=s[:start]+body+s[end:]
        elif name=='adapt':
            first='    call detail_enter(DP_WT_DRIVER)'
            s=routine(s,'WT_after_step',first,1)
        elif block:
            s=routine(s,'prepare_block_multistage_tendency','    call detail_enter(DP_CAPTURE)\n    call begin_block_scalar_divergence_capture',2)
        else:
            s=routine(s,'trend_ml','    call detail_enter(DP_SHARED)',2)
        p.write_text(s)

    # Prove that removing only observer additions restores every original byte.
    # No numerical statement, branch or expression may be changed by this tool.
    for name,before in original.items():
        after=(root/name).read_text()
        if name=='Makefile':
            after=replace_once(after,'SRC = workflow_profile_195.f90 '+chr(92)+'\n      ','SRC = ')
        else:
            after=''.join(line for line in after.splitlines(True) if not re.match(r'^\s*(?:use workflow_profile_195_mod|call wf_\w+)\b',line,re.I))
        if after!=before:raise ValueError('Observer transformation changes original source: '+name)

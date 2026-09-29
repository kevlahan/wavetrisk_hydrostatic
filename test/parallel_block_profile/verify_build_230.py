"""Exercise the real three-build orchestrator with explicitly nonexecuted compiler stubs."""
import json,tarfile,tempfile,subprocess
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import build_capability_230 as b

def verify(kit):
    with tempfile.TemporaryDirectory(prefix='stage230-build-workflow-') as tmp:
        root=Path(tmp);ref=root/'reference';ref.mkdir()
        with tarfile.open(kit/'stage230-source.tar.gz') as ar:ar.extractall(ref,filter='data')
        commands=[];real_run=subprocess.run
        def run(command,**kwargs):
            if command[0]=='mpirun':return SimpleNamespace(stdout='test MPI identity; no execution')
            if command[0]=='bash':
                obj=ref/'src/physics/simple_physics/phyparam/obj';obj.mkdir(exist_ok=True)
                for p in obj.glob('*.o'):p.unlink()
                for i in range(23):(obj/f'physics{i}.o').write_text('test object '+str(i))
            elif command[0]=='make':
                commands.append(command)
                opts={s.split('=',1)[0]:s.split('=',1)[1] for s in command if '=' in s}
                dst=kwargs['cwd']/opts['BIN_DIR'];dst.mkdir();(dst/'climate').write_text('nonexecuted test binary '+opts['BIN_DIR'])
                kwargs['stdout'].write('compiler stub for orchestration only\n')
            elif command[:2]==['git','apply']:return real_run(command,**kwargs)
            else:raise AssertionError(command)
            return subprocess.CompletedProcess(command,0)
        with patch.object(b.subprocess,'run',side_effect=run),patch.object(b.subprocess,'check_output',return_value='test compiler identity'):
            b.build(kit,root)
        record=json.loads((root/'capability-230-builds.json').read_text())
        assert record['status']=='passed' and set(record['builds'])==set(b.KEYS) and len(commands)==3
        assert len(record['physics_objects'])==23
        for key,v in record['builds'].items():
            kind=key.split('-',1)[1]
            expected='on'
            assert 'STAGE203='+expected in v['command'] and v['stage203']==expected
            assert 'FP_CONTRACT=off' in v['command'] and 'STAGE202=on' in v['command']
            assert v['observer']==kind.endswith('_observer')
            assert 'DEBUG='+('check' if kind=='checked' else 'false') in v['command']
    print('PASS: actual three-build workflow, selectors, source hashes and shared physics; compiler execution stubbed.')
if __name__=='__main__':
    import sys
    verify(Path(sys.argv[1]).resolve())

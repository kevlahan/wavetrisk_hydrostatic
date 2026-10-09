"""Preserve Slurm elapsed time even when accounting is delayed or unavailable."""
from datetime import datetime, timezone
import subprocess


def parse_accounting(stdout, job):
    rows=[]
    for line in stdout.splitlines():
        fields=line.strip().split('|')
        if not fields or fields[0]!=str(job):continue
        if len(fields)!=6 or not fields[4].isdigit():raise ValueError('Malformed Slurm job accounting')
        rows.append(dict(job_id=fields[0],state=fields[1],start=fields[2],end=fields[3],
                         elapsed_seconds=int(fields[4]),elapsed=fields[5]))
    if len(rows)>1:raise ValueError('Duplicate Slurm job accounting')
    return rows[0] if rows else None


def capture(job):
    command=['sacct','-j',str(job),'-X','-n','-P','-o','JobIDRaw,State,Start,End,ElapsedRaw,Elapsed']
    try:
        result=subprocess.run(command,capture_output=True,text=True,timeout=30)
        code,stdout,stderr=result.returncode,result.stdout,result.stderr
    except (OSError,subprocess.TimeoutExpired) as error:
        code,stdout,stderr=-1,'',str(error)
    return dict(command=command,captured_utc=datetime.now(timezone.utc).isoformat(),
                returncode=code,stdout=stdout,stderr=stderr,
                job=parse_accounting(stdout,job) if code==0 else None)

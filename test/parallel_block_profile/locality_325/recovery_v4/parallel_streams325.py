"""CLI adapter for bounded optimized capture shard construction."""
from pathlib import Path
import json,sys
from recover325 import save,require
from parallel_replay325 import cache_shards325,tasks_from_results325

def main(kit,original,scratch,output):
 anchors=json.loads((Path(kit)/'prior-native-input-anchors.json').read_text())
 def progress(kind,arm,index,count):print(f'{kind} {arm} index={index} streams={count}',flush=True)
 r=cache_shards325(tasks_from_results325(original),Path(scratch),workers=16,original_root=original,anchors=anchors,progress=progress)
 r['passed']=True;r['actual_capture_reads']=r.pop('capture_reads');save(output,r)
if __name__=='__main__':
 require(len(sys.argv)==5,'Usage: parallel_streams325.py KIT ORIGINAL SCRATCH OUTPUT');main(*sys.argv[1:])

#!/bin/bash
#grep -A 1 'd abs' --no-group-separator gcm.log | grep -v 'd abs'  > gcm.output
grep -A 1 'd abs' --no-group-separator gcm.log | grep -v 'd abs' | head > gcm.output
diff -s gcm.output gcm.output.$1 && diff gcm.output gcm.output.$1

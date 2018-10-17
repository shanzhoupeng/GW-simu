import os
import sys
for i in range(0,5):
    a="pkurun-cnlong 1 20 python clustergwtest.py %s %s %s %s"%(i,sys.argv[1],sys.argv[2],5)
    os.system(a)
    os.system('sleep 1')

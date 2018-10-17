import os
import sys
import numpy as np
def listnpy(path, list_name):
    for file in os.listdir(path):  
        file_path = os.path.join(path, file)  
        if os.path.isdir(file_path):  
            listnpy(file_path, list_name)  
        else: 
            if os.path.splitext(file_path)[1]=='.npy':
                list_name.append(file_path)
pathname='/home/chenxian_pkuhpc/lustre1/test'
lm=[]
listnpy(pathname,lm)
g=np.array(lm)
data=np.load(g[0])
for i in g[1:g.size]:
    a=np.load(i)
    data=np.concatenate((data,a))
dmax=data[:,2].max()
dmin=data[:,2].min()
filename='%s-%s'%(dmin,dmax)
np.save(filename,data)
print "done!"

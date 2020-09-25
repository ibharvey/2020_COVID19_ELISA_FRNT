# coding: utf-8

import sys

t1 = sys.argv[1]
t2 = sys.argv[2]
from os import listdir
files = sorted([x for x in listdir('.') if t1 in x])
files2 = sorted([x for x in listdir('.') if t2 in x])

import numpy as np
for i in range(len(files)):
    f1 = [x.split(',') for x in open(files[i],'r').readlines()]
    f2 = [x.split(',') for x in open(files2[i],'r').readlines()]
    fw = open(files[i][:-4] + '_avg.csv','w+')
    fw.write(','.join(f1[0]))
    for j in range(1,len(f1)):
        temp = [str(np.mean([int(f1[j][k]),int(f2[j][k])])) for k in range(1,len(f1[j]))]
        temp = [x.replace('\n','') for x in temp]
        fw.write('{0},{1}\n'.format(f1[j][0],','.join(temp)))
    fw.close()
    

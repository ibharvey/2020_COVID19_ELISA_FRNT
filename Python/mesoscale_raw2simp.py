#! /bin/python3

import sys

out_files = ['320','1280','10240']

for o in range(len(out_files)):
    dilution_series = {}
    sample_names = []
    frl = open(sys.argv[1],'r').readlines()
    for a in range(3,93,9):
        f1 = frl[a:a+9]
        antigen = f1[0].split(',')[1]
        dilution_series[antigen] = []
        for i in range(4):
            for f in f1[1:]:
                temp = f.split(',')[0+4*i:4+4*i]
                if a == 3:
                    sample_names.append(temp[0])
                #Pick dilution
                dilution_series[antigen].append(temp[o+1])
        
    fw = open('{0}_simp_{1}.csv'.format(sys.argv[1][:-4],out_files[o]),'w+')
    fw.write(",")
    fw.write('{0}\n'.format(','.join(dilution_series.keys())))
    for i in range(len(sample_names)):
        fw.write('{0},'.format(sample_names[i]))
        sample_set = []
        for k in dilution_series.keys():
            sample_set.append(dilution_series[k][i])
        sample_set = [x.replace('\n','') for x in sample_set]
        fw.write('{0}\n'.format(','.join(sample_set)))
    
    fw.close()

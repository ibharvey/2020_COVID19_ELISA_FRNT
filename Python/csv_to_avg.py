# coding: utf-8

import numpy as np

frl = open('20200420_ELISA_SARS2_Uchcago_16_samples.csv','r').readlines() #replicate 1
brl = open('20200424_ELISA_SARS2_UChicago_16_samples_reformat.csv','r').readlines() #replicate 2

names = []
dt = [[] for x in range(12)]
temp_dt = [[] for x in range(12)]
bemp_dt = [[] for x in range(12)]
counter = 0
bcounter = 0
boolDNR = False
boolDNR2 = False
for i in range(len(frl)):
    line = frl[i]
    if counter == 0:
        a_name = line.replace('\n','').split(',')[0]
        if a_name == 'DNR1':
            if boolDNR:
                break
            else:
                boolDNR = True
        names.append(a_name)
        counter += 1
    elif counter == 9:
        for j in range(12):
            dt[j].append(temp_dt[j])
        temp_dt = [[] for x in range(12)]
        counter = 0
    else:
        temp = line.replace('\n','').split(',')
        for j in range(1,13):
            temp_dt[j-1].append(temp[j])
        counter += 1
    bine = brl[i]
    if bcounter == 0:
        b_name = bine.replace('\n','').split(',')[0]
        if b_name == 'DNR1':
            if boolDNR2:
                break
            else:
                boolDNR2 = True
        bcounter += 1
    elif bcounter == 9:
        for j in range(12):
            dt[j].append(bemp_dt[j])
        bemp_dt = [[] for x in range(12)]
        bcounter = 0
    else:
        temp = bine.replace('\n','').split(',')
        for j in range(1,13):
            bemp_dt[j-1].append(temp[j])
        bcounter += 1
        
for j in range(12):
    dt[j].append(temp_dt[j])
    
for j in range(12):
    dt[j].append(bemp_dt[j])
    
bemp_dt = [[] for x in range(12)]
temp_dt = [[] for x in range(12)]
counter = 0
bcounter = 0
for i in range(len(frl2)):
    line = frl2[i]
    if counter == 0:
        a_name = line.replace('\n','').split(',')[0]
        names.append(a_name)
        counter += 1
    elif counter == 9:
        for j in range(12):
            dt[j].append(temp_dt[j])
        temp_dt = [[] for x in range(12)]
        counter = 0
    else:
        temp = line.replace('\n','').split(',')
        for j in range(1,13):
            temp_dt[j-1].append(temp[j])
        counter += 1
    bine = brl2[i]
    if bcounter == 0:
        b_name = bine.replace('\n','').split(',')[0]
        bcounter += 1
    elif bcounter == 9:
        for j in range(12):
            dt[j].append(bemp_dt[j])
        bemp_dt = [[] for x in range(12)]
        bcounter = 0
    else:
        temp = bine.replace('\n','').split(',')
        for j in range(1,13):
            bemp_dt[j-1].append(temp[j])
        bcounter += 1
        
                
for j in range(12):
    dt[j].append(temp_dt[j])
    
for j in range(12):
    dt[j].append(bemp_dt[j])
    
dtn = [np.array(x) for x in dt]
dtni = [np.transpose(x) for x in dtn]
fw = open('ELISA_round_1_and_2.csv','w+')
for f in dtni:
    fw.write('\n')
    for line in f:
        fw.write(','.join(line))
        fw.write('\n')
        
fw.close()

frl = open('20200424_ELISA_Chicago_round_1_and_2.csv','r').readlines()
fw = open('20200424_ELISA_Chicago_round_1_and_2_avg.csv','w+')                                              
for line in frl:                                          
    temp = line.replace('OVRFLW','4').split(',')
    temp2 = [np.mean([float(temp[i]),float(temp[i+1])]) for i in range(0,len(temp)-1,2)]
    fw.write(','.join(map(str,temp2)))
    fw.write('\n')   

fw.close()
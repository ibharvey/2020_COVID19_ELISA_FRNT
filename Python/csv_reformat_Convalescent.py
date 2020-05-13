# coding: utf-8

import numpy as np
import re

d = "/home/ian/Dropbox/Fremont/Bioinformatics/ELISA/AliConvalescent/"

frl = open(d+'20200503_ELISA_80_convalescent_samples.csv','r').readlines() #replicate 1
brl = open(d+'20200507_ELISA_80_convalescent_samples2.csv','r').readlines() #replicate 2

names = []
dt = [[] for x in range(12)]
temp_dt = [[] for x in range(12)]
bemp_dt = [[] for x in range(12)]
counter = 0
bcounter = 0
for i in range(len(frl)):
    line = frl[i]
    if counter == 0:
        a_name = line.replace('\n','').split(',')[0]
        start,end = map(int,re.findall("\d+\-\d+",a_name)[0].split('-'))
        new_names = ["ConvSera_{0}".format(x) for x in range(start,end+1)]
        names += new_names
        counter += 1
    elif counter == 9:
        counter = 0
    elif counter % 2 == 0:
        temp = line.replace('\n','').split(',')
        for j in range(1,13):
            temp_dt[j-1].append(temp[j])
            dt[j-1].append(temp_dt[j-1])
        temp_dt = [[] for x in range(12)]
        counter += 1
    else:
        temp = line.replace('\n','').split(',')
        for j in range(1,13):
            temp_dt[j-1].append(temp[j])
        counter += 1
    bine = brl[i]
    if bcounter == 0:
        b_name = bine.replace('\n','').split(',')[0]
        bcounter += 1
    elif bcounter == 9:
        bcounter = 0
    elif bcounter % 2 == 0:
        bemp = bine.replace('\n','').split(',')
        for j in range(1,13):
            bemp_dt[j-1].append(bemp[j])
            dt[j-1].append(bemp_dt[j-1])
        bemp_dt = [[] for x in range(12)]
        bcounter += 1
    else:
        bemp = bine.replace('\n','').split(',')
        for j in range(1,13):
            bemp_dt[j-1].append(bemp[j])
        bcounter += 1
        
'''for j in range(12):
    dt[j].append(temp_dt[j])
    
for j in range(12):
    dt[j].append(bemp_dt[j])'''
    
'''bemp_dt = [[] for x in range(12)]
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
    dt[j].append(bemp_dt[j])'''
    
dtn = [np.array(x) for x in dt]
dtni = [np.transpose(x) for x in dtn]
fw = open(d+'20200513_ELISA_Convalescent_2_Replicates.csv','w+')
for f in dtni:
    fw.write('\n')
    for line in f:
        fw.write(','.join(line))
        fw.write('\n')
        
fw.close()

# Average the two replicates
frl = open(d+'20200513_ELISA_Convalescent_2_Replicates.csv','r').readlines()
fw = open(d+'20200513_ELISA_Convalescent_2_Replicates_avg.csv','w+')                                              
for line in frl:                                          
    temp = line.replace('OVRFLW','4').split(',')
    temp2 = [np.mean([float(temp[i]),float(temp[i+1])]) for i in range(0,len(temp)-1,2)]
    fw.write(','.join(map(str,temp2)))
    fw.write('\n')   

fw.close()

# Split each dataset by dilution factor
frl = open(d+'20200513_ELISA_Convalescent_2_Replicates_avg.csv','r').readlines()
frl = frl[1:]
fw = [None for f in range(2)]
dilutions = [320,1280]
for i in range(len(fw)):                                                                   
    fw[i] = open(d+"Conv_" + str(dilutions[i]) + '.csv','w+')

dat = [[] for x in range(2)]
counter = 0                                                           
for i in range(len(frl)):                                                                  
    if counter == 2:                                           
        counter = 0                                       
    else:                                                                             
        dat[counter].append(frl[i].replace('\n','').split(','))                            
        counter += 1

ndat = [np.array(x) for x in dat]
ndati = [np.transpose(x) for x in ndat]

for f in range(2):
    fw[f].write('BSA, HA, MERS Trimer, SARS1 Trimer, SARS2 Trimer, SARS1 mRBD, SARS2 mRBD, SARS1 bRBD, SARS2 bRBD, SARS2 NP-RNA, SARS2 NP-FL, SARS2 ORF8\n')

for f in range(len(ndati)):                                           
    for line in ndati[f]:                                                                  
        fw[f].write(','.join(line))                            
        fw[f].write('\n')                                 
    fw[f].write('\n')

for f in range(2):                                                    
    fw[f].close()
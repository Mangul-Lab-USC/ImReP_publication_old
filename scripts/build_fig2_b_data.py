import pandas as pd
import os
import numpy as np
from math import log as ln
# need to rerun mixcr for '009-0231', for now wi exclude this 

samples=['009-0192',
'009-0184',
'009-0148',
'009-0171',
'009-0203',
'009-0186',
'009-0174',
'009-0109',
'009-0249',
'009-0112',
'009-0122',
'009-0202',
'009-0103']

dict={}
newDF = pd.DataFrame()



for s in samples:
    dict[s]={}
    
for s in samples:
    dict[s]['true']=newDF
    dict[s]['imrep']=newDF
    dict[s]['mixcr']=newDF
    dict[s]['true.imrep']=newDF
    dict[s]['true.mixcr']=newDF
    dict[s]['imrep.mixcr']=newDF
    

directory_true='../raw_data/BCR-Seq/'




for root, dirs, files in os.walk(directory_true):
    for file in files:
        if file.endswith('.cdr3'):
            for s in samples:
                if s in file:
                    dict[s]['true']=pd.read_csv(directory_true+file)


directory_imrep='../raw_data/imrep/'



for root, dirs, files in os.walk(directory_imrep):
    for file in files:
        if file.endswith('.csv'):
            for s in samples:
                if s in file:
                    dict[s]['imrep']=pd.read_csv(directory_imrep+file)
                    dict[s]['imrep'] = dict[s]['imrep'].rename(columns={'relative.frequency': 'FREQ.imrep'})

directory_mixcr='../raw_data/mixcr/'



for root, dirs, files in os.walk(directory_mixcr):
    for file in files:
        if file.endswith('.clean.cdr3'):
            for s in samples:
                if s in file:
                    dict[s]['mixcr']=pd.read_csv(directory_mixcr+file)
                    dict[s]['mixcr'] = dict[s]['mixcr'].rename(columns={'FREQ': 'FREQ.mixcr'})


for s in samples:

    if dict[s]['imrep'].size==0:
        dict[s]['true.imrep']=pd.DataFrame()
    if dict[s]['mixcr'].size==0:
        dict[s]['true.mixcr']=pd.DataFrame()
        
    if dict[s]['imrep'].size!=0 and dict[s]['mixcr'].size!=0:
    
        dict[s]['true.imrep']=pd.merge( dict[s]['true'],dict[s]['imrep'], on='CDR3') #imrep and true
        dict[s]['true.mixcr']=pd.merge( dict[s]['true'],dict[s]['mixcr'], on='CDR3') #mixcr and true
        dict[s]['imrep.mixcr']=pd.merge( dict[s]['imrep'],dict[s]['mixcr'], on='CDR3') # mixcr and imrep
    

file=open('../raw_data/BCR.Seq.validation.csv',"w")
file.write("ID,n_true,n_imrep_true,n_mixcr_true,max_true,min_true,sum_imrep_true,max_imrep_true,min_imrep_true,sum_mixcr_true,max_mixcr_true,min_mixcr_true")
file.write("\n")

for s in samples:
    
    
    
    sum_imrep_true='0'
    max_imrep_true='0'
    min_imrep_true='0'
    sum_mixcr_true='0'
    max_mixcr_true='0'
    min_mixcr_true='0'
    
    sdi_imrep='0'
    sdi_mixcr='0'
    
    ID=s
    n_true=str(dict[s]['true']['FREQ'].size)
    
    if dict[s]['true.imrep'].size==0:
        n_imrep_true='0' 
    else:
        n_imrep_true=str(dict[s]['true.imrep']['FREQ.imrep'].size)
        
    if dict[s]['true.mixcr'].size==0:
        n_mixcr_true='0' 
    else:
        n_mixcr_true=str(dict[s]['true.mixcr']['FREQ.mixcr'].size) 
        
        
    max_true=str(dict[s]['true']['FREQ'].max()*100)
    min_true=str(dict[s]['true']['FREQ'].min()*100)
    
    
    #sum max min
    if dict[s]['true.imrep'].size==0:
        sum_imrep_true='0'
        max_imrep_true='0'
        min_imrep_true='0'
    else:
        sum_imrep_true=str(dict[s]['true.imrep']['FREQ'].sum()*100)
        max_imrep_true=str(dict[s]['true.imrep']['FREQ'].max()*100)
        min_imrep_true=str(dict[s]['true.imrep']['FREQ'].min()*100)
        
    if dict[s]['true.mixcr'].size==0:
        sum_mixcr_true='0'
        max_mixcr_true='0'
        min_mixcr_true='0'
    else:
        sum_mixcr_true=str(dict[s]['true.mixcr']['FREQ'].sum()*100)
        max_mixcr_true=str(dict[s]['true.mixcr']['FREQ'].max()*100)
        min_mixcr_true=str(dict[s]['true.imrep']['FREQ'].min()*100)
        
    
    
    #print sum_imrep,sum_mixcr
    
    file.write(ID+","+n_true+","+n_imrep_true+","+n_mixcr_true+","+max_true+","+min_true+","+sum_imrep_true+","+max_imrep_true+","+min_imrep_true+","+sum_mixcr_true+","+max_mixcr_true+","+min_mixcr_true)
    file.write("\n")

file.close()

def p(n, N):
    """ Relative abundance """
    if n is  0:
        return 0
    else:
        return (float(n)/N) * ln(float(n)/N)

def sdi(data):
    N = sum(data)
    if len(data)==1:
        return 0.0
    return -sum(p(n, N) for n in data if n is not 0)


bigdata_imrep = pd.DataFrame()

for s in samples:
    bigdata_imrep = bigdata_imrep.append(dict[s]['true.imrep'], ignore_index=True)


bigdata_imrep.to_csv(path_or_buf='../summary_data/Figure2c_data.csv', index=False)

bigdata_mixcr = pd.DataFrame()

for s in samples:
    bigdata_mixcr = bigdata_mixcr.append(dict[s]['true.mixcr'], ignore_index=True)

bigdata_mixcr.to_csv(path_or_buf='../summary_data/Figure2d_data.csv', index=False)

file=open('../summary_data/Figure2b_data.csv',"w")
file2=open('../summary_data/portion.captured.items.csv',"w")

file.write("th,imrep_portion_frequency,mixcr_portion_frequency\n")
file2.write("th,imrep_portion_items,mixcr_portion_items\n")

previous=0

for th in np.linspace(0.1,0,2000):
    
    
    
    k_imrep=0
    k_mixcr=0
    k_true=0
    
    s_imrep=0
    s_mixcr=0
    s_true=0
    for s in samples:
        if dict[s]['true.imrep'].size !=0:
            freq_local_imrep=dict[s]['true.imrep']['FREQ']
            for i in freq_local_imrep:
                if i>=th:
                    k_imrep+=1.0
                    s_imrep+=i
        if dict[s]['true.mixcr'].size !=0:
            freq_local_mixcr=dict[s]['true.mixcr']['FREQ']
            for i in freq_local_mixcr:
                if i>=th:
                    k_mixcr+=1.0
                    s_mixcr+=i
        if dict[s]['true'].size !=0:
            freq_local_true=dict[s]['true']['FREQ']
            for i in freq_local_true:
                if i>=th:
                    k_true+=1.0
                    s_true+=i


    imrep_portion=str((k_imrep/k_true))
    mixcr_portion=str((k_mixcr/k_true))
    file2.write(str(th)+","+imrep_portion+","+mixcr_portion)
    file2.write("\n")
    
    imrep_portion2=str((s_imrep/s_true))
    mixcr_portion2=str((s_mixcr/s_true))
    file.write(str(th)+","+imrep_portion2+","+mixcr_portion2)
    file.write("\n")
    
    
file.close()
file2.close()

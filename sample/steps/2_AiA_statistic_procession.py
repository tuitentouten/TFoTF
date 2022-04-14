# -*- coding: utf-8 -*-
"""
Created on Thu Jul 15 19:16:10 2021

@author: wangf
"""


import pandas as pd
import numpy as np
import os

os.chdir('F:/OneDrive/1_工作和学习/manuscripts/TF预测/raw data/fig5/scripts/sample/sources')

gene_selected = 'STAT1'

df = pd.read_csv('out/AiAcancer_' + gene_selected + '.csv')

source_df = pd.read_csv('expression/CHOL_expression.csv')#读取基因名列表
source_df = source_df.set_index('sample')

source = source_df.index.to_list()


ENSG_list = df.ID.to_list()

gene_list = []

cc = 0
for ENSG in source:
    print(gene_selected)
    cc += 1
    print(cc/len(source))
    if ENSG in ENSG_list:
        gene_list.append(ENSG)


df = df.drop(df[df.p_value == 'False'].index)

df = df.apply(pd.to_numeric,errors='ignore')

df = df.drop(df[df.p_value > 0.05].index)


NAME = []
ID = []

MAXTCGA = []
MAXR = []

P_MAXTCGA = []
P_MAXR = []
N_MAXTCGA = []
N_MAXR = []

MINTCGA = []
MINR = []

MEANR = []
MEDIANR = []
SDR = []

ALLS = []
SIGNIS = []

P_SIGNIS = []
N_SIGNIS = []

counter = 0

for gene in gene_list:
    
    counter+=1
    print(gene_selected)
    print('%.2f%%' % ((counter/len(gene_list) * 100)))

    try:
        

        df_s = df[df['ID'] == gene]
        
        df_s = df_s.loc[df_s['r_value'] != 'False']
        
        df_p = df_s[df_s['r_value'] > 0]
        df_n = df_s[df_s['r_value'] < 0]
        
        P_SIGNIS_count = len(df_s[df_s['r_value'] >= 0.4])#
        N_SIGNIS_count = len(df_s[df_s['r_value'] <= -0.4])#
        
        p_l = df_p.r_value.to_list()
        n_l = df_n.r_value.to_list()
        
        pTCGA = df_p.TCGA.to_list()
        nTCGA = df_n.TCGA.to_list()
        
        pTCGA_tup = tuple(pTCGA)
        nTCGA_tup = tuple(nTCGA)
        
        p_tup = tuple(p_l)
        n_tup = tuple(n_l)
        
        if len(p_l) != 0:
            max_p = max(p_tup)#
            max_pNO = p_tup.index(max_p)
            max_pTCGA = pTCGA[max_pNO]#
        else:
            max_p = 'False'
            max_pTCGA = 'False'
        
        if len(n_l) != 0:
            max_n = min(n_tup)#
            max_nNO = n_tup.index(max_n)
            max_nTCGA = nTCGA[max_nNO]#
        else:
            max_n = 'False'
            max_nTCGA = 'False'
        
    #获取所需的数据
        r_list = df_s.r_value.tolist()
        
        absr_list = []
        for r in r_list:
            absr = abs(r)
            absr_list.append(absr)
            
        signi_list = []
        for absr in absr_list:
            if absr >= 0.4:
                signi_list.append(absr)
                
        
        
        TCGA_list = df_s.TCGA.tolist()
        
        absr_tup = tuple(absr_list)
        TCGA_tup = tuple(TCGA_list)
        
        
        all_samples = len(TCGA_tup)#
        signi_samples = len(signi_list)#
            
        max_s = max(absr_tup)#
        max_No = absr_tup.index(max(absr_tup))
        max_TCGA = TCGA_tup[max_No]#
        
        min_s = min(absr_tup)#
        min_No = absr_tup.index(min(absr_tup))
        min_TCGA = TCGA_tup[min_No]#
        
        mean_s = np.mean(absr_tup)#
        
        median_s = np.median(absr_tup)#
        
        SD_s = np.std(absr_tup)#
        
        name_s = df_s.gene.tolist()[0]#
        
        NAME.append(name_s)
        ID.append(gene)
        
        MAXTCGA.append(max_TCGA)
        MAXR.append(max_s)
        
        P_MAXTCGA.append(max_pTCGA)
        P_MAXR.append(max_p)
        
        N_MAXTCGA.append(max_nTCGA)
        N_MAXR.append(max_n)
        
        MINTCGA.append(min_TCGA)
        MINR.append(min_s)
        
        MEANR.append(mean_s)
        MEDIANR.append(median_s)
        SDR.append(SD_s)
        
        ALLS.append(all_samples)
        SIGNIS.append(signi_samples)
        
        P_SIGNIS.append(P_SIGNIS_count)
        N_SIGNIS.append(N_SIGNIS_count)
        
        
    
    except:
        pass
    
    dt = {'Name':NAME,'ID':ID,
          'Max_TCGA':MAXTCGA,'Max_R':MAXR,
          
          'Positive_MAX_TCGA':P_MAXTCGA,'Positive_MAX_R':P_MAXR,
          'Negative_MAX_TCGA':N_MAXTCGA,'Negative_MAX_R':N_MAXR,
          
          'Min_TCGA':MINTCGA,'Min_R':MINR,
          'Mean_R':MEANR,'Median_R':MEDIANR,'SD_R':SDR,
          
          'all_samples':ALLS,'Significant_samples':SIGNIS,
          
          'Positive_sig_samples':P_SIGNIS,'Negative_sig_samples':N_SIGNIS}        

        
dataF = pd.DataFrame(dt)

dataF.to_csv('out/AiAcancer_statat_' + gene_selected + '.csv')
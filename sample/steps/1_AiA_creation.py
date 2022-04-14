# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 02:25:26 2021

@author: wangf
"""

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
import pandas as pd
import scipy.stats as st
import math
import os

os.chdir('F:/OneDrive/1_工作和学习/manuscripts/TF预测/raw data/fig5/scripts/sample/sources')
notation = pd.read_csv('ENSGs_collection_single.csv')
notation_ID_first=notation.set_index(['Gene stable ID'])  
notation_name_first=notation.set_index(['Gene name'])

CHOL  = pd.read_csv('expression/CHOL_expression.csv')
DLBC  = pd.read_csv('expression/DLBC_expression.csv')
USC  = pd.read_csv('expression/USC_expression.csv')
#...
#load other cancer datasets one by one. These data were NOT provided in this sample
    
#----
def shan(lt, item):
    count = lt.count(item)
    if count == 0:   
        return
    for i in range(count):
        lt.remove(item)
    return lt

def logize(data):
    lenth=len(data)
    for i in range(lenth):
        data[i]=math.log(float(data[i]))
    return data  
#----Logarithmic

def pick(data):
    if data[0:4]=='ENSG':
        data=notation_ID_first.at[data,'Gene name']
    else:
        data=notation_name_first.at[data,'Gene stable ID'] 
    return data
#----translation

target_gene = 'STAT1'

file_li = []
for file in os.listdir('expression'):
    file_li.append(file)

AC = CHOL.set_index('sample')
ENSG_list = AC.index.to_list()

counter = 0

for gene in ENSG_list[:10]:
      
    counter+=1
    
    print('%.2f%%' % ((counter/len(ENSG_list) * 100)))
    
    try:
        gene_sym = pick(gene)
    except:
        gene_sym = 'False'
    
    ENSG_col = []
    gene_col = []
    TCGA_col = []
    Co_col = []
    p_col = []
    len_col = []
    
    
    y_name_li = []
    
    y_name_li.append(gene)
    
    for file_name in file_li:
        
        
            TCGA_ID = file_name.split('_')[0]
            
            code = 'GeneEX = '+ TCGA_ID
            
            exec(code)
            
            x_name = target_gene
            
            x_ID = x_name#x axis gene is waht
            
            EX2 = GeneEX.set_index(['sample'])#ENSG id as index
    
            xvalue_list = EX2.reindex([pick(x_ID)]).values.tolist()[0]
            
            xvalue_tuple = tuple(xvalue_list)
            
            for i in y_name_li:
                
                xvalue_list=list(xvalue_tuple)
                
                
                try:
            
                    yi=i

                    yvalue=EX2.reindex([yi])
            
                    yvalue_list=yvalue.values.tolist()
                           
                    yvalue_list=yvalue_list[0]

                    for j in range(0,len(xvalue_list)):#remove 0
                        if type(xvalue_list[j])==str:
                            xvalue_list[j]=0
                            yvalue_list[j]=0
                    
                        if type(yvalue_list[j])==str:
                            xvalue_list[j]=0
                            yvalue_list[j]=0           
                    
                    for j in range(0,len(xvalue_list)):#remove 0
                        if xvalue_list[j]==0:
                            xvalue_list[j]='b'
                            yvalue_list[j]='b'
                    
                        if yvalue_list[j]==0:
                            xvalue_list[j]='b'
                            yvalue_list[j]='b'
            
                    shan(xvalue_list,'b')    
                    shan(yvalue_list,'b')
                    
                   
                    if len(xvalue_list)>40:

                        slope, intercept, r_value, p_value, std_err = st.linregress(xvalue_list, yvalue_list)

                        
                        ENSG_col.append(gene)
                        gene_col.append(gene_sym)
                        TCGA_col.append(TCGA_ID)
                        Co_col.append(r_value)
                        p_col.append(p_value)
                        len_col.append(len(xvalue_list))
                        
                    else:
                        ENSG_col.append(gene)
                        gene_col.append(gene_sym)
                        TCGA_col.append(TCGA_ID)
                        Co_col.append('False')
                        p_col.append('False')
                        len_col.append(len(xvalue_list))                           
                            
                except:
                    ENSG_col.append(gene)
                    gene_col.append(gene_sym)
                    TCGA_col.append(TCGA_ID)
                    Co_col.append('False')
                    p_col.append('False')
                    len_col.append('False')
                    
    d={'gene':gene_col,'ID':ENSG_col,'TCGA':TCGA_col,'r_value':Co_col,'p_value':p_col,'sample_size':len_col}
    
    df=pd.DataFrame(d)
    
    print(df)
    
    if counter == 1:
        
        df.to_csv('out/AiAcancer_' + target_gene + '.csv',mode='a', header=True)
        
    else:
        df.to_csv('out/AiAcancer_' + target_gene + '.csv',mode='a', header=False)

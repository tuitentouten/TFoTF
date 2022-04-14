# -*- coding: utf-8 -*-
"""
Created on Sat Jan 15 15:34:15 2022

@author: wangf
"""

import pandas as pd
import os

os.chdir('F:/OneDrive/1_工作和学习/manuscripts/TF预测/raw data/fig5/scripts/sample/sources')

def caculateMAX(score_li):
    
    score_li.sort(reverse=True) 
    
    total = sum(score_li)
    
    max1 = max(score_li)
    
    if len(score_li) >= 3:
        ssscore_li = score_li[:3]
        max3 = sum(ssscore_li)
    if len(score_li) < 3:
        max3 = -1
        
    return total,max1,max3
   

     

TF = 'STAT1'

a = pd.read_csv('out/'+TF+'_pro5000.csv')

genelist=pd.read_csv('promoter_forTCGA_RNAseq.csv')#读取promoter整理文档



#sa = a[a['gene'] == 'ENSG00000234958']





genelist=genelist.set_index('gene_ID')#用gene名作为index
 
genelist = genelist.dropna()

genes = genelist.index.to_list()#获取gene名为列表

gene_li = []
names_li = []
total_li = []
max1_li = []
max3_li = []


con = 0
alln = len(genes)
for i in genes:
    
    con+=1
    print(i)
    print('%.3f%%' %((con/alln)*100))
    
    try:
        
        sle = a[a['gene'] == i]
        
        n = sle.name.to_list()

        name = n[0]#
        
        score_li = sle.value.to_list()
        
        total,max1,max3 = caculateMAX(score_li)
        
        
        
        gene_li.append(i)
        
        names_li.append(name)
        
        total_li.append(total)
        
        max1_li.append(max1)
        
        max3_li.append(max3)
        
    except:
        
        gene_li.append(i)
        
        names_li.append('somethw')
        
        total_li.append(-1)
        
        max1_li.append(-1)
        
        max3_li.append(-1)
        
d = {'gene':gene_li,'name':names_li,'total':total_li,'max1':max1_li,'max3':max3_li}
df = pd.DataFrame(d)

df.to_csv('out/'+TF+'_pro5000_caculated.csv')


    
        
        


         

         

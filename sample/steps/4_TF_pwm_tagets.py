# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 18:23:01 2022

@author: wangf
"""

import pandas as pd
import time
import os

os.chdir('F:/OneDrive/1_工作和学习/manuscripts/TF预测/raw data/fig5/scripts/sample/sources')

timeID = time.strftime('%Y%m%d_%H%M%S')

target = 'STAT1'

pwm=pd.read_csv('pwm_'+target+'.csv')#读取pwm打分矩阵

genelist=pd.read_csv('out/promoter210703.csv')#读取promoter整理文档


genelist=genelist.set_index('gene_ID')#用gene名作为index
 
genelist = genelist.dropna()

gene=genelist.index.to_list()#获取gene名为列表

motif_len=pwm.shape[1] - 1

output_file = 'out/'+target+ '_pro5000.csv'
output_file_error = 'out/rst_'+target+'_'+ timeID +'_error.csv'

a={'gene':[],'name':[],'start':[],'sequence':[],'value':[],'chain':[],'promoter':[]}
df=pd.DataFrame(a)
df.to_csv(output_file,mode='a', header=True)
df.to_csv(output_file_error,mode='a', header=True)

#gene = ['ENSG00000264462']
n = 0
for j in gene:
    n+=1
    alln = len(gene)
    print('%.3f%%' %((n/alln)*100))
    
    try:
        
    
        promoter=genelist.loc[j,'promoter']#获取gene j的promoter区序列
        
        if promoter == 'error':
            print('have error')
            pass   
        
        else:
            
            genename = genelist.loc[j,'gene_name']
            num=len(promoter)-motif_len+1
            #print(num)
            k=0
            #print('0')
            
            while k<num:
                #print('1')
                k+=1
                test=promoter[k:(k+motif_len)]
                #print(k)
             
                    
                if 'N' in test:
                    print('have N')
                    pass
    
                
                else:
                    #print('2')
                    p=1
                    value=0
                    for i in test:
                        
                        y=p
                        p+=1
                        if i=='A':
                            x=pwm.iloc[0,y]
                        if i=='C':
                            x=pwm.iloc[1,y]
                        if i=='G':
                            x=pwm.iloc[2,y]
                        if i=='T':
                            x=pwm.iloc[3,y]
        
                        #print(x)
                        value+=x
                        #print(value)
                        #print('3')
                    if value > 0:   
                        #print('4')
                        #print(target,'  ',genename,'  ',j,'  ',k,'  ',test,'  ',value,'  +')
                            
                        a={'gene':j,'name':genename,'start':k,'sequence':test,'value':value,'chain':'+','promoter':promoter}
                        df=pd.DataFrame(a,index=[0])
                        
                        df.to_csv(output_file,mode='a', header=False)
            print(j+'completed')
            
        
            
    except:
        #print('5')
        a={'gene':j,'name':'genename','start':-1,'sequence':'ATCG','value':-1,'chain':'+','promoter':promoter}
        df=pd.DataFrame(a,index=[0])
        
        df.to_csv(output_file_error,mode='a', header=False)
        print(j+'failed')
        




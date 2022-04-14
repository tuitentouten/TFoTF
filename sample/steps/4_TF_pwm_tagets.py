# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 18:23:01 2022

@author: wangf
"""

import pandas as pd
import time
import os

os.chdir('F:/sample/sources')### set workdir, for example, 'F:/sample/sources'

timeID = time.strftime('%Y%m%d_%H%M%S')

target = 'STAT1'

pwm=pd.read_csv('pwm_'+target+'.csv')

genelist=pd.read_csv('out/promoter210703.csv')


genelist=genelist.set_index('gene_ID')
 
genelist = genelist.dropna()

gene=genelist.index.to_list()

motif_len=pwm.shape[1] - 1

output_file = 'out/'+target+ '_pro5000.csv'
output_file_error = 'out/rst_'+target+'_'+ timeID +'_error.csv'

a={'gene':[],'name':[],'start':[],'sequence':[],'value':[],'chain':[],'promoter':[]}
df=pd.DataFrame(a)
df.to_csv(output_file,mode='a', header=True)
df.to_csv(output_file_error,mode='a', header=True)


n = 0
for j in gene:
    n+=1
    alln = len(gene)
    print('%.3f%%' %((n/alln)*100))
    
    try:
        
    
        promoter=genelist.loc[j,'promoter']
        
        if promoter == 'error':
            print('have error')
            pass   
        
        else:
            
            genename = genelist.loc[j,'gene_name']
            num=len(promoter)-motif_len+1

            k=0

            while k<num:

                k+=1
                test=promoter[k:(k+motif_len)]

             
                    
                if 'N' in test:
                    print('have N')
                    pass
    
                
                else:

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
        

                        value+=x

                    if value > 0:   

                            
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
        




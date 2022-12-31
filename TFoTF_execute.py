# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 02:25:26 2021

@author: wangf
"""



# step 1: Set workdir, for example, 'F:/sample/sources'.
set_dir = 'O:/TFoTF_main/sample/sources'

# step 2: Set target TF. 
# You must have the corresponding pwm data file, like './sample/sources/pwm_STAT1.csv'. Use uppercase letters.
target_gene = 'CREB1' 

# step 3: Set cut-offs to determine the significance.
pvalue_cut = 0.05 # from 0 to 1

# step4: Run this script. 





















#DO NOT modify the following codes if you do not know their meaning.
import pandas as pd
import scipy.stats as st
import math
import os
import numpy as np
import time
os.chdir(set_dir) 
#------------------------
try:
    for i in os.listdir('temporary files/'):
        print(i)
        txt_path=os.path.join('temporary files/',i)
        os.remove(txt_path)
except:
    pass
#------------------------
peasonr_cut = 0.4 
#load other cancer datasets one by one. These data were NOT provided in this sample
notation = pd.read_csv('ENSGs_information.csv')
notation_ID_first=notation.set_index(['Gene stable ID'])  
notation_name_first=notation.set_index(['Gene name'])

print('step(1/4) loading files...')
ACC  = pd.read_csv('expression/ACC_expression.csv')
BLCA = pd.read_csv('expression/BLCA_expression.csv')
BRCA = pd.read_csv('expression/BRCA_expression.csv')
CESC = pd.read_csv('expression/CESC_expression.csv')
CHOL = pd.read_csv('expression/CHOL_expression.csv')
COAD = pd.read_csv('expression/COAD_expression.csv')
DLBC = pd.read_csv('expression/DLBC_expression.csv')
ESCA = pd.read_csv('expression/ESCA_expression.csv')
GBM  = pd.read_csv('expression/GBM_expression.csv')
HNSC = pd.read_csv('expression/HNSC_expression.csv')
KICH = pd.read_csv('expression/KICH_expression.csv')
KIRC = pd.read_csv('expression/KIRC_expression.csv')
KIRP = pd.read_csv('expression/KIRP_expression.csv')
LAML = pd.read_csv('expression/LAML_expression.csv')
LGG  = pd.read_csv('expression/LGG_expression.csv')
LIHC = pd.read_csv('expression/LIHC_expression.csv')
LUAD = pd.read_csv('expression/LUAD_expression.csv')
LUSC = pd.read_csv('expression/LUSC_expression.csv')
MESO = pd.read_csv('expression/MESO_expression.csv')
OV   = pd.read_csv('expression/OV_expression.csv')
PAAD = pd.read_csv('expression/PAAD_expression.csv')
PCPG = pd.read_csv('expression/PCPG_expression.csv')
PRAD = pd.read_csv('expression/PRAD_expression.csv')
READ = pd.read_csv('expression/READ_expression.csv')
SARC = pd.read_csv('expression/SARC_expression.csv')
SKCM = pd.read_csv('expression/SKCM_expression.csv')
STAD = pd.read_csv('expression/STAD_expression.csv')
TGCT = pd.read_csv('expression/TGCT_expression.csv')
THCA = pd.read_csv('expression/THCA_expression.csv')
THYM = pd.read_csv('expression/THYM_expression.csv')
UCEC = pd.read_csv('expression/UCEC_expression.csv')
UVM  = pd.read_csv('expression/UVM_expression.csv')
USC  = pd.read_csv('expression/USC_expression.csv')
#...    
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
file_li = []
for file in os.listdir('expression'):
    file_li.append(file)

AC = CHOL.set_index('gene')
ENSG_list = AC.index.to_list()
counter = 0

for gene in ENSG_list: 
    counter+=1
    print('step(2/4) Calculation in progress'+'  '+'%.2f%%' % ((counter/len(ENSG_list) * 100)))
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
        
        EX2 = GeneEX.set_index(['gene'])#ENSG id as index

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
                yvalue_list=logize(yvalue_list)
                xvalue_list=logize(xvalue_list)
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
    if counter == 1:
        df.to_csv('temporary files/AiAcancer_' + target_gene + '.csv',mode='a', header=True)
    else:
        df.to_csv('temporary files/AiAcancer_' + target_gene + '.csv',mode='a', header=False)

#------------------------
print(' ')
df = pd.read_csv('temporary files/AiAcancer_' + target_gene + '.csv')

source_df = pd.read_csv('expression/CHOL_expression.csv')#读取基因名列表
source_df = source_df.set_index('gene')

source = source_df.index.to_list()
ENSG_list = df.ID.to_list()

gene_list = []

cc = 0
print('data organization in process')
for ENSG in source:
    cc += 1
    if ENSG in ENSG_list:
        gene_list.append(ENSG)
df = df.drop(df[df.p_value == 'False'].index)

df = df.apply(pd.to_numeric,errors='ignore')

df = df.drop(df[df.p_value > pvalue_cut].index)
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
    print('step(3/4) data organization in process'+'  '+'%.2f%%' % ((counter/len(gene_list) * 100)))

    try:
        df_s = df[df['ID'] == gene]
        
        df_s = df_s.loc[df_s['r_value'] != 'False']
        
        df_p = df_s[df_s['r_value'] > 0]
        df_n = df_s[df_s['r_value'] < 0]
        
        P_SIGNIS_count = len(df_s[df_s['r_value'] >= peasonr_cut])#
        N_SIGNIS_count = len(df_s[df_s['r_value'] <= -peasonr_cut])#
        
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
            if absr >= peasonr_cut:
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

dataF.to_csv('out/rvalue_file_' + target_gene + '.csv')


#----------------------------
timeID = time.strftime('%Y%m%d_%H%M%S')
pwm=pd.read_csv('pwm_'+target_gene+'.csv')
genelist=pd.read_csv('promoter_forTCGA_RNAseq.csv')
genelist=genelist.set_index('gene_ID')
genelist = genelist.dropna()
gene=genelist.index.to_list()
motif_len=pwm.shape[1] - 1

output_file = 'temporary files/'+target_gene+ '_pro5000.csv'
output_file_error = 'temporary files/rst_'+target_gene+'_'+ timeID +'_error.csv'

a={'gene':[],'name':[],'start':[],'sequence':[],'value':[],'chain':[],'promoter':[]}
df=pd.DataFrame(a)
df.to_csv(output_file,mode='a', header=True)
df.to_csv(output_file_error,mode='a', header=True)
n = 0
for j in gene:
    print('\npwm calculation in process')
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
            print(j+' completed')
    except:
        #print('5')
        a={'gene':j,'name':'genename','start':-1,'sequence':'ATCG','value':-1,'chain':'+','promoter':promoter}
        df=pd.DataFrame(a,index=[0])
        
        df.to_csv(output_file_error,mode='a', header=False)
        print(j+' failed')
        
#----------------------------       
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
   
a = pd.read_csv('temporary files/'+target_gene+'_pro5000.csv')
genelist=pd.read_csv('promoter_forTCGA_RNAseq.csv')
genelist=genelist.set_index('gene_ID')
genelist = genelist.dropna()
genes = genelist.index.to_list()
gene_li = []
names_li = []
total_li = []
max1_li = []
max3_li = []
con = 0
alln = len(genes)
for i in genes:
    con+=1
    print('step(4/4)' + i)
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
df.to_csv('out/pwm_file_'+target_gene+'.csv')

print('TFoTF execute: completed')
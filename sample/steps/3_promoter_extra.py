# -*- coding: utf-8 -*-
"""
Created on Fri Jul  2 20:04:26 2021

@author: wangf
"""


from Bio import SeqIO
import  pandas as pd
import os

os.chdir('F:/OneDrive/1_工作和学习/manuscripts/TF预测/raw data/fig5/scripts/sample/sources')

data0 = pd.read_csv('expression/CHOL_expression.csv')
ENSG_li = data0['sample'].to_list()

data1 = pd.read_csv('ENSGs_translation.csv').set_index('ID')

data2 = pd.read_csv('translation210702.csv')
data2 = data2.set_index('Gene_ID')

data3 = data2.set_index('Transcript_ID')

upstream = 5000


def get_sequence(chromosome,start,end):
    
    in_file='chrs_hg38/{c}.fa'.format(c = 'chr'+ chromosome) #genome file
    handle_in=open(in_file,"rU")
    Seq_ref=SeqIO.parse(handle_in,"fasta")
    scaffold_name='chr'+ chromosome # scaffold/chromosome name in the  fasta header, without ">" 
    pos_start=int(start+1) # start of the target sequence
    pos_end=int(end) # end of the target sequence, start and end coordinates can be identical to get a single nucleotide position
    
    Extracted_seq=[]                                # create blank list to contain the 1.5kb sequences to be extracted
              # parse the genome fasta file
    seq_sp=[]
    for record in Seq_ref:				# iterate over each scaffold of the reference genome
    	if record.id==scaffold_name:		# if the reference record has the same name as the target scaffold/chromosome
    		#print (">"+scaffold_name+":"+str(pos_start)+":"+str(pos_end))	# append the scaffold name into the final sequence file
    		seq_sp=list(str(record.seq))		# turn the sequence of the scaffold into a list containing single-letter elements
    		Extracted_seq.append(''.join(seq_sp[pos_start-1:pos_end]))
    for i in Extracted_seq:
    		print (start)
    handle_in.close()
    return(i)

counter = 0
#----------------------------
for i in ENSG_li[:10]:
    #print(i)
    
    ENSG = []
    name = []
    chro = []
    TSSR = []
    stra = []
    upst = []
    prom = []
    
    
    try:
        i_name = data1.at[i,'name']
        #i = 'ENSG00000204713'
        #获取主要转录本的ID------------------------------------------------
        length = data2.at[i,'length'].tolist()#某ENSG的所有转录本长度
        
        typeC = data2.at[i,'type']#某ENSG的所有转录本类型
        
        transC = data2.at[i,'Transcript_ID']#某ENSG的所有转录本ID
        
        if type(typeC) == str:#如果这个基因只有一个转录本，则直接组成list
            
            gtype = []
            gtrans = []
            
            gtype.append(typeC)
            gtrans.append(transC)
            
        else:#如果这个基因有多个转录本，则tolist
                
            gtype = typeC.tolist()#某ENSG所有转录本类型
            
            gtrans = transC.tolist()#某ENSG所以转录本ID
        
        tempd = pd.DataFrame({'l':length,'gt':gtype,'gtr':gtrans})#上述三list编成表格
        
        tempd = tempd.sort_values("l",ascending=False).reset_index(drop=True)#上表内容按照转录本长度从大到小排序
        
        #tempd_gt = tempd.set_index('gt')#上表用转录本类型作为索引
        
        if 'protein_coding' in gtype:#如果该基因存在蛋白编码转录本，则选择最长的蛋白编码转录本作为主要转录本
                      
            gtrp = tempd[tempd['gt'] == 'protein_coding'].reset_index()#_gt.at['protein_coding','gtr']#选出已按照长度从大到小排列的所有蛋白编码转录本的转录本ID
            
            major_trans = gtrp.at[0,'gtr']
            
        else:
            
            major_trans = tempd.at[0,'gtr']#如果该基因不存在蛋白编码转录本，则直接选择长度排名第一的转录本作为主要转录本
            
        TSS = int(data3.at[major_trans,'TSS'])
        chromosome = str(data3.at[major_trans,'chr'])
        strand = str(data3.at[major_trans,'Strand'])
        
        
        #-----------------------获得启动子区（TSS上游2000bp）
        try:
        
            if strand == '-1':
                
                x = get_sequence(chromosome,TSS,TSS+upstream)
                
                x = x.upper()
                
                x = list(x)
                
                x.reverse()
                
                xr = []
                
                for base in x:
                    
                    if base == 'A':
                        xr.append('T')
                    if base == 'T':
                        xr.append('A')
                    if base == 'C':
                        xr.append('G')
                    if base == 'G':
                        xr.append('C')
                    if base == 'N':
                        xr.append('N')
                
                x = ''.join(xr)
                
            else:
                
                x = get_sequence(chromosome,TSS-upstream-1,TSS-1).upper()
            
            check = x[1]
            
        except:
            
            x = 'error'
            
        
        ENSG.append(i)
        name.append(i_name)
        chro.append(chromosome)
        TSSR.append(TSS)
        stra.append(strand)
        upst.append(upstream)
        prom.append(x)
        
    except:
        
        x = 'error'
        ENSG.append(i)
        name.append('error')
        chro.append('error')
        TSSR.append('error')
        stra.append('error')
        upst.append('error')
        prom.append('error')
    
    a = {'gene_ID':ENSG,'gene_name':name,'chr':chro,'TSS':TSSR,'strand':stra,'upstream':upst,'promoter':prom}
    
    df = pd.DataFrame(a)
    
    if i == ENSG_li[0]:
        df.to_csv('out\promoter210703.csv',mode='a', header=True)
    else:
        df.to_csv('out\promoter210703.csv',mode='a', header=False)
        
    print('   ' + str(counter))
    counter+=1
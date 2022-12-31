# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 17:09:47 2022

@author: wangf
"""




# set workdir, for example, 'F:/sample/sources/out'
set_dir = 'O:/TFoTF_main/sample/sources/out'

# Set the pearson correlation R score cut-off, Recommend > 0.6
rvalue_cutoff = 0.6

# Set the type of pwm score, you can choose between 'max1', 'max3', 'total'
pwm_type = 'max1'

# Set the pwm score cut-off, Recommended > 8 for max1
pwm_cutoff = 8

# The previously calculated processed data
rvalue_file = 'rvalue_file_CREB1.csv'
pwm_file = 'pwm_file_CREB1.csv'

# name the result file, and the file will be created in your workdir
savefile_to = 'predict_result_CREB1.csv'















#DO NOT modify the following codes if you do not know their meaning.
import pandas as pd
import os
os.chdir(set_dir)
p = pd.read_csv(pwm_file)
r = pd.read_csv(rvalue_file)

ps = p[p[pwm_type] > pwm_cutoff]
rs = r[r['Max_R'] > rvalue_cutoff]

lr = rs.ID.to_list()
lp = ps.gene.to_list()

tmp = [val for val in lr if val in lp]

pps = ps[ps['gene'].isin(tmp)]
rrs = rs[rs['ID'].isin(tmp)]

rrrs = rrs[['Name','ID','Max_R']].set_index('Name',drop=True)
ppps = pps[['name','total','max1','max3']].set_index('name',drop=True)

ss = pd.concat([rrrs,ppps],axis=1)

ss.to_csv(savefile_to)
pa = os.path.abspath (savefile_to)
print('The prediction results, i.e. the list of target genes, are saved in the path: \'{}\''.format(pa))





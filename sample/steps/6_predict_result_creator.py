# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 17:09:47 2022

@author: wangf
"""

import pandas as pd
import os

os.chdir('F:/OneDrive/1_工作和学习/manuscripts/TF预测/raw data/fig5/scripts/sample/sources')

p = pd.read_csv('out/STAT1_pro5000_caculated.csv')

r = pd.read_csv('out/AiAcancer_statat_STAT1.csv')

ps = p[p['max1'] > 8]

rs = r[r['Max_R'] > 0.7]

li = rs.ID.to_list()

s = ps[ps['gene'].isin(li)]

s.to_csv('out/predict_result.csv')
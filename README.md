# TFoTF
1_ Download and uncompress TFoTF-main. Copy 'TFoTF-main/sample' to a certain directory, for example, 'F:/sample'.

2_ Download and uncompress 'sources' from 

3_ Execute the scripts in the 'sample/steps/' directory in numerical order by file name. Before running each script, make sure you have set the correct workdir.


3_1 1_AiA_creation.py: Calculate the expression correlation of a given TF with other genes, based on the gene expression data.

3_2 2_AiA_statistic_procession.py: Collation and statistics of the pan-cancer expression correlation data obtained from 1_ calculation.

3_3 3_promoter_extra.py: Get the promoter sequence from Genome Reference hg38. 

3_4 4_TF_pwm_tagets.py: Calculate the PWM score of given TF with other genes promoter sequence.

3_5 5_TF_pwm_tagets_caculator.py: Collation and statistics of the pan-cancer PWM score obtained from 4_ calculation.

3_6 6_predict_result_creator.py: Set personalized cut-offs and obtain predict results.

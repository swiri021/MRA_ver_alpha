from __init__ import dl,cal

import matplotlib as mlab
mlab.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns; sns.set_style("whitegrid")

expr_df, clinical_df = dl.gdx_loading(set_index_expr='selectedEntrezIds', set_index_clinical='celfile_name')
clinical_df = clinical_df.loc[clinical_df.index.notnull()]

ext_id, ext_expr = dl.data_extract(expr=expr_df, clinical=clinical_df, extract_columns='race')
bi_dict = {1:'Black', 3:'White'}
ext_id = map(bi_dict.get,ext_id)

df1 = ext_expr[ext_id.index('Black')]
df2 = ext_expr[ext_id.index('White')]

# 0:TFID, 1: df1 result, 2: df2 result
# 0:es, 1:nes, 2:pval
tfs = cal.GSEA_exec(df1, df2, tf_set='data/PDI_181211.csv', perm_type='phenotype', weighted_score_type=0, nperm=100, plot=True)
print tfs

"""
df1_filtered = []
df2_filtered = []
rw1 = open('df1_result.csv','w')
rw2 = open('df2_result.csv','w')

for tf_enriched in tfs:
	if tf_enriched[1][2]<0.05:
		#df1_filtered.append([tf_enriched[0], tf_enriched[1]])
		rw1.write(str(tf_enriched[0])+','+str(tf_enriched[1][0])+','+str(tf_enriched[1][1])+','+str(tf_enriched[1][2])+'\n')

	if tf_enriched[2][2]<0.05:
		#df2_filtered.append([tf_enriched[0], tf_enriched[2]])
		rw2.write(str(tf_enriched[0])+','+str(tf_enriched[2][0])+','+str(tf_enriched[2][1])+','+str(tf_enriched[2][2])+'\n')

rw1.close()
rw2.close()
"""

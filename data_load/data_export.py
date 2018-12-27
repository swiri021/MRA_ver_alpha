from multiprocessing import Process, Queue
import numpy as np

def gene_set_zscore( arr1, gene_set=[] ,sample_status="multiple", form='pandas'):

	def cal_z(temp_df, inter):
		selected_adj = temp_df.loc[inter].values.astype(float)
		total_adj = temp_df.values.astype(float)

		diff_mean = np.mean(selected_adj) - np.mean(total_adj)
		result = diff_mean*np.sqrt(len(selected_adj))/np.std(total_adj,ddof=1)
		return result

	ft_dat_index = gene_set
	arr1_index = arr1.index.tolist()
	inter = list(set(arr1_index).intersection(ft_dat_index))

	if sample_status=="single":
		zscore = [cal_z(arr1, inter)]

	elif sample_status=="multiple":
		diff_mean = arr1.loc[inter].mean(axis=0).subtract(arr1.mean(axis=0))
		len_norm = arr1.std(ddof=1, axis=0).apply(lambda x: np.sqrt(len(inter))/x)
		zscore = diff_mean*len_norm
		zscore = zscore.to_frame()
		zscore.columns = ['Zscore']

		if form=='pandas':
			#zscore = pd.DataFrame(data=zscore, columns=['Zscore'], index=arr1.columns.tolist())
			return zscore
		elif form=='array':
			return zscore['Zscore'].values.tolist()
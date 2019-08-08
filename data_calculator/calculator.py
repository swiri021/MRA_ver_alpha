from __future__ import division
from scipy import stats
import pandas as pd
import numpy as np
import random

import matplotlib as mlab
mlab.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
from multiprocessing import Process, Queue




class GSEA_calculator:

	def plot(self, rank_d, NES, pval, RES, hit_ind, filen = 'test1.png'):
		plt.clf()

		gs = plt.GridSpec(12,1)
		fig = plt.figure(figsize=(5,6))

		#RES = [float("{:.3f}".format(x)) for x in RES.tolist()]
		#print RES

		rank_d['ES'] = RES.tolist()
		rank_d['hit'] = 0
		rank_d['hit'].iloc[hit_ind]=1

		ax1 = fig.add_subplot(gs[0:8])
		ax1 = rank_d['ES'].plot(kind='line', ylim=(-1,1))
		ax1.axhline(linewidth=.8, color='k')
		ax1.text(.8, .8, "NES : {:.3f}".format(NES), transform=ax1.transAxes)
		ax1.text(.8, .9, "P-value : {:.3f}".format(pval), transform=ax1.transAxes)

		ax2 = fig.add_subplot(gs[8], sharex=ax1)
		trans2 = transforms.blended_transform_factory(ax2.transData, ax2.transAxes)
		ax2.vlines(hit_ind, 0, 1,linewidth=.5,transform=trans2)
		ax2.spines['bottom'].set_visible(False)
		ax2.get_yaxis().set_visible(False)
		#ax2.tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off', right='off')
		ax2.grid(False)

		ax3 = fig.add_subplot(gs[9], sharex=ax1)
		im_matrix = np.tile(rank_d['S2N'].values.tolist(), (2,1))
		ax3.imshow(im_matrix, aspect='auto', cmap=plt.cm.seismic, interpolation='none') # cm.coolwarm
		ax3.spines['bottom'].set_visible(False)
		ax3.get_yaxis().set_visible(False)
		#ax3.tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off', right='off')
		ax3.grid(False)

		ax4 = fig.add_subplot(gs[10:13], sharex=ax1)
		ax4.fill_between(range(len(rank_d['S2N'].values.tolist())), y1=rank_d['S2N'].values.tolist(), y2=0, color='#C9D3DB')
		ax4.set_xlabel("Rank in Ordered Dataset", fontsize=14)
		ax4.get_yaxis().set_visible(False)
		#ax4.tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off', right='off')
		ax4.grid(False)

		fig.savefig(filen)


	def GSEA_exec(self, df1, df2, tf_set, weighted_score_type=1, nperm=1000, perm_type='gene_set', rs=np.random.RandomState(), single=False, scale=False, plot=False, plot_path='/home/ubuntu/codes/MRA_GSEA_proj/example_result/gsea_plot/'):
		assert perm_type=='gene_set' or perm_type=='phenotype', "Please, clarify permutation type(gene_set or phenoteyp)"

		# 0 : TF ID, 1 : df1 result, 2: df2 result
		tf_result = TF_process().run(df1, df2, data_location=tf_set)
		rank_d1 = self.ranking_metric(df1, df2, nperm=nperm, perm_type=perm_type, combined=True)
		rank_d2 = self.ranking_metric(df2, df1, nperm=nperm, perm_type=perm_type, combined=True)

		gsea_result = []

		for tf in tf_result:

			if perm_type=='gene_set':
				es1, esnull1, hit_ind1, RES1, leading_edge1, portion_le1 = self.enrichment_score(gene_list=rank_d1.index.tolist(), correl_vector=rank_d1.values.tolist(), gene_set=tf[1], weighted_score_type=weighted_score_type, nperm=nperm, rs=rs, single=single, scale=scale)
				pval1 = self.gsea_pval([es1], [esnull1])
				nes1 = self.normalize(es1,esnull1)
				#es2, esnull2, hit_ind2, RES2, leading_edge2, portion_le2 = self.enrichment_score(gene_list=rank_d2.index.tolist(), correl_vector=rank_d2.values.tolist(), gene_set=tf[2], weighted_score_type=weighted_score_type, nperm=nperm, rs=rs, single=single, scale=scale)
				#pval2 = self.gsea_pval([es2], [esnull2])
				#nes2 = self.normalize(es2,esnull2)

			elif perm_type=='phenotype':
				es1, esnull1, hit_ind1, RES1, leading_edge1, portion_le1 = self.enrichment_score_phenotype(ranking_metric=rank_d1, gene_set=tf[1], weighted_score_type=weighted_score_type, nperm=nperm, rs=rs, single=single, scale=scale)
				pval1 = self.gsea_pval([es1], [esnull1])
				nes1 = self.normalize(es1,esnull1)
				#es2, esnull2, hit_ind2, RES2,leading_edge2, portion_le2 = self.enrichment_score_phenotype(ranking_metric=rank_d2, gene_set=tf[2], weighted_score_type=weighted_score_type, nperm=nperm, rs=rs, single=single, scale=scale)
				#pval2 = self.gsea_pval([es2], [esnull2])
				#nes2 = self.normalize(es2,esnull2)

			if plot==True:
				#print rank_d1
				if perm_type=='phenotype':
					rank_d1_p = rank_d1[-1]
					#rank_d2_p = rank_d2[-1]

				self.plot(rank_d1_p, nes1, pval1, RES1, hit_ind1, filen=plot_path+tf[0]+'_1side.png')
				#self.plot(rank_d2_p, nes2, pval2, RES2, hit_ind2, filen=plot_path+tf[0]+'_2side.png')

			#gsea_result.append([tf[0], [es1, nes1, pval1], [es2,nes2,pval2]])
			gsea_result.append([tf[0], es1, nes1, pval1])

		return gsea_result


	def ranking_metric(self, df1, df2, nperm, perm_type='gene_set',combined=True):
		def s2n(df1, df2):
			df1_mean = df1.mean(axis=1)
			df2_mean = df2.mean(axis=1)

			df1_std = df1.std(axis=1)
			df2_std = df2.std(axis=1)

			sub_mean = df1_mean-df2_mean
			sum_std = df1_std+df2_std

			result = sub_mean/sum_std

			return result

		if perm_type=='gene_set':
			result = s2n(df1, df2)
			if combined==True:
				result1 = result
				result1.name='S2N'
				result1 = result1.to_frame()
				result1.index = result1.index.map(lambda x:x+'+')

				result2 = result
				result2.name='S2N'
				result2 = result2.to_frame()
				result2 = result2.apply(lambda x:x*-1)
				result2.index = result2.index.map(lambda x:x+'-')

				combined_result = result1.append(result2)
				combined_result = combined_result.sort_values(by=['S2N'],ascending=False)

			else:
				result1 = result.sort_values(ascending=False)
				result1.name='S2N'
				result1 = result1.to_frame()
				combined_result = result1

		elif perm_type=='phenotype':
			df1_len, df2_len = len(df1.columns.tolist()), len(df2.columns.tolist())
			tdf = pd.concat([df1, df2], axis=1)

			pdf = np.random.choice(len(tdf.columns.tolist()), (nperm,len(tdf.columns.tolist())))
			pdf = [[list(x)[:df1_len], list(x)[df1_len:]]for x in pdf]
			pdf = pdf+[[range(len(tdf.columns.tolist()))[:df1_len], range(len(tdf.columns.tolist()))[df1_len:]]]
			s2n_arr = [s2n(tdf[tdf.columns[x[0]]], tdf[tdf.columns[x[1]]]) for x in pdf]

			permutated_combined = []
			for x in s2n_arr:
				result = x
				if combined==True:
					result1 = result
					result1.name='S2N'
					result1 = result1.to_frame()
					result1.index = result1.index.map(lambda x:x+'+')

					result2 = result
					result2.name='S2N'
					result2 = result2.to_frame()
					result2 = result2.apply(lambda x:x*-1)
					result2.index = result2.index.map(lambda x:x+'-')

					combined_result = result1.append(result2)
					combined_result = combined_result.sort_values(by=['S2N'],ascending=False)
					permutated_combined.append(combined_result)

				else:
					result1 = result.sort_values(ascending=False)
					result1.name='S2N'
					result1 = result1.to_frame()
					combined_result = result1
					permutated_combined.append(combined_result)

			combined_result = permutated_combined

		return combined_result


	####### FROM GSEAPY, It must be belong to GSEApy, and its copyright
	def enrichment_score_phenotype(self, ranking_metric, gene_set, weighted_score_type=1, nperm=1000, perm_type='gene_set', rs=np.random.RandomState(), single=False, scale=False):
		"""This is the most important function of GSEApy. It has the same algorithm with GSEA and ssGSEA.
		:param gene_list:       The ordered gene list gene_name_list, rank_metric.index.values
		:param gene_set:        gene_sets in gmt file, please used gsea_gmt_parser to get gene_set.
		:param weighted_score_type:  It's same with gsea's weighted_score method. weighting by the correlation
								is a very reasonable choice that allows significant gene sets with less than perfect coherence.
								options: 0(classic),1,1.5,2. default:1. if one is interested in penalizing sets for lack of
								coherence or to discover sets with any type of nonrandom distribution of tags, a value p < 1
								might be appropriate. On the other hand, if one uses sets with large number of genes and only
								a small subset of those is expected to be coherent, then one could consider using p > 1.
								Our recommendation is to use p = 1 and use other settings only if you are very experienced
								with the method and its behavior.
		:param correl_vector:   A vector with the correlations (e.g. signal to noise scores) corresponding to the genes in
								the gene list. Or rankings, rank_metric.values
		:param nperm:           Only used this parameter when computing esnull for statistical testing. set the esnull value
								equal to the permutation number.
		:param rs:              Random state for initialize gene list shuffling. Default: np.random.RandomState(seed=None)
		:return:
		 ES: Enrichment score (real number between -1 and +1)
		 ESNULL: Enrichment score calculated from random permutation.
		 Hits_Indices: index of a gene in gene_list, if gene included in gene_set.
		 RES: Numerical vector containing the running enrichment score for all locations in the gene list .
		"""
		correl_vector_array = []
		tag_indicator_array = []

		gene_list = ranking_metric[len(ranking_metric)-1].index.tolist()

		for x in ranking_metric:
			temp_gl=x.index.tolist()
			correl_vector=x.values.tolist()

			N = len(gene_list)
			# Test whether each element of a 1-D array is also present in a second array
			# It's more intuitived here than orginal enrichment_score source code.
			# use .astype to covert bool to intergers
			tag_indicator = np.in1d(temp_gl, gene_set, assume_unique=True).astype(int)  # notice that the sign is 0 (no tag) or 1 (tag)
			if weighted_score_type == 0 :
				correl_vector = np.repeat(1, N)
			else:
				correl_vector = np.abs(correl_vector)**weighted_score_type

			# get indices of tag_indicator

			# if used for compute esnull, set esnull equal to permutation number, e.g. 1000
			# else just compute enrichment scores
			# set axis to 1, because we have 2 dimentional array
			correl_vector_array.append(correl_vector)
			tag_indicator_array.append(tag_indicator)

		correl_vector = np.array(correl_vector_array)
		tag_indicator = np.array(tag_indicator_array)

		axis = 1
		hit_ind = np.flatnonzero(tag_indicator[-1]).tolist()

		Nhint = tag_indicator.sum(axis=axis, keepdims=True)
		sum_correl_tag = np.sum(correl_vector*tag_indicator, axis=axis, keepdims=True)
		# compute ES score, the code below is identical to gsea enrichment_score method.
		no_tag_indicator = 1 - tag_indicator
		Nmiss =  N - Nhint
		norm_tag =  1.0/sum_correl_tag
		norm_no_tag = 1.0/Nmiss

		RES = np.cumsum(tag_indicator * correl_vector * norm_tag - no_tag_indicator * norm_no_tag, axis=axis)

		if scale: RES = RES / N
		if single:
			es_vec = RES.sum(axis=axis)
		else:
			max_ES, min_ES =  RES.max(axis=axis), RES.min(axis=axis)
			es_vec = np.where(np.abs(max_ES) > np.abs(min_ES), max_ES, min_ES)
		# extract values
		es, esnull, RES = es_vec[-1], es_vec[:-1], RES[-1,:]

		if es<0:

			#print list(RES).index(min(RES))
			tag_range = tag_indicator[-1]
			tag_range = tag_range[:list(RES).index(min(RES))]
			indices = np.nonzero(tag_range)
			leading_edge = np.array(gene_list)[indices]

		else:

			#print list(RES).index(min(RES))
			tag_range = tag_indicator[-1]
			if max(RES)==np.nan:
				leading_edge = np.nan

			else:
				tag_range = tag_range[:list(RES).index(max(RES))]
				indices = np.nonzero(tag_range)
				leading_edge = np.array(gene_list)[indices]

		return es, esnull, hit_ind, RES, leading_edge, float(len(leading_edge))/float(len(gene_set))

	def enrichment_score(self, gene_list, correl_vector, gene_set, weighted_score_type=1, nperm=1000, rs=np.random.RandomState(), single=False, scale=False):
		"""This is the most important function of GSEApy. It has the same algorithm with GSEA and ssGSEA.
		:param gene_list:       The ordered gene list gene_name_list, rank_metric.index.values
		:param gene_set:        gene_sets in gmt file, please used gsea_gmt_parser to get gene_set.
		:param weighted_score_type:  It's same with gsea's weighted_score method. weighting by the correlation
								is a very reasonable choice that allows significant gene sets with less than perfect coherence.
								options: 0(classic),1,1.5,2. default:1. if one is interested in penalizing sets for lack of
								coherence or to discover sets with any type of nonrandom distribution of tags, a value p < 1
								might be appropriate. On the other hand, if one uses sets with large number of genes and only
								a small subset of those is expected to be coherent, then one could consider using p > 1.
								Our recommendation is to use p = 1 and use other settings only if you are very experienced
								with the method and its behavior.
		:param correl_vector:   A vector with the correlations (e.g. signal to noise scores) corresponding to the genes in
								the gene list. Or rankings, rank_metric.values
		:param nperm:           Only used this parameter when computing esnull for statistical testing. set the esnull value
								equal to the permutation number.
		:param rs:              Random state for initialize gene list shuffling. Default: np.random.RandomState(seed=None)
		:return:
		 ES: Enrichment score (real number between -1 and +1)
		 ESNULL: Enrichment score calculated from random permutation.
		 Hits_Indices: index of a gene in gene_list, if gene included in gene_set.
		 RES: Numerical vector containing the running enrichment score for all locations in the gene list .
		"""

		N = len(gene_list)
		# Test whether each element of a 1-D array is also present in a second array
		# It's more intuitived here than orginal enrichment_score source code.
		# use .astype to covert bool to intergers
		tag_indicator = np.in1d(gene_list, gene_set, assume_unique=True).astype(int)  # notice that the sign is 0 (no tag) or 1 (tag)
		if weighted_score_type == 0 :
			correl_vector = np.repeat(1, N)
		else:
			correl_vector = np.abs(correl_vector)**weighted_score_type

		# get indices of tag_indicator
		hit_ind = np.flatnonzero(tag_indicator).tolist()
		# if used for compute esnull, set esnull equal to permutation number, e.g. 1000
		# else just compute enrichment scores
		# set axis to 1, because we have 2 dimentional array
		axis = 1
		tag_indicator = np.tile(tag_indicator, (nperm+1,1))
		correl_vector = np.tile(correl_vector,(nperm+1,1))
		# gene list permutation
		for i in range(nperm): rs.shuffle(tag_indicator[i])
		# np.apply_along_axis(rs.shuffle, 1, tag_indicator)

		Nhint = tag_indicator.sum(axis=axis, keepdims=True)
		sum_correl_tag = np.sum(correl_vector*tag_indicator, axis=axis, keepdims=True)
		# compute ES score, the code below is identical to gsea enrichment_score method.
		no_tag_indicator = 1 - tag_indicator
		Nmiss =  N - Nhint
		norm_tag =  1.0/sum_correl_tag
		norm_no_tag = 1.0/Nmiss

		RES = np.cumsum(tag_indicator * correl_vector * norm_tag - no_tag_indicator * norm_no_tag, axis=axis)

		if scale: RES = RES / N
		if single:
			es_vec = RES.sum(axis=axis)
		else:
			max_ES, min_ES =  RES.max(axis=axis), RES.min(axis=axis)
			es_vec = np.where(np.abs(max_ES) > np.abs(min_ES), max_ES, min_ES)
		# extract values
		es, esnull, RES = es_vec[-1], es_vec[:-1], RES[-1,:]

		if es<0:

			#print list(RES).index(min(RES))
			tag_range = tag_indicator[-1]
			tag_range = tag_range[:list(RES).index(min(RES))]
			indices = np.nonzero(tag_range)
			leading_edge = np.array(gene_list)[indices]

		else:

			#print list(RES).index(min(RES))
			tag_range = tag_indicator[-1]
			tag_range = tag_range[:list(RES).index(max(RES))]
			indices = np.nonzero(tag_range)
			leading_edge = np.array(gene_list)[indices]

		return es, esnull, hit_ind, RES, leading_edge, float(len(leading_edge))/float(len(gene_set))


	def gsea_pval(self, es, esnull):
		"""Compute nominal p-value.
		From article (PNAS):
		estimate nominal p-value for S from esnull by using the positive
		or negative portion of the distribution corresponding to the sign
		of the observed ES(S).
		"""

		# to speed up, using numpy function to compute pval in parallel.
		es = np.array(es)
		esnull = np.array(esnull)

		# try:
		condlist = [ es < 0, es >=0]
		choicelist = [np.sum(esnull < es.reshape(len(es),1), axis=1)/ np.sum(esnull < 0, axis=1),
					  np.sum(esnull >= es.reshape(len(es),1), axis=1)/ np.sum(esnull >= 0, axis=1)]
		pval = np.select(condlist, choicelist)

		return pval[0]

	def normalize(self, es, esnull):
		"""normalize the ES(S,pi) and the observed ES(S), separately rescaling
			the positive and negative scores by dividing the mean of the ES(S,pi).
		"""

		try:
			if es == 0:
				return 0.0
			if es >= 0:
				meanPos = np.mean([a for a in esnull if a >= 0])
				#print es, meanPos
				return es/meanPos
			else:
				meanNeg = np.mean([a for a in esnull if a < 0])
				#print es, meanNeg
				return -es/meanNeg
		except:
			return 0.0 #return if according mean value is uncalculable
	####### FROM GSEAPY, It must be belong to GSEApy, and its copyright

class TF_process:

	def run(self, df1, df2, data_location):
		gene_set = pd.read_csv(data_location, index_col=0)
		gene_set.index = gene_set.index.astype(str)

		#gene_set = gene_set.loc['190']####For testing

		tf_list = list(set(gene_set.index.tolist()))
		tf_list = list(set(tf_list).intersection(set(df1.index.tolist())))
		tf_list = [x for x in tf_list if isinstance(gene_set.loc[x], pd.DataFrame)]

		tf_result = []

		for tf in tf_list:
			test_gene_set = gene_set.loc[tf]
			change_type = lambda x:[str(y) for y in x]


			test_gene_set = change_type(test_gene_set['Target ID'].values.tolist())
			test_gene_set = list(set(test_gene_set).intersection(set(df1.index.tolist())))
			if tf in test_gene_set:
				test_gene_set.pop(test_gene_set.index(tf))

			result1, result2 = self.TF_corr(df1, df2, test_gene_set, TF=tf)
			tf_result.append([tf, result1, result2])

		return tf_result


	def TF_corr(self, df1, df2, gene_set, TF):

		def corr_test(test, gs, tf):
			corr_target = []

			for x in gs:

				r,p = stats.spearmanr(test.loc[x].values.tolist(), test.loc[tf].values.tolist())
				if p<0.05 and r>0:
					corr_target.append(x+'+')
				elif p<0.05 and r<0:
					corr_target.append(x+'-')

			return corr_target

		test1 = df1.loc[gene_set+[TF]]
		test2 = df2.loc[gene_set+[TF]]

		test1 = test1.subtract(test1.median(axis=1), axis=0).dropna()
		test2 = test2.subtract(test2.median(axis=1), axis=0).dropna()
		test1_result = corr_test(test1, gene_set, TF)
		test2_result = corr_test(test2, gene_set, TF)
		return test1_result, test2_result


#class FET_calculator:

# Master Regulator Analysis by GSEA

Reference : gseapy (link : https://gseapy.readthedocs.io/en/latest/, copyright : Zhuoqing Fang Revision)
# Example :
```Python
from data_calculator.calculator import GSEA_calculator
import pandas as pd

##Group1 and Group2 have to be DataFrame format in pandas, and index has to be EntrezGene ID
group1 = pd.read_csv('Your_dataset1.csv', index_col=0)
group2 = pd.read_csv('Your_dataset2.csv', index_col=0)

##Index type = str
group1.index = group1.index.astype(int).astype(str)
group2.index = group2.index.astype(int).astype(str)

##Call class
cal = GSEA_calculator();

#tfs = cal.GSEA_exec(group1, group2, tf_set='data/PDI_181211.csv', perm_type='phenotype', weighted_score_type=0, nperm=100, plot=True)
tfs = cal.GSEA_exec(group1, group2, tf_set='data/PDI_181211.csv', perm_type='genotype', weighted_score_type=0, nperm=100, plot=True)

##MRA result
## Result order =
## EntrezID of MR candidate, Enrichment score(ES), Normalizaed ES, pvalue
print tfs
```

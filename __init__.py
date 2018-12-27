import pandas as pd
from data_load.GDX_data_loading import DATA_LOAD
from data_calculator.calculator import GSEA_calculator

dl = DATA_LOAD();
cal = GSEA_calculator();

#gene_set = pd.read_csv('test_input/test_InterF_alpha.csv', index_col=0)






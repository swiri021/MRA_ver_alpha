from data_load import DATA_PATH, EXPRESSION_DATA, CLINICAL_DATA, pd

class DATA_LOAD:

	def gdx_loading(self, set_index_expr='', set_index_clinical=''):
		if set_index_expr=='':
			expr_df = pd.read_csv(DATA_PATH+EXPRESSION_DATA, engine='c', low_memory=False, index_col=0)
			expr_df = self.gdx_data_preprocessing(expr_df)

		else:
			expr_df = pd.read_csv(DATA_PATH+EXPRESSION_DATA, engine='c', low_memory=False, index_col=0)
			expr_df = expr_df.set_index(set_index_expr)
			expr_df = self.gdx_data_preprocessing(expr_df)

		if set_index_clinical=='':
			clinical_df = pd.read_csv(DATA_PATH+CLINICAL_DATA, index_col=0)

		else:
			clinical_df = pd.read_csv(DATA_PATH+CLINICAL_DATA)
			clinical_df = clinical_df.set_index(set_index_clinical)
			clinical_df = clinical_df[clinical_df.index.notnull()]

		return expr_df, clinical_df

	def gdx_data_preprocessing(self, df):
		input = df
		input = input.loc[input.index.notnull()]
		input.index = input.index.astype(int)
		input.index = input.index.astype(str)
		input = input.groupby(input.index).max()
		return input

	def data_extract(self, expr, clinical, extract_columns=''):
		assert isinstance(expr, pd.DataFrame), "data_extract function needs DataFrame for input."
		assert len(expr)!=0 and len(clinical)!=0, "Cannot Extract GDX, because one of dataframe input is empty."
		assert extract_columns=='' or extract_columns in clinical.columns.tolist(), "Cannot Extract GDX, please clarify extract_point."
		assert len(set(expr.columns.tolist()).intersection(set(clinical.index.tolist())))!=0, "Nothing matched between expr and clinical."

		extract_template = clinical[extract_columns].unique().tolist()
		extracted_ids = [clinical[clinical[extract_columns]==i].index.tolist() for i in extract_template]
		extracted_expr = [expr[i] for i in extracted_ids]

		return extract_template, extracted_expr




import numpy as np
import pandas as pd

def adjust_frequency(df, filter_few=False, if_adj=False):
	count = df['count']
        
	if filter_few or not if_adj:
		total_sum = count.sum()
		freq = count / total_sum            
	else:
		freq = df['freq']

	chain = df['chain']            

	if if_adj:
		count_adj = df['count_adj'].astype("Int64")
		total_sum = count_adj.sum()
		freq_adj = count_adj / total_sum
		chainAdjSum = df.groupby('chain', sort=False)['count_adj'].transform('sum')
		freq_adj_chain = np.nan_to_num(count_adj / chainAdjSum)

	#reorder columns
	for col in ['count', 'freq', 'chain', 'count_adj', 'countAdj', 'freqAdj', 'freqAdjChain']:
		if col in df.columns:
			df.drop(columns=[col], inplace=True)

	if if_adj:
		df.insert(0, 'freqAdjChain', freq_adj_chain)
		df.insert(0, 'freqAdj', freq_adj)
	df.insert(0, 'freq', freq)
	df.insert(0, 'chain', chain)
	if if_adj:
		df.insert(0, 'countAdj', count_adj)
	df.insert(0, 'count', count)
	df = df.sort_values(by=["count"], ascending=False).reset_index(drop=True)

	return df
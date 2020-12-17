#!/usr/bin/env python
import sys
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
from scr.modules import adjust_frequency

def GlobalAlignment(v, w):
		sig = 1

		s = [[0 for i in range(len(w)+1)] for j in range(len(v)+1)]
		backtrack = [[0 for i in range(len(w)+1)] for j in range(len(v)+1)]
		
		for i in range(1, len(w)+1):
			s[0][i] = -i*sig

		for j in range(1, len(v)+1):
			s[j][0] = -j*sig

		for i in range(1, len(v)+1):
			for j in range(1, len(w)+1):
				mu = 0 if v[i-1] == w[j-1] else -1
				score_list = s[i-1][j] - sig, s[i][j-1] - sig, s[i-1][j-1] + mu
				s[i][j] = max(score_list)
				backtrack[i][j] = score_list.index(s[i][j])

		a, b = len(v), len(w)
		score = s[a][b]

		while a*b != 0:
			if backtrack[a][b] == 0:
				a -= 1
				w = w[:b] + "-"+ w[b:]
			elif backtrack[a][b] == 1:
				b -= 1
				v = v[:a] + "-"+ v[a:]
			else:
				a -= 1
				b -= 1

		while a != 0:
			w = w[:0] + "-"+ w[0:]
			a -= 1
		while b != 0:
			v = v[:0] + "-"+ v[0:]
			b -= 1

		return v, w
	
	
def count_indels(align1, align2):
	#Учитываются только инсерции и делеции, без замен
	#Если есть хоть 1 замена, возвращается inf
	#Делается на основе выравниваний
	e = 0
	for i, j in zip(align1, align2):
		if i != j:
			if (i == "-") or (j == "-"):
				e += 1
			else:
				return np.inf
	return e


def check_subclones(df, indel_thr, seq_err):
	#return an array, where index - subclone, value - real clone
	cdr3_cor = np.full(len(df), None)
	cdr3_cor_unique = set()

	for i in tqdm(reversed(range(len(df))), total=len(df)):
		#If this clone doesn't contain another subclones
		if i not in cdr3_cor_unique:
			s1 = df.iloc[i]["cdr3nt"]
			c1 = df.iloc[i]["count"]
			for x in range(0, i):

				s2 = df.iloc[x]["cdr3nt"]
				c2 = df.iloc[x]["count"]

				if c1 / (c1 + c2) <= seq_err:
					a1, a2 = GlobalAlignment(s1, s2)
					indels = count_indels(a1, a2)

					if indels <= indel_thr:
						cdr3_cor[i] = x
						cdr3_cor_unique.add(x)
						break
				else:
					break
	return cdr3_cor


def filter_subclones(df, indel_thr, seq_err):
	
	df = df.sort_values(by=["count"], ascending=False).reset_index(drop=True)
	cdr3_cor = check_subclones(df, indel_thr, seq_err)
	df_filt = df.copy()

	for i, v in enumerate(cdr3_cor):
		if v is not None:
			df_filt.at[v, "count"] += df_filt.at[i, "count"]
			df_filt = df_filt.drop(index=i)
	df_filt = df_filt.sort_values(by=["count"], ascending=False).reset_index(drop=True)
	df_filt = adjust_frequency(df_filt, filter_few=0, if_adj=False)

	n_filt = len(df) - len(df_filt)
	if n_filt == 1:
		print(f'{n_filt} clonotype was filtered')
	else:
		print(f'{n_filt} clonotypes were filtered')

	return df_filt


def main(**kwargs):
	#### Argument parsers
	args = kwargs.get('args', None)
	if args is None:
		parser = argparse.ArgumentParser(description='Filter clonotypes by indel and sequencing error thresholds.', formatter_class=argparse.RawTextHelpFormatter)

		parser.add_argument("-i", "--input", help = "Path of input VDJtools table", type=str, metavar='<input>', required=True)
		parser.add_argument("-o", "--output",  help = "Path of output VDJtools table", type=str, metavar='<output>', required=True)
		parser.add_argument("-se", "--seq_error", help = "Probable error of sequencing (default=0.01)", default=0.01, type=float, metavar='<float>')
		parser.add_argument("-id", "--indel",  help = "Maximal amount of indels to concidering CDR3s to be identical (default=1)", default=1, type=int, metavar='<int>')
		args = parser.parse_args()
		
		if len(sys.argv)==1:
			parser.print_help(sys.stderr)
			sys.exit(1) 		
	####
	#main part
	df = pd.read_csv(args.input, delimiter="\t")	
	df_filt = filter_subclones(df, indel_thr=args.indel, seq_err=args.seq_error)
	df_filt.to_csv(args.output, sep="\t", index=False)

	print("done")
	
	
if __name__ == "__main__":
	main()
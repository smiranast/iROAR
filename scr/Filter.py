#!/usr/bin/env python
import os
import sys
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
from scr.modules import *


class FilterSubclones:
    
    def __init__(self, indel_thr, seq_err, filter_type, iroar_path):
        self.indel_thr = indel_thr
        self.seq_err = seq_err
        self.iroar_path = iroar_path
        self.filter_types = filter_type
        
    @staticmethod
    def _DL_distance(s1, s2):
        d = dict()
        n, m = len(s1), len(s2)
        for i in range(-1,n+1):
            d[(i,-1)] = i+1
        for j in range(-1,m+1):
            d[(-1,j)] = j+1

        for i in range(n):
            for j in range(m):
                cost = 0 if s1[i] == s2[j] else 1

        for i in range(n):
            for j in range(m):
                d[(i,j)] = min(d[(i-1,j)] + 1, d[(i,j-1)] + 1, d[(i-1,j-1)] + cost)

        return d[n-1,m-1]


    def _compare_clonal_df(self, df1, df2, cdr3_cor, cdr3_cor_unique, filter_type):
        mess = "out-of-frame in in-frame" if filter_type == "OinI" else "in-frame in out-of-frame" if filter_type == "IinO" else "all"
        print(f'Filter {mess} clones')
        for i in tqdm(list(df2.index)[::-1]):
            #If this clone doesn't contain another subclones
            if i not in cdr3_cor_unique:
                s1 = df2.loc[i, "cdr3nt"]
                c1 = df2.loc[i, "count"]
                for x in [x for x in list(df1.index) if x < i]:

                    s2 = df1.loc[x, "cdr3nt"]
                    c2 = df1.loc[x, "count"]

                    if c1 / (c1 + c2) <= self.seq_err:
                        dl = self._DL_distance(s1, s2)

                        if dl <= self.indel_thr:
                            cdr3_cor[i] = x
                            cdr3_cor_unique.add(x)
                            break
                    else:
                        break

        return cdr3_cor, cdr3_cor_unique
        

    def _check_subclones(self, df):
        #return an array, where index - subclone, value - real clone
        cdr3_cor = np.full(len(df), None)
        cdr3_cor_unique = set()

        #Compare all clones vs all clones
        if "all" in self.filter_types:
            cdr3_cor, cdr3_cor_unique = self._compare_clonal_df(df, df, cdr3_cor, cdr3_cor_unique, "all")
            return cdr3_cor

        gene_names = gene_names_parser(self.iroar_path, df["chain"].unique())
        df_inframe = gene_names.filter_outframes(df, outframes=False, save_index=False)
        df_outframe = gene_names.filter_outframes(df, outframes=True, save_index=False)

        #Compare outframe clones within inframe clones
        if "OinI" in self.filter_types:
            cdr3_cor, cdr3_cor_unique = self._compare_clonal_df(df_inframe, df_outframe, cdr3_cor, cdr3_cor_unique, "OinI")

        #Compare inframe clones within outframe clones
        if "IinO" in self.filter_types:
            cdr3_cor, cdr3_cor_unique = self._compare_clonal_df(df_outframe, df_inframe, cdr3_cor, cdr3_cor_unique, "IinO")

        df.drop(columns=["ind"], inplace=True)
        return cdr3_cor


    def filter_subclones(self, df):

        df = df.sort_values(by=["count"], ascending=False).reset_index(drop=True)
        df["chain"] = get_chain(df)
        cdr3_cor = self._check_subclones(df)
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
        parser.add_argument("-ft", "--filter_type",  help = "Which frame groups are compared during the filtering (default=all)", choices=['IinO', 'OinI', 'all'], default="all", metavar='<list>', nargs='+')
        args = parser.parse_args()
        
        if len(sys.argv)==1:
            parser.print_help(sys.stderr)
            sys.exit(1)         
    ####
    #main part
    iroar_path = os.path.dirname(os.path.realpath(__file__))
    iroar_path = "/".join(iroar_path.split("/")[:-1])

    filter_type = args.filter_type
    
    df = pd.read_csv(args.input, delimiter="\t")   
    Filter = FilterSubclones(indel_thr=args.indel, seq_err=args.seq_error, filter_type=args.filter_type, iroar_path=iroar_path)
    df_filt = Filter.filter_subclones(df)
    df_filt.drop(columns=["chain"], inplace=True)
    df_filt.to_csv(args.output, sep="\t", index=False)

    print("done")
    
    
if __name__ == "__main__":
    main()
#!/usr/bin/env python3

import sys
import argparse
import numpy as np
import pandas as pd

class OARbias:
    
    def __init__(self, seed=0, method='vjmplex', noise=False):
        self.seed = seed
        self.m = method
        self.noise = False

    def oarbias_gene(self, df, gene_type):
        # Add multiplex bias, a particular OAR for each gene
        gene_list = df[gene_type].unique()
        gene_dict = dict(zip(gene_list, logbias(len(gene_list), self.seed)))
        df['count'] = df.apply(lambda x: x['count'] * gene_dict.get(x[gene_type]), axis=1)
        return df
        
    def add_bias(self, df):
        df_adj = self.oarbias_gene(df, 'v')
        df_adj = df_adj if self.m == 'vmplex' else self.oarbias_gene(df_adj, 'j')

        if self.noise:
            #Random clone-secific bias from normal distribution
            df_adj["count"] = (df_adj["count"]*normnoise(len(df_adj), seed=self.seed)).astype(int)
        
        df_adj["count"] = (df_adj["count"]).astype(int)
        df_adj.loc[df_adj["count"] == 0, "count"] = 1
        df_adj["freq"] = (df_adj["count"].astype(float) / df_adj["count"].sum()).apply(self.formatfreqs)
        df_adj.sort_values(by="count", ascending=False, inplace=True)
        return df_adj
    
    @staticmethod 
    def formatfreqs(freq):
        # Returns freq column in the save format as in original VDJtools table
        if freq > 0.001:
            freq = format(freq, ".16f")
            return freq.rstrip("0")
        else:
            freq = "{:.16e}".format(freq).split("e")
            return freq[0].rstrip("0") + "E" + freq[1]

        
def logbias(n=1, seed=0):
    np.random.seed(seed)
    a = 2**np.random.normal(1, 0.21, n) - 1
    ind = sorted(np.random.choice(np.arange(0,n), n // 2, replace=False))
    a[ind] **= -1
    return a

def normnoise(n=1, v=0.025, seed=0):
    np.random.seed(seed)
    return np.random.normal(1, v, n)
            
    
def main(**kwargs):
    #### Argument parsers
    args = kwargs.get('args', None)
    if args is None:
        parser = argparse.ArgumentParser(description='Introduce an artificial OAR bias in VDJtools repertoire table.', formatter_class=argparse.RawTextHelpFormatter)
        parser.add_argument("-i", "--input", help = "Path of input VDJtools table", type=str, metavar='<input>', required=True)
        parser.add_argument("-o", "--output",  help = "Path of biased VDJtools table ", type=str, metavar='<output>', required=True)
        parser.add_argument("--seed",  help = "Randomization seed (default=0)", type=int, metavar='<int>', default=0)
        parser.add_argument("-m", "--method", help = "Method of sequencing library preparation (default=vjmplex)", choices=['vjmplex', 'vmplex'], default='vjmplex')
        parser.add_argument("--noise", help = "Add random noise in bias for each clone", action='store_true')
        
        args = parser.parse_args()
        
        if len(sys.argv)==1:
            parser.print_help(sys.stderr)
            sys.exit(1)       
    ####
    
    #main part
    df = pd.read_csv(args.input, sep="\t")
    oarbias = OARbias(seed=args.seed, method=args.method, noise=args.noise)
    df_adj = oarbias.add_bias(df)
    df_adj.to_csv(args.output, index=False, sep="\t")
    
if __name__ == "__main__":
    main()

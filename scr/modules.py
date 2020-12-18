import json
import numpy as np
import pandas as pd


class gene_names_parser:
    
    def __init__(self, iroar_path, chains):
        self.iroar_path = iroar_path
        self.chains = chains
        
    def get_gene_list(self, region, func):
        #Functional genes - "F" only
        #Unfunctional - "ORF" and pseudogenes ("P")
        gene_file = f"{self.iroar_path}/aux/germline_{region}"
        with open(gene_file, "r") as file:
            if func == "all":
                content = {k: [i["name"] for i in v] for k,v in json.load(file).items() if k in self.chains}
            elif func == "y":
                content = {k :[i["name"] for i in v if i["type"] == "F"]
                           for k,v in json.load(file).items() if k in self.chains}
            elif func == "n":
                content = {k :[i["name"] for i in v if i["type"] == "ORF" or i["type"] == "P"]
                           for k,v in json.load(file).items() if k in self.chains}
            return content
     
    def filter_outframes(self, df, outframes, save_index=True):
        genes_notfunc_v = sum(self.get_gene_list("V", func="n").values(), [])
        genes_notfunc_j = sum(self.get_gene_list("J", func="n").values(), [])
        outframe_genes_mask = (df['v'].isin(genes_notfunc_v) | df['j'].isin(genes_notfunc_j))
        uncomplete_mask = (df["v"].str.contains("\+") |
                           df["j"].str.contains("\+") | 
                           df["v"].str.contains("Intron") |
                           df["j"].str.contains("KDE"))
        outframe_cdr3_mask = (df["cdr3aa"].str.contains("_") | 
                            df["cdr3aa"].str.contains("\*") |
                            ~df["cdr3aa"].str.startswith("C") |
                            ~(df["cdr3aa"].str.endswith("W") | df["cdr3aa"].str.endswith("F")))
            
        df["ind"] = list(df.index)
        if outframes:
            df = df[uncomplete_mask | outframe_cdr3_mask | outframe_genes_mask ]
        else:
            df = df[~(uncomplete_mask | outframe_cdr3_mask | outframe_genes_mask)]
        if not save_index:
            df.index = df["ind"]
            del df.index.name
        df.drop(columns=["ind"], inplace=True)
        return df
        

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


#Add chain column because standard VDJtools output format doesn't have it
def get_chain(df): 
    return df['j'].astype(str).str[:3]
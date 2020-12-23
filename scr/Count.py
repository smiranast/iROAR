#!/usr/bin/env python

import json
import re
import os
import sys
import argparse
import requests
import random
import string
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import ceil, isnan
from datetime import datetime
from scr.modules import *
from scr.Merge import mean_oar_counter
from scr.Filter import FilterSubclones
pd.options.mode.chained_assignment = None
np.seterr(divide='ignore', invalid='ignore') #Disable numpy warnings (on OAR_module.clones_statistics)


class Bad_segment_counter:
    def __init__(self):
        self.bad_list = []
    
    def add(self, segment):
        self.bad_list.append(segment)
    
    def check(self, segment):
        if segment in self.bad_list:
            return True
        else:
            return False

        
class OAR_counter:
    
    def __init__(self, iroar_path, outframe, own_table, no_outliers, oar_lib_v, oar_lib_j, n_iter, err, filter_few, upper_only, min_oof_clones, short, prefilt, collision_filter):
        
        self.outframe = outframe
        self.own_table = own_table
        self.no_outliers = no_outliers
        self.oar_lib_v = oar_lib_v
        self.oar_lib_j = oar_lib_j
        self.err = err
        self.n_iter = n_iter
        self.iroar_path = iroar_path
        self.filter_few = filter_few
        self.upper_only = upper_only
        self.min_oof_clones = min_oof_clones
        self.freq_clones_v = None
        self.freq_clones_j = None
        self.genes_v = None
        self.genes_j = None
        self.clonal_freq_warning = True
        self.short = short
        self.prefilt=prefilt
        self.collision_filter=collision_filter
            
            
    #Add "-D" or "-A" to the end of TRA/TRD genes to distinguish normal and hybrid
    def _process_hyrid_v(self, df):
        mask = df["v"].str.contains('TRAV') | df["v"].str.contains('TRDV')
        df.loc[mask, "v"] = df[mask].apply(lambda x: "/".join([i + "-" + x["j"][2] for i in x["v"].split("/")]), axis=1)
        return df


    #filter outliers by count in each chain using z-score
    def filter_outliers(self, df):
        df_new = pd.DataFrame()
        """
        Save only clones which:
        (EITHER zscore is less than 3 -> only upper border cause we need to filter out supposed tumor clones
        OR has zscore = NaN -> only one clone within chain)
        AND has more than 10 reads -> zscore test works incorrectly with small ranges
        """
        mask_zscore = (df["count"] - df["count"].mean())/df["count"].std()
        mask_zscore = mask_zscore.fillna(0)
        mask_threshold = df["count"] < 10
        df_new = pd.concat([df_new, df.loc[(mask_zscore  <= 3) | (mask_threshold)]], axis=0)
        df_bad = pd.concat([df_new, df.loc[(mask_zscore  > 3) & (~mask_threshold)]], axis=0)
        #warning in case there are unique V/J genes in outlier clones that don't meet in normal ones
        all_bad = list(df_bad["v"].unique()) + list(df_bad["j"].unique())
        all_new = list(df_new["v"].unique()) + list(df_new["j"].unique())
        uniq_bad = [i for i in set(all_bad) if i not in set(all_new)]
        if uniq_bad:
            print("WARNING! Following genes meet only in outlier clones: {}".format(" ,".join([str(i) for i in uniq_bad])), file = sys.stderr)   
        return df_new
    
    
    """
    #previous function, which calculates for each chain separately
    def filter_outliers_old(self, df):
        df_new = pd.DataFrame()
        for chain in df["chain"].unique():
            df_chain = df.loc[df["chain"] == chain]
            mask_zscore = (df_chain["count"] - df_chain["count"].mean())/df_chain["count"].std()
            mask_zscore = mask_zscore.fillna(0)
            mask_threshold = df_chain["count"] < 10
            df_new = pd.concat([df_new, df_chain.loc[(mask_zscore  <= 3) | (mask_threshold)]], axis=0)
            df_bad = pd.concat([df_new, df_chain.loc[(mask_zscore  > 3) & (~mask_threshold)]], axis=0)
            #warning in case there are unique V/J genes in outlier clones that don't meet in normal ones
            all_bad = list(df_bad["v"].unique()) + list(df_bad["j"].unique())
            all_new = list(df_new["v"].unique()) + list(df_new["j"].unique())
            uniq_bad = [i for i in set(all_bad) if i not in set(all_new)]
            if uniq_bad:
                print("WARNING! Following genes meet only in outlier clones: {}".format(" ,".join([str(i) for i in uniq_bad])), file = sys.stderr)
        return df_new
    """
        
    def _check_clones(self, df, gene_dict, gene_type):
        #get the precalculated numbers of clones for each gene
        freq_clones_precalc_raw, freq_clones_precalc = dict(), dict()
        
        freq_clones_precalc_raw = pd.read_csv(f"{self.iroar_path}/aux/clonal_freq_{gene_type}.csv",
                                              sep=",", index_col="genes")
        freq_clones_precalc_raw = freq_clones_precalc_raw.to_dict(orient='index')
        
        #form freq_clones_precalc_raw without stdev
        for k, v in freq_clones_precalc_raw.items():
            freq_clones_precalc[k] = v["AVG"]
            
        #get numbers of clones from the input dataset
        freq_clones = {chain: dict() for chain in self.chains}
        
        if_good = True

        for chain in self.chains:
            total_clones_chain = df.loc[df["chain"] == chain].shape[0]
            if total_clones_chain == 0:
                for gene in gene_dict[chain]:
                    freq_clones[gene] = 0
            else:
                #use pseudocounts to get rid from freq = 0
                #despite situations when there is no clones for the chain
                for gene in gene_dict[chain]:
                    freq_clones[gene] = df.loc[df[gene_type.lower()].str.contains(gene), ["count", gene_type.lower()]].shape[0]
                    freq_clones[gene] += 1
                    freq_clones[gene] /= total_clones_chain
                    popul_data = freq_clones_precalc_raw.get(gene, None)
                    if popul_data: #if population frequence for this gene calculated
                        #three sigma rule is used
                        sigma = 3*popul_data["STDEV"]
                        if freq_clones[gene] > popul_data["AVG"] + sigma or freq_clones[gene] < popul_data["AVG"] - sigma:
                            if_good = False

        #if all frequencies are similar to populational ones
        return freq_clones, freq_clones_precalc, if_good
    

    def _clones_statistics(self, df, gene_list, outframe, gene_type):
        gene_reads, gene_clones, chains = [], [], []
        gene_X, gene_Y, gene_OAR, gene_oof_n = [], [], [], []
        freq_clones = self.freq_clones_v if gene_type == "V" else self.freq_clones_j
        
        total_reads = {chain: df.loc[df["chain"] == chain]["count"].sum() for chain in
                       ["IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG"]}

        for gene in gene_list:
            df_gene = df.loc[df[gene_type.lower()].str.contains(gene), ["count", gene_type.lower()]]
            X = freq_clones[gene]
            if df_gene.empty:
                reads, clones = 0, 0
            else:
                reads = df_gene["count"].astype(float).sum()
                reads = df_gene.apply(lambda x: x["count"] / (x[gene_type.lower()].count("/") + 1), axis=1).sum()
                clones = df_gene.shape[0]
            if X == 0 and gene[-2:] != "-D" and gene[-2:] != "-A":
                chains.append(gene[:3])
                Y = 0
            else:
                if gene[-2:] == "-D":
                    chains.append("TRD")
                    Y = reads / total_reads["TRD"] if X != 0 else 0
                elif gene[-2:] == "-A":
                    chains.append("TRA")
                    Y = reads / total_reads["TRA"] if X != 0 else 0
                else:
                    chains.append(gene[:3])
                    Y = reads / total_reads[gene[:3]]

            gene_reads.append(int(reads))
            gene_clones.append(int(clones))
            
            OAR = Y / X if X != 0 and Y != 0 else np.nan
            if (not (isnan(OAR))
                and ((outframe and int(clones) < self.min_oof_clones) # If this gene has less out-of-frame than N
                or (self.upper_only and OAR < 1) #Change OAR coefficient to not less than 1
                or (self.v_only and gene_type == "j"))): # Ignore J-genes OAR for the further calculation
                    OAR = 1.
                
            gene_OAR.append(OAR)
            
            if not isnan(OAR) and outframe and int(clones) >= self.min_oof_clones:
                gene_oof_n.append(int(clones))
            else:
                gene_oof_n.append(0)                      
                  
        return {k: {"chain": v1, "reads": v2, "clones": v3, "OAR": v4, "out-of-frame": v5, "sample": None}
                for k, v1,v2,v3,v4,v5 in zip(gene_list, chains, gene_reads, gene_clones, gene_OAR, gene_oof_n)}


    def count_OAR(self, df, outframe, verbosity):
        gene_names = gene_names_parser(self.iroar_path, self.chains)
        if outframe:
            df = gene_names.filter_outframes(df, outframes=True)
            
        if not (self.freq_clones_v and self.freq_clones_j):
            self.freq_clones_v, freq_clones_precalc_v, if_good_v = self._check_clones(df, self.genes_dict_v, "V")
            self.freq_clones_j, freq_clones_precalc_j, if_good_j = self._check_clones(df, self.genes_dict_j, "J")

            if (not if_good_v or not if_good_j) and self.clonal_freq_warning:
                print("Clonal frequencies of the input data differs sufficiently to populational!")
                yn = input("Use populational frequencies instead? yes/no\n").lower().strip("\n")
                while not (yn == "yes" or yn == "y" or yn == "no" or yn == "n"):
                    yn = input("Print 'yes' or 'no' ('y' or 'n')\n").lower().strip().strip("\n")
                if yn == "yes" or yn == "y":
                    self.freq_clones_v = freq_clones_precalc_v
                    self.freq_clones_j = freq_clones_precalc_j
                self.clonal_freq_warning = False
        
        V_stat_dict = self._clones_statistics(df, self.genes_v, outframe, "V")
        J_stat_dict = self._clones_statistics(df, self.genes_j, outframe, "J")
 
        s = "out-of-frame" if outframe else "all"
    
        for v in V_stat_dict.values():
            if np.isnan(v["OAR"]) == False:
                v["sample"] = s
        for v in J_stat_dict.values():
            if np.isnan(v["OAR"]) == False:
                v["sample"] = s

        return V_stat_dict, J_stat_dict


    @staticmethod        
    def clone_recount(row, V_stat, J_stat, bsc, verbosity):
        
        def count_aver_OAR(reg_stat, region):
            sum_OAR = 0
            for i in row[region.lower()].split("/"):
                if i in reg_stat.keys() and not isnan(reg_stat[i]["OAR"]):
                    sum_OAR += reg_stat[i]["OAR"]
                else:
                    if not bsc.check(i):
                        if verbosity:
                            print(f"Warning! There is no {i} segment in OAR dictionary", file = sys.stderr)
                        bsc.add(i)
                    sum_OAR += 1.
                    
            aver_OAR = sum_OAR / (row[region.lower()].count("/") + 1)
            return aver_OAR
        
        try:
            aver_V_OAR = count_aver_OAR(V_stat, "V")
            aver_J_OAR = count_aver_OAR(J_stat, "J")
            
            return row["count"] / (aver_V_OAR * aver_J_OAR)
        except:
            print("Something went wrong with row \"{}\"".format(" ".join([str(i) for i in row.tolist()])), file = sys.stderr)
            raise
    
    
    def _cloneCount_adjust_pair(self, tmp_file, bsc, recalc, verbosity):
        with open(tmp_file, 'rb') as f:
            df = pickle.load(f)
        
        #Calculates V_dict and J_dict
        if self.outframe:
            V_dict, J_dict =  self.count_OAR(df, outframe=True, verbosity=verbosity)
            V_dict_a, J_dict_a = self.count_OAR(df, outframe=False, verbosity=verbosity)
            V_dict.update({k:v for k,v in V_dict_a.items() if np.isnan(V_dict[k]["OAR"]) or V_dict[k]["OAR"] == 0})
            J_dict.update({k:v for k,v in J_dict_a.items() if np.isnan(J_dict[k]["OAR"]) or J_dict[k]["OAR"] == 0})
        else:
            V_dict, J_dict = self.count_OAR(df, outframe=False, verbosity=verbosity)
        
        if recalc:            
            #adjust previuos count with new V_dict and J_dict for the next iteration. Use simple variant of dict
            df["count"] = df.apply(self.clone_recount, args=(simple_dict(V_dict), simple_dict(J_dict), bsc, verbosity), axis=1)

            with open(tmp_file, 'wb') as f:
                pickle.dump(df, f)
            
        return V_dict, J_dict
                
        
    def _calc_iter_error(self, V_dict, J_dict, i_iter):
        v_OARs = np.array([v["OAR"] for v in V_dict.values()])
        j_OARs = np.array([v["OAR"] for v in J_dict.values()])
        
        v_OARs = v_OARs[~np.isnan(v_OARs)]

        V_error = abs(v_OARs - np.ones(v_OARs.shape[0]))
        V_error_max = 0 if len(V_error) == 0 else V_error.max()
        if not self.v_only:
            j_OARs = j_OARs[~np.isnan(j_OARs)]
            J_error = abs(j_OARs - np.ones(j_OARs.shape[0]))
            J_error_max = 0 if len(J_error) == 0 else J_error.max()
            error = std_one(np.append(v_OARs, j_OARs))
            verb_mess = f'Iteration: {i_iter}\tStd error: {round(error, 7)}\tMax abs V error: {round(V_error_max, 7)}\tMax abs J error: {round(J_error_max, 7)}'
        else:
            J_error = np.zeros(1)
            J_error_max = 0 if len(J_error) == 0 else J_error.max()
            error = std_one(v_OARs)
            verb_mess = f'Iteration: {i_iter}\tStd error: {round(error, 7)}\tMax abs V error: {round(V_error_max, 7)}'
        return error, verb_mess
        
    
    def cloneCount_adjust(self, input_list, output_dir, chains, v_only, verbosity):
        """
        V_dict(_a), J_dict(_a) - dictionaries of OAR, containing meta info, for current iteration
        self.V_dict_meta, self.J_dict_meta - dictionaries for all iterations with meta info. "OAR_log" - logged, "OAR" - final
        "simple dictionary" - dictionary, used for clone_recount. Format: {<gene_name, str>: <OAR, float>}
        """
        
        tmp_dir = output_dir + '/tmp'
        
        self.chains = chains
        self.v_only = v_only
        bsc = Bad_segment_counter()
        
        input_list_processed = [tmp_dir + "/" + os.path.splitext(os.path.basename(file))[0]+ ".pr" for file in input_list]
        input_list_filtered = [tmp_dir + "/" + os.path.splitext(os.path.basename(file))[0]+ ".flt" for file in input_list]
        
        
        if not (self.genes_v and self.genes_j):
            gene_names = gene_names_parser(self.iroar_path, self.chains)
            self.genes_dict_v = gene_names.get_gene_list("V", func="all")
            self.genes_dict_j = gene_names.get_gene_list("J", func="all")
            self.genes_v = sum(self.genes_dict_v.values(), [])
            self.genes_j = sum(self.genes_dict_j.values(), [])
        
        self.V_dict_meta = {k: {"chain": None, "reads": [], "clones": [], 
                                "OAR_log": [], "out-of-frame": None,
                                "sample": None} for k in self.genes_v}
        self.J_dict_meta = {k: {"chain": None, "reads": [], "clones": [], 
                                "OAR_log": [], "out-of-frame": None,
                                "sample": None} for k in self.genes_j}
        
        for file, processed_file, filtered_file in zip(input_list, input_list_processed, input_list_filtered):
            df = pd.read_csv(file, delimiter="\t")
            df["chain"] = get_chain(df)
            df = df[df['chain'].isin(self.chains)]

            df = self._process_hyrid_v(df)
            with open(processed_file, 'wb') as f:
                pickle.dump(df, f)
                
            df_f = self.collision_filter.filter_subclones(df) if self.prefilt else df
            df_f = self.filter_outliers(df_f) if self.no_outliers else df_f
            with open(filtered_file, 'wb') as f:
                pickle.dump(df_f, f)
            
            #get clone statistics for initial clone counts
            V_dict, J_dict = self._cloneCount_adjust_pair(filtered_file, bsc, recalc=False, verbosity=verbosity)
            
            if verbosity:
                print(f'Calculating primary reads/clones statistics for {os.path.basename(file)}')

            for k, v in self.V_dict_meta.items():
                v["chain"] = V_dict[k]["chain"]
                v["reads"].append(V_dict[k]["reads"])
                v["clones"].append(V_dict[k]["clones"])
                          
            for k, v in self.J_dict_meta.items():
                v["chain"] = J_dict[k]["chain"]
                v["reads"].append(J_dict[k]["reads"])
                v["clones"].append(J_dict[k]["clones"])
        
        i_iter = 0 if self.own_table else 1
        error = 0              
        
        while (i_iter <= 1) or (i_iter <= self.n_iter) and (error > self.err):
            """
            Iteration 1: recalculate using input libraries and save adjusted tables to temporary files           
            Of note: all files must have suffix "V_OAR_stat" or "J_OAR_stat!!"   
            """
            if (i_iter == 1 and not self.own_table) or (i_iter == 0 and self.own_table):
                if not self.own_table:
                    if self.v_only:
                        mean_oar = mean_oar_counter(search_dir=None,
                                                    file_list=self.oar_lib_v,
                                                    min_outframe=self.min_oof_clones)
                        mean_j = {k:{"OAR": 1., "out-of-frame": 0} for k in self.genes_j}
                        self.oar_lib_j = [None for _ in range(len(self.oar_lib_v))]
                    else:
                        mean_oar = mean_oar_counter(search_dir=None,
                                                    file_list=self.oar_lib_v + self.oar_lib_j,
                                                    min_outframe=self.min_oof_clones)
                        mean_j = mean_oar.count_mean("J")
                        mean_j = {k: ({"OAR": 1., "out-of-frame": 0} if (self.upper_only and v < 1) or
                                      np.isnan(v["OAR"]) else v)
                                      for k,v in mean_j.items() if k in self.genes_j}

                    mean_v = mean_oar.count_mean("V")
                    mean_v = {k: ({"OAR": 1., "out-of-frame": 0} if (self.upper_only and v < 1) or
                                  np.isnan(v["OAR"]) else v)
                                  for k,v in mean_v.items() if k in self.genes_v}
                    
                    for k, v in mean_v.items():
                        self.V_dict_meta[k]["OAR_log"] = [v["OAR"]]
                        
                    for k, v in mean_j.items():
                        self.J_dict_meta[k]["OAR_log"] = [v["OAR"]]
                    
                else:
                    mean_v = {k:{"OAR": 1., "out-of-frame": 0} for k in self.genes_v}
                    mean_j = {k:{"OAR": 1., "out-of-frame": 0} for k in self.genes_j}
                    
                tmp_df_files, tmp_V_dict_files, tmp_J_dict_files = [], [], []
                    
                for i, (oar_lib_v, oar_lib_j) in enumerate(zip(self.oar_lib_v, self.oar_lib_j)):
                    if self.v_only and verbosity and not self.own_table:
                        print(f"Loading library {oar_lib_v}")
                    elif verbosity and not self.own_table:
                        print(f"Loading libraries {oar_lib_v} and {oar_lib_j}")
                            
                    for filtered_file in input_list_filtered:
                        #Create temp binary file for each pair of libraries
                        tmp_file = tmp_dir + "/" + ''.join([random.choice(string.ascii_lowercase) for _ in range(10)]) + str(i)
                        tmp_df_files.append(tmp_file)
                        tmp_V_dict_files.append(tmp_file + ".V_OAR_stat.json")
                        tmp_J_dict_files.append(tmp_file + ".J_OAR_stat.json")

                        if self.v_only and not self.own_table:
                            with open(oar_lib_v) as f1:
                                V_dict = json.load(f1)
                                J_dict = {k:1. for k in j_list}
                        elif not self.own_table:
                            with open(oar_lib_v) as f1, open(oar_lib_j) as f2:
                                V_dict, J_dict = json.load(f1), json.load(f2)
                        else:
                            V_dict = {k:{"OAR": 1., "out-of-frame": 0} for k in self.genes_v}
                            J_dict = {k:{"OAR": 1., "out-of-frame": 0} for k in self.genes_j}

                        with open(filtered_file, 'rb') as f:
                            df_f = pickle.load(f)
                            
                        df_f["count"] = df_f.apply(self.clone_recount, args=(V_dict, J_dict, bsc, verbosity), axis=1)

                        with open(tmp_file, 'wb') as f:
                            pickle.dump(df_f, f)
                
            if (i_iter != 1 and not self.own_table) or (i_iter != 0 and self.own_table):
                for tmp_file, tmp_V_dict_file, tmp_J_dict_file in zip(tmp_df_files, tmp_V_dict_files, tmp_J_dict_files):
                    V_dict, J_dict = self._cloneCount_adjust_pair(tmp_file, bsc, recalc=True, verbosity=verbosity)

                    with open(tmp_V_dict_file, 'w') as f1, open(tmp_J_dict_file, 'w') as f2:
                        json.dump(simple_dict(V_dict), f1)
                        json.dump(simple_dict(J_dict), f2)
                                
                mean_oar = mean_oar_counter(search_dir=None,
                                            file_list=tmp_V_dict_files + tmp_J_dict_files,
                                            min_outframe=0)
                mean_v = mean_oar.count_mean("V")                    
                mean_j = mean_oar.count_mean("J")
                
                for k, v in self.V_dict_meta.items():
                    if i_iter == 2:
                        v["out-of-frame"] = V_dict[k]["out-of-frame"]
                        v["sample"] = V_dict[k]["sample"]
                    v["OAR_log"].append(mean_v[k]["OAR"])

                for k, v in self.J_dict_meta.items():
                    if i_iter == 2:
                        v["out-of-frame"] = J_dict[k]["out-of-frame"]
                        v["sample"] = J_dict[k]["sample"]
                    v["OAR_log"].append(mean_j[k]["OAR"])   
                
            if self.own_table and i_iter == 0:
                i_iter += 1
                continue
                
            error, verb_mess = self._calc_iter_error(mean_v, mean_j, i_iter)
            
            i_iter += 1
            if verbosity:
                print(verb_mess)
                    
                    
        #remove created temp files
        for file in input_list_filtered + tmp_df_files + tmp_V_dict_files + tmp_J_dict_files:
            if os.path.isfile(file):
                os.remove(file)
                 
        #After all iterations count total OAR
        for k, v in self.V_dict_meta.items():
            v["OAR"] = np.prod(np.array(v["OAR_log"]))
        for k, v in self.J_dict_meta.items():
            v["OAR"] = np.prod(np.array(v["OAR_log"]))
        #Finally, create simple dictionaries for recalculation
        V_dict, J_dict = simple_dict(self.V_dict_meta), simple_dict(self.J_dict_meta)
        
        for file, processed_file in zip(input_list, input_list_processed):
            with open(processed_file, 'rb') as f:
                df = pickle.load(f)
                
            #Recount using ready simple dict
            df["count_adj"] = df.apply(self.clone_recount, args=(V_dict, J_dict, bsc, verbosity), axis=1)

            #Normalize to new summary readcount
            df.loc[:,'count_adj'] *= df["count"].sum() / df.loc[~np.isnan(df["count_adj"]), "count_adj"].sum()

            #Filter clones which counts before ceiling are less than X
            if self.filter_few:
                df = df.loc[df['count_adj'] > self.filter_few]

            df.loc[~np.isnan(df["count_adj"]), "count_adj"] = df.loc[~np.isnan(df["count_adj"]), "count_adj"].apply(ceil)

            #remove -D and -A suffixes
            df["v"] = df["v"].apply(lambda x: "/".join(i.strip('-D') for i in x.split("/")))
            df["v"] = df["v"].apply(lambda x: "/".join(i.strip('-A') for i in x.split("/")))
            
            #recalculate frequency
            df = adjust_frequency(df, self.filter_few, if_adj=True)
    
            #adjust table to VDJtools output format
            if self.short:
                df = to_vdjtools_format(df)
            #write adjusted table
            prefix = os.path.splitext(os.path.basename(file))[0]
            output_file = output_dir + "/" + prefix + ".OAR.txt"
            df.to_csv(output_file, sep="\t", index=False)
            
            if os.path.isfile(processed_file):
                os.remove(processed_file)
    
    
def std_one(x):
    #input - np.array
    #Calculated standard deviation where mean = 1
    std = np.sqrt(np.power(x - np.ones(x.shape[0]), 2).mean())
    std = 0 if np.isnan(std) else std
    return std
 

def simple_dict(data):
    return {k: {"OAR": v["OAR"], "out-of-frame": v["out-of-frame"]} for k, v in data.items()}            


def download_germline(region, iroar_path, verbosity):
    gene_file = f"{iroar_path}/aux/germline_{region}"
    chains = ["IG/IGH", "IG/IGK", "IG/IGL", "TR/TRA", "TR/TRB", "TR/TRD", "TR/TRG"]
    gene_names_total_tmp = {k.split("/")[1]: [] for k in chains}
    gene_names_total = {k.split("/")[1]: [] for k in chains}
    for chain in chains:
        if verbosity:
            print(f"Downloading {chain} {region} region dataset")
        url = f"http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/{chain}{region}.fasta"
        content = requests.get(url).content.decode("utf-8").split("\n")
        #Get gene names and types (functional, ORF or pseudogene)
        gene_names = []
        gene_names_tmp = sorted(list(set([(i.split("|")[1].split("*")[0].replace("/", ""),
                                       i.split("|")[3].replace("(", "").replace(")", "").replace("[", "").replace("]", ""))
                                      for i in content if i and i[0] == ">"])))
        
        #Get rid of dublicates of functional, ORF and pseudogenes alleles the same time
        for gene_name in gene_names_tmp:
            if gene_name[1] == "F":
                gene_names.append(gene_name)
            elif not (gene_name[0], "F") in gene_names_tmp:
                if gene_name[1] == "P" and not (gene_name[0], "ORF") in gene_names_tmp or gene_name[1] == "ORF":
                    gene_names.append(gene_name)
        
        if (chain == "TR/TRA" or chain == "TR/TRD") and region == "V":
            gene_names_total_tmp["TRA"] += [(i[0] + "-A", i[1]) for i in gene_names]
            gene_names_total_tmp["TRD"] += [(i[0] + "-D", i[1]) for i in gene_names]
        else:
            if chain == "IG/IGK" and region == "V":
                gene_names.append(("IGK-Intron", "P"))
            if chain == "IG/IGK" and region == "J":
                gene_names.append(("IGK-KDE", "P"))
            gene_names_total_tmp[chain.split("/")[1]] = gene_names
        
    #change to JSON format
    for chain in chains:
        gene_names_total[chain.split("/")[1]] = [{"name":i[0], "type":i[1]}
                                                 for i in gene_names_total_tmp [chain.split("/")[1]]]
    with open(gene_file, 'w') as file:
        json.dump(gene_names_total, file, indent=4)

           
def vdjtools_format_check(file):
    df = pd.read_csv(file, delimiter="\t")
    #check "count" column separately as it can be both int and float
    if not "count" in df.columns:
        return False
    elif not all((df["count"].map(type) == int).array + (df["count"].map(type) == float).array):
        return False
        
    columns = {"freq": float, "cdr3nt": str, "cdr3aa": str, "v": str, "d": str, "j":str}
    for column, exp_type in columns.items():
        if not column in df.columns:
            return False
        elif not all(df[column].map(type) == exp_type):
            return False
        
    return True
    
    
def to_vdjtools_format(df):
    df["count"] = df["countAdj"]
    df["freq"] = df["freqAdj"]
    df.drop(columns=['countAdj', 'freqAdj', 'freqAdjChain', 'chain'], inplace=True)
    df = df.sort_values(by=["count"], ascending=False).reset_index(drop=True)
    return df
    
    
def json_df_converter(json_data):
    genes, chain, reads, clones, oar, sample = [], [], [], [], [], []
    for k, v in json_data.items():
        genes.append(k)
        chain.append(v["chain"])
        reads.append(v["reads"])
        clones.append(v["clones"])
        oar.append(v["OAR"])
        sample.append(v["sample"])
    return pd.DataFrame({"genes": genes, "chain": chain, "reads": reads, "clones": clones, "OAR": oar, "sample": sample})
    
    
def draw_plot(std_dict, iter_range, threshold, gene_type, path):
    fig, ax = plt.subplots(figsize=(12, 8))
    for chain, errors in std_dict.items():
        plt.plot(iter_range, errors, label=chain)
    plt.plot(iter_range, np.ones(max(iter_range))*threshold, "--", label="Threshold")
    plt.xlabel("Iterations", fontsize=20)
    plt.ylabel("STD", fontsize=20)
    
    xtics_names = np.append(np.arange(1, max(iter_range), max(iter_range) // 5 or 1), np.array([max(iter_range)]))
    plt.xticks(ticks=xtics_names, labels=xtics_names, fontsize=14)
    plt.yticks(fontsize=14)
    plt.title(f"Standard deviation of {gene_type} genes OAR for each iteration", fontsize=20)
    plt.grid()
    plt.legend(fontsize=14)
    plt.savefig(path, format="png")

    
def make_report(args, tables, reports_dir, prefix, region_dict, region, threshold):
    pd.set_option('display.max_colwidth', -1)
    #Info about input arguments
    df_args = pd.DataFrame.from_dict(args.__dict__, orient='index')
    df_args.drop(index=['func'], inplace=True)
    df_args.loc['voar',0] = '\n'.join(map(str, df_args.loc['voar',0]))
    df_args.loc['joar',0] = '\n'.join(map(str, df_args.loc['joar',0]))
    df_args.columns=[f'Input arguments for {prefix}']
    report_args = df_args.to_html(justify="center").replace("\\n","<br>")
    
    #Info about input data
    input_file_base = [os.path.basename(i) for i in tables]
    input_dir = args.input
    df_info = pd.DataFrame(input_file_base, index=list(np.arange(1, len(tables) + 1)))
    df_info.columns=['List of input VDJtools tables']
    report_info = df_info.to_html(justify="center")
    
    #Standard statistics tables
    df_stat = json_df_converter(region_dict)
    df_stat = df_stat.replace(np.nan, '', regex=True)
    
    #split "reads" and "clones" columns
    reads_cols = [f"reads_{i+1}" for i in range(len(tables))]
    clones_cols = [f"clones_{i+1}" for i in range(len(tables))]
    df_stat[reads_cols] = pd.DataFrame(df_stat["reads"].tolist(), index=df_stat.index)
    df_stat[clones_cols] = pd.DataFrame(df_stat["clones"].tolist(), index=df_stat.index)
    df_stat = df_stat[["genes", "chain"] + reads_cols + clones_cols + ["OAR", "sample"]]
    
    df_stat.columns=pd.MultiIndex.from_product([[f'General statistics information for {region} region'],df_stat.columns])
    report_stat = df_stat.to_html(index=False)
    
    #filter NaN
    region_dict = {k:v for k, v in region_dict.items() if ~np.isnan(v["OAR"])}
    
    # Report for final value
    chains = sorted(list(set([v["chain"] for v in region_dict.values()])))
    
    std_tot = [std_one(np.array([v["OAR"] for v in region_dict.values() if v['chain'] == chain])) for chain in chains]
    std_tot = [round(i, 5) for i in std_tot]

    df_std_tot = pd.DataFrame({k:v for k, v in zip(chains, std_tot)}, index=["std"])
    df_std_tot.columns=pd.MultiIndex.from_product([[f'Standard deviation of final OAR values of {region} region'],df_std_tot.columns])
    report_std_tot = df_std_tot.to_html(index=False)

    # Report for each iteration
    iter_std = dict()
    for chain in chains:
        iter_std[chain] = []
        meta_chain = {k: v for k, v in region_dict.items() if v["chain"] == chain}
        for step in np.array([v["OAR_log"] for v in meta_chain.values()]).T:
            iter_std[chain].append(std_one(step))

    iter_range = range(1, len(list(iter_std.values())[0]) + 1)
    df_std_iter = pd.DataFrame(iter_std, index=iter_range)
    df_std_iter.columns=pd.MultiIndex.from_product([[f'Standard deviation of {region} genes OAR for each iteration'],df_std_iter.columns])
    report_std_iter = df_std_iter.to_html()

    #merge all previous HTML subreports to one file
    report_total = "<br><br>".join([report_args, report_info, report_stat, report_std_tot, report_std_iter])
    
    # show report_std_iter data as plot
    if len(list(iter_std.values())[0]) > 1:
        draw_plot(iter_std, iter_range, threshold, region, path=f'{reports_dir}/{prefix}_plot_{region}.png')
        report_total += f"<br><br><img src='{prefix}_plot_{region}.png'>"
    
    #write report to file
    with open(f'{reports_dir}/report_{prefix}_{region}.html', "w") as f:
        f.write(report_total)
        
    
def write_stat_jsons(prefix, region_dict, region):
    with open(prefix + f".{region}_OAR_stat.json", 'w') as f:
        json.dump(simple_dict(region_dict), f, indent=4, allow_nan=True)

    
def main(**kwargs):
    
    ####
    #Argument parsers
    args = kwargs.get('args', None)
    if args is None:
        parser = argparse.ArgumentParser(description='Adjust clonal counts with OAR statistics.', formatter_class=argparse.RawTextHelpFormatter)
        parser.add_argument("input" , help = "Input directory containing VDJtools tables", type=str, metavar='<input>')
        parser.add_argument("output",  help = "Output directory path", type=str, metavar='<output>')
        parser.add_argument("--long", help = "Do not overwrite standard VDJtools columns instead of adding new ones\n(default=False)", action='store_false')
        parser.add_argument("-c", "--chains", help = "List of chains to analyse, sepatated by comma", type=str, metavar='<chain1>,<chain2>...<chainN>', default="IGH,IGK,IGL,TRA,TRB,TRD,TRG")
        parser.add_argument("-z", "--outliers", help = "Calculate OAR on data after filtering outliers with Z-test\n(if owntable = True; default=False)", action='store_true')
        parser.add_argument("-f" ,"--all_frame", help = "Calculate OAR using all clones, not only out-of-frame", action='store_true')
        parser.add_argument("-min_outframe", help = "Minimal out-of-frame clones threshold for OAR calculation.\nIf a segment is present within less clones than the specified value\nOAR is equal to 1(if outframe = True)", default=0, type=int, metavar='<int>')
        parser.add_argument("-filter_few", help = "Don't show clones with counts less than N before ceiling (default=show all)", default=0, type=float, metavar='<float>')
        parser.add_argument("-u", "--upper_only", help = "Adjust only counts which have OAR > 1", action='store_true')
        parser.add_argument("-m", "--method", help = "Method of sequencing library preparation (default=vjmplex)", choices=['vjmplex', 'vmplex'], default='vjmplex')
        parser.add_argument("--voar", help = "Library of OAR stats for V region in JSON format", type=str, default=[], metavar='<str>', nargs='+')
        parser.add_argument("--joar", help = "Library of OAR stats for J region in JSON format\nIf multiple V and J libs are used the number of --vlib files\nmust be equal to --jlib files!", type=str, default=[], metavar='<str>', nargs='+')
        parser.add_argument("-r" ,"--report", help = "Create report", action='store_true')
        parser.add_argument("-wj" ,"--writejson", help = "Write OAR statistics in json format (only OAR)", action='store_true')
        parser.add_argument("-iter", help = "Maximal number of iteration\n(if owntable = True; default=infinite)", type=int, default=np.inf, metavar='<int>')
        parser.add_argument("-err", help = "Maximal absolute deviation\n(if owntable = True; default=0.1)", type=float, default=0.1, metavar='<float>')
        parser.add_argument("--filter", help = "Prefilter clonotypes by indels in CDR3 nucleotide sequences\n(similar to stand-alone Filter command, default=False)", action='store_true')
        parser.add_argument("-se", "--seq_error", help = "Probable error of sequencing (default=0.01)", default=0.01, type=float, metavar='<float>')
        parser.add_argument("-id", "--indel", help = "Maximal amount of indels to concidering CDR3s to be identical (default=1)", default=1, type=int, metavar='<int>')
        parser.add_argument("-ft", "--filter_type",  help = "Which frame groups are compared during the filtering (default=all)", choices=['IinO', 'OinI', 'all'], default="all", metavar='<list>', nargs='+')
        parser.add_argument("-if", "--iterative_filter",  help = "Apply itterative collision merge (default for --filter_type=all)", action='store_true')   
        parser.add_argument("-v", "--verbosity", help = "Print messages to stdout, not only warnings and errors\n(default=False)", action='store_true')
        args = parser.parse_args()
        
        if len(sys.argv)==1:
            parser.print_help(sys.stderr)
            sys.exit(1)      
    ####
    prefix = f'iROAR_run_{datetime.now().strftime("%Y%m%d%H%M%S")}'
    input_dir = args.input.rstrip("/")
    output_dir = args.output.rstrip("/")
    verbosity = args.verbosity
    outframe = not args.all_frame
    oar_lib_v = args.voar
    oar_lib_j = args.joar
    v_only = True if args.method == "vmplex" else False
    chains = args.chains.split(",")
    
    iroar_path = os.path.dirname(os.path.realpath(__file__))
    iroar_path = "/".join(iroar_path.split("/")[:-1])
    
    #Instead of "-t" key for the program
    owntable = False if oar_lib_v or oar_lib_j else True
    
    if not oar_lib_v and oar_lib_j:
        print('Error! Vlib isn\'t specified', file = sys.stderr)  
        sys.exit(1)  
    elif oar_lib_v and not oar_lib_j and not v_only:
        print('Error! Jlib isn\'t specified', file = sys.stderr)
        sys.exit(1)
        
    if not owntable and not v_only and len(oar_lib_v) != len(oar_lib_j):
        print('Error! Different numbers of Vlibs and Jlibs specified', file = sys.stderr)
        sys.exit(1)
        
    if owntable:
        oar_lib_v, oar_lib_j = [None], [None]
    
    if not os.path.exists(iroar_path + '/aux'):
        os.makedirs(iroar_path + '/aux')
    
    if not os.path.exists(output_dir + '/tmp'):
        os.makedirs(output_dir + '/tmp')
    
    if not os.path.isfile(iroar_path + "/aux/germline_V"):
        download_germline("V", iroar_path, verbosity)
    if not os.path.isfile(iroar_path + "/aux/germline_J"):
        download_germline("J", iroar_path, verbosity)
        
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    #check for vdjtools tables the input directory
    vdjtools_tables = []
    for file in [f for f in os.listdir(input_dir) if os.path.isfile(os.path.join(input_dir, f))]:
        if vdjtools_format_check(input_dir + "/" + file):
            vdjtools_tables.append(input_dir + "/" + file)
            if verbosity:
                print(f'VDJtools table {file} is found!')
                
    if len(vdjtools_tables) == 0:
        print("No VDJtools tabled is found. Exiting")
        sys.exit(1)
    
    #main part
    collision_filter =  FilterSubclones(args.indel,
                              args.seq_error,
                              args.filter_type,
                              args.iterative_filter,
                              iroar_path)
            
    recounter = OAR_counter(iroar_path=iroar_path,
                            outframe=outframe,
                            own_table=owntable,
                            no_outliers=args.outliers,
                            oar_lib_v=oar_lib_v,
                            oar_lib_j=oar_lib_j,
                            n_iter=args.iter,
                            err=args.err,
                            filter_few=args.filter_few,
                            upper_only=args.upper_only,
                            min_oof_clones=args.min_outframe,
                            short=args.long,
                            prefilt=args.filter,
                            collision_filter=collision_filter)                    
    
    recounter.cloneCount_adjust(vdjtools_tables, output_dir, chains, v_only, verbosity)
    
    #write statistics files
    if args.report:
        reports_dir = output_dir + "/" + "iROAR_reports"
        if not os.path.exists(reports_dir):
            os.makedirs(reports_dir)
        make_report(args, vdjtools_tables, reports_dir, prefix, recounter.V_dict_meta, "V", args.err)
        if not v_only:
            make_report(args, vdjtools_tables, reports_dir, prefix, recounter.J_dict_meta, "J", args.err)
    if args.writejson: #
        write_stat_jsons(output_dir + "/" + prefix, recounter.V_dict_meta, "V")
        if not v_only:
            write_stat_jsons(output_dir + "/" + prefix, recounter.J_dict_meta, "J")
    
    if verbosity:
        print("done")


if __name__ == "__main__":
    main()
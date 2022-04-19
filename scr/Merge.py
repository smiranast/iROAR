#!/usr/bin/env python3

import argparse
import json
import re
import os
import sys
import pandas as pd
import numpy as np

class mean_oar_counter:
    
    #create lists of all V and J libraries
    def __init__(self, search_dir, file_list, min_outframe):
        self.V_stat_files = []
        self.J_stat_files = []
        self.min_outframe = min_outframe
        
        all_files = []            
        if search_dir:
            all_files += [search_dir + "/" + i for i in os.listdir(search_dir)]
        if file_list:
            all_files += file_list
          
        #remove possible duplicates and divide to V_stat and J_stat files
        for file in list(set(all_files)):
            if "J_OAR_stat" in file:
                self.J_stat_files.append(file)
            if "V_OAR_stat" in file:
                self.V_stat_files.append(file)
        
    def count_mean(self, region):

        def file_parser(file_list):
            for file in file_list:
                yield file
        
        stat_files = file_parser(self.V_stat_files) if region == "V" else file_parser(self.J_stat_files)
        #initiate null df
        df_total = pd.DataFrame()
        df_total.index.name = "genes"

        for file in stat_files:
            with open(file) as f:
                df = json.load(f)

                df_OAR = {k: (v["OAR"] if ~np.isnan(v["OAR"]) and v["out-of-frame"] >= self.min_outframe else 1)
                      for k,v in df.items()}
                df_OAR = pd.Series(df_OAR).rename('OAR')

                df_oof = {k: v["out-of-frame"] if ~np.isnan(v["out-of-frame"]) and ~np.isnan(v["OAR"]) else 0
                      for k,v in df.items()}
                df_oof = pd.Series(df_oof).rename('out-of-frame')

            #df = pd.read_json(file, typ='series', dtype=False).rename('OAR')
            df_total = df_total.merge(df_OAR, how='outer', left_index=True, right_index=True)
            df_total = df_total.merge(df_oof, how='outer', left_index=True, right_index=True)
        
        df_OAR = df_total.filter(regex='^OAR',axis=1).mean(axis=1, numeric_only=True).to_dict()
        df_oof = df_total.filter(regex='^out-of-frame',axis=1).mean(axis=1, numeric_only=True).to_dict()

        return {k: {"OAR": OAR, "out-of-frame": oof} for (k, OAR), oof in zip(df_OAR.items(), df_oof.values())}

    
def write_mean_json(prefix, region, mean_json):
    with open(f'{prefix}.{region}_OAR_mean.json', 'w') as f:
        json.dump(mean_json, f, indent=4)

            
def main(**kwargs):
    #### Argument parsers
    args = kwargs.get('args', None)
    if args is None:
        parser = argparse.ArgumentParser(description='Merge multiple OAR statistics files to one JSON library.', formatter_class=argparse.RawTextHelpFormatter, epilog="At least one type of input data list (separate files -f or directory -p)  must be given!")

        parser.add_argument("-f", "--files", help = "List of full paths of input files separated with space", type=str, metavar='<input>', nargs='+', default=None)
        parser.add_argument("-p", "--path",  help = "Full path of directory containing file for merging", type=str, metavar='<input>', default=None)
        parser.add_argument("-min_outframe", help = "Minimal out-of-frame clones threshold for mean OAR calculation.\nIf a segment in this JSON file is present within less clones\nthan the specified value OAR is equal to 1", default=0, type=int, metavar='<int>')
        parser.add_argument("-o", "--output",  help = "Prefix of merged library", type=str, metavar='<output>', required=True)
        args = parser.parse_args()
        
        if len(sys.argv)==1:
            parser.print_help(sys.stderr)
            sys.exit(1) 
            
    if not args.path and not args.files:
        parser.error("No files/directories is specified. Leaving")
        sys.exit(1)
    ####
    
    #main part
    mean_oar = mean_oar_counter(args.path, args.files, args.min_outframe)
    mean_v = mean_oar.count_mean("V")
    if mean_v:
        write_mean_json(args.output, "V",mean_v)
        
    mean_j = mean_oar.count_mean("J")
    if mean_j:
        write_mean_json(args.output, "J",mean_j)
    
    print("done")
    
    
if __name__ == "__main__":
    main()

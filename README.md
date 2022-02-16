# iROAR (Immune Repertoire Over Amplification Removal)
iROAR is a python toolkit for multiplex PCR-specific bias evaluation and correction in BCR and TCR repertoires.

[The use of non-functional clonotypes as a natural spike-in for multiplex PCR bias correction in immune receptor repertoire profiling](https://www.biorxiv.org/content/10.1101/2021.03.24.436794v1.full): biorxiv, (2020), Alexander Komkov et al

## Table of content
* [Installation](#installation)
  *  [Recommended system configuration](#recommended-system-configuration)
* [Input and output formats](#input-and-output-formats)
* [OAR statistics file](#oar-statistics-file)
* [iROAR report](#iroar-report)
* [Examples](#examples)
* [Usage](#usage)
  * [Count](#count)
  * [Merge](#merge)
  * [Filter](#filter)
  * [Bias](#bias)
* [License](#license)

## Installation
1. Install the proper python3 version and all listed packages
2. Download iROAR repository
  * via <em>Code > Download ZIP</em> and unpack the archieve, OR
  * ```git clone https://github.com/smiranast/iROAR.git ```
4. Add iROAR directory to PATH in `.bash_profile`
5. Add the runfile execution permission
```shell
chmod +x .../iROAR/iroar
```
4. On the first run germline V/D/J segments of will be downloaded if aux/germline_* files are absent


## Recommended system configuration
* Linux operation system (not tested on macOS), 2 CPU, 8GB RAM
* python 3.7.3
* matplotlib 3.0.3
* numpy 1.16.2
* pandas 0.24.2
* requests 2.21.0
* tqdm 4.43.0
* scipy  1.3.1

## Input and output formats
iROAR uses TCR and BCR clonotypes tables in the default [VDJtools](https://github.com/mikessh/vdjtools) format. 

As an input it accepts full directory path for **Count** and **Merge** commands and full path of separate VDJtools table for **Merge** and **Filter** commands.
In the former case iROAR checks the presense of appropriate RepSeq tables in the given directory.

After each table processing, in output files the tool overwrites <em>count</em> and <em>freq</em> columns with new values (by default) or adds new columns
<em>countAdj</em>, <em>freqAdj</em> and <em>freqAdjChain</em> (with `--long` argument)

## OAR statistics file
Running **Count** command with `-wj` allows to export OAR statistics, used for the recalculation of the input tables, into <em>.json</em> file. These files are named as `iROAR_run_XXXXXXXXXXXXXX.V_OAR_stat.json` and `iROAR_run_XXXXXXXXXXXXXX.J_OAR_stat.json` (the latter is created only if `--method` is set to ` vjmplex`), where XXXXXXXXXXXXXX is date and time of the run in format yyy/mm/dd/hh/mm/ss.

### The structure of OAR_stat.json files:
```
{
  gene1: {
          "OAR": OAR of this gene used for the count adjustment (<int>),
          "out-of-frame": Number of out-of-frame clones used for the calculation (<int>)
          },
   gene2: {
          ...
          },
   ...
  }
```

## iROAR report
The detailed report is saved in `iROAR_reports/` directory as `report_iROAR_run_XXXXXXXXXXXXXX.html` file (naming format is the same to [OAR statistics file](#oar-statistics-file)). It contains the following tables:
**Basic Information**
* Input arguments for iROAR_run_XXXXXXXXXXXXXX;
* List of input VDJtools tables
* Reads/clones ratio for out-of-frame clones
**For V/J genes**
* General statistics information about the region (for each gene):
  * <em>chain</em>
  * For each repertoire
    * <em>reads_N</em> - number of reads
    * <em>clones_N</em> - number of clones
    * <em>OAR_before_N</em> - OAR value before calculation
    * <em>OAR_after_N</em> - OAR value after calculation
  * <em>OAR</em> - OAR value used to adjust clonal count
  * <em>OAR norm</em> - normalized OAR value over the mean
  * <em>sampleN</em> - which clones were used for the calculation ("out-of-frame" or "all")
* Standard deviation of final OAR values of the region
* Standard deviation of V genes OAR for each iteration


## Examples
For illustration purposes iROAR is supported by an example clonotype table <em>test_dataset.txt</em> in <em>example/</em> folder.
It represents a TCR repertoire table extracted from 5'RACE library (ENA run <a href="https://www.ebi.ac.uk/ena/browser/view/ERR2869430">ERR2869430</a>,
see [this paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7346890/) for more details), downsampled to 100.000 reads.

Besides, the subdirectory <em>biased/</em> contains multiple versions of this repertoire with an artificial bias obtained by **iroar Bias** command:
one-side multiplex PCR (<em>vmplex_simulated</em>) and two-side multiplex PCR (<em>vjmplex_simulated</em>) 
1. Fix repertoires with VMPlex libraries, apply only one iteration</b>
```shell
iroar Count -iter 1 example/vmplex_simulated results/
```

2. Fix repertoires with VJMPlex libraries, until max absolute deviation of OAR from 1 will reach 0.05, write a report on statistics
```shell
iroar Count -m vjmplex example/vjmplex_simulated results/
```

3. The same as 1), but add recalculated counts as new columns and export OAR statisics to json file
```shell
iroar Count -iter 1 --long -wj example/vmplex_simulated results/
```

4. Adjust counts of VMPlex repertoires with OAR statistics, calculated in 3)
```shell
iroar Count --voar results/iROAR_run_XXXXXXXXXXXXXX.V_OAR_stat.json \
    --joar results/iROAR_run_XXXXXXXXXXXXXX.J_OAR_stat.json example/vmplex_simulated results/
```

5. The same as 4), but without the further iteration
```shell
iroar Count -iter 1 --voar results/iROAR_run_XXXXXXXXXXXXXX.V_OAR_stat.json \
    --joar results/iROAR_run_XXXXXXXXXXXXXX.J_OAR_stat.json example/vmplex_simulated results/
```

## Usage
### Count
`Count` command evaluates Over Aplification Rate for input VDJtools-formated tables.
```
iroar Count [options] <input> <output>
      
<input>: Input directory containing VDJtools tables
<output>: Output directory path

**Optional arguments:**
--long: Do not overwrite standard VDJtools columns instead of adding new ones (default=False)
-c <str>, --chains <str>: List of chains to analyse, sepatated by comma (default=IGH,IGK,IGL,TRA,TRB,TRD,TRG)
-z <int>, --outliers <int>: Do not include top-N clones for each gene in OAR calculation (default=1)
-zd <int>, --outliers_depth <int>: Minilal number of clones for each gene, where outliers filtration is applied (default=10)
-f ,--all_frame: Calculate OAR using all clones, not only out-of-frame
-min_outframe <int>: Minimal out-of-frame clones threshold for OAR calculation.
                   If a segment is present within less clones than the specified value. 
                   OAR is equal to 1 (if outframe = True), default=0
-filter_few <float>: Don't show clones with counts less than N before ceiling (default - show all)
-u, --upper_only: Adjust only counts which have OAR > 1
-m <str>, --method <str>: Method of sequencing library preparation. Possible values: 'vjmplex', 'vmplex' (default='vjmplex')
--voar <str>: Library of OAR stats for V region in JSON format
--joar <str>: Library of OAR stats for J region in JSON format.
            If multiple V and J libs are used the number of --voar files must be equal to --joarlib files
-r <str>, --report <str>: Save the report into file
-wj ,--writejson: Write OAR statistics in json format (only OAR)
-iter <int>: Maximal number of iteration (if owntable = True; default=infinite)
-err <float>: Maximal absolute deviation (if owntable = True; default=0.1)
-mt <float>, --momentum: OAR momentum on each step. Possible values: [0.01;1] (default=1.0)
--filter: Prefilter clonotypes by indels in CDR3 nucleotide sequences (similar to stand-alone Filter command, default=False)
-se <float>, --seq_error <float>: Probable error of sequencing (default=0.01)
-id <int>, --indel <int>: Maximal amount of indels to concidering CDR3s to be identical (default=1)
-ft <str>, --filter_type <str>: Which frame groups are compared during the filtering.
                              Possible values: 'IinO', 'OinI', 'all' (default=all)
-if, --iterative_filter: Apply itterative collision merge (default if --filter_type=all)
-v, --verbosity: Print messages to stdout, not only warnings and errors (default=False)
```
### Merge
`Merge` multiple OAR statistics files to one JSON library
```
iroar Merge [options] -f <input>|<input1>,<input2>...<inputN> -o <output>

OR

iroar Merge [options] -p <input> -o <output>

  -f <input1>,<input2>...<inputN>, --files <input1>,<input2>...<inputN>:
                              List of full paths of input files separated with space
  -p <input1>,<input2>...<inputN>, --path <input1>,<input2>...<inputN>:
                              Full path of directory containing file for merging         
  At least one type of input data list (-f or -p) must be given!

  -o <output>, --output <output>: Prefix of merged library.

**Optional arguments:**

  -min_outframe <int>: Minimal out-of-frame clones threshold for mean OAR calculation.
                               If a segment in this JSON file is present within less clones
                               than the specified value OAR is equal to 1. (default=0)
```
### Filter
`Filter` clonotypes by indel and sequencing error thresholds.
```
iroar Filter [options] -i <input> -o <output>

  -i <input>, --input <input>: Path of input VDJtools table
  -o <output>, --output <output>: Path of output VDJtools table

**Optional arguments:**    
  -se <float>, --seq_error <float>: Probable error of sequencing (default=0.01)
  -id <int>, --indel <int>: Maximal amount of indels to concidering CDR3s to be identical (default=1, type=int)
  -ft <str>, --filter_type <str>: Which frame groups should be compared during the filtering.  Possible values: 'IinO', 'OinI', 'all'.
  \                               If multiple, a list separated by comma (default=all).
  -if, --iterative_filter,  help = Apply itterative collision merge (default for --filter_type=all)
```
### Bias
`Bias` command introduces an artificial OAR bias in VDJtools repertoire table (for test purposes).
```
iroar Bias [options] <input> <output>
      
  <input>: Path of input VDJtools table
  <output>: Path of biased VDJtools table

**Optional arguments:**
  --seed <int>: Randomization seed (default=0)
  -m <str>, --method <str>: Method of sequencing library preparation. Possible values: 'vjmplex', 'vmplex' (default='vjmplex')
  --noise: Add random noise in bias for each clone (default=False)
```

## License
Copyright (c) 2020, Alexander Komkov, Anastasia Smirnova, Ilgar Mamedov
(here and after addressed as Inventors)

All Rights Reserved

<p>Permission to use, copy, modify and distribute any part of this program for
educational, research and non-profit purposes, by non-profit institutions
only, without fee, and without a written agreement is hereby granted,
provided that the above copyright notice, this paragraph and the following
three paragraphs appear in all copies.

Those desiring to incorporate this work into commercial products or use for
commercial purposes should contact Alexander Komkov using the following email address:
[alexandrkomkov@yandex.ru](mailto:alexandrkomkov@yandex.ru)

IN NO EVENT SHALL THE INVENTORS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN IF THE INVENTORS HAS BEEN
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE INVENTORS HAS
NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
MODIFICATIONS. THE INVENTORS MAKES NO REPRESENTATIONS AND EXTENDS NO
WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A
PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY
PATENT, TRADEMARK OR OTHER RIGHTS.
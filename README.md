# NucHMM: A Quantitative modeling of nucleosome organization identifies functional nucleosome states.  

## Introduction

Nucleosome organization, often described as its positioning, spacing and regularity, is the interplay among nucleosome and nucleosome-binding factors such as DNA-binding factors, histone chaperones, and ATP-dependent chromatin remodelers. To address the lacking power of determining the combinational effects of the different influencing factors on nucleosome organization, we presented **NucHMM** for identifying functional nucleosome states. NucHMM integrates a hidden Markov model (HMM), and estimated nucleosome regularity, spacing as well as positioning information, to identify nucleosome states associated with cell type-specific combinatorial histone marks.

## Recent Changes for NucHMM (version 1.0)

*N/A*

## Workflow
<img src="https://github.com/KunFang93/NucHMM/blob/master/workflow/NucHMM_workflow.png" width="900">

## Installation

The program is tested on *Linux 3.10.0/CentOS 7*.

The installation currently requires the [**miniconda**](https://docs.conda.io/en/latest/miniconda.html). 

**Install miniconda:**

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

**Download NucHMM:**

```
git clone https://github.com/KunFang93/NucHMM.git
```


**Create the conda environment with dependents and Activate the environment:**

```
conda env create -n nuchmm -f <path to NucHMM>/scripts/env.yaml
conda activate nuchmm
```

**Install NucHMM:**

```
cd <path to the NucHMM>/scripts
pip install --editable .
cd NucHMM_Cplus/
make
```

## Quick Start
*__Note:__ The context in <> should be replace by user according to real data*

### **Step1: ChIP-seq peaks and MNase-seq nucleosome locations files preparation.**  
 

__Skip__ this step if you alreay have the histone marks' peak files and iNPS derived nucleosome files (remove the header, and please keep the all like_wig files). You can check those files' format in example_files folder. We __recommend__ use your own favored pipeline.  

If the input all files in fastq format
```
# not include MNase-seq fq
NucHMM nuchmm-prep --fastq -p 20 -ifql fqs_list.txt -qc -bip <Full Path to bowtie index>/<idx-basename>

# include MNase-seq fq
NucHMM nuchmm-prep --fastq -p 20 -ifql fqs_list.txt -qc -bip <Full Path to bowtie index>/<idx-basename> -inps <Full Path to NucHMM>/scripts/iNPS_V1.2.2.py
```
If the input all files in bam format
```
# not include MNase-seq bam
NucHMM nuchmm-prep --bam -p 20 -ibl bams_list.txt

# include MNase-seq bam
NucHMM nuchmm-prep --bam -p 20 -ibl bams_list.txt -inps <Full Path to NucHMM>/scripts/iNPS_V1.2.2.py
```
The output peak files will locate at **peakcalling_result/ folder**, and the output nucleosome location files will locate at **nuc_calling_result/ folder**.

### **Step2: NucHMM initialization.** 


First manually create 
```
# Check the those files format in exmaple_files folder.
1. <celltype>_histone_marks.txt that contains all histone mark peak files;
2. histonelists_list.txt contains all <celltype>_histone_marks.txt file;
3. nucposfiles_list.txt contains all <celltype>_nucleosome_locations.bed.
```

Then, use command
```
NucHMM nuchmm-init -iplf histonelists_list.txt -nucf nucposfile_list.txt -gf <Full Path to NucHMM>/annotation/genebody_anno_hg19.txt -rmf
```

The default output file will be named <celltype>_<# histone mark>.precomp in the current directory. User can specify the ouptut name by -ofl parameter.  
  
### **Step3: NucHMM training.**


First manually create.  
```
precompfiles_list.txt that contains all precomp files result from nuchmm-init. 
(check the file format in example_files/histone_marks.txt)
```

Then, 

```
NucHMM --hmm-directory <Full Path to NucHMM>/scripts/NucHMM_Cplus/bin/ nuchmm-train -refg <Full Path to NucHMM>/annotation/hg19.chrom.sizes.txt -pl precompfile_list.txt -numh <number of histone marks. e.g. 8> -rmf
```

Note that training process could take few hours up to a day depends on how much data you input. We normally use nohup command to run the this nuchmm-train command in the background by

```
nohup NucHMM --hmm-directory <Full Path to NucHMM>/scripts/NucHMM_Cplus/bin/ nuchmm-train -refg <Full Path to NucHMM>/annotation/hg19.chrom.sizes.txt -pl precompfile_list.txt -numh <number of histone marks. e.g. 8> -rmf &> train.log &
```

The default output file will be name HMM_<# histone marks>.rawhmm in the current directory. User can specify the output name by -ohmm parameter.

### **Step4: Matrix visualization and background state detection.**


Manually create 
```
histone_marks.txt file that contains all histone marks. 
(check the file format in example_files/histone_marks.txt)  
```

Then, visualize the Transition and Mark-state matrix first.
```
NucHMM matrix-visualize -rhf HMM_<# histone marks>.rawhmm -hlf histone_marks.txt
```

The default output files will be two png files. 

```
View the Mark_state matrix and define no histone mark HMM states as background state. 
(e.g. in example_files/Mark_state.png, 8 and 12 is background states)
```

### **Step5: NucHMM screen initialization.**


Note the order of -ct parameter should be same with the cell type order in histonelists_list.txt
```
NucHMM nuchmm-screen-init -rhf HMM_<# histone marks>.rawhmm -hlf histone_marks.txt \
-bg <background state1> -bg <background state2> .. -bg <background stateN> \
-gf <Full Path to NucHMM>/annotation/genebody_anno_hg19.txt \
-ct <celltype1> -ct <celltype2> .. -ct <celltypeN> -ptm -rmf
```

The default output files:   
  Total_distribution_up.png, Total_distribution_gbdown.png, Total_distribution.png  
  nuchmm_screen_init_result_files.txt  
  states_genomic_location.txt.

### **Step6: NucHMM screen.**

Firstly,
```
Double-Check states_genomic_location.txt and Total_distribution.png.
Make sure the predicted genomic location in states_genomic_location.txt fits your empirical knowledge.
If not, you can manually edit the states_genomic_location.txt file.
```

Then,
```
NucHMM nuchmm-screen -gf <Full Path to NucHMM>/annotation/genebody_anno_hg19.txt \
-lwfl like_wigs_list.txt -nucf nucposfiles_list.txt -bg <background state1> .. -bg <background stateN> \
-sn <Total number of states, e.g. Mark_state.png show total 13 states> -pm -rmf
```

The default output files:
```
<celltypes>_gl_an_resp_pos_filt.bed
functional_nucleosome_state_post.txt
Nuc_pos_vio.png, Regularity_score_filt.png, State_Spectral_density_filt.png, Array_distribution_filt.png
```

The final function nucleosome states features are store in functional_nucleosome_state_post.txt.  
The final kept nucleosomes are in < celltypes >_gl_an_resp_pos_filt.bed.   

Enjoy :blush:!

## Command Line

* ### Usage

```
NucHMM [OPTIONS] COMMAND [ARGS]

Options:
  --version               Show the version and exit
  --hmm-directory PATH    the path of the NucHMM_Cplus/bin folder
  --help                  Show this message and exit

Commands:
  nuchmm-prep             Prepare ChIP-seq peaks files or nucleosome locations files from fastq/bam files.
  nuchmm-init             Assign histone marks to nucleosomes and create precomp bins for nuchmm-train.
  nuchmm-train            Use Hidden Markov Model(HMM) and Viterbi Algorithm to decode HMM states for each nucleosomes.
  nuchmm-screen-init      Initializae the screen step by providing sorted state files and the suggested genomic locations of HMM states.
  nuchmm-screen           Filter nucleosomes by genomic location, array number, nucleosome regularity, spacing and positioning.
  matrix-visualize        Visualize the Transition and Mark-state matrix.
```
 * ### Options
   * `--version`:  
   Show NucHMM current version and exit.
   
   * `--hmm-directory PATH`:  
   Used in nuchmm-train for input the path of the NucHMM_Cplus/bin folder.
   
   * `--help`:  
   Show help message.
   
 * ### Commands
 
   * #### nuchmm-prep
   
     NucHMM provides basic ChIP-seq and MNase-seq pipeline to handle the fastq/bam files. If you have other favored ChIP-seq pipeline, we recommend use your favored way to process the raw fastq or bam files. However, for MNase-seq, we only accept result from iNPS currently.
     * `--fastq/--bam`:    
     Indicate the input data type
     
     * `--inputfqslist(-ifql)`:  
     Used with --fastq command. Input the fastqs list file. Input format: 'FILE READ-LEN SEQ-TYPE(chip/mnase) PEAK-TYPE(narrow?broad/none)'. Check fqs_list.txt in example_files folder for the detailed format.
     
     * `--inputbamslist(-ibl)`:  
     Used with --bam command. Input the bams list. Input format: for each line 'FILE PE?SE SEQ-TYPE(chip/mnase) PEAK-TYPE(narrow?broad/none)'. Check bams_list.txt in example_files folder for the detailed format.
     
     * `--qualitycontrol(-qc)`:  
     Used with --fastq command. Flag for whether run QC trim step for fastq files.
     
     * `--bowtieindexpath(-bip)`:  
     Required if input format is fastq, the path and basename of bowtie index (bowtie -x)
     
     * `--bowtie2indexpath(-b2ip)`:  
     Required if input format is fastq, the path and basename of bowtie2 index (bowtie2 -x)
     
     * `--inpspath(-inps)`:  
     Required if input has MNase-seq, the path of iNPS.py. For example, /data/NucHMM/scripts/iNPS_V1.2.2.py.
     
     * `--threads(-p)`:  
     Number of threads
     
   * #### nuchmm-init
  
     Assign histone marks to nucleosomes and create precomp bins for nuchmm-train. All .precomp files name are writing to Precompfiles_list.txt.
   
     * `--inputpeakslistfiles(-iplf)`[**Required**]:  
     Input TF-peaks_files list. Check the histonelists_list.txt in example_files folder for the detailed format.
     
     * `--nucpositionfiles(-nucf)`[**Required**]:  
     Input the nucleosome position files list. Check the nucposfiles_list.txt in example_files folder for the detailed format.
     
     * `--genefile(-gf)`[**Required**]:  
     Input the interested genes list file. Check the genebody_anno_hg19.txt in annotation_files for the detailed format.
     
     * `--outputfilelist(-ofl)`:  
     Specify the output files path and name.   
     Format:   
     < Path >/< celltype1 >.precomp   
     < Path >/< celltype2 >.precomp   
     ..  
     < Path >/< celltypeN >.precomp
     
     * `--refgenome(-refg)`:  
     The reference genome file contains the length of each chromosomes. Default: built-in hg19.
     
     * `--intersect_cutoff(-ic)`:  
     The intersect threshold to assign histone mark peaks to the nucleosomes. Default is 0.3. For example, for a 150 bp nucleosome core region, if a histone mark peak have 45 (150 * 0.3) bp overlap with this nucleosome. Then the program will consider this nucleosome has this histone mark.
     
     * `--gap(-g)`:  
     Flag of giving a mark to inter-nucleosome region (gaps between nucleosomes). 
     
     * `--upboundary(-up)`:  
     Upstream boundary of TSS for selecting the training region. Default:100000. The larger the number, the more computational resources (memory) and training time needed.
     
     * `--downboundary(-down)`:    
     Downstream boundary of TTS for selecting the training region. Default:10000. The larger the number, the more computational resources (memory) and training time needed.
     
     * `--removetmpfile(-rmf)`: 
     Flag of removing all temporary files.
     
   * #### nuchmm-train
   
     Use Hidden Markov Model(HMM) and Viterbi Algorithm to decode HMM states for each nucleosomes. 
     
     * `--refgenome(-refg)`[**Required**]:  
     Input the reference genome file that contains the length of each chromosomes. Check the hg19.chrom.sizes.txt in annotation_files folder.
     
     * `--precomp_list(-pl)`[**Required**]:  
     Input precomp files list resulted from nuchmm-init. Check the precompfiles_list.txt in example_files folder.
     
     * `--numhistonemarks(-numh)`[**Required**]:  
     The number of histone marks used. For example, if user use H3K4me1/K4me3/K27ac in nuchmm-init, then input 3. 
     
     * `--numStates(-nums)`:  
     The number of states for first round hmm training. Default: 20.
     
     * `-b`:  
     Calculate and report BIC after every iteration. WARNING: This parameter will severely slow the algorithm!
     
     * `--bicfile(-obic)`:  
     Specify the path and name of the BIC score file. If -b not specified, only report BIC score of the first and the last iteration.
     
     * `-i`:  
     Maximum number of iterations for first round hmm training. Default: 300.
     
     * `--num/--perc`:  
     The crietrion to remove the redundant states in the first round hmm training. `--num` uses the specific number as cutoff while `--perc` uses a specific percentage as cutoff. Default use is --perc (0.005). If specify --num, then default number is 350. For example, --num -nn 300 function as: if a HMM state has less than 300 nucleosomes in the first round hmm training, then this HMM state will be removed from the tranisition and emission matrix which are used as the initial matrix in the second round HMM. --perc -pn 0.01 function as: if a HMM state has less than 0.01% of the total number of nucleosomes, then this HMM will be removed from the tranisition and emission matrix which are used as the initial matrix in the second round HMM.
     
     * `--nucnum(-nn)`:  
     Used with --num to specify the number cutoff.
     
     * `--percentage(-pn)`:
     Used with --perc to specify the percentage cutoff.
     
     * `--outrawhmmfile(-ohmm)`:  
     Specify the path and name of the output rawhmm file.
     
     * `--outstats(-os)`:  
     Output the nucleosome counts for each HMM states in each cell type.
     
     * `--removetmpfile(-rmf)`:  
     Remove temporary files.
   
   * #### nuchmm-screen-init
     
     Initializae the screen step by providing sorted state files and the suggested genomic locations of HMM states.
     
     * `--rawhmmfile(-rhf)`[**Required**]:  
     Input the secondr.rawhmm file resulted from nuchmm-train.
     
     * `--histonelistfile(-hlf)`[**Required**]:  
     Input histone_marks.txt file that contains all histone marks used in nuchmm-init. Check the histone_marks.txt for the detail format.
     
     * `--genesfile(-gf)`[**Required**]:  
     Input the interested genes list file. Check the genebody_anno_hg19.txt in annotation_files for the detailed format.
     
     * `--celltypes(-ct)`[**Required**]:   
     Input all cell types name used in nuchmm-init. For example, if user train model with data from MCF7, H1 and IMR90, user should input -ct MCF7 -ct H1 -ct IMR90. The order of the -ct should be the same with the celltype order in nucposfiles_list.txt.
     
     * `--bgstate(-bg)`[**Required if exists background states**]:  
     Specify background states that identified from the Mark-state matrix. If there are multiple background states, for example HMM state1 and HMM state2 are background states, then use -bg 1 -bg 2. 
     
     * `--statesfilelist(-sfl)`:    
     Specify the path and name of the output < celltype > _states_secondr.bed. The order of the < celltype > should be same with -ct command. The command will detect the states_secondr.bed in current directory automatically. If user want to run command in the background, this parameter is required.  
     format:  
     < celltype1 > _states_secondr.bed  
     ..  
     < celltypeN > _states_secondr.bed  
     
     * `--outputsfilelist(-ofl)`:  
     Specify the path and name of the output < celltype > _output_secondr.bed. The order of the < celltype > should be same with -ct command. The command will detect the output_secondr.bed in current directory automatically. If user want to run command in the background, this parameter is required.  
     format:  
     < celltype1 > _output_secondr.bed  
     ..  
     < celltypeN > _output_secondr.bed  
     
     * `--plotcellmark(-pcm)`:    
     Plot the cell type specific HMM states distribution. 
     
     * `--plottotalmark(-ptm)`:  
     Plot the all cell averaged HMM states distribution. 
     
     * `--overwrite(-ow)`:  
     Overwrite all existing files. Useful when user need to re-run command in the same directory.
     
     * `--winmethod/--refmethod`:  
     Pick a method to identify unique HMM state for the nucleosomes that have multi-states. `--winmethod` uses window method to identify the unique state; `--refmethod` uses output_secondr.bed as reference to identify the unique state. Default method: `--refmethod`. See the *pick_one* function in NucHMM_Feature_screen for the detailed inplementation.
     
     * `--winsize(-ws)`:  
     Used with `--winmethod`, Default value: 5. Window method use the adjcent nucleosome states in the selected window to identify the unique state for the target nucleosome. For example, if a nucleosome has both HMM state1(associated with H3K4me1) and HMM state2(associated with H3K4me3) labelled, we then select the upstream < 5 > nucleosomes and downstream < 5 > nucleosomes of this target nucleosome to find which HMM state is most likely appeared under this environment. See the *pick_one* function in NucHMM_Feature_screen for the detailed inplementation.
     
     * `--updistal(-uD)`:  
     Specify upstream boundary for plotting the distribution, cannot exceed than `--upboundary` used in nuchmm-init. Default: 100000.
     
     * `--uppromoter(-uProm)`:  
     Specify upstream boundary of the promoter region. Used in identifying the genomic location of HMM states. Default: 1000.
     
     * `--downbound(-db)`: 
     Specify downstream boundary for plotting the distribution, cannot exceed than `--downboundary` used in nuchmm-init. Default: 10000.
     
     * `--rescalelength(-rl)`:  
     Specify the rescaled gene body length. Default: 30000.
     
     * `--refhg19/--refhg38`:  
     Specify the reference genome. Default: built-in hg19.
     
     * `--outfile(-of)`:  
     Specify the path and name of the identified genomic location file.
     
     * `--removetmpfile(-rmf)`:  
     Remove temporary files.
     
   * #### nuchmm-screen
   
   * #### matrix-visualize





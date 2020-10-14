#!/usr/bin/env python
#--coding:utf-8 --


#################################################################################################
#################################################################################################
########                                                                                 ########
########    Integrated HMM Methods to identify functional nucleosome states              ########
########                                                                                 ########
########    Author:  Kun Fang                                                            ########
########                                                                                 ########
########                                                                                 ########
########    Working Environment:  Python3                                                ########
########                                                                                 ########
########    Date:      2020-10-13                                                        ########
########                                                                                 ########
########                                                                                 ########
#################################################################################################
#################################################################################################


import re
import sys
import click
import subprocess
import numpy as np
import pandas as pd
import multiprocessing
from collections import defaultdict
from NucHMM_prep import QC_step, Mapping_step, Peak_Nuc_calling_step
from NucHMM_precompile import transcript_factor_assign, segment_gene, precompile_step
from NucHMM_hmm import hmm, hmm_calling, modifyhmm, hmm_second
from NucHMM_states_coverage import create_state_coverage_files
from NucHMM_Feature_screen import select_unique_state, genomic_loc_finder, genomic_loc_filter, array_num_filter
from NucHMM_Feature_screen import nuc_positioning_filter,regularity_spacing_filter,statem_to_histm
from NucHMM_utilities import get_time, filter_info_print, count_file_rows, file_check, query_yes_no, input_new_name, which
from NucHMM_common_data import hg19,hg38,val2chr,info2strand,spe_colors
from NucHMM_Visualization import plot_state_coverage,HMM_matrix_visualization,plot_violin_nuc_pos
from NucHMM_Load_Write_files import load_genome_size_file, write_binnum_file, write_precomp_file,load_assignedtf_file, load_fq_list
from NucHMM_Load_Write_files import load_histonefile, load_nuc_position_file, load_rawhmm, load_fq_list, load_bam_list, write_prep_resultlist


# class for setting up subcommand order.
# Inspired by https://stackoverflow.com/questions/47972638/how-can-i-define-the-order-of-click-sub-commands-in-help
class SpecialHelpOrder(click.Group):

    def __init__(self, *args, **kwargs):
        self.help_priorities = {}
        super(SpecialHelpOrder, self).__init__(*args, **kwargs)

    def get_help(self, ctx):
        self.list_commands = self.list_commands_for_help
        return super(SpecialHelpOrder, self).get_help(ctx)

    def list_commands_for_help(self, ctx):
        """reorder the list of commands when listing the help"""
        commands = super(SpecialHelpOrder, self).list_commands(ctx)
        return (c[1] for c in sorted(
            (self.help_priorities.get(command, 1), command)
            for command in commands))

    def command(self, *args, **kwargs):
        """Behaves the same as `click.Group.command()` except capture
        a priority for listing command names in help.
        """
        help_priority = kwargs.pop('help_priority', 1)
        help_priorities = self.help_priorities

        def decorator(f):
            cmd = super(SpecialHelpOrder, self).command(*args, **kwargs)(f)
            help_priorities[cmd.name] = help_priority
            return cmd

        return decorator

required_programs = ["bedops","bedtools","deeptools","macs2","epic2","fastqc","samtools","trim_galore","bowtie","bowtie2"]
for prog in required_programs:
    p = which(prog)
    # print("Checking if %s is available" %(prog))
    if p is None:
        print(" - Cannot find %s! Please install it Or Some functions will be impaired!" %(prog))
        if not query_yes_no("Keep running without installation?"):
            exit(1)
        else:
            continue
    else:
        continue
        #print("  - %s is available" %(prog))


class Config(object):

    def __init__(self):
        self.verbose = False

pass_config = click.make_pass_decorator(Config, ensure=True)
@click.group(cls=SpecialHelpOrder)
@click.version_option(version=1.1)
@click.option('--hmm-directory', type=click.Path(),help='the path of NucHMM_Cplus bin')
@pass_config
def cli(config, hmm_directory):
    if hmm_directory is None:
        hmm_directory = '.'
    config.hmm_directory = hmm_directory

@cli.command(help_priority=1)
@click.option('--fastq/--bam', default=False, required = True, help = 'Indicate input data type')
@click.option('--inputfqslist', '-ifql', type=click.Path(exists=True),
              help='Input fastq files list, used with --fastq command. '
                   'Input format: for each line \'FILE READ-LEN SEQ-TYPE(chip/mnase) PEAK-TYPE(narrow?broad/none)\'. '
                   'For exmaple, if inputfile is paired-end narrowpeak ChIP-seq data, the format should be \'F_R1.fq F_R2.fq Read-length chip narrow\'. '
                   'If inputfile is single-end broad ChIP-seq data, the format should be \'F.fq Read-length chip broad\'.'
                   'If inputfile is single-end MNase-seq data, the format should be \'F.fq Read-length mnase none\'.'
                   'If inputfile is paired-end MNase-seq data, the format should be \'F_R1.fq F_R2.fq Read-length mnase none\'.')
@click.option('--inputbamslist','-ibl', type=click.Path(exists=True),
              help='Input bam files list, used with --bam command. '
                   'Input format: for each line \'FILE PE?SE SEQ-TYPE(chip/mnase) PEAK-TYPE(narrow?broad/none)\'. '
                   'For example, if inputfile is paired-end narrowpeak ChIP-seq data, the format should be \'F.bam PE chip narrow\'.'
                   'If inputfile is single-end MNase-seq bam file, the format should be \'F.bam SE mnase none\'.')
@click.option('--qualitycontrol', '-qc', is_flag=True, help = 'whether run QC trim step for fastq file, use with --fastq.')
@click.option('--bowtieindexpath', '-bip', type=click.Path(), help = 'Required if input format is fastq, the path and basename of bowtie index (bowtie -x)')
@click.option('--bowtie2indexpath', '-b2ip', type=click.Path(), help = 'Required if input format is fastq, the path and basename of bowtie2 index (bowtie2 -x)')
@click.option('--inpspath','-inps',type=click.Path(), help = 'Required if input has MNase-seq, the path of iNPS.py. '
                                                             'For example, /data/NucHMM/scripts/iNPS_V1.2.2.py.')
@click.option('--threads','-p',default=5, help = 'Number of threads')
def NucHMM_prep(fastq,inputfqslist,inputbamslist,fastqc,bowtieindexpath,bowtie2indexpath,inpspath,threads):
    '''Prepare files for nuchmm-init from fastq/bam files.
    NucHMM only provide basic fastq and peak calling pipeline. We recommend user use their favorite pipeline to process
    their fastq files and peaking the peaks. But we only accept iNPS result of MNase-seq data currently '''
    # MAPQ cutoff
    mapq_cutoff = 20
    if fastq:
        # elements in read_len_list already in int type
        fq_info_dict = load_fq_list(inputfqslist)
        if fastqc:
            # add filt_file_list into fq_info_dict
            fq_info_dict = QC_step(fq_info_dict,threads)
            # add aligned_file_list into fq_info_dict
            fq_info_dict = Mapping_step(fq_info_dict,bowtieindexpath,bowtie2indexpath,threads,mapq_cutoff)
        else:
            fq_info_dict = Mapping_step(fq_info_dict,bowtieindexpath,bowtie2indexpath,threads,mapq_cutoff)

        fq_info_dict = Peak_Nuc_calling_step(fq_info_dict,inpspath)
        write_prep_resultlist(fq_info_dict, 'NucHMM_prep_result_files' + get_time() + '.txt')
    else:
        bam_info_dict = load_bam_list(inputbamslist)
        bam_info_dict = Peak_Nuc_calling_step(bam_info_dict,inpspath)
        write_prep_resultlist(bam_info_dict, 'NucHMM_prep_result_files' + get_time() + '.txt')


@cli.command(help_priority=2)
@click.option('--inputpeakslistfiles', '-iplf',type=click.Path(exists=True), required = True,
              help='Input TF-peaks_files list')
@click.option('--nucpositionfiles', '-nucf', type=click.Path(exists=True), required=True,
              help='The list of nucleosome position files')
@click.option('--genefile', '-gf', type=click.Path(exists=True), required=True, help='Interested genes\' file')
@click.option('--outputfilelist', '-ofl', type=click.Path(), help='Specify the output files path and name.' )
@click.option('--refgenome', '-refg', type=click.Path(), help='The reference genome file contains the length of each chromosomes.')
@click.option('--intersect_cutoff', '-ic',default=0.3, help='The intersect threshold to assign histone mark peaks to the nucleosomes.')
@click.option('--gap', '-g', is_flag=True, help='Flag of giving a mark to inter-nucleosome region (gaps between nucleosomes)')
@click.option('--upboundary', '-up',default=100000, help='Upstream boundary of TSS for selecting the training region.')
@click.option('--downboundary', '-down',default=10000, help='Downstream boundary of TTS for selecting the training region.')
@click.option('--removetmpfile','-rmf',is_flag=True,help='Flag of removing all temporary files.')
def NucHMM_init(inputpeakslistfiles, nucpositionfiles, intersect_cutoff, gap, genefile, upboundary, downboundary,
                outputfilelist, refgenome, removetmpfile):
    '''Assign histone marks to nucleosomes and create precomp bins for nuchmm-train. All .precomp files name
    are writing to Precomps_list.txt'''
    print("\nNucHMM init start!")
    if refgenome is None:
        refgene = hg19
    else:
        refgene = load_genome_size_file(refgenome)

    # notice the celltypes should have the same order with peaksfile_list and nucposfile list
    celltype_list = []
    peaksfile_list = load_histonefile(inputpeakslistfiles)
    nucposfile_list = load_histonefile(nucpositionfiles)
    outputfile_list = []

    # create histone_marks.txt
    print("Createing histone_marks.txt")
    if file_check('histone_marks.txt'):
        print("histone_marks.txt exists, skip and please manually check it.")
    else:
        with open('histone_marks','w') as his_file:
            marks_list = peaksfile_list[0]
            for mark_file in marks_list:
                # if user following the standard naming rule, the histone marks should behind celltype
                mark_name = mark_file.split('/')[-1].split('_')[1]
                if str.lower(mark_name[0])=='h':
                    his_file.write(mark_name+'\n')
                else:
                    print("Not standard file name, need manually create histone_marks.txt")
        his_file.close()

    # assign peaks mark to nucleosome
    peaks_num = 0
    for idx,nucposfile in enumerate(nucposfile_list):
        celltype = nucposfile.split('/')[-1].split('_')[0]
        celltype_list.append(celltype)
        print("\nAssign Peaks for %s" % celltype)
        inputpeakslist = peaksfile_list[idx]
        nucposition = nucposfile
        peaks_list = load_histonefile(inputpeakslist)
        peaks_num = len(peaks_list)
        outputfile = celltype + '_' + str(peaks_num) + '_assign.bed'
        outputfile_list.append(outputfile)
        if file_check(outputfile):
            query_mark = query_yes_no(outputfile + ' exists, do you want to overwrite it?')
            if query_mark:
                transcript_factor_assign(intersect_cutoff, gap, inputpeakslist, outputfile, nucposition)
            else:
                print("Use the exist %s." % outputfile)
        else:
            transcript_factor_assign(intersect_cutoff, gap, inputpeakslist, outputfile, nucposition)

    # create gene segement file for later precompile step
    print("\nCreating gene segment file..")
    segoutfile = 'gene_seg_' + str(upboundary//1000) + 'k_' + str(downboundary//1000) + 'k.txt'
    if file_check(segoutfile):
        query_mark2 = query_yes_no(segoutfile + ' exists, do you want to overwrite it?')
        if query_mark2:
            segment_gene(genefile, segoutfile, upboundary, downboundary, refgene)
        else:
            print("Use the exist %s." % segoutfile)
    else:
        segment_gene(genefile, segoutfile, upboundary, downboundary, refgene)

    # precompile step
    print("\nPrecompiling the files")
    assignedfilelist = "Assigned_files.txt"
    f = open(assignedfilelist,'w')
    for file in outputfile_list:
        f.write(file+'\n')
    f.close()

    if outputfilelist is None:
        prec_outputfile_list = 'Precompfiles_list.txt'
        f = open(prec_outputfile_list,'w')
        for celltype in celltype_list:
            f.write(celltype + '_' +str(peaks_num) + '.precomp'+'\n')
        f.close()
    else:
        prec_outputfile_list = outputfilelist

    prec_skip_mark = True
    for file in load_histonefile(prec_outputfile_list):
        prec_skip_mark = prec_skip_mark and file_check(file)

    if prec_skip_mark:
        query_mark3 = query_yes_no('Precomp files exist, do you want to overwrite them?')
        if query_mark3:
            precompile_step(assignedfilelist, segoutfile, prec_outputfile_list, refgene)
        else:
            print("Use the existed .precomp files.")
    else:
        precompile_step(assignedfilelist, segoutfile, prec_outputfile_list, refgene)

    print("\nWriting precomp files list to %s." % prec_outputfile_list)

    if removetmpfile:
        print("\nRemoving temporary files..")
        for file in outputfile_list:
            subprocess.call("rm " + file, shell=True)
        subprocess.call("rm " + segoutfile, shell=True)
        subprocess.call("rm " + assignedfilelist, shell=True)

    print('\nNucHMM init Finish!')

@cli.command(help_priority=3)
@click.option('--refgenome', '-refg', type=click.Path(exists=True), required=True, help='Input the reference genome file '
                                                                                        'that contains the length of each chromosomes.')
@click.option('--precomp_list', '-pl', type=click.Path(exists=True), required=True, help='Input precomp files list resulted from nuchmm-init.')
@click.option('numhist','--numhistonemarks', '-numh',type=int, required=True,help='The number of histone marks used in nuchmm-init')
@click.option('numstates', '--numStates', '-nums',default=20, help='The number of states for first round hmm training. Default: 20.')
@click.option('-b', is_flag=True,
              help='Calculate and report BIC after every iteration. WARNING: This parameter'
                   'will slow the algorithm a lot!')
@click.option('--bicfile', '-obic', type=click.Path(),help='BIC score file, if -b not specified, only report first and '
                                                           'last BIC score')
@click.option('--initmatrix', '-im', type=click.Path(exists=True), help='Use a pre-trained matrix as initial matrix.')
@click.option('-i', default=300, help='Maximum number of iterations for first round hmm training')

@click.option('--num/--perc',default=False,help='Two methods to remove the redundant states, recommended method is '
                                               'percentage method')
@click.option('--nucnum', '-nn', default=350, help='The minimum number of nucleosomes to keep the states')
@click.option('--percentage', '-pc', default=0.005, help='The minimum percentage of nucleosomes to keep the states,'
                                                 'default 0.005 (0.5%)')
@click.option('--outrawhmmfile', '-ohmm', type=click.Path(), help='Output file\'s path and name')
@click.option('--outstats', '-os', type=click.Path(), help='The nucleosomes counts for each states in each cell type')
@click.option('--removetmpfile','-rmf',is_flag=True,help='Remove temporary files')
@pass_config
def NucHMM_train(config, numhist, refgenome, b, i, numstates, precomp_list, num, nucnum, percentage, outstats,
                    outrawhmmfile, bicfile, initmatrix, removetmpfile):
    '''Use Hidden Markov Model and Viterbi algorithm to decode HMM state for each nucleosomes'''

    # check if NucHMM-learn etc in directory
    command_list = ['NucHMM-learn', 'NucHMM-output_results']
    for command in command_list:
        file_abs_path = config.hmm_directory + '/' + command
        if file_check(file_abs_path):
            continue
        else:
            print(file_abs_path + " not exists!")
            exit(1)

    # Read precompiles list file to a list
    precompilefiles_list = load_histonefile(precomp_list)

    # create chromosome-bin file
    binoutfile = 'binnum.txt'
    if file_check(binoutfile):
        query_mark0 = query_yes_no(binoutfile + ' exists, do you want to overwrite it?')
        if query_mark0:
            write_binnum_file(precompilefiles_list[0],binoutfile)
        else:
            print("Use exist %s." % binoutfile)
    else:
        print("Creating %s." % binoutfile)
        write_binnum_file(precompilefiles_list[0],binoutfile)

    # if there is a initial matrix
    if initmatrix is not None:
        trans_matrix_modi,emit_matrix_modi,numstate_modi,numoutput_modi = load_rawhmm(initmatrix)
        if outrawhmmfile is None:
            out_rawhmm = 'HMM_'+str(numstate_modi)+'.rawhmm'
        else:
            out_rawhmm = outrawhmmfile

        if bicfile is None:
            bic_file = 'BIC' + get_time() +'.txt'
        else:
            bic_file = bicfile

        hmm_second(config, refgenome, out_rawhmm, b, 200, numstate_modi, numoutput_modi, precompilefiles_list, bic_file,
                   initmatrix, binoutfile)

        # HMM calling
        print("HMM Viterbi calling")
        statefiles_list_second = []
        for precompfile in precompilefiles_list:
            celltype = precompfile.split('/')[-1].split('_')[0]
            print("Viterbi calling state and output file for %s" % celltype)
            state_outfile = celltype + '_states_i_secondr.bed'
            statefiles_list_second.append(state_outfile)
            output_outfile = celltype + '_output_i_secondr.bed'
            if file_check(state_outfile) and file_check(output_outfile):
                query_mark2 = query_yes_no("HMM model Viterbi calling %s file exist, do you want to overwrite it?"%celltype)
                if query_mark2:
                    hmm_calling(config, refgenome, state_outfile, output_outfile, precompfile, out_rawhmm, binoutfile)
                else:
                    print("Use the exists calling files")
            else:
                hmm_calling(config, refgenome, state_outfile, output_outfile, precompfile, out_rawhmm, binoutfile)
    else:
        # number of outputs is 2 ** number of histone modifications
        numoutputs = 2 ** numhist

        # First round hmm training
        rawhmm_file = 'HMM_' + str(numhist) +'.rawhmm'
        bic_file_1 = 'bicfile_1' + get_time() + '.txt'
        if file_check(rawhmm_file):
            query_mark = query_yes_no(rawhmm_file + ' exists, do you want to overwrite it?')
            if query_mark:
                print("First round HMM training..")
                hmm(config, refgenome, rawhmm_file, b, i, numstates, numoutputs, precompilefiles_list, bic_file_1, binoutfile)
            else:
                if query_yes_no('Do your want to write to a new name?'):
                    rawhmm_file = input_new_name('New name:')
                    print("First round HMM training..")
                    hmm(config, refgenome, rawhmm_file, b, i, numstates, numoutputs, precompilefiles_list, bic_file_1, binoutfile)
                else:
                    print("Use the exists %s." % rawhmm_file)
        else:
            print("First round HMM training..")
            hmm(config, refgenome, rawhmm_file, b, i, numstates, numoutputs, precompilefiles_list, bic_file_1, binoutfile)


        # modify hmm according the number of nucleosomes
        statefiles_list = []
        outputfiles_list = []
        for precompfile in precompilefiles_list:
            celltype = precompfile.split('/')[-1].split('_')[0]
            print("Viterbi calling state and output file for %s" % celltype)
            state_outfile = celltype + '_states_firstr'  + '.bed'
            statefiles_list.append(state_outfile)
            output_outfile = celltype + '_output_firstr'  + '.bed'
            outputfiles_list.append(output_outfile)
            if file_check(state_outfile) and file_check(output_outfile):
                query_mark2 = query_yes_no('First HMM model Viterbi calling files exist, do you want to overwrite them?')
                if query_mark2:
                    hmm_calling(config, refgenome, state_outfile, output_outfile, precompfile, rawhmm_file, binoutfile)
                else:
                    print("Use the exists calling files")
            else:
                hmm_calling(config, refgenome, state_outfile, output_outfile, precompfile, rawhmm_file, binoutfile)

        print("Modifying rawhmm file")
        modi_rawhmm = 'HMM_' + str(numhist) +'_modi.rawhmm'
        if file_check(modi_rawhmm):
            query_mark3 = query_yes_no('First HMM model Viterbi calling file exists, do you want to overwrite it?')
            if query_mark3:
                modifyhmm(statefiles_list,rawhmm_file,modi_rawhmm,num,nucnum,percentage,outstats)
            else:
                print("Use the exist %s." % modi_rawhmm)
        else:
            modifyhmm(statefiles_list,rawhmm_file,modi_rawhmm,num,nucnum,percentage,outstats)

        if bicfile is None:
            bic_file = 'BIC_second_run' + get_time() +'.txt'
        else:
            bic_file = bicfile

        trans_matrix_modi,emit_matrix_modi,numstate_modi,numoutput_modi = load_rawhmm(modi_rawhmm)
        # second round HMM training, 200 iterations are enough for converging
        print("Second round HMM training..")
        if outrawhmmfile is None:
            out_rawhmm = 'HMM_'+str(numhist)+'_secondr.rawhmm'
        else:
            out_rawhmm = outrawhmmfile
        hmm_second(config, refgenome, out_rawhmm, b, 200, numstate_modi, numoutput_modi, precompilefiles_list, bic_file,
                   modi_rawhmm, binoutfile)

        # second round HMM calling
        print("Seconda round HMM Viterbi calling")
        statefiles_list_second = []
        for precompfile in precompilefiles_list:
            celltype = precompfile.split('/')[-1].split('_')[0]
            print("Viterbi calling state and output file for %s" % celltype)
            state_outfile = celltype + '_states_secondr.bed'
            statefiles_list_second.append(state_outfile)
            output_outfile = celltype + '_output_secondr.bed'
            if file_check(state_outfile) and file_check(output_outfile):
                query_mark2 = query_yes_no("First HMM model Viterbi calling %s file exist, do you want to overwrite it?"%celltype)
                if query_mark2:
                    hmm_calling(config, refgenome, state_outfile, output_outfile, precompfile, out_rawhmm, binoutfile)
                else:
                    print("Use the exists calling files")
            else:
                hmm_calling(config, refgenome, state_outfile, output_outfile, precompfile, out_rawhmm,binoutfile)

        print('States file: ' + ','.join(statefiles_list_second))
        if removetmpfile:
            subprocess.call("rm " + binoutfile, shell=True)
            subprocess.call("rm " + rawhmm_file, shell=True)
            subprocess.call("rm " + bic_file_1, shell=True)
            for file in statefiles_list:
                subprocess.call("rm " + file,shell=True)
            for file in outputfiles_list:
                subprocess.call("rm " + file, shell=True)
            subprocess.call("rm " + modi_rawhmm, shell=True)
        print('Nuchmm training Finish!')

@cli.command(help_priority=4)
@click.option('--rawhmmfile','-rhf',type=click.Path(exists=True),required=True,help='Input the secondr.rawhmm file resulted from nuchmm-train.')
@click.option('--histonelistfile', '-hlf', type=click.Path(exists=True), required=True,
              help='Input histone_marks.txt file that contains all histone marks used in nuchmm-init.')
@click.option('--genesfile','-gf',type=click.Path(exists=True),required=True,help='Input the interested genes list file.')
@click.option('--celltypes', '-ct',multiple=True, required=True,help=' Input all cell types name used in nuchmm-init. '
                                                                     'Multiple celltype input as -ct MCF7 -ct H1')
@click.option('--bgstate', '-bg',multiple=True, help='Specify background states that identified from the Mark-state matrix.'
                                                     'Multiple states input as -bg 1 -bg 2')
@click.option('--statesfilelist','-sfl',type=click.Path(),help='The file contains the path and name of the decoding '
                                                               'secondr states files, which should be in same order with -ct command')
@click.option('--outputsfilelist','-ofl',type=click.Path(),help='The file contains the path and name of the decoding '
                                                               'secondr output files, which should be in same order with -ct command,'
                                                                'and use with statesfilelist command')
@click.option('--plotcellmark','-pcm',is_flag=True,help='Plot the cell type specific HMM states distribution.')
@click.option('--plottotalmark','-ptm',is_flag=True,help='Plot the all cell averaged HMM states distribution.')
@click.option('--overwrite','-ow',default=True,help='Overwrite all existing files. Useful when user need to re-run command in the same directory.')
@click.option('--winmethod/--refmethod',default=False,help='Pick a method to identify unique HMM state for the nucleosomes that have multi-states. '
                                                           '`--winmethod` uses window method to identify the unique state; '
                                                           '`--refmethod` uses output_secondr.bed as reference to identify the unique state.')
@click.option('--winsize', '-ws', default=5, help='Window size for selecting the unique state. Default: 5.')
@click.option('--updistal','-uD',default=100000,help='Specify upstream boundary for plotting the distribution, cannot exceed than upboundary parameter used in nuchmm-init')
@click.option('--upproximal','-uProx',default=5000,help='Specify upstream boundary of the proximal region. '
                                                        'Used in identifying the genomic location of HMM states. Default: 5000.')
@click.option('--uppromoter','-uProm',default=1000,help='Specify upstream boundary of the promoter region. '
                                                        'Used in identifying the genomic location of HMM states. Default: 1000.')
@click.option('--downbound','-db',default=10000,help='Specify downstream boundary for plotting the distribution, '
                                                     'cannot exceed than `--downboundary` used in nuchmm-init. Default: 10000.')
@click.option('--rescalelength','-rl',default=30000,help='Specify the rescaled gene body length. Default: 30000.')
@click.option('--markthreshold','-mt',default=0.25,help='Specify the mark threshold for state-mark matrix. Default: 0.25.')
@click.option('--refhg19/--refhg38',default=True,help='Specify the reference genome. Default: built-in hg19.')
@click.option('--outfile','-of',type=click.Path(),help='Specify the path and name of the identified genomic location file.')
@click.option('--removetmpfile','-rmf',is_flag=True,help='Remove temporary files')
def NucHMM_screen_init(rawhmmfile,histonelistfile,bgstate,genesfile,celltypes, statesfilelist, outputsfilelist, winmethod,
                       overwrite, winsize, updistal, upproximal, uppromoter, downbound,rescalelength, markthreshold, refhg19,
                       plotcellmark, plottotalmark, outfile,removetmpfile):
    '''Initialize the screen step by providing sorted state file and the suggested genomic locations of HMM states'''
    background_state = list(bgstate)
    celltypes_list = list(celltypes)
    const_interval = 100
    samplepoints = rescalelength//const_interval
    filt_statefiles_list = []
    statefile_path = './'
    trans_matrix,emit_matrix,numstate,numoutput = load_rawhmm(rawhmmfile)
    state_mark2histone_mark,state4write = statem_to_histm(histonelistfile,numstate,numoutput,emit_matrix, markthreshold,background_state)
    # output mark combination for each state
    print('State Corresponding Histone marks')
    for state in sorted(state4write):
        if str(state) in background_state:
            continue
        else:
            print('S'+str(state) + ' ' +','.join(state4write[state]))
            
    if statesfilelist is None and outputsfilelist is None:
        for celltype in celltypes_list:
            outstatefile = celltype + '_states_srt_uniq.bed'
            if file_check(outstatefile):
                query_mark = query_yes_no("%s file detected, use this one?" % outstatefile)
                if query_mark:
                    print("Use exist %s" % outstatefile)
                    filt_statefiles_list.append(outstatefile)
                else:
                    print("Overwrite %s" % outstatefile)
                    filt_statefiles_list.append(outstatefile)
                    statefile_name = statefile_path + '/' + celltype + '_states_secondr.bed'
                    if file_check(statefile_name):
                        query_mark2 = query_yes_no("%s state file detected, use this one?" % celltype)
                        if query_mark2:
                            select_unique_state(state_mark2histone_mark, statefile_name, None, outstatefile, background_state,
                                                winmethod, winsize, celltype)
                        else:
                            statefile_name = input_new_name('Please input ' + celltype + ' state file from second round hmm: ' )
                            statefile_path = '/'.join((statefile_name.split('/')[:-1]))
                            if file_check(statefile_name):
                                select_unique_state(state_mark2histone_mark, statefile_name, None, outstatefile, background_state,
                                                    winmethod, winsize, celltype)
                            else:
                                print(celltype + 'secondr round file not exists!')
                                exit(1)
                    else:
                        statefile_name = input_new_name('Please input ' + celltype + ' state file from second round hmm: ' )
                        statefile_path = '/'.join((statefile_name.split('/')[:-1]))
                        if file_check(statefile_name):
                            select_unique_state(state_mark2histone_mark, statefile_name, None, outstatefile, background_state,
                                                winmethod, winsize, celltype)
                        else:
                            print(celltype + 'secondr round file not exists!')
                            exit(1)
            else:
                filt_statefiles_list.append(outstatefile)
                statefile_name = statefile_path + '/' + celltype + '_states_secondr.bed'
                if file_check(statefile_name):
                    query_mark2 = query_yes_no("%s state file detected, use this one?" % celltype)
                    if query_mark2:
                        select_unique_state(state_mark2histone_mark, statefile_name, None, outstatefile, background_state,
                                            winmethod, winsize, celltype)
                    else:
                        statefile_name = input_new_name('Please input ' + celltype + ' state file from second round hmm: ' )
                        statefile_path = '/'.join((statefile_name.split('/')[:-1]))
                        if file_check(statefile_name):
                            select_unique_state(state_mark2histone_mark, statefile_name, None, outstatefile, background_state,
                                                winmethod, winsize, celltype)
                        else:
                            print(celltype + 'secondr round file not exists!')
                            exit(1)
                else:
                    statefile_name = input_new_name('Please input ' + celltype + ' state file from second round hmm: ' )
                    statefile_path = '/'.join((statefile_name.split('/')[:-1]))
                    if file_check(statefile_name):
                        select_unique_state(state_mark2histone_mark, statefile_name, None, outstatefile, background_state,
                                            winmethod, winsize, celltype)
                    else:
                        print(celltype + 'secondr round file not exists!')
                        exit(1)
    else:
        # load secondr files
        statesfiles_list = load_histonefile(statesfilelist)
        outputsfiles_list = load_histonefile(outputsfilelist)
        for idx,statefile_name in enumerate(statesfiles_list):
            outputfile = outputsfiles_list[idx]
            celltype = celltypes_list[idx]
            outstatefile = celltype + '_states_srt_uniq.bed'
            select_unique_state(state_mark2histone_mark, statefile_name, outputfile, outstatefile, background_state,
                                winmethod, winsize, celltype)
            filt_statefiles_list.append(outstatefile)

    if refhg19:
        refgene = hg19
    else:
        refgene = hg38

    if outfile is None:
        outputfile = 'states_genomic_location.txt'
        if file_check(outputfile):
            if overwrite:
                query_mark3 = True
            else:
                query_mark3 = query_yes_no("%s exists! Do you want ot overwrite it?" % outputfile)
            if query_mark3:
                print('overwriting')
                genomic_loc_finder(genesfile,filt_statefiles_list,rawhmmfile,histonelistfile,updistal,upproximal,uppromoter,
                                   downbound,bgstate,outputfile,samplepoints,rescalelength,refgene,plotcellmark,plottotalmark,
                                   markthreshold,removetmpfile)
            else:
                print("Please use existed %s" % outputfile)
        else:
            genomic_loc_finder(genesfile,filt_statefiles_list,rawhmmfile,histonelistfile,updistal,upproximal,uppromoter,
                               downbound,bgstate,outputfile,samplepoints,rescalelength,refgene,plotcellmark,plottotalmark,
                               markthreshold,removetmpfile)
    else:
        outputfile = outfile
        if file_check(outputfile):
            if overwrite:
                query_mark3 = True
            else:
                query_mark3 = query_yes_no("%s exists! Do you want ot overwrite it?" % outputfile)
            if query_mark3:
                genomic_loc_finder(genesfile,filt_statefiles_list,rawhmmfile,histonelistfile,updistal,upproximal,uppromoter,
                                   downbound,bgstate,outputfile,samplepoints,rescalelength,refgene,plotcellmark,plottotalmark,
                                   markthreshold,removetmpfile)
            else:
                print("Please use existed %s" % outputfile)
        else:
            genomic_loc_finder(genesfile,filt_statefiles_list,rawhmmfile,histonelistfile,updistal,upproximal,uppromoter,
                               downbound,bgstate,outputfile,samplepoints,rescalelength,refgene,plotcellmark,plottotalmark,
                               markthreshold,removetmpfile)

    # write sorted and unique state file list for nuchmm-screen use
    f_name = 'nuchmm_screen_init_result_files.txt'
    if file_check(f_name):
        if overwrite:
            query_mark4 = True
        else:
            query_mark4 = query_yes_no("Do you want to overwrite %s" % f_name)
        if query_mark4:
            f_file = open(f_name, 'w')
            for file in filt_statefiles_list:
                f_file.write(file+'\n')
            f_file.close()
        else:
            print("Use the exist %s" % f_name)
    else:
        f_file = open(f_name, 'w')
        for file in filt_statefiles_list:
            f_file.write(file+'\n')
        f_file.close()


@cli.command(help_priority=5)
@click.option('--genesfile','-gf',type=click.Path(exists=True),required=True,help='Input the interested genes list file. Check the genebody_anno_hg19.txt'
                                                                                  'in annotation_files for the detailed format.')
@click.option('--like_wig_fileslist','-lwfl',type=click.Path(exists=True),required=True,
              help='Files that contains the like_wig list files name and path, the order of it should be the same with statesfile list')
@click.option('--nucpositionfiles','-nucf',type=click.Path(exists=True),required=True,
              help='Input the nucleosome position files list')
@click.option('--statesnumber', '-sn',required=True, type=int,help='Input the total number of HMM states (Include background states.)')
@click.option('--bgstate', '-bg',multiple=True,help=' Specify background states that identified from the Mark-state matrix.')
@click.option('--inputfileslist','-ifl',type=click.Path(),help='Input the list of files resulted from nuchmm-screen-init. '
                                                               'Check nuchmm_screen_init_result_files.txt in example_files folder for the detailed information. '
                                                               'The program will automatically look for the nuchmm_screen_init_result_files.txt in current directory.')
@click.option('--statesregionfile','-srf',type=click.Path(),help='a file that contain the gene region for states')
@click.option('--naregstates','-nrs',multiple=True,help='NA regularity score states (deprecated)')
@click.option('--updistal','-uD',default=100000,help='up boundary of the Distal')
@click.option('--upproximal','-uProx',default=5000,help='up boundary of the Proximal')
@click.option('--uppromoter','-uProm',default=1000,help='up boundary of the Promoter')
@click.option('--cutoffdist','-cfd',default=350,
              help='Specify the maximum distance between two nucleosomes in a nucleosome array. '
                   'If the distance between two nucleosome larger than this cutoff distance, '
                   'the program will consider these two nucleosomes are in different nucleosome arrays. Default: 350.')
@click.option('--downratio','-dr',default=5,help='Specify the bottom quantile cutoff in array-num and nucleosoem positioning filtering process. Default: 5. ')
@click.option('--upratio','-ur',default=95,help='Specify the top quantile cutoff in array-num and nucleosome positioning filtering process. Default: 95.')
@click.option('--arraydown','-adown',default=1000,help='Specify the range of the nucleosome array for calculating nucleosoem regularity and spacing. '
                                                       'Currently only accept 1000 and 2000 as input. Default: 1000.')
@click.option('--rankcoef','-rc',default=10,help='Specify the cofficient of the regularity rank for filtering unmatched nucleosome. '
                                                 'The larger the coefficient, the looser the filter. Default: 10.')
@click.option('--refhg19/--refhg38',default=True,help='Specify the reference genome. Default: built-in hg19.')
@click.option('--writeinfo','-wi',is_flag=True,help='Write files contain the detail array-num, nucleosome regularity, spacing and positioning information.')
@click.option('--plotmark','-pm',is_flag=True,help='Plot state-array distribution, bar-plot of regularity score,'
                                                   'spectral-density plot and nuc-pos violin plot')
@click.option('--max/--mean',default=True,help='Select the method to calculate the regularity score')
@click.option('--removetmpfile','-rmf',is_flag=True,help='remove temporary files')
def NucHMM_screen(genesfile,like_wig_fileslist,nucpositionfiles,bgstate,statesnumber,inputfileslist,statesregionfile,
                  naregstates,updistal,upproximal,uppromoter,cutoffdist,downratio,upratio,arraydown,rankcoef,refhg19,
                  writeinfo,plotmark,max,removetmpfile):
    '''Filter nucleosomes by genomic location, array number, nucleosome regularity, spacing and positioning'''
    if refhg19:
        refgene = hg19
    else:
        refgene = hg38

    background_state = list(bgstate)
    statesfilelist_name = 'nuchmm_screen_init_result_files.txt'
    if inputfileslist is None:
        if file_check(statesfilelist_name):
            query_mark1 = query_yes_no("Detect %s exists! Do you want to use it?" % statesfilelist_name)
            if query_mark1:
                print("Use the exist %s" % statesfilelist_name)
            else:
                statesfilelist_name = input_new_name("Please input the absolute path and name of the sorted and unique state files list: ")
        else:
            statesfilelist_name = input_new_name("Please input the absolute path and name of the sorted and unique state files list: ")
    else:
        statesfilelist_name = inputfileslist

    stateregion_name = 'states_genomic_location.txt'
    if statesregionfile is None:
        if file_check(stateregion_name):
            query_mark2 = query_yes_no("Detect %s exists! Do you want to use it?" % stateregion_name)
            if query_mark2:
                print("Use the exist %s" % stateregion_name)
            else:
                stateregion_name = input_new_name("Please input the path and name of the predicted region file: ")
        else:
            stateregion_name = input_new_name("Please input the path and name of the predicted region file: ")
    else:
        stateregion_name = statesregionfile

    statesfile_list = load_histonefile(statesfilelist_name)
    stateregionfile = stateregion_name
    array_info_dict = defaultdict(list)
    row_index = []
    cell_types = []
    print("\nFiltering by genomic location and array number..")
    for nucstatefile in statesfile_list:
        celltype = nucstatefile.split('/')[-1].split('_')[0]
        print('\n')
        print(celltype)
        cell_types.append(celltype)
        # Genomic location filtering
        filter_file_list,count_nuc,count_gl_nuc = genomic_loc_filter(nucstatefile,background_state,stateregionfile,refgene,
                                                                     celltype,genesfile,updistal,upproximal,uppromoter,
                                                                     removetmpfile)
        filter_info_print('Genomic location',count_nuc,count_gl_nuc)

        # array number filtering
        showarraylengthdistribute = False
        count_an_nuc = 0
        for file in filter_file_list:
            name_element = file.split('_')
            gl_an_filtfile = name_element[0] + '_state_' + name_element[2] + '_gl_an_filt.bed'
            print('S'+name_element[2])
            down_num, up_num, weighted_num = array_num_filter(file,gl_an_filtfile,downratio,upratio,cutoffdist,
                                                              showarraylengthdistribute)
            count_an_nuc += count_file_rows(gl_an_filtfile)
            array_info_dict[celltype].append(str(weighted_num)+' ('+str(down_num)+'-'+str(up_num)+')')

        row_index = ['S'+file.split('_')[2] for file in filter_file_list]

        filter_info_print('Array number',count_gl_nuc,count_an_nuc)

    # calculate the average number of array
    weighted_aver = [0] * len(row_index)
    down_aver = [0] * len(row_index)
    up_aver = [0] * len(row_index)
    for key in array_info_dict:
        for idx, array_info in enumerate(array_info_dict[key]):
            weighted_aver[idx] += round(float(array_info.split()[0]),2)
            down_aver[idx] += int(array_info.split()[1].split('-')[0].split('(')[1])
            up_aver[idx] += int(array_info.split()[1].split('-')[1].split(')')[0])
    weighted_aver = np.around(np.array(weighted_aver)/float(len(array_info_dict)),decimals=2)
    down_aver = np.around(np.array(down_aver)/float(len(array_info_dict)),decimals=2)
    up_aver = np.around(np.array(up_aver)/float(len(array_info_dict)),decimals=2)

    # less_3_nuc_state stores the potential NAstates
    less_3_nuc_state = {}
    for idx,state in enumerate(row_index):
        if weighted_aver[idx] < 3:
            less_3_nuc_state[state] = weighted_aver
        array_info_dict['Ave. No'].append(str(weighted_aver[idx])+' ('+str(down_aver[idx])+'-'+str(up_aver[idx])+')')

    # identify NA state
    # if naregstates is None:
    #     min_weight_aver = min(less_3_nuc_state.values())
    #     NA_nuc_state = [key[1] for key in less_3_nuc_state if less_3_nuc_state[key] == min_weight_aver]
    # else:
    #     NA_nuc_state = list(naregstates)
    if naregstates is not None:
        NA_nuc_state = list(naregstates)
    if writeinfo:
        df_array = pd.DataFrame(array_info_dict,index=row_index)
        cols = cell_types[:]
        cols.append('Ave. No')
        df_array = df_array[cols]
        out_filename = 'Array_number_information_' + get_time() +'.txt'
        df_array.to_csv(out_filename,sep='\t')

    # regularity and spacing filtering section
    print("\nFiltering by regularity score and spacing..")
    arrayup = 0
    hzup = 5
    hzdown = 7
    gl_an_filt_file_path = '.'
    like_wig_files_list = load_histonefile(like_wig_fileslist)
    numcpu = multiprocessing.cpu_count()//2
    if numcpu > statesnumber:
        numcpu = statesnumber
    '''Filter the nucleosomes by their regularity and spacing'''
    out_suffix = '_gl_an_resp_filt.bed'
    # the output file name format is celltype + state_* + out_suffix
    nuc_count,nuc_filt_count,resp_info_dict = regularity_spacing_filter(gl_an_filt_file_path,cell_types,like_wig_files_list,
                                                                        background_state,numcpu,statesnumber,arrayup,arraydown,hzup,
                                                                        hzdown,rankcoef,plotmark,max,writeinfo,cutoffdist,out_suffix,
                                                                        removetmpfile,NA_nuc_state)

    # the NAstates gl_an_filt file directly become gl_an_resp_filt file
    # for celltype in cell_types:
    #     for state in NA_nuc_state:
    #         subprocess.call("cp " + celltype + '_state_' + state + '_gl_an_filt.bed' + " " +
    #                         celltype + '_state_' + state + out_suffix)

    filter_info_print('regularity and spacing',nuc_count,nuc_filt_count)
    # combine gl_an_resp_filt files
    gl_an_resp_filt_files = []
    for celltype in cell_types:
        combined_file_name = celltype + '_total' + out_suffix
        if file_check(combined_file_name):
            query_mark3 = query_yes_no("%s is already exists, do you want to overwrite it?")
            if not query_mark3:
                gl_an_resp_filt_files.append(combined_file_name)
                print("Use exist %s." % combined_file_name)
            else:
                gl_an_resp_filt_files.append(combined_file_name)
                tmp_file = 'tmp_combine.bed'
                print("Combining the regualrity score and spacing filtered files")
                subprocess.call("cat " + celltype + '_state_*' + out_suffix + " > " + tmp_file, shell=True)
                print("Sorting the combined file")
                subprocess.call("sort -k1,1V -k2,2n -k3,3n " + tmp_file + " > " + combined_file_name, shell=True)
                subprocess.call("rm " + tmp_file, shell=True)
        else:
            gl_an_resp_filt_files.append(combined_file_name)
            tmp_file = 'tmp_combine.bed'
            print("Combining the regualrity score and spacing filtered files")
            subprocess.call("cat " + celltype + '_state_*' + out_suffix + " > " + tmp_file, shell=True)
            print("Sorting the combined file")
            subprocess.call("sort -k1,1V -k2,2n -k3,3n " + tmp_file + " > " + combined_file_name, shell=True)
            subprocess.call("rm " + tmp_file, shell=True)

    if removetmpfile:
        for celltype in cell_types:
            subprocess.call("rm " + celltype + '_state_*' + out_suffix, shell=True)
            subprocess.call("rm " + celltype + '_state_*' + '_gl_filt.bed', shell=True)
            subprocess.call("rm " + celltype + '_state_*' + '_gl_an_filt.bed', shell=True)
        subprocess.call("rm " + '*bin300*', shell=True)

    # nucleosome positioning filtering section
    print("\nFiltering by nucleosome positioning..")
    k = 1.5
    states_list = []
    df_row = []
    for state in range(1, statesnumber + 1):
        if str(state) not in background_state:
            states_list.append(str(state))
            df_row.append('S'+str(state))

    # the files in inputfilelist and nucdetailfilelist should be in same order
    inputfile_list = gl_an_resp_filt_files
    nucdetailfile_list = load_histonefile(nucpositionfiles)

    outfile_list = []
    for file in inputfile_list:
        celltype = file.split('/')[-1].split('_')[0]
        outfile_list.append(celltype+'_gl_an_resp_pos_filt.bed')

    info_wrm_dict = {}
    info_worm_dict = {}
    celltype_list = []
    up_ratio = float(upratio)/100
    down_ratio = float(downratio)/100
    showinfo = False
    states_pos_total = {state:[] for state in states_list}
    for idx,file in enumerate(inputfile_list):
        nuc_detail_file = nucdetailfile_list[idx]
        out_file = outfile_list[idx]
        celltype = file.split('/')[-1].split('_')[0]
        celltype_list.append(celltype)
        print(celltype)
        nuc_count,nuc_filt_count,pos_mean_dict,states_pos = nuc_positioning_filter(file,out_file,states_list,
                                                                        nuc_detail_file,k,up_ratio,down_ratio,
                                                                        showinfo)
        for state in states_pos:
            states_pos_total[state] += states_pos[state]

        info_wrm_dict[celltype] = pos_mean_dict['Mean_wrm']
        info_worm_dict[celltype] = pos_mean_dict['Mean_worm']
        filter_info_print('nucleosome positioning',nuc_count,nuc_filt_count)

    if plotmark:
        figname = "Nuc_pos_vio_" + get_time() + ".png"
        plot_violin_nuc_pos(states_pos_total,states_list,spe_colors, figname)

    # Calculate the average positioning for each state
    aver_wrm = np.zeros(len(df_row))
    aver_worm = np.zeros(len(df_row))
    for cell_type in info_wrm_dict:
        aver_wrm += np.array(info_wrm_dict[cell_type])
        aver_worm += np.array(info_worm_dict[cell_type])
    aver_wrm = aver_wrm/len(info_wrm_dict)
    aver_worm = aver_worm/len(info_worm_dict)
    info_wrm_dict['Ave. Pos'] = aver_wrm.tolist()
    info_worm_dict['Ave. Pos'] = aver_worm.tolist()

    if writeinfo:
        df_wrm = pd.DataFrame(info_wrm_dict,index=df_row)
        df_worm = pd.DataFrame(info_worm_dict,index=df_row)
        celltype_list.append('Ave. Pos')
        cols = celltype_list[:]
        df_wrm = df_wrm[cols].round(2)
        df_worm = df_worm[cols].round(2)
        out_filename_wrm = 'Nucleosome_positioning_information_wrm_' + get_time() +'.txt'
        out_filename_worm = 'Nucleosome_positioning_information_worm_' + get_time() +'.txt'
        df_wrm.to_csv(out_filename_wrm,sep='\t')
        df_worm.to_csv(out_filename_worm,sep='\t')

    print("Final filtered files are: " + '\t'.join(outfile_list))
    genomic_loc = []
    for line in load_histonefile(stateregion_name):
        genomic_loc.append('/'.join(line.strip().split()[1:]))
    final_table_name_post = 'functional_nucleosome_state_post.txt'
    final_info_dict_post = {}
    final_cols = ['Genomic Location','Ave. No. of Nucs','Ave. Spacing','Regularity score', 'Degree of Positioning']
    final_info_dict_post['Ave. No. of Nucs'] = array_info_dict['Ave. No']
    final_info_dict_post['Genomic Location'] = genomic_loc
    final_info_dict_post['Ave. Spacing'] = np.round(resp_info_dict['Post.Ave. Spacing'],2)
    final_info_dict_post['Regularity score'] = np.round(resp_info_dict['Post.Regularity Score'],2)
    final_info_dict_post['Degree of Positioning'] = np.round(info_wrm_dict['Ave. Pos'],2)
    try:
        df_final_info_post = pd.DataFrame(final_info_dict_post,index=df_row)
    except ValueError:
        print("Please check the -sn parameter")
        exit(1)
    df_final_info_post = df_final_info_post[final_cols]
    df_final_info_post.to_csv(final_table_name_post,sep='\t')
    print("Final nucleosome feature tables is %s" % (final_table_name_post))

    if writeinfo:
        final_table_name_pre = 'functional_nucleosome_state_pre.txt'
        final_info_dict_pre = {}
        final_info_dict_pre['Ave. No. of Nucs'] = array_info_dict['Ave. No']
        final_info_dict_pre['Genomic Location'] = genomic_loc
        final_info_dict_pre['Ave. Spacing'] = np.round(resp_info_dict['Pre.Ave. Spacing'],2)
        final_info_dict_pre['Regularity score'] = np.round(resp_info_dict['Pre.Regularity Score'],2)
        final_info_dict_pre['Degree of Positioning'] = np.round(info_worm_dict['Ave. Pos'],2)
        df_final_info_pre = pd.DataFrame(final_info_dict_pre,index=df_row)
        df_final_info_pre = df_final_info_pre[final_cols]
        df_final_info_pre.to_csv(final_table_name_pre,sep='\t')
        print("Final nucleosome feature tables are %s and %s" % (final_table_name_post,final_table_name_pre))
    print("Screen Finish!")

@cli.command(help_priority=6)
@click.option('--rawhmmfile', '-rhf',type=click.Path(exists=True), required=True, help='Input the secondr.rawhmm file resulted from nuchmm-train.')
@click.option('--histonelistfile', '-hlf', type=click.Path(exists=True), required=True,
              help='Input histone_marks.txt file that contains all histone marks used in nuchmm-init.')
@click.option('--matrixcolor','-mc',default=2, help='Specify the color palette of the matrix. 0 is red-white; '
                                                    '1 is red-yellow(YlOrRd) and 2 is red-blue(coolwarm). Default: 2.')
@click.option('--markthreshold','-mt',default=0.25,help='Specify the threshold to show the probability in the matrix. Default: 0.25.')
@click.option('--transmat','-tmat',type=click.Path(),
              help='Specify the path and name of the transition probability matrix, otherwise will automatically save to trans.< current time >.png')
@click.option('--markstatemat','-msmat',type=click.Path(),
              help=' Specify the path and name of the mark-state matrix, otherwise will automatically save to mark_state.< current time >.png')
@click.option('--emitmat','-emat',type=click.Path(),help='Specify the path and name of the emission probability matrx. Default: None.')
@click.option('mark','--M',is_flag=True,help='Divide the element of emission probability matrix by the number of the mark in that observation. '
                                             'Very strong assumption. Not recommend use. ')
@click.option('--writematrix', is_flag=True, help='Write mark_state matrix and transistion probability matrix in txt files.')
def matrix_visualize(rawhmmfile,histonelistfile,transmat,markstatemat,mark,writematrix,emitmat,matrixcolor,markthreshold):
    '''Visualize the Transition and Mark-state matrix.'''
    if markthreshold == 0.25:
        print("Use default mark threshold 0.25, -mt command is used to specify mark threshold.")
    HMM_matrix_visualization(rawhmmfile,histonelistfile,transmat,markstatemat,mark,writematrix,matrixcolor,emitmat,markthreshold)

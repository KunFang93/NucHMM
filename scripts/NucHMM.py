#!/usr/bin/env python
#--coding:utf-8 --


#################################################################################################
#################################################################################################
########                                                                                 ########
########    Integrated HMM Methods to identify nucleosome states                         ########
########                                                                                 ########
########    Author:  Kun Fang                                                            ########
########                                                                                 ########
########                                                                                 ########
########    Working Environment:  Python3                                                ########
########                                                                                 ########
########    Date:      2020-07-08                                                        ########
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
from NucHMM_Feature_screen import nuc_positioning_filter,regularity_spacing_filter
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
@click.group(cls=SpecialHelpOrder,invoke_without_command=True)
@click.version_option(version=1.0)
@click.option('--verbose', is_flag=True)
@click.option('--hmm-directory', type=click.Path(),help='the path of NucHMM_Model bin')
@pass_config
def cli(config, verbose, hmm_directory):
    config.verbose = verbose
    if hmm_directory is None:
        hmm_directory = '.'
    config.hmm_directory = hmm_directory

@cli.command(help_priority=1)
@click.option('--fastq/--bam', default=False, required = True, help = 'Input data type')
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
@click.option('--fastqc', '-qc', is_flag=True, help = 'whether run QC trim step for fastq file, use with --fastq.')
@click.option('--bowtieindexpath', '-bip', type=click.Path(exists=True), help = 'Required if input format is fastq, the path of bowtie index fold path')
@click.option('--bowtie2indexpath', '-b2ip', type=click.Path(exists=True), help = 'Required if input format is fastq, the path of bowtie2 index fold path')
@click.option('--inpspath','-inps',type=click.Path(), help = 'The path of iNPS.py. For example, /data/iNPS/iNPS.py use with --mnase')
@click.option('--threads','-p',default=5, help = 'number of threads')
def NucHMM_prep(fastq,inputfqslist,inputbamslist,fastqc,bowtieindexpath,bowtie2indexpath,inpspath,threads):
    '''Prepare files for nuchmm-init from fastq/bam files.
    NucHMM only provide basic fastq and peak calling pipeline. So we strongly recommend user use their favorite pipeline to process
    their fastq files and peaking the peaks. MNase-seq data should use this function'''
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
              help='The input TF peaks\' files list')
@click.option('--nucpositionfiles', '-nucf', type=click.Path(exists=True), required=True,
              help='The list of nucleosome position files')
@click.option('--genefile', '-gf', type=click.Path(exists=True), required=True, help='Interested genes\' file')
@click.option('--intersect_cutoff', '-ic',default=0.3, help='The intersect threshold to determine whether nucleosome have this peak')
@click.option('--gap', '-g', is_flag=True, help='give mark to the gaps between nucleosomes')
@click.option('--upboundary', '-up',default=100000, help='Upstream boundary of gene segmentation')
@click.option('--downboundary', '-down',default=10000, help='Downstream boundary of gene segmentation')
@click.option('--outputfilelist', '-outf', type=click.Path(), help='Output file\'s path and name' )
@click.option('--refgenome', '-refg', type=click.Path(), help='The reference genome file, default is hg19')
@click.option('--removetmpfile','-rmf',is_flag=True,help='Remove temporary files')
def NucHMM_init(inputpeakslistfiles, nucpositionfiles, intersect_cutoff, gap, genefile, upboundary, downboundary,
                outputfilelist, refgenome, removetmpfile):
    '''Assign histone modification peaks to nucleosomes and create .precomp files for hmm training. All .precomp files name
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
        prec_outputfile_list = 'Precomp_outlist.txt'
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
@click.option('--refgenome', '-refg', type=click.Path(exists=True), required=True, help='The reference genome file')
@click.option('--precomp_list', '-pl', type=click.Path(exists=True), required=True, help='Precomp file\'s list from nuchmm-init')
@click.option('numhist','--numhistonemarks', '-numh',type=int, required=True,help='Number of histone marks used in nuchmm-init')
@click.option('numstates', '--numStates', '-nums',default=20, help='Number of states for first round hmm training')
@click.option('-b', is_flag=True,
              help='Calculate and report BIC after every iteration, output file name of BIC score will be need. WARNING: This parameter'
                   'will slow the algorithm a lot!')
@click.option('-i', default=300, help='Maximum number of iterations for first round hmm training')

@click.option('--num/--perc',default=False,help='Two methods to remove the redundant states, recommended method is '
                                               'percentage method')
@click.option('--nucnum', '-nn', default=350, help='The minimum number of nucleosomes to keep the states')
@click.option('--percentage', '-pn', default=0.005, help='The minimum percentage of nucleosomes to keep the states,'
                                                 'default 0.005 (0.5%)')
@click.option('--outrawhmmfile', '-ohmm', type=click.Path(), help='Output file\'s path and name')
@click.option('--outstats', '-os', type=click.Path(), help='The nucleosomes counts for each states in each cell type')
@click.option('--bicfile', '-obic', type=click.Path(),help='BIC score file, if -b not specified, only report first and '
                                                           'last BIC score')
@click.option('--removetmpfile','-rmf',is_flag=True,help='Remove temporary files')
@pass_config
def NucHMM_train(config, numhist, refgenome, b, i, numstates, precomp_list, num, nucnum, percentage, outstats,
                    outrawhmmfile, bicfile, removetmpfile):
    '''Use hidden markov model to acquire initial states for each nucleosomes'''

    # check if NucHMM-learn etc in directory

    # number of outputs is 2 ** number of histone modifications
    numoutputs = 2 ** numhist

    # Read precompiles list file to a list
    precompilefiles_list = load_histonefile(precomp_list)

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
@click.option('--rawhmmfile','-rhf',type=click.Path(exists=True),required=True,help='rawhmm file from hmm second round.')
@click.option('--histonelistfile', '-hlf', type=click.Path(exists=True), required=True,
              help='The order of histones in histone file must be same with the inputlistfile for tfassign')
@click.option('--bgstate', '-bg',multiple=True, help='background states, multiple states input as -bg 1 -bg 2')
@click.option('--genesfile','-gf',type=click.Path(exists=True),required=True,help='Gene lists, format:chr   start   end strand  '
                                                                          'gene_symbol(unique).')
@click.option('--celltypes', '-ct',multiple=True, required=True,help='All cell types name, multiple name input as -ct MCF7 -ct H1')
@click.option('--statesfilelist','-sfl',type=click.Path(),help='the file contains the path and name of the decoding '
                                                               'secondr states files, which should be in same order with -ct command')
@click.option('--outputsfilelist','-ofl',type=click.Path(),help='the file contains the path and name of the decoding '
                                                               'secondr output files, which should be in same order with -ct command,'
                                                                'and use with statesfilelist command')
@click.option('--overwrite','-ow',default=True,help='overwrite all existing files')
@click.option('--winmethod/--refmethod',default=False,help='winmethod is using window method to identify unique state '
                                                          'and refmethod is using output.bed file as reference to '
                                                          'identify unique state')
@click.option('--winsize', '-ws', default=5, help='window size for selecting the unique state')
@click.option('--updistal','-uD',default=100000,help='up boundary of the Distal, cannot exceed than upboundary parameter used in nuchmm-init')
@click.option('--upproximal','-uProx',default=5000,help='up boundary of the Proximal')
@click.option('--uppromoter','-uProm',default=1000,help='up boundary of the distribution')
@click.option('--downbound','-db',default=10000,help='down boundary of the distribution,cannot exceed than downboundary parameter used in nuchmm-init')
@click.option('--rescalelength','-rl',default=10000,help='the length of rescaled gene body')
@click.option('--refhg19/--refhg38',default=True,help='select reference genome')
@click.option('--plotcellmark','-pcm',is_flag=True,help='plot the cell-specific distribution')
@click.option('--plottotalmark','-ptm',is_flag=True,help='plot the distribution of all cell types')
@click.option('--outfile','-of',type=click.Path(),help='The path and name of the identified genomic location file')
@click.option('--removetmpfile','-rmf',is_flag=True,help='remove temporary files')
def NucHMM_screen_init(rawhmmfile,histonelistfile,bgstate,genesfile,celltypes, statesfilelist, outputsfilelist, winmethod,
                       overwrite, winsize, updistal, upproximal,uppromoter, downbound,rescalelength, refhg19,
                       plotcellmark, plottotalmark, outfile,removetmpfile):
    '''Initialize the screen step by providing sorted and unqiue state file and the suggested genomic locations of each states'''
    background_state = list(bgstate)
    celltypes_list = list(celltypes)
    const_interval = 100
    samplepoints = rescalelength//const_interval
    filt_statefiles_list = []
    if statesfilelist is None and outputsfilelist is None:
        for celltype in celltypes_list:
            outstatefile = celltype + '_states_srt_uniq.bed'
            if file_check(outstatefile):
                query_mark = query_yes_no("%s file detected, use this one?" % outstatefile)
                if query_mark:
                    print("Use exist %s" % outstatefile)
                    filt_statefiles_list.append(outstatefile)
                else:
                    pirnt("Overwrite %s" % outstatefile)
                    filt_statefiles_list.append(outstatefile)
                    statefile_name = celltype + '_states_secondr.bed'
                    if file_check(statefile_name):
                        query_mark2 = query_yes_no("%s state file detected, use this one?" % celltype)
                        if query_mark2:
                            select_unique_state(rawhmmfile, histonelistfile, statefile_name, None, outstatefile, background_state,
                                                winmethod, winsize, celltype)
                        else:
                            statefile_name = input_new_name('Please input ' + celltype + ' state file from second round hmm: ' )
                            if file_check(statefile_name):
                                select_unique_state(rawhmmfile, histonelistfile, statefile_name, None, outstatefile, background_state,
                                                    winmethod, winsize, celltype)
                            else:
                                print(celltype + 'secondr round file not exists!')
                                exit(1)
                    else:
                        statefile_name = input_new_name('Please input ' + celltype + ' state file from second round hmm: ' )
                        if file_check(statefile_name):
                            select_unique_state(rawhmmfile, histonelistfile, statefile_name, None, outstatefile, background_state,
                                                winmethod, winsize, celltype)
                        else:
                            print(celltype + 'secondr round file not exists!')
                            exit(1)
            else:
                filt_statefiles_list.append(outstatefile)
                statefile_name = celltype + '_states_secondr.bed'
                if file_check(statefile_name):
                    query_mark2 = query_yes_no("%s state file detected, use this one?" % celltype)
                    if query_mark2:
                        select_unique_state(rawhmmfile, histonelistfile, statefile_name, None, outstatefile, background_state,
                                            winmethod, winsize, celltype)
                    else:
                        statefile_name = input_new_name('Please input ' + celltype + ' state file from second round hmm: ' )
                        if file_check(statefile_name):
                            select_unique_state(rawhmmfile, histonelistfile, statefile_name, None, outstatefile, background_state,
                                                winmethod, winsize, celltype)
                        else:
                            print(celltype + 'secondr round file not exists!')
                            exit(1)
                else:
                    statefile_name = input_new_name('Please input ' + celltype + ' state file from second round hmm: ' )
                    if file_check(statefile_name):
                        select_unique_state(rawhmmfile, histonelistfile, statefile_name, None, outstatefile, background_state,
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
            select_unique_state(rawhmmfile, histonelistfile, statefile_name, outputfile, outstatefile, background_state,
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
                                   downbound,bgstate,outputfile,samplepoints,rescalelength,refgene,plotcellmark,plottotalmark,removetmpfile)
            else:
                print("Please use existed %s" % outputfile)
        else:
            genomic_loc_finder(genesfile,filt_statefiles_list,rawhmmfile,histonelistfile,updistal,upproximal,uppromoter,
                               downbound,bgstate,outputfile,samplepoints,rescalelength,refgene,plotcellmark,plottotalmark,removetmpfile)
    else:
        outputfile = outfile
        if file_check(outputfile):
            if overwrite:
                query_mark3 = True
            else:
                query_mark3 = query_yes_no("%s exists! Do you want ot overwrite it?" % outputfile)
            if query_mark3:
                genomic_loc_finder(genesfile,filt_statefiles_list,rawhmmfile,histonelistfile,updistal,upproximal,uppromoter,
                                   downbound,bgstate,outputfile,samplepoints,rescalelength,refgene,plotcellmark,plottotalmark,removetmpfile)
            else:
                print("Please use existed %s" % outputfile)
        else:
            genomic_loc_finder(genesfile,filt_statefiles_list,rawhmmfile,histonelistfile,updistal,upproximal,uppromoter,
                               downbound,bgstate,outputfile,samplepoints,rescalelength,refgene,plotcellmark,plottotalmark,removetmpfile)

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
@click.option('--genesfile','-gf',type=click.Path(exists=True),required=True,help='Gene lists, format:chr   start   end strand  '
                                                                          'gene_symbol(unique).')
@click.option('--like_wig_fileslist','-lwfl',type=click.Path(exists=True),required=True,
              help='Files that contains the like_wig list files name and path, the order of it should be the same with statesfile list')
@click.option('--nucdetailfilelist','-ndfl',type=click.Path(exists=True),required=True,
              help='A file contains all the nucleosome detail information files')
@click.option('--bgstate', '-bg',multiple=True, required=True,help='background states, multiple states input as -bg 1 -bg 2')
@click.option('--statesnumber', '-sn',required=True, type=int,help='Total number of states')
@click.option('--inputfileslist','-ifl',type=click.Path(),help='a list of files resulted from nuchmm-screen-init')
@click.option('--statesregionfile','-srf',type=click.Path(),help='a file that contain the gene region for states')
@click.option('--naregstates','-nrs',multiple=True,help='NA regularity score states (deprecated)')
@click.option('--updistal','-uD',default=100000,help='up boundary of the Distal')
@click.option('--upproximal','-uProx',default=5000,help='up boundary of the Proximal')
@click.option('--uppromoter','-uProm',default=1000,help='up boundary of the Promoter')
@click.option('--cutoffdist','-cfd',default=350,help='the max distance between nucleosomes to be considered as array')
@click.option('--downratio','-drl',default=5,help='bottom cutoff for filtering process')
@click.option('--upratio','-url',default=95,help='top cutoff in array filtering process')
@click.option('--arraydown','-adown',default=1000,help='The end point of the array, right now only can choose between 1000 and 2000')
@click.option('--rankcoef','-rc',default=10,help='Rank cofficient that use for screening the nucleosome')
@click.option('--refhg19/--refhg38',default=True,help='select reference genome')
@click.option('--writeinfo','-wi',is_flag=True,help='Write the state features information')
@click.option('--plotmark','-pm',is_flag=True,help='Plot state-array distribution, bar-plot of regularity score,'
                                                   'spectral-density plot and nuc-pos violin plot')
@click.option('--max/--mean',default=True,help='Select the method to calculate the regularity score')
@click.option('--removetmpfile','-rmf',is_flag=True,help='remove temporary files')
def NucHMM_screen(genesfile,like_wig_fileslist,nucdetailfilelist,bgstate,statesnumber,inputfileslist,statesregionfile,
                  naregstates,updistal,upproximal,uppromoter,cutoffdist,downratio,upratio,arraydown,rankcoef,refhg19,
                  writeinfo,plotmark,max,removetmpfile):
    '''Filter nucleosoems by genomic location, array number, regularity score, spacing and nucleosome positioning'''
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
    nucdetailfile_list = load_histonefile(nucdetailfilelist)

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
@click.option('--rawhmmfile', '-rh',type=click.Path(exists=True), required=True, help='rawhmm file for visualization')
@click.option('--histonelistfile', '-hf', type=click.Path(exists=True), required=True,
              help='The order of histones in histone file must be same with the inputlistfile for tfassign')
@click.option('--transmat','-tmat',type=click.Path(),
              help='path and name of transition probability matrix figure, default save to trans.jpg')
@click.option('--markstatemat','-msmat',type=click.Path(),
              help='path and name of mark-state matrix figure, derived from emission probability matrix, '
                   'default save to mark_state.jpg')
@click.option('mark','--M',is_flag=True,help='divide probability by the number or mark, '
                                             'very strong assumption, carefully use')
@click.option('--writematrix', is_flag=True, help="Write mark_state_matrix.txt and trans_matrix.txt for some following "
                                                  "analysis")
def hmm_visualize(rawhmmfile,histonelistfile,transmat,markstatemat,mark,writematrix):
    '''Visualization of Transition matrix and mark-state matrix'''
    HMM_matrix_visualization(rawhmmfile,histonelistfile,transmat,markstatemat,mark,writematrix)


# @cli.command(help_priority=6)
# @click.option('--rawhmmfile','-rhf',type=click.Path(exists=True),required=True,help='rawhmm file from hmm second round.')
# @click.option('--genesfile','-gf',type=click.Path(exists=True),required=True,help='Gene lists, format:chr   start   end strand  '
#                                                                           'gene_symbol(unique).')
# @click.option('--stateslistfile','-sf',type=click.Path(exists=True),required=True,help='State.merged.srt.bed list '
#                                                                              'file to plot distribution.')
# @click.option('--histonelistfile', '-hlf', type=click.Path(exists=True), required=True,
#               help='The order of histones in histone file must be same with the inputlistfile for tfassign')
# @click.option('--stateregionfile','-srf',type=click.Path(exists=True),required=True,
#               help='state region file which is the result from state_loc_finder')
# @click.option('--bgstate', '-bg',multiple=True, required=True,help='background states, multiple states input as -bg 1 -bg 2')
# @click.option('--updistal','-uD',default=100000,help='up boundary of the Distal')
# @click.option('--upproximal','-uProx',default=5000,help='up boundary of the Proximal')
# @click.option('--uppromoter','-uProm',default=1000,help='up boundary of the Promoter')
# @click.option('--cutoffdist','-cfd',default=350,help='the max distance between nucleosomes to be considered as array')
# @click.option('--downratiolimit','-drl',default=25,help='down ratio cutoff in array filtering process, '
#                                                         'e.g 25 mean remove array with length below 25% array length distribution')
# @click.option('--upratiolimit','-url',default=90,help='up ratio cutoff in array filtering process, '
#                                                       'e.g 90 mean remove array with length up 90% array length distribution')
# @click.option('--refhg19/--refhg38',default=True,help='select reference genome')
# @click.option('--writearrayinfo','-wai',is_flag=True,help='write the array number information')
# @click.option('--writefilename','-wfn',type=click.Path(), help = 'the path and name of array information file, use with'
#                                                                  'writearrayinfo')
# @click.option('--showarraylengthdistribute','-showald',is_flag=True, help='show the array length distribution during filtering process')
# @click.option('--removetmpfile','-rmf',is_flag=True,help='remove temporary files')
# def location_array_filter(genesfile,stateslistfile,rawhmmfile,stateregionfile,bgstate,updistal,upproximal,uppromoter,cutoffdist,
#                          downratiolimit,upratiolimit,refhg19,removetmpfile,writearrayinfo,writefilename,
#                          showarraylengthdistribute):
#     '''
#     Filter nucleosomes by the genomic location and array number.
#     '''
#     if refhg19:
#         refgene = hg19
#     else:
#         refgene = hg38
#
#     background_state = list(bgstate)
#     statesfile_list = load_histonefile(stateslistfile)
#     array_info_dict = defaultdict(list)
#     row_index = []
#     cell_types = []
#     for nucstatefile in statesfile_list:
#         celltype = nucstatefile.split('/')[-1].split('_')[0]
#         print(celltype)
#         cell_types.append(celltype)
#         # Genomic location filtering
#
#         filter_file_list,count_nuc,count_gl_nuc = genomic_loc_finder(genesfile,statesfile_list,rawhmmfile,histonelistfile,
#                                                                      updistal,upproximal,uppromoter,downbound,bgstate,
#                                                                      outputfile,samplepoints,rescalelength,refgene,
#                                                                      plotcellmark,plottotalmark,removetmpfile)
#         filter_info_print('Genomic location',count_nuc,count_gl_nuc)
#
#         # array number filtering
#         count_an_nuc = 0
#         for file in filter_file_list:
#             name_element = file.split('_')
#             gl_an_filtfile = name_element[0] + '_state_' + name_element[2] + '_gl_an_filt.bed'
#             down_num, up_num, weighted_num = array_num_filter(file,gl_an_filtfile,downratiolimit,upratiolimit,cutoffdist,
#                                                               showarraylengthdistribute)
#             count_an_nuc += count_file_rows(gl_an_filtfile)
#             array_info_dict[celltype].append(str(weighted_num)+' ('+str(down_num)+'-'+str(up_num)+')')
#
#         row_index = ['S'+file.split('_')[2] for file in filter_file_list]
#
#         filter_info_print('Array number',count_gl_nuc,count_an_nuc)
#
#     # calculate the average number of array
#     weighted_aver = [0] * len(row_index)
#     down_aver = [0] * len(row_index)
#     up_aver = [0] * len(row_index)
#     for key in array_info_dict:
#         for idx, array_info in enumerate(array_info_dict[key]):
#             weighted_aver[idx] += int(array_info.split()[0])
#             down_aver[idx] += int(array_info.split()[1].split('-')[0].split('(')[1])
#             up_aver[idx] += int(array_info.split()[1].split('-')[1].split(')')[0])
#     weighted_aver = np.around(np.array(weighted_aver)/float(len(array_info_dict)),decimals=2)
#     down_aver = np.around(np.array(down_aver)/float(len(array_info_dict)),decimals=2)
#     up_aver = np.around(np.array(up_aver)/float(len(array_info_dict)),decimals=2)
#
#     for idx,state in enumerate(row_index):
#         array_info_dict['Ave. No'].append(str(weighted_aver[idx])+' ('+str(down_aver[idx])+'-'+str(up_aver[idx])+')')
#
#     if writearrayinfo:
#         df_array = pd.DataFrame(array_info_dict,index=row_index)
#         cell_types.append('Ave. No')
#         cols = cell_types[:]
#         df_array = df_array[cols]
#         if writefilename is None:
#             out_filename = 'Array_number_information_' + get_time() +'.txt'
#         else:
#             out_filename = writefilename
#
#         df_array.to_csv(out_filename,sep='\t')


@cli.command(help_priority=7)
@click.option('--inputfilespath','-ifp',type=click.Path(),required=True,help='A path contain all input files')
@click.option('--celltypes', '-ct',multiple=True, required=True,help='All cell types name, multiple name input as -ct MCF7 -bg H1')
@click.option('--like_wig_fileslist','-lwfl',type=click.Path(exists=True),required=True,help='Files that contains the like_wig '
                                                                                             'list files name and path, the order of it'
                                                                                             'should be the same with -ct command.')
@click.option('--bgstate', '-bg',multiple=True, required=True,help='Background states, multiple states input as -bg 1 -bg 2')
@click.option('--statesnumber', '-sn',required=True, type=int,help='Total number of states')
@click.option('--writeinfo','-wi',is_flag=True,help='Write the regularity score and spacing information')
@click.option('--numcpu','-ncpu',default=4,help='The number of cpu use in coverage calculation')
@click.option('--arrayup','-aup',default=0,help='The start point of the array, not recommand change it ')
@click.option('--arraydown','-adown',default=2000,help='The end point of the array, right now only can choose between 1000 and 2000')
@click.option('--hzup','-hzu',default=5,help='up-boundary of the window for calculating regularity score, 5hz corresponding 200bp peak')
@click.option('--hzdown','-hzd',default=7,help='down-boundary of the window for calculating regularity score, 5hz corresponding 200bp peak')
@click.option('--rankcoef','-rc',default=10,help='Rank cofficient that use for screening the nucleosome')
@click.option('--plotmark','-pm',is_flag=True,help='Plot state-array distribution, bar-plot of regularity score, and '
                                                   'spectral-density plot')
@click.option('--max/--mean',default=True,help='Select the method to calculate the regularity score')
@click.option('--cutoffdist','-cfd',default=350,help='The max distance between nucleosomes to be considered as array')
@click.option('--removetmpfile','-rmf',is_flag=True,help='Remove temporary files')
def reg_spacing_filter(inputfilespath,celltypes,like_wig_fileslist,bgstate,numcpu,statesnumber,arrayup,arraydown,hzup,
                       hzdown,rankcoef,plotmark,max,writeinfo,cutoffdist,removetmpfile):
    '''Use regularity and spacing information to filter the nucleosomes'''
    gl_an_filt_file_path = inputfilespath
    celltype_list = list(celltypes)
    background_state = list(bgstate)
    like_wig_files_list = load_histonefile(like_wig_fileslist)
    out_suffix = input_new_name('Please input filtered file suffix')
    '''Filter the nucleosomes by their regularity and spacing'''
    nuc_count,nuc_filt_count,resp_info_dict = regularity_spacing_filter(gl_an_filt_file_path,celltype_list,like_wig_files_list,
                                                         background_state,numcpu,statesnumber,arrayup,arraydown,hzup,
                                                         hzdown,rankcoef,plotmark,max,writeinfo,cutoffdist,out_suffix,removetmpfile)
    filter_info_print('regularity and spacing',nuc_count,nuc_filt_count)


@cli.command(help_priority=8)
@click.option('--inputfilelist','-ifl',type=click.Path(exists=True),required=True,help='A file contains all files needed '
                                                                                       'to be filter')
@click.option('--nucdetailfilelist','-ndfl',type=click.Path(exists=True),required=True,help='A file contains all the '
                                                                                            'nucleosome detail inforoamtion'
                                                                                            'files')
@click.option('--statesnumber', '-sn',required=True, type=int,help='total number of states')
@click.option('--bgstate', '-bg',multiple=True, required=True,help='background states, multiple states input as -bg 1 -bg 2')
@click.option('-k',default=1.5,help='Coefficient for calculating outliers')
@click.option('--upratio','-upr',default=0.95,help='range [0-1], the higher the number keeps more nucleosoems')
@click.option('--downratio','-downr',default=0.05,help='range [0-1], the lower the number keeps more nucleosoems')
@click.option('--outputfilelist','-ofl',type=click.Path(), help = 'A file contain the path and name of all fitlered files')
@click.option('--showinfo','-show',is_flag=True, help='show the positioning information during filtering process')
@click.option('--writeposinfo','-wpi',is_flag=True,help='write the nucleosome positioning information')
def nuc_pos_filter(inputfilelist,nucdetailfilelist,statesnumber,bgstate,k,upratio,downratio,outputfilelist,showinfo,
                   writeposinfo):
    '''Use nucleosome positioning information to filtering nucleosomes'''
    background_states = list(bgstate)
    states_list = []
    df_row = []
    for state in range(1, statesnumber + 1):
        if str(state) not in background_states:
            states_list.append(str(state))
            df_row.append('S'+str(state))

    # the files in inputfilelist and nucdetailfilelist should be in same order
    inputfile_list = load_histonefile(inputfilelist)
    nucdetailfile_list = load_histonefile(nucdetailfilelist)

    if outputfilelist is None:
        outfile_list = []
        for file in inputfile_list:
            celltype = file.split('/')[-1].split('_')[0]
            outfile_list.append(celltype+'_positioning_filtered_' + get_time() + '.bed')
    else:
        outfile_list = outputfilelist

    info_wrm_dict = {}
    info_worm_dict = {}
    celltype_list = []
    states_pos_total = {state:[] for state in states_list}
    for idx,file in enumerate(inputfile_list):
        nuc_detail_file = nucdetailfile_list[idx]
        out_file = outfile_list[idx]
        celltype = file.split('/')[-1].split('_')[0]
        celltype_list.append(celltype)
        print(celltype)
        nuc_count,nuc_filt_count,pos_mean_dict,states_pos = nuc_positioning_filter(file,out_file,states_list,
                                                                        nuc_detail_file,k,upratio,downratio,
                                                                        showinfo)
        for state in states_pos:
            states_pos_total[state] += states_pos[state]

        info_wrm_dict[celltype] = pos_mean_dict['Mean_wrm']
        info_worm_dict[celltype] = pos_mean_dict['Mean_worm']
        filter_info_print('nucleosome positioning',nuc_count,nuc_filt_count)

    if writeposinfo:
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


    if writeposinfo:
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

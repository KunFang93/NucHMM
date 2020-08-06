#!/usr/bin/env python
#--coding:utf-8 --

import re
import sys
import operator
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict, Counter
from scipy import signal, interpolate
from NucHMM_common_data import spe_colors, Histone_location, hg19, hg38
from NucHMM_Visualization import emit2mark, plot_distribution
from NucHMM_utilities import Chr2Num, file_check, query_yes_no, boundary_check, bedops_diff, get_time, count_file_rows
from NucHMM_utilities import outlier_threshold, input_new_name
from NucHMM_Load_Write_files import load_rawhmm, load_histonefile, load_states_file, load_output_file
from NucHMM_Load_Write_files import load_nuc_detail_file, load_genome_size_file, write_multi_statesfile_from_dict
from NucHMM_states_coverage import create_state_coverage_files, file2dict


# sort-unique-state block start
def binary_count(bin_string):
    """calculate how many 1 in the bin string. eg 0b10101 has 3 count"""
    count = 0
    for bits in bin_string:
        try:
            count += int(bits)
        except ValueError:
            continue
    return count

def bitwise_relevant(states_list, state_mark2histone_mark, env_states_list_mark):
    sum_count = 0
    for states in states_list:
        tmp_list = [state_mark2histone_mark[states]] * len(env_states_list_mark)
        # AND bitwise two list
        tmp_result = [a&b for a,b in zip(tmp_list,env_states_list_mark)]
        # transfer to bit_string list
        binary_tmp_result = [bin(i) for i in tmp_result]
        # use binary_count_transfer
        bicount_result = [binary_count(i) for i in binary_tmp_result]
        # only save the largest one
        try:
            if sum(bicount_result) > sum_count:
                sum_count = sum(bicount_result)
                final_picked_states = states
            else:
                continue
        except TypeError:
            print(bicount_result,sum(bicount_result),sum_count,key,diction[key])
            exit(1)
    return final_picked_states

def reference_relevant(key, states_list,state_mark2histone_mark,ref_diction):
    ref_mark = ref_diction[key]
    convert_list = []
    mark_skip = False
    picked_state = '0'
    for state in states_list:
        current_histm = state_mark2histone_mark[state]
        if current_histm == ref_mark:
            picked_state = state
            mark_skip = True
        bit_result = ref_mark&current_histm
        convert_bit_result = binary_count(bin(bit_result))
        convert_list.append(convert_bit_result)

    if not mark_skip:
        m = max(convert_list)
        m_index = [i for i, j in enumerate(convert_list) if j == m]
        # print(states_list,ref_mark, convert_list,m_index)
        # if len(m_index) == 1:
        #     picked_state = states_list[m_index[0]]
        # else:
        #     # arbitrarily pick the first one
        #     picked_state = states_list[m_index[0]]
        picked_state = states_list[m_index[0]]

    return picked_state

def pick_one(idx, key, sorted_key, diction, pick_one_window, state_mark2histone_mark, winmethod, ref_diction):
    '''pick up the state according to its environment with some window size
    :param idx: the index of current nucleosome
    :param key: the coordinate of current nucleosome
    :param sorted_key: the sorted coordinates of nucleosome ['chr1_1_2','chr2_1-2']
    :param diction: the dictionary which has key (nucleosome coordinate) and value (states)
    :return: picked states for current nucleosome
    '''
    list_len = len(sorted_key)

    # select all nucleosome in the window size
    if idx - pick_one_window < 0:
        key_env = sorted_key[slice(0,idx+pick_one_window+1)]
    elif (idx + 1 + pick_one_window) > list_len:
        key_env = sorted_key[slice(idx-pick_one_window,list_len)]
    else:
        key_env = sorted_key[slice(idx-pick_one_window,idx+pick_one_window+1)]

    env_states_list = []
    for keys in key_env:
        if keys != key:
            env_states_list += diction[keys]
    # pay attention to dont use shallow copy here
    states_list = diction[key][:]

    try:
        env_states_list_mark = [state_mark2histone_mark[states] for states in env_states_list]
    except KeyError:
        print('Check the background state')
        for keys in key_env:
            print(keys, diction[keys])
        exit(1)

    #filter other state. Rules: 1) if one is in the env while others not, pick this one
    # 2) if all states are in/not in the env, pick the most relevant one
    exist_mark = 0
    for state in states_list:
        if state in env_states_list:
            exist_mark += 1
            continue
        else:
            states_list.remove(state)

    if len(states_list) == 0:
        # both of them are not exist in the env
        states_list = diction[key][:]
        if winmethod:
            final_picked_states = bitwise_relevant(states_list,state_mark2histone_mark,env_states_list_mark)
        else:
            final_picked_states = reference_relevant(key, states_list,state_mark2histone_mark,ref_diction)

    elif len(states_list) == 1:
        # only one of state exist in the env
        final_picked_states = states_list[0]
    else:
        # all of states are exist in the env
        # all exist, outnumbered one win, if equal, follow bitwise relevant function'''
        tmp = dict(Counter(env_states_list))
        # check if there is equal situation and key is belong to states_list
        tmp_max1 = max(tmp.items(), key=operator.itemgetter(1,0))[0]
        tmp_max2 = max(tmp.items(), key=operator.itemgetter(1,1))[0]
        if winmethod:
            if tmp_max1 == tmp_max2 and tmp_max1 in states_list:
                final_picked_states = max(tmp.items(), key=operator.itemgetter(1))[0]
            else:
                final_picked_states = bitwise_relevant(states_list,state_mark2histone_mark,env_states_list_mark)
        else:
            final_picked_states = reference_relevant(key, states_list,state_mark2histone_mark,ref_diction)

    return final_picked_states

def statem_to_histm(histonelistfile,numstate,numoutput,emit_matrix):
    '''
    Generate statem2histm dictionary
    :param rawhmmfile:
    :param histonelistfile:
    :return: statem2histm: state_mark to histone mark dictionary e.g {'1':0,'2':24,'3':8,'4':17,'5':2,'6':65,
    '7':16,'8':5,'9':6,'10':32,'11':21,'12':64}. key is state ('1':state1) and value is histone mark value
    state4write: a dictionary with key: state and value: the histone marks list
    '''
    # 2 mean only those marks that are more than max_prob/2 can be consider as functional mark for this state
    mark_threshold = 2
    histone_list = load_histonefile(histonelistfile)

    # hist2value format {'H3K4me3':1,'H3K4me1':2...}
    # idx2hist format {0:''H3K4me3,1:'H3K3me1,...}
    hist2value = {hist:2**idx for idx, hist in enumerate(histone_list)}

    idx2hist = {idx:hist for idx,hist in enumerate(histone_list)}

    out = emit2mark(emit_matrix,numstate,numoutput,False)

    bar_max = out.max()
    # filtout_array format is np.array([[False, True],[True, False]])
    statem2histm = {}
    state4write = defaultdict(list)
    #
    filtout_array = out >=bar_max/mark_threshold

    for idy, state in enumerate(filtout_array):
        sum = 0
        for idx,value in enumerate(state):
            if value:
                sum += hist2value[idx2hist[idx]]
                state4write[idy+1].append(idx2hist[idx])
        statem2histm[str(idy+1)] = sum
    return statem2histm,state4write


def select_unique_state(rawhmmfile, histonelistfile, statefile, outputfile, filtstatefile, background_state, winmethod,
                        winsize, celltype):
    '''
    select the unique state for some multi-state nucleosome.
    Idea here is pick up the state according to its environment with some window size.
    Rules:
    1) if one is in the env while others not, pick this one.
    2) if all states are in/not in the env, pick the most relevant one
    :param statefile: states.bed file from hmm-second and hmm-calling subcommand
    :param filtstatefile: the output file path and name. This file is sorted and filtered. Each nucleosome in this file
    only have one state.
    :param background_state: list of background state. e.g ['1']
    :param winmethod: True/False
    :param winsize: window size for picking the unique state
    :param statem2hism: see statem_to_histm function
    :return: None
    '''
    trans_matrix,emit_matrix,numstate,numoutput = load_rawhmm(rawhmmfile)
    state_mark2histone_mark,state4write = statem_to_histm(histonelistfile,numstate,numoutput,emit_matrix)

    # print(state_mark2histone_mark)
    # output mark combination for each state
    print('State Corresponding Histone marks')
    for state in sorted(state4write):
        print('S'+str(state) + ' ' +','.join(state4write[state]))

    # window size for picking the unique state
    pick_one_window = winsize
    diction = load_states_file(statefile,background_state)
    # load ref output file if necessary
    ref_diction = None
    if not winmethod:
        if outputfile is None:
            reffile = celltype + '_output_secondr'  + '.bed'
            if file_check(reffile):
                query_mark = query_yes_no("Detected outputfile exists. Do you want to use this %s" % reffile)
                if query_mark:
                    ref_diction = load_output_file(reffile)
                else:
                    reffile = input_new_name('The name and path of the output reference file:')
                    ref_diction = load_output_file(reffile)
            else:
                reffile = input_new_name('The name and path of the output reference file:')
                ref_diction = load_output_file(reffile)
        else:
            ref_diction = load_output_file(outputfile)

    with open(filtstatefile, 'w') as outputfile:
        # sort the key first by chromosome, then start
        print("Sorting Step")
        sorted_key = sorted(diction, key=lambda x: (Chr2Num(x.partition('_')[0]),
                                              int(x.partition('_')[2].partition('_')[0])))
        out_line = 0
        print("Write merged file")
        for idx, keys in enumerate(sorted_key):
            K = keys.split('_')
            klen = len(diction[keys])
            out_line += 1
            sys.stdout.write('\rLine: ' + str(out_line))
            if klen != 1:
                # pick the most possible state according to the environment
                state = pick_one(idx,keys,sorted_key,diction,pick_one_window,state_mark2histone_mark,winmethod,ref_diction)
                diction[keys] = [state]
                outputfile.write(K[0] + '\t' + K[1] + '\t' + K[2] + '\t' + state + '\n')
            else:
                outputfile.write(K[0] + '\t' + K[1] + '\t' + K[2] + '\t' + diction[keys][0] + '\n')
        print('\n')
        print("Finish!")
    outputfile.close()

# sort-unique-state block end

# genomic_location_finder start
def create_genomic_seg_diction(upbound_distal,upbound_proximal,upbound_promoter,gene_rescale_length,sample_rate):
    '''
     Create a dictionary with key: Region name (e.g. Distal, Proximal) and
    the index of those region in States_x dictionary.
    :param upbound: Up boundary of the distal region
    :param upbound_distal: the distance upstream from TSS. e.g 100
    :param upbound_proximal:
    :param upbound_promoter:
    :param gene_rescale_length: rescaled gene body length
    :param sample_rate: unit bp/per point e.g 100 means 100bp we consider it as a point
    :return: a dictionary with key: Region name (e.g. Distal, Proximal) and
    the index of those region in States_x dictionary.
    '''
    # correspond to slice index usage slice(0,100) give point for index 0 to index 99
    distal_index1 = 0
    distal_index2 = (upbound_distal-upbound_proximal)//sample_rate
    proximal_index1 = distal_index2
    proximal_index2 = proximal_index1 + (upbound_proximal - upbound_promoter)//sample_rate
    promoter_index1 = proximal_index2
    promoter_index2 = promoter_index1 + upbound_promoter//sample_rate
    genebody_index1 = promoter_index2
    genebody_index2 = genebody_index1 + gene_rescale_length//sample_rate
    genebody5_index1 = genebody_index1
    genebody5_index2 = genebody5_index1 + (gene_rescale_length//sample_rate)//2
    genebody3_index1 = genebody5_index2
    genebody3_index2 = genebody_index2
    genomic_region_index = {'Distal':[distal_index1,distal_index2],'Proximal':[proximal_index1,proximal_index2],
                            'Promoter':[promoter_index1,promoter_index2],'Genebody':[genebody_index1,genebody_index2],
                            '5\'-Genebody':[genebody5_index1,genebody5_index2],
                            '3\'-Genebody':[genebody3_index1,genebody3_index2]}
    return genomic_region_index


def assign_rank_to_state(States_x_total,genome_region_index,region_name):
    '''
    Assign rank to each state in specific region
    :param States_x_total: see identify_genomic_location function
    :param genome_region_index: see create_genomic_seg_diction
    :param region_name: region name e.g Distal
    :return: a dictionary with key: state, value: rank. Notice that rank 1 means this state has the most coverage count
    '''
    # initial region dictionary
    region_dict = {state:0 for state in States_x_total}
    region_rank_dict = {state:0 for state in States_x_total}
    for state in States_x_total:
        bins_region = States_x_total[state][slice(genome_region_index[region_name][0],
                                                  genome_region_index[region_name][1])]
        # use sum as criterion
        region_dict[state] = sum(bins_region)
    # rank the states according to the sum of the coverage
    state_rank = 1
    for key, val in sorted(region_dict.items(),reverse=True, key=lambda item:item[1]):
        region_rank_dict[key] = state_rank
        state_rank += 1

    return region_rank_dict


def intra_state_comparison(state,States_x_total,genome_region_index):
    '''
    If state doesn't stand out from inter-state comparison, we perform this intra-state comparison to check the location
    it mainly distributed
    :param state: the state to check
    :param States_x_total: the dictionary with key:state and value the coverage count
    :param genome_region_index: the dictionary contains the region's index in the States_x_total values
    :return:
    '''
    num_region_return = 3
    region_sum = {}
    for region in genome_region_index:
        bins_region = States_x_total[state][slice(genome_region_index[region][0],genome_region_index[region][1])]
        # use mean as criterion because the regions have different length
        region_sum[region] = np.mean(bins_region)

    rank = 1
    region_result = []
    for key,val in sorted(region_sum.items(),reverse=True, key=lambda item:item[1]):
        if rank <= num_region_return:
            region_result.append(key+'_'+str(rank))
        else:
            break
        rank += 1

    return region_result


def identify_genomic_location(States_x_total,state4write,genome_region_index,Histone_location):
    '''
    Identify genomic location according to states distribution and empirical histone marks' location
    :param States_x_total: the dictionary contains the coverage counts of each states. Key: states; Value: coverage count.
    :param state4write: the dictionary contains the histone mark
    :return: a dictionary with key: state; and value: the list of locations. For a single state, the locations should be
    less than two but at least have one
    '''
    # The judgement comes from two aspect: 1) the distribution of the state in different region and
    # 2) empirical histone marks' location

    # for the 1): the one that is out-numbered (use the sum of the coverage as the criterion)
    # in the region will consider main distributed in this region
    # initialize region_coverage_sum_dict
    # allregion_rank_dict: key is region, value is a dictionary with key:state, value: rank
    allregion_rank_dict = {}
    for region in genome_region_index:
        allregion_rank_dict[region] = assign_rank_to_state(States_x_total,genome_region_index,region)

    rank_identfied_dict = {state:[] for state in States_x_total}
    for region in genome_region_index:
        # first round assignment
        for state in allregion_rank_dict[region]:
            # choose top 3 state as the enriched state for this region
            if allregion_rank_dict[region][state] <= 3:
                rank_identfied_dict[state].append(region+'_'+str(allregion_rank_dict[region][state]))


    # for those states that don't stand out by inter-state comparison, perform intra-state comparision
    for state in rank_identfied_dict:
        if len(rank_identfied_dict[state]) == 0:
            intra_comp_result = intra_state_comparison(state,States_x_total,genome_region_index)
            rank_identfied_dict[state] = intra_comp_result

    # empirical knowledge
    empirical_dict = {state:[] for state in States_x_total}
    for state in state4write:
        for histonem in state4write[state]:
            try:
                empirical_dict[state] += Histone_location[histonem]
            except KeyError:
                print(state4write)
                print(Histone_location)
                print(empirical_dict)
                print("Detect background state that enriched with histone modifications, please the re-input "
                      "correct background state.")
                exit(1)


    # combine rank information and empirical knowledge
    final_identified_dict = {state:[] for state in States_x_total}
    state_need_ref = []
    for state in rank_identfied_dict:
        # empirical mark. If all regions in rank_identified_dict not in empirical_dict, add empirical dict
        empirical_mark = 0
        for region_rank in rank_identfied_dict[state]:
            # if exist at both dictionary
            if region_rank.split('_')[0] in empirical_dict[state]:
                print('1',region_rank)
                final_identified_dict[state].append(region_rank)
                empirical_mark += 1
            else:
                # rank top 2 can be added
                if int(region_rank.split('_')[1]) <= 2:
                    print('2',region_rank)
                    final_identified_dict[state].append(region_rank)
        if empirical_mark == 0:
            state_need_ref.append(state)

    for state in state_need_ref:
        for region in empirical_dict[state]:
            if region not in final_identified_dict[state]:
                # avoid overlap region
                print('3',region)
                final_identified_dict[state].append(region)

    # final modifications for the dictionary
    print(final_identified_dict)
    max_keep_regions = 3
    for state in final_identified_dict:
        # the max number of a region should be 3, keep the top1 rank regions
        if len(final_identified_dict[state]) > max_keep_regions:
            for region_rank in final_identified_dict[state]:
                tmp_rank = int(region_rank.split('_')[1])
                if tmp_rank != 1:
                    final_identified_dict[state].remove(region_rank)

        for idx,region in enumerate(final_identified_dict[state]):
            final_identified_dict[state][idx] = region.split('_')[0]

        if '5\'-Genebody' in final_identified_dict[state] and '3\'-Genebody' in final_identified_dict[state]:
            final_identified_dict[state].remove('5\'-Genebody')
            final_identified_dict[state].remove('3\'-Genebody')
            if 'Genebody' not in final_identified_dict[state]:
                final_identified_dict[state].append('Genebody')
        elif '5\'-Genebody' in final_identified_dict[state] and 'Genebody' in final_identified_dict[state]:
            final_identified_dict[state].remove('Genebody')
            final_identified_dict[state].remove('5\'-Genebody')
            final_identified_dict[state].append('Genebody(5\'-Genebody)')
        elif '3\'-Genebody' in final_identified_dict[state] and 'Genebody' in final_identified_dict[state]:
            final_identified_dict[state].remove('Genebody')
            final_identified_dict[state].remove('3\'-Genebody')
            final_identified_dict[state].append('Genebody(3\'-Genebody)')
        else:
            continue

    return final_identified_dict


def sum_diction_plot(fileslist,interval_other_area,samplepoints,up_boundary,down_boundary,background_state,
                states,totalgenenum,States_x_total,plot_mark,celltype,spe_colors):
    '''
    Transfer mapped files to dictionary and add it to State_x_total
    :param fileslist: mapped files list
    :rest params: see genomic_loc_finder
    :return: State_x_total with key:state, and value: the accumulated coverage count
    '''
    States_x = file2dict(fileslist,interval_other_area,samplepoints,up_boundary,down_boundary,
                         background_state,states,totalgenenum)
    if plot_mark:
        plot_distribution(celltype,States_x,interval_other_area,up_boundary,down_boundary,samplepoints,spe_colors)

    for state in States_x:
        States_x_total[state] += States_x[state]
    return States_x_total

def output_plot_text(genesfile,file,mappedfiles,samplepoints,refgene,interval_other_area,up_boundary,down_boundary,
                     step1,step2,totalgenenum,celltype,spe_colors,states,background_state,States_x_total,
                     plot_mark,createfile_mark):
    '''
    Just to make the code concise. it contains two functions: 1) create the files for plotting and identifying if it is
    necessary. 2): accumulated the coverage count to State_x_total
    :param step1: see create_state_coverage_files step1
    :param step2: see create_state_coverage_files step2
    :param totalgenenum:
    :return:
    '''
    if createfile_mark:
        create_state_coverage_files(genesfile,file,samplepoints,refgene,interval_other_area,up_boundary,down_boundary,
                                    step1,step2)
    States_x_total_out = sum_diction_plot(mappedfiles,interval_other_area,samplepoints,up_boundary,down_boundary,
                                     background_state,states,totalgenenum,States_x_total,plot_mark,celltype,spe_colors)

    return States_x_total_out


def genomic_loc_finder(genesfile,filelist,rawhmmfile,histonelistfile,upDistal,up_proximal,up_promoter,
                       down_boundary,bgstate,outputfile,samplepoints,rescalelength,refgene,plot_mark,
                       plot_total_mark,removetmpfile):
    '''
    The main body of state-loc-finder
    :param genesfile: the file that contains all genes we use for selecting region
    :param stateslistfile: the merged.srt.unique.states.bed files list
    :param rawhmmfile: .rawhmm file from hmm-second
    :param histonelistfile: histone marks file contain histone marks with the same order of it hmm traninng peak file
    :param upDistal: up boundary of the Distal region (ref. from TSS)
    :param up_proximal: up boundary of the proximal region (ref. from TSS)
    :param up_promoter: up boundary of the promoter region (ref. from TSS)
    :param down_boundary: down boundary of TTS (still positive int)
    :param bgstate: the background state
    :param statenumber: the total number of states including background state
    :param outputfile: the path and name of the output state-region file
    :param samplepoints: sample_rate
    :param rescalelength: rescaled length of the gene
    :param refgene: reference gene size
    :param plot_mark: whether to plot the distribution pic for each cell type
    :param plot_total_mark: whether to plot the distribution pic for the summed coverage count for each state
    :param removetmpfile: whether to remove the middle files
    '''
    rescale_gene_length = rescalelength
    background_state = list(bgstate)
    '''sample rate 100 means for every 100 bp length, we sample 1 sample point'''
    sample_rate = rescalelength//samplepoints
    interval_other_area = sample_rate
    total_bins = (upDistal+down_boundary+rescale_gene_length)//sample_rate
    '''Count total gene number'''
    totalgenenum = count_file_rows(genesfile)

    trans_matrix,emit_matrix,numstate,numoutput = load_rawhmm(rawhmmfile)

    '''Create list of states'''
    states = []
    for i in range(numstate):
        if str(i+1) in background_state:
            continue
        else:
            states.append(i+1)

    States_x_total = {state:np.zeros(total_bins) for state in states}
    for file in filelist:
        '''file here is states.merged.srt.bed'''
        celltype = file.split('/')[-1].split('_')[0]
        print('\n')
        print(celltype)
        mapped_upbed = celltype+'_TSS_up_states.bed'
        mapped_genebody = celltype +'_genebody_states.bed'
        mapped_downbed = celltype + '_TTS_down_states.bed'
        # check if post-step2 file exists
        if file_check(mapped_upbed) and file_check(mapped_genebody) and file_check(mapped_downbed):
            mappedfiles = [mapped_upbed,mapped_genebody,mapped_downbed]
            States_x_total = output_plot_text(genesfile,file,mappedfiles,samplepoints,refgene,interval_other_area,upDistal,
                                     down_boundary,False,True,totalgenenum,celltype,spe_colors,states,background_state,
                                     States_x_total,plot_mark,False)
        else:
            query_mark = query_yes_no("Do you have "+mapped_upbed+'/'+mapped_genebody+'/'+mapped_downbed
                                      +'/'+" in other directory?")
            if query_mark:
                if sys.version_info > (3,0):
                    new_path = input( "Please input the directory(Folder path): ")
                else:
                    new_path = raw_input( "Please input the directory(Folder path): ")
                new_mapped_upbed = new_path + '/' + mapped_upbed
                new_mapped_genebody = new_path +  '/' + mapped_genebody
                new_mapped_downbed = new_path + '/' + mapped_downbed
                mappedfiles = [new_mapped_upbed,new_mapped_genebody,new_mapped_downbed]
                States_x_total = output_plot_text(genesfile,file,mappedfiles,samplepoints,refgene,interval_other_area,
                                                  upDistal,down_boundary,False,True,totalgenenum,celltype,spe_colors,
                                                  states,background_state,States_x_total,plot_mark,False)
            else:
                mappedfiles = [mapped_upbed,mapped_genebody,mapped_downbed]
                if file_check('./TSS_up_seg.srt.filt.bed') and file_check('./genebody_seg.srt.bed') and \
                        file_check('./TTS_down_seg.srt.filt.bed'):
                    print('Seg.srt.bed exists. Start from Step2')
                    '''start from step2'''
                    print(' Get mapped-state file...')
                    # createfile mark is True, step1 is False and step2 is True
                    States_x_total = output_plot_text(genesfile,file,mappedfiles,samplepoints,refgene,interval_other_area,
                                                      upDistal,down_boundary,False,True,totalgenenum,celltype,spe_colors,
                                                      states,background_state,States_x_total,plot_mark,True)
                else:
                    '''start from step1'''
                    print('Start from Step1')
                    # createfile mark is True, step1 is True and step2 is True
                    States_x_total = output_plot_text(genesfile,file,mappedfiles,samplepoints,refgene,interval_other_area,
                                                      upDistal,down_boundary,True,True,totalgenenum,celltype,spe_colors,
                                                      states,background_state,States_x_total,plot_mark,True)


    '''Loading mark-state matrix'''
    print('\n')
    print('Predicting genomic location for state...')
    state_mark2histone_mark,state4write = statem_to_histm(histonelistfile,numstate,numoutput,emit_matrix)
    genome_region_index = create_genomic_seg_diction(upDistal,up_proximal,up_promoter, rescale_gene_length, sample_rate)
    state_region_dict = identify_genomic_location(States_x_total,state4write,genome_region_index,Histone_location)
    # because the value in States_x is already the frequency, we need to reform it
    reform_total = np.zeros(((upDistal + down_boundary + rescale_gene_length) // sample_rate))
    for state in States_x_total:
        reform_total += States_x_total[state]
    for state in States_x_total:
        States_x_total[state] = States_x_total[state] / reform_total * 100
    if plot_total_mark == True:
        plot_distribution('Total',States_x_total,interval_other_area,upDistal,down_boundary,samplepoints,spe_colors)

    print('Writing results...')
    with open(outputfile,'w') as out_file:
        for state in states:
            out_file.write('S'+str(state)+'\t'+'\t'.join(state_region_dict[state])+'\n')
    out_file.close()
    print('Finish predicting!')

    if removetmpfile:
        print('Remove tmp files')
        subprocess.call("rm ./genebody_seg.srt.bed",shell=True)
        subprocess.call("rm ./TSS_up.bed",shell=True)
        subprocess.call("rm ./TSS_up_seg.srt.filt.bed",shell=True)
        subprocess.call("rm ./TTS_down_seg.srt.filt.bed",shell=True)


# genomic location finder end

# genomic location filter start
def create_ref_filter_files(genesfile,up_distal,up_proximal,up_promoter,refgene):
    '''
    create three basic files for following analysis
    :param genesfile:
    :param up_distal:
    :param up_proximal:
    :param up_promoter:
    :return:
    '''
    # step1 select the region
    distalfile = 'TSS_'+str(up_distal//1000)+'k_'+str(up_proximal//1000)+'k.bed'
    proximalfile = 'TSS_'+str(up_proximal//1000)+'k_'+str(up_promoter//1000)+'k.bed'
    promoterfile = 'TSS_'+str(up_promoter//1000)+'k_'+str(0)+'k.bed'
    print('Writing temporary region files..')
    with open(genesfile,'r') as gene_file, open(distalfile,'w') as distal_file, open(proximalfile,'w') as proximal_file,\
            open(promoterfile,'w') as promoter_file:
        for line in gene_file:
            line_info = line.strip().split()
            line_chr = line_info[0]
            line_start = int(line_info[1])
            line_end = int(line_info[2])
            line_strand = line_info[3]
            line_name = line_info[4]
            if line_strand == '+':
                distal_start = boundary_check(line_start - up_distal,line_chr, refgene)
                distal_end = boundary_check(line_start - up_proximal,line_chr, refgene)
                proximal_start = boundary_check(line_start - up_proximal, line_chr, refgene)
                proximal_end = boundary_check(line_start - up_promoter, line_chr, refgene)
                promoter_start = line_start - up_promoter
                promoter_end = line_start
                out_distal = [line_chr, str(distal_start), str(distal_end), line_strand, line_name]
                out_proximal = [line_chr, str(proximal_start), str(proximal_end), line_strand, line_name]
                out_promoter = [line_chr, str(promoter_start), str(promoter_end), line_strand, line_name]
                distal_file.write('\t'.join(out_distal)+'\n')
                proximal_file.write('\t'.join(out_proximal)+'\n')
                promoter_file.write('\t'.join(out_promoter)+'\n')
            else:
                distal_start = boundary_check(line_end + up_proximal,line_chr,refgene)
                distal_end = boundary_check(line_end + up_distal, line_chr, refgene)
                proximal_start = boundary_check(line_end + up_promoter, line_chr, refgene)
                proximal_end = boundary_check(line_end + up_proximal, line_chr, refgene)
                promoter_start = line_end
                promoter_end = line_end + up_promoter
                out_distal = [line_chr, str(distal_start), str(distal_end), line_strand, line_name]
                out_proximal = [line_chr, str(proximal_start), str(proximal_end), line_strand, line_name]
                out_promoter = [line_chr, str(promoter_start), str(promoter_end), line_strand, line_name]
                distal_file.write('\t'.join(out_distal)+'\n')
                proximal_file.write('\t'.join(out_proximal)+'\n')
                promoter_file.write('\t'.join(out_promoter)+'\n')
    gene_file.close()
    distal_file.close()
    proximal_file.close()
    promoter_file.close()

    # step2 filter
    distalfiltfile = 'TSS_'+str(up_distal//1000)+'k_'+str(up_proximal//1000)+'k.filt.bed'
    proximalfiltfile = 'TSS_'+str(up_proximal//1000)+'k_'+str(up_promoter//1000)+'k.filt.bed'
    promoterfiltfile = 'TSS_'+str(up_promoter//1000)+'k_'+str(0)+'k.filt.bed'
    # start from promoter filter file, order is important
    print('Filtering temporary region files..')
    promoter_filt_files = [genesfile]
    bedops_diff(promoterfile,promoter_filt_files,promoterfiltfile)
    proximal_filt_files = [genesfile,promoterfiltfile]
    bedops_diff(proximalfile,proximal_filt_files,proximalfiltfile)
    distal_filt_files = [genesfile,promoterfiltfile,proximalfiltfile]
    bedops_diff(distalfile,distal_filt_files,distalfiltfile)

    # remove unfiltered files
    subprocess.call("rm " + distalfile, shell=True)
    subprocess.call("rm " + proximalfile, shell=True)
    subprocess.call("rm " + promoterfile, shell=True)

# main function of genomic location filter
def genomic_loc_filter(nucstatefile,background_state,state_regions_file,refgene,celltype,genesfile,up_distal,up_proximal,up_promoter,
                       rmtmpfile):
    '''
    Use to filter nucleosome state file according the identified genomic location of each state
    :param nucstatefile: merged.srt.unique.state.bed
    :param bgstate: background_state
    :param state_regions_file: state region file which is the result of state_loc_finder
    :param refgene: reference genome size
    :param celltype: cell type
    :param genesfile: the txt file contain all genes we used
    :param up_distal: upstream boundary of the distal region (from TSS)
    :param up_proximal: upstream boundary of the proximal region (from TSS)
    :param up_promoter: upstream boundary of the promoter region (from TSS)
    :param rmtmpfile: True/False whether to remove temporary file
    :return: a list of file names which are the genomic location filtered file
    '''

    # step1: divide the merged.srt.unique.state.bed to state-specific.bed
    nucstate_dict = defaultdict(list)
    print('Reading nuclsoeome state file..')
    count_nuc = 0
    with open(nucstatefile,'r') as nucstate_file:
        for line in nucstate_file:
            count_nuc += 1
            line_info = line.strip().split()
            line_chr = line_info[0]
            line_start = line_info[1]
            line_end = line_info[2]
            # line_state should be same type with input background state
            line_state = line_info[3]
            value = '_'.join([line_chr,line_start,line_end])
            # the value list is ordered, background states have been removed
            if line_state not in background_state:
                nucstate_dict[line_state].append(value)
    nucstate_file.close()
    print('Writing separate state file..')
    # make the order from small number to large number, import for later on array information writing
    state_files_name = [celltype+'_state_'+ state + '.bed' for state in sorted(nucstate_dict.keys(), key=lambda x:int(x))]
    for file in state_files_name:
        with open(file, 'w') as out_file:
            current_state = file.split('.')[0].split('_')[-1]
            for nuc in nucstate_dict[current_state]:
                out_line = nuc.split('_')
                out_file.write('\t'.join(out_line)+'\t'+current_state+'\n')
        out_file.close()

    # step2: create basic region file
    create_ref_filter_files(genesfile,up_distal,up_proximal,up_promoter,refgene)

    # step3: create state-specific filter file
    # filtered file name should be same with the file name in create_ref_filter_files
    distalfiltfile = 'TSS_'+str(up_distal//1000)+'k_'+str(up_proximal//1000)+'k.filt.bed'
    proximalfiltfile = 'TSS_'+str(up_proximal//1000)+'k_'+str(up_promoter//1000)+'k.filt.bed'
    promoterfiltfile = 'TSS_'+str(up_promoter//1000)+'k_'+str(0)+'k.filt.bed'
    # in order to make the format consistent
    genefile4filt = 'genebody_for_filter.bed'
    subprocess.call("awk -v OFS='\t' '{print $1,$2,$3}' " + genesfile + " > " + genefile4filt, shell=True)
    print('Loading state-genomic-region file information..')
    files2merge = defaultdict(list)
    with open(state_regions_file,'r') as input_file:
        for line in input_file:
            line_info = line.strip().split()
            # line_state is a string, [1:] because 10,11,12
            line_state = line_info[0][1:]
            line_region = line_info[1:]
            for region in line_region:
                if 'Genebody' == region.split('(')[0] or 'Genebody' in region.split('-'):
                    files2merge[line_state].append(genefile4filt)
                elif region == 'Distal':
                    files2merge[line_state].append(distalfiltfile)
                elif region == 'Proximal':
                    files2merge[line_state].append(proximalfiltfile)
                elif region == 'Promoter':
                    files2merge[line_state].append(promoterfiltfile)
    input_file.close()

    print('Create state-specific filter region file..')
    filter_region_file_list = []
    for state in files2merge:
        filter_region_file = state + '_for_filter.bed'
        filter_region_file_list.append(filter_region_file)
        subprocess.call("cat " + ' '.join(files2merge[state]) + " > " + filter_region_file,shell=True)

    print('Genomic location filtering..')
    filtered_files_list = []
    count_gl_nuc = 0
    for state_file in state_files_name:
        current_state = state_file.split('.')[0].split('_')[-1]
        filter_region_file = current_state + '_for_filter.bed'
        filtered_file = celltype + '_state_' + current_state + '_gl_filt.bed'
        filtered_files_list.append(filtered_file)
        subprocess.call('bedtools intersect -wa -a ' + state_file + ' -b ' + filter_region_file + ' -u > ' + filtered_file,
                        shell=True)
        count_gl_nuc += count_file_rows(filtered_file)

    if rmtmpfile:
        for state_file in state_files_name:
            subprocess.call("rm "+ state_file,shell=True)
        for file in filter_region_file_list:
            subprocess.call("rm " + file, shell=True)

    return filtered_files_list,count_nuc,count_gl_nuc

# genomic location filter end

# array number filter start
def find_up_down_array_num(diction_array_num,down_ratio_limit,up_ratio_limit,show_distribute_mark):
    '''
    find the up and down range of array length of the state
    :param diction_array_num: the dictionary store all array information
    :param down_ratio_limit:
    :param up_ratio_limit:
    :return: up and down cutoff fot filtering
    '''
    total_array = 0
    weighted_total_array = 0
    weighted_num = 0
    down_number = 0
    upper_number = 0
    diction_stats = {}
    # count for total number of the array with certain number of nucleosomes
    for array_num in diction_array_num:
        total_array += len(diction_array_num[array_num])
        diction_stats[array_num] = len(diction_array_num[array_num])
    # calculate the upper and down cutoff
    sum_ratio = 0
    for array_num in sorted(diction_array_num.keys()):
        current_ratio = round(len(diction_array_num[array_num])/(total_array*1.0),3)*100
        sum_ratio += current_ratio
        if sum_ratio >= down_ratio_limit and down_number == 0:
            down_number = array_num
        if sum_ratio >= up_ratio_limit and upper_number == 0:
            upper_number = array_num
        if show_distribute_mark:
            current_ratio = round(current_ratio,3)
            sum_ratio = round(sum_ratio,3)
            print("array length is %d, current ratio is %f, and sum ratio is %f " % (array_num,current_ratio,sum_ratio))
    for array_num in sorted(diction_array_num.keys()):
        if array_num >= down_number:
            weighted_num += array_num * len(diction_array_num[array_num])
            weighted_total_array += len(diction_array_num[array_num])
    weighted_num /= weighted_total_array
    # print("Weighted average array number: %d" % weighted_num)
    # print("Array number range bewteen %d to %d" % (down_number,upper_number))
    return down_number,upper_number, weighted_num

def array_num_filter(gl_filtfile,gl_an_filtfile,down_ratio_limit,up_ratio_limit,cutoff_distance,show_distribute_mark):
    with open(gl_filtfile,'r') as input_file:
        # acquire stats information of this state
        # diction_array_num: key: array number and value: a dictionary with key: the count number of this certain length array,
        # and value: the chr, start, end information of this array. for example diction_array_num[3][3] represent
        # the 3rd(in second []) found array with length 3(in the first [])
        # compare to original codes, we don't consider the nuc_number here for minimizing the input file we need
        diction_array_num = {}
        count = 0
        L_last_chr = 0
        L_last_nuc_num = 0
        num_array_count = {}
        current_array = []
        for line in input_file:
            L = line.strip().split()
            L_chr = L[0]
            L_start = L[1]
            L_end = L[2]
            L_state = L[3]
            L_dyad = (int(L_start) + int(L_end))//2
            L_info = L_chr + '_' + L_start + '_' + L_end + '_' + L_state
            # initialize last parameters
            if count == 0:
                L_last_chr = L_chr
                L_last_dyad = L_dyad
                current_array.append(L_info)
                count += 1
                continue
            if L_last_chr == L_chr:
                if L_dyad - L_last_dyad <= cutoff_distance:
                    current_array.append(L_info)
                    L_last_dyad = L_dyad
                else:
                    array_len = len(current_array)
                    if array_len in num_array_count:
                        num_array_count[array_len] += 1
                    else:
                        num_array_count[array_len] = 0

                    if array_len in diction_array_num:
                        diction_array_num[array_len][num_array_count[array_len]] = current_array
                    else:
                        diction_array_num[array_len] = {}
                        diction_array_num[array_len][num_array_count[array_len]] = current_array
                    current_array = []
                    current_array.append(L_info)
                    L_last_dyad = L_dyad
            else:
                array_len = len(current_array)
                if array_len in num_array_count:
                    num_array_count[array_len] += 1
                else:
                    num_array_count[array_len] = 0

                if array_len in diction_array_num:
                    diction_array_num[array_len][num_array_count[array_len]] = current_array
                else:
                    diction_array_num[array_len] = {}
                    diction_array_num[array_len][num_array_count[array_len]] = current_array
                current_array = []
                current_array.append(L_info)
                L_last_chr = L_chr
    input_file.close()

    down_number,upper_number,weighted_number = find_up_down_array_num(diction_array_num,down_ratio_limit,up_ratio_limit,
                                                                      show_distribute_mark)
    if sys.version_info <(3,0):
        keys = diction_array_num.keys()
        print('Filtering the outlier arrays..')
        for array_num in keys:
            if array_num < down_number or array_num > upper_number:
                diction_array_num.pop(array_num)
    else:
        print('Filtering the outlier arrays..')
        # print("down number:",down_number,"upper number:",upper_number)
        for array_num in list(diction_array_num):
            if array_num < down_number or array_num > upper_number:
                diction_array_num.pop(array_num)

    tmpfile = 'tmp' + get_time() + '.txt'
    print('Writing the results..')
    with open(tmpfile,'w') as out_file:
        for array_num in diction_array_num:
            for array_sub_number in diction_array_num[array_num]:
                for info in diction_array_num[array_num][array_sub_number]:
                    nuc_list = info.split('_')
                    out_file.write(('\t').join(nuc_list)+'\n')
    out_file.close()

    subprocess.call("sort -k1,1V -k2,2n -k3,3n " + tmpfile + " > " + gl_an_filtfile, shell=True)
    print("Surviving Nucleosomes: " + str(count_file_rows(gl_an_filtfile)))
    subprocess.call("rm "+ tmpfile, shell=True)

    return down_number, upper_number, weighted_number
# array number filter end

# nucleosome positioning filter start
def nuc_pos_state_outlier(nucs_coor,states_list,positioning_score_list,nuc_states,up_ratio,down_ratio,k,showinfo):
    '''
    Screen out some nucleosome outliers according to the positioning
    :param positoning_socre_list:
    :param states_list:
    :return:
    '''

    # separate the states
    states_pos_score = defaultdict(list)
    for idx,state in enumerate(nuc_states):
        states_pos_score[state].append(positioning_score_list[idx])

    # find up and down cutoff for the state, pos_mean_dict is for writing mean of the positioning score of each state.
    # pos_mean_dict value is in same order with states_list
    pos_mean_dict = defaultdict(list)
    states_cutoff = defaultdict(list)
    for state in states_list:
        down_cut = np.quantile(states_pos_score[state],down_ratio)
        up_cut = np.quantile(states_pos_score[state],up_ratio)
        if showinfo:
            print("S%s poisitioning score range bewteen %3f to %3f" % (state, down_cut, up_cut))

        interval_qr = up_cut - down_cut
        up_limit = interval_qr * k + up_cut
        down_limit = down_cut - interval_qr * k
        states_cutoff[state].append(down_limit)
        states_cutoff[state].append(up_limit)

        state_worm_pos_mean = np.mean(states_pos_score[state])
        tmp_score = [value for value in states_pos_score[state] if value <= up_limit or value >= down_limit]
        state_wrm_pos_mean = np.mean(tmp_score)
        pos_mean_dict['Mean_wrm'].append(state_wrm_pos_mean)
        pos_mean_dict['Mean_worm'].append(state_worm_pos_mean)

        if showinfo:
            print("S%s mean without remove outlier %3f" % (state,state_worm_pos_mean))
            print("S%s mean with remove outlier %3f" % (state,state_wrm_pos_mean))

    # filtering step
    new_pos_score_list = []
    new_nuc_coor = []
    new_nuc_states = []
    for idx,pos_score in enumerate(positioning_score_list):
        current_state = nuc_states[idx]
        if pos_score > states_cutoff[current_state][1] or pos_score < states_cutoff[current_state][0]:
            continue
        else:
            new_pos_score_list.append(pos_score)
            new_nuc_coor.append(nucs_coor[idx])
            new_nuc_states.append(nuc_states[idx])

    return new_pos_score_list,new_nuc_coor,new_nuc_states,pos_mean_dict

def nuc_positioning_filter(gl_an_resp_filtfile,gl_an_pos_filtfile,states_list,nuc_detail_file,k,up_ratio,down_ratio,showinfo):

    # value: list[0]:Nucleosome_index; list[1]:Width_between_inflection; list[2]:Peak_height; list[3]:Area_under_curve
    # list[4]:Physical_property; list[5]:-log10(Pvalue_of_peak); list[6]:-log10(Pvalue_of_valley)
    nuc_detail_dict = load_nuc_detail_file(nuc_detail_file)

    nucs_width = []
    nucs_height = []
    nucs_area = []
    nucs_coef = []
    nucs_pvalpeak = []
    nucs_pvalvalley = []
    nuc_states = []
    nucs_coor = []
    # read the post-array-filtered file
    nuc_count = 0
    print('Read input file..')
    with open(gl_an_resp_filtfile, 'r') as input_file:
        for line in input_file:
            line_info = line.strip().split()
            line_chr = line_info[0]
            line_start = line_info[1]
            line_end = line_info[2]
            line_state = line_info[3]
            nuc_states.append(line_state)
            nucs_coor.append(line_chr+'_'+line_start+'_'+line_end)
            nuc_count += 1
            sys.stdout.write('\rRead Line:'+str(nuc_count))
            try:
                nuc_info = nuc_detail_dict[line_chr + '_' + line_start + '_' + line_end]
            except IndexError:
                print('Please check the order of the file in inputfileslist and nucdetailfilelist')
                exit(1)
            line_width = float(nuc_info[1])
            line_height = float(nuc_info[2])
            line_area = float(nuc_info[3])
            # actually not use it in the equation
            if nuc_info[4] == "MainPeak":
                line_coef = 1
            else:
                line_coef = 0.5
            line_ppeak = round(float(nuc_info[5]), 2)
            line_pvalley = round(float(nuc_info[6]), 2)
            nucs_width.append(line_width)
            nucs_height.append(line_height)
            nucs_area.append(line_area)
            nucs_coef.append(line_coef)
            nucs_pvalpeak.append(line_ppeak)
            nucs_pvalvalley.append(line_pvalley)
    input_file.close()

    print('\nIdentify the parameters\' outliers..')
    # find k times iqr value for normalizing parameters
    width_up_limit = outlier_threshold(nucs_width,k)
    height_up_limit = outlier_threshold(nucs_height, k)
    area_up_limit = outlier_threshold(nucs_area, k)
    pp_up_limit = outlier_threshold(nucs_pvalpeak, k)
    pv_up_limit = outlier_threshold(nucs_pvalvalley, k)

    # pvalley is not accurate when nucleosome at edge of an array of nucleosomes, this is important!
    nucs_pvalvalley = np.array(nucs_pvalvalley)
    nucs_pvalvalley[nucs_pvalvalley>pv_up_limit] = np.mean(nucs_pvalvalley)
    nucs_pvalvalley = nucs_pvalvalley.tolist()
    pv_up_limit = outlier_threshold(nucs_pvalvalley, k)

    # calculate positioning score by equation: height*log2(PP*PV+1)*sqrt(area)/(width**3)
    print('\nCalculating the positioning score..')
    print('Positioning score equation is (height+log2(PP*PV+1)+sqrt(area))/(width*3).')

    # transfer to nd-array for calculation and normalize it
    width_norm = np.array(nucs_width)/width_up_limit
    height_norm = np.array(nucs_height)/height_up_limit
    area_norm = np.array(nucs_area)/area_up_limit
    pvalpeak_norm = np.array(nucs_pvalpeak)/pp_up_limit
    pvalvalley_norm = np.array(nucs_pvalvalley)/pv_up_limit

    # calculate the positioning scores
    positioning_score = (height_norm + np.log2(pvalpeak_norm * pvalvalley_norm + 1) + area_norm)/(width_norm*3)

    # normalize positioning score to 20
    score_up_limit = outlier_threshold(positioning_score, k)
    positioning_score = (positioning_score/score_up_limit)*20
    print("Calculating Finish!")

    # filtering outlier nucleosomes
    print("Filtering outlier nucleosomes..")

    new_pos_score_list,new_nuc_coor,new_nuc_states,pos_mean_dict = nuc_pos_state_outlier(nucs_coor,states_list,
                                                                                         positioning_score,
                                                                                         nuc_states,up_ratio,
                                                                                         down_ratio,k,showinfo)


    print("Writing Output file..")
    nuc_filt_count = 0
    states_pos = defaultdict(list)
    with open(gl_an_pos_filtfile,'w') as output_file:
        for idx, coor in enumerate(new_nuc_coor):
            nuc_filt_count += 1
            sys.stdout.write('\rWrite Line:'+str(nuc_filt_count))
            coor_info = coor.split('_')
            line_chr = coor_info[0]
            line_start = coor_info[1]
            line_end = coor_info[2]
            line_state = new_nuc_states[idx]
            line_pos_score = round(new_pos_score_list[idx],3)
            states_pos[line_state].append(line_pos_score)
            out_line = [line_chr,line_start,line_end,line_state,str(line_pos_score)]
            output_file.write('\t'.join(out_line)+'\n')
    output_file.close()
    print('\n')
    return nuc_count,nuc_filt_count,pos_mean_dict,states_pos
# nucleosome positioning filter end

# nucleosome regularity and spacing filter start
def select_array_region(gl_an_filt_file,seleted_region_file,up_boundary,down_boundary):
    '''
    transfer individual nucleosomes file to nucleosomes array region file
    :param gl_an_filt_file: output file of the location_array_filter
    :param seleted_region_file: output file
    :return:
    '''
    with open(gl_an_filt_file,'r') as input_file, open(seleted_region_file, 'w') as out_file:
        count = 0
        nuc_count = 0
        for line in input_file:
            L = line.strip().split()
            L_chr = L[0]
            L_start = int(L[1])
            if count == 0:
                L_start_new = L_start + up_boundary
                L_end_new = L_start + down_boundary
                L_last_chr = L_chr
                nuc_count = 1
            else:
                if L_start > L_end_new and L_last_chr == L_chr:
                    out_file.write(L_chr+'\t'+str(L_start_new)+'\t'+str(L_end_new)+'\t'+str(nuc_count)+'\n')
                    L_start_new = L_start + up_boundary
                    L_end_new = L_start + down_boundary
                    L_last_chr = L_chr
                    nuc_count = 1
                else:
                    if L_last_chr != L_chr:
                        L_start_new = L_start + up_boundary
                        L_end_new = L_start + down_boundary
                        L_last_chr = L_chr
                        nuc_count = 1
                    else:
                        nuc_count += 1
                        continue
            count += 1
    input_file.close()
    out_file.close()


def cal_spacing(smoothed_y):
    '''
    function for calculating the spacing
    :param smoothed_y:
    :return:
    '''
    smooth_index = smoothed_y.index.values
    peaks,_ = signal.find_peaks(smoothed_y,distance=13,prominence=0.005)
    peaks = peaks + smooth_index[0]
    spacing = np.diff(peaks)
    sum = 0
    count = 0
    if len(spacing) <= 4:
        for idx,space_value in enumerate(spacing):
            if space_value < 25 :
                sum += space_value
                count += 1
            else:
                continue
    else:
        for idx,space_value in enumerate(spacing[:-1]):
            if space_value < 25 :
                sum += space_value
                count += 1
            else:
                continue
    current_ave_spacing = float(sum)/count
    return peaks, current_ave_spacing, spacing

def create_bin300_file(celltype,outfile_name,state,numcpu,like_wig_listfile,suffix):
    print("Create temporary bin.sh file..")
    tmp_binfile = "tmp_bin300.sh"
    f = open(tmp_binfile,'w')
    # get the celltype
    f.write("celltype="+celltype+'\n')
    # get the chromosome name from like_wig file name
    f.write("v=${1%.*}"+'\n')
    f.write("chr=${v##*_}"+'\n')
    # awk inside bedtools command is to create a file with this kind of format chr1 start end gaussian_smoothed_score
    # tail command it to skip the header. bedtools command find the intersect part between like_wig file and
    # gl_an_filt_file
    f.write("bedtools intersect -wa -a <( awk -v var=\"$chr\" 'BEGIN{OFS=\"\t\"}{print var,$1,$1,$3}' "
            "\"$1\"|tail -n +4) -b " + outfile_name + "|awk '{print $0,\"\t\",(NR-1)%301}'|awk '{A[$5]+=$4}END"
            "{for(key in A){print key\"\t\"A[key]}}' | sort -k1,1n - > \"$celltype\"_state_" + state + "_bin300_\"$chr\""+suffix+".txt")
    f.close()

    subprocess.call("cat " + like_wig_listfile + "| xargs -P " + str(numcpu) + " -n 1 bash " + tmp_binfile, shell=True)
    subprocess.call("rm " + tmp_binfile, shell=True)


def create_file_for_regs_spacing(gl_an_filt_file_path,inputfile_suffix,outputfile_suffix,coverfile_suffix,bin300_suffix,
                                 celltypelist,like_wig_files,numcpu,states_list,rmtmpfile):
     # from formal experience, the average array number in 99% case will less than 10.
    # In this case, we think 200*10 = 2000 bp range is big enough for following anaylsis
    region_up = -1000
    region_down = 2000
    remember_mark1 = True
    remember_mark2 = True
    for idx,celltype in enumerate(celltypelist):
        like_wig_listfile = like_wig_files[idx]
        for state in states_list:
            print("Writing %s state %s array region.." % (celltype,state))
            inputfile_name = gl_an_filt_file_path + '/' + celltype + '_state_' + state + inputfile_suffix
            print("Input file is %s" % inputfile_name)
            if not file_check(inputfile_name):
                inputfile_name = gl_an_filt_file_path + celltype + '_state_' + state + inputfile_suffix
                if not file_check(inputfile_name):
                    print('Please check if ' + inputfile_name + ' exist!')
                    exit(1)

            outfile_name = celltype + '_state_' + state + outputfile_suffix
            print("Output file is %s" % outfile_name)
            if file_check(outfile_name):
                if remember_mark1:
                    query_mark = query_yes_no(outfile_name + ' already exists, do you want to overwrite it?')
                    if query_mark:
                        select_array_region(inputfile_name,outfile_name,region_up,region_down)
                    else:
                        query_mark2 = query_yes_no('Save the choose for later state_region?')
                        if query_mark2:
                            remember_mark1 = False
                        print("Use the old %s " % outfile_name)
            else:
                select_array_region(inputfile_name,outfile_name,region_up,region_down)

            # check if all celltype_state_*_bin300_total.txt exist
            if not file_check(celltype + '_state_' + state + coverfile_suffix):
                bin300_chr_exist_mark = True
                # check if bin300 file exists, this is time consuming step
                for chrom in hg19:
                    check_file_name = celltype + '_state_' + state + '_bin300_' + chrom + bin300_suffix +'.txt'
                    bin300_chr_exist_mark = bin300_chr_exist_mark and file_check(check_file_name)

                if bin300_chr_exist_mark:
                    if remember_mark2:
                        query_mark = query_yes_no(celltype + '_state_' + state + '_bin300_chr*'+bin300_suffix+'.txt'
                                                  + ' files already exist, do you want to overwrite it?')
                        if query_mark:
                            create_bin300_file(celltype,outfile_name,state,numcpu,like_wig_listfile,bin300_suffix)
                        else:
                            query_mark3 = query_yes_no('Save the choose for later state_bin300 files?')
                            if query_mark3:
                                remember_mark2 = False
                            print("Use the exist %s" % celltype + '_state_' + state + '_bin300_chr*'+bin300_suffix+'.txt')
                else:
                    create_bin300_file(celltype,outfile_name,state,numcpu,like_wig_listfile,bin300_suffix)

                # sum all chromosome coverage for a state
                tmp_cal_total_file = "total_bin300_cal"+bin300_suffix+".sh"
                if file_check(tmp_cal_total_file):
                    print("Calculating coverage..")
                    subprocess.call("bash " + tmp_cal_total_file + " state_" + state + " " + celltype, shell=True)
                    # remove temporary file
                    if rmtmpfile:
                        subprocess.call("rm " + outfile_name, shell=True)
                        subprocess.call("rm " + celltype + "_state_" + state + "_bin300_chr*"+bin300_suffix+".txt", shell = True)
                else:
                    print("Writing %s" % tmp_cal_total_file)
                    f2 = open(tmp_cal_total_file,'w')
                    f2.write("paste  $2_\"$1\"_bin300_chr*"+bin300_suffix+".txt > $2_\"$1\"_total_raw"+bin300_suffix+".txt" + '\n')
                    f2.write("awk '{ for (i=1;i<=NF;i+=2) $i=\"\" }1' $2_\"$1\"_total_raw"+bin300_suffix+".txt |awk "
                             "'{print NR,\"\t\",$0}'|awk '{for(i=t=0;i<NF;) t+=$++i; $0=t}1'|awk '{print NR,\"\t\",$0}' > "
                             "$2_\"$1\"" + coverfile_suffix + '\n')
                    f2.write("rm $2_\"$1\"_total_raw"+bin300_suffix+".txt")
                    f2.close()
                    print("Writing Finish!")
                    print("Calculating coverage..")
                    subprocess.call("bash " + tmp_cal_total_file + " state_" + state + " " + celltype, shell=True)
                    if rmtmpfile:
                        subprocess.call("rm " + outfile_name, shell=True)
                        subprocess.call("rm " + celltype + "_state_" + state + "_bin300_chr*"+bin300_suffix+".txt", shell = True)
                print('\n')
            else:
                print("Use the old %s" % celltype + '_state_' + state + coverfile_suffix)

def cal_regularity_spacing(states_list,celltypelist,array_up,array_down,hz_up,hz_down,plotmark,regularity_method,
                           coverfile_suffix,plot_suffix):

    if regularity_method:
        reg_method = 'Max'
    else:
        reg_method = 'Mean'

    States = ['state_'+i for i in states_list]
    original2bin= {-1000:0, 0:98, 1000:198, 2000:298}
    all_info = {celltype:{} for celltype in celltypelist}

    for celltype in celltypelist:
        for state in States:
            # read coverage
            print("Reading file: " + celltype + '_' + state + coverfile_suffix)
            if sys.version_info < (3,0):
                tmp_data = pd.read_csv(celltype + '_' + state + coverfile_suffix, sep='\t',
                                                                header=None)[1]
            else:
                # remain modified
                tmp_data = pd.read_csv(celltype + '_' + state + coverfile_suffix, sep='\t',
                                                                header=None)[1]

            # normalize and re-center the coverage distribution
            # tmp_data = tmp_data/max(tmp_data) - np.mean(tmp_data/max(tmp_data))
            tmp_data = tmp_data/max(tmp_data)
            all_info[celltype][state] =  tmp_data


    # select the array region. range from -1000 to 2000
    original_up = array_up
    original_down = array_down
    # up 6 down 7 means only choose 6 Hz, 6-8 mean choose 6,7, 6~166bp 5~200bp'''
    Hz_range_up = hz_up
    Hz_range_down = hz_down

    region_up = original2bin[original_up]
    region_down = original2bin[original_down]
    selected_length = region_down - region_up

    selected_region_sum = {state:np.zeros(selected_length) for state in States}
    for celltype, cover_info in all_info.items():
        for state in States:
            selected_region_sum[state] += cover_info[state][region_up:region_down]


    ave_spacing_dict = {}
    detail_spacing_dict = {}
    legend = ['S'+state.split('_')[1] for state in States]
    if plotmark:
        plt.figure()
    # calculate the nucleosome spacing and plot the array distribution if needed
    print('Calculating nucleosome spacing..')
    for state in States:
        peaks,current_spacing,spacings = cal_spacing(selected_region_sum[state])
        ave_spacing_dict[state] = current_spacing * 10
        detail_spacing_dict[state] = spacings
        if plotmark:
            plt.plot(selected_region_sum[state],color=spe_colors[state.split('_')[1]])
            plt.scatter(peaks,selected_region_sum[state][peaks],marker="^")

    if plotmark:
        if region_up != 0:
            ticks_a = [ticks for ticks in range(region_up,region_down+100,100)]
            ticks_b = [str(ticks) for ticks in range(original_up, original_down+1000,1000)]
            plt.xticks(ticks_a, ticks_b)
        else:
            ticks_a = [ticks for ticks in range(98,region_down+100,100)]
            ticks_b = [str(ticks) for ticks in range(original_up, original_down+1000,1000)]
            ticks_a.insert(0,0)
            plt.xticks(ticks_a,ticks_b)

        plt.xlabel('bp')
        plt.ylabel('Normalized coverage')
        plt.legend(legend,loc='upper right')
        fig_name = 'Array_distribution_' + plot_suffix+ get_time() + '.png'
        plt.savefig(fig_name,dpi=300)
        plt.close()

    # calculate regularity score based on total distribution'''
    # frequence 4 to 10 range correspond 250 to 100 bp. In most case, nucleosomes spacing are in this range.
    frequency_range_lower = 4
    frequency_range_higher = 10
    fs = 1000
    regularity_score = []
    regularity_score_dict = {}

    if plotmark:
        bar_color = []
        plt.figure()
    for state in States:
        if len(detail_spacing_dict[state])<=3:
            data = selected_region_sum[state][:100]
            x = np.arange(0,100,1)
            x_new = np.arange(0,100,0.1)
        else:
            data = selected_region_sum[state]
            x = np.arange(0,selected_length,1)
            x_new = np.arange(0,selected_length,0.1)
        # data = selected_region_sum[state]
        tck = interpolate.splrep(x,data)
        y_bspline = interpolate.splev(x_new,tck)
        f, Pxx_den = signal.welch(y_bspline, fs, nperseg=1000)
        f2bp = 1000/f[frequency_range_lower:frequency_range_higher]

        if plotmark:
            plt.semilogy(f2bp,Pxx_den[frequency_range_lower:frequency_range_higher],color=spe_colors[state.split('_')[1]])
            bar_color.append(spe_colors[state.split('_')[1]])

        if reg_method == 'Max' or reg_method == 'max':
            regularity_score.append(max(Pxx_den[Hz_range_up:Hz_range_down]))
        else:
            regularity_score.append(np.mean(Pxx_den[Hz_range_up:Hz_range_down]))
    if plotmark:
        plt.legend(legend,loc='upper left')
        plt.xlabel('Periode (bp)')
        plt.ylabel('Spectral density')
        fig_name = 'State_Spectral_density_' + plot_suffix + get_time() + '.png'
        plt.savefig(fig_name,dpi=300)
        plt.close()

    # normalize regularity score to 1-21
    scalar_value = int(1/np.min(regularity_score))
    regularity_score = np.array(regularity_score)*scalar_value
    regularity_score = np.log(regularity_score)
    regularity_score /= np.max(regularity_score)
    regularity_score *= 20
    regularity_score += 1

    for idx,state in enumerate(States):
        regularity_score_dict[state] = regularity_score[idx]

    if plotmark:
        plt.figure()
        plt.bar(legend,np.array(regularity_score),color=bar_color)
        plt.ylabel('Regularity score')
        fig_name = 'Regularity_score_' + plot_suffix + get_time() + '.png'
        plt.savefig(fig_name,dpi=300)
        plt.close()

    return ave_spacing_dict, detail_spacing_dict, regularity_score_dict


def regularity_spacing_filter(gl_an_filt_file_path,celltypelist,like_wig_files_list,background_state,numcpu,statesnum,array_up,
                              array_down,hz_up,hz_down,rank_coefficient,plotmark,regularity_method,writeinfo,cutoff_dist,out_suffix,
                              rmtmpfile,NAstates):
    '''

    :param gl_an_filt_file_path: the path contains all genomic
    :param celltypelist:
    :param like_wig_filelist:
    :param bgstate:
    :param numcpu:
    :param statesnum:
    :param rmtmpfile:
    :return:
    '''

    # state in states_list is string type
    states_list = []
    for i in range(1,statesnum+1):
        if str(i) not in background_state:
            states_list.append(str(i))

    # remove background state and state with average number of nucleosomes in array less than 3
    # states_list = []
    # for i in range(1,statesnum+1):
    #     if str(i) not in background_state and str(i) not in NAstates:
    #         states_list.append(str(i))

    # step1 create the coverage files for anaylsis
    create_file_for_regs_spacing(gl_an_filt_file_path,'_gl_an_filt.bed','_gl_an_region.bed','_bin300_total.txt','',celltypelist,
                                 like_wig_files_list,numcpu,states_list,rmtmpfile)

    # step2 calculate the average regularity score and the average nucleosome spacing for each state
    ave_spacing_dict, detail_spacing_dict, regularity_score_dict = cal_regularity_spacing(states_list,celltypelist,
                                                                                          array_up,array_down,hz_up,hz_down,
                                                                                          plotmark,regularity_method,
                                                                                          '_bin300_total.txt','')

    # step3 filter the abnormal nucleosomes in array based on the spacing and regularity we calculated from step2
    # step3.1 rank the states by their regualrity score, the state with the smallest regularity score is rank1
    state_regularity_rank = {}
    rank = 1
    for key in sorted(regularity_score_dict, key=regularity_score_dict.get):
        state_regularity_rank[key] = rank
        rank += 1

    # step3.2 filtering according to the following rule the spacing between the nucleosomes
    # will be nucleosome spacing (result from step2 )  +/- interval * (5 + regularity rank * rank_coefficient) bp (inspired by paper :
    # "Fuzziness and Noise in nucleosomal architecture")
    nuc_count = 0
    nuc_cell_filtered_dict = {celltype:defaultdict(list) for celltype in celltypelist}
    for celltype in celltypelist:
        for state in states_list:
            tobefilter_file = celltype + '_state_' + state + '_gl_an_filt.bed'
            base_spacing = float(ave_spacing_dict['state_'+state])
            line_last_dyad = 0
            line_count = 0
            with open(tobefilter_file,'r') as input_file:
                for line in input_file:
                    nuc_count += 1
                    line_count += 1
                    line_info = line.strip().split()
                    line_chr = line_info[0]
                    line_start = int(line_info[1])
                    line_end = int(line_info[2])
                    line_dyad = (line_start + line_end) / 2
                    line_state = line_info[3]
                    key = line_chr + '_' + str(line_start) + '_' + str(line_end)
                    if line_count == 1:
                        line_last_dyad = line_dyad
                        array_count = 1
                        line_last_chr = line_chr
                        current_array_key = [key]
                        current_array_state = [line_state]
                        continue
                    current_spacing = line_dyad - line_last_dyad
                    if current_spacing <= cutoff_dist and line_last_chr == line_chr:
                        ref_spacing_low = base_spacing - array_count * ( 5 + state_regularity_rank['state_'+state] * rank_coefficient)
                        ref_spacing_high = base_spacing + array_count * (5 + state_regularity_rank['state_'+state] * rank_coefficient)
                        array_count += 1
                        line_last_dyad = line_dyad
                        line_last_chr = line_chr
                        if current_spacing < ref_spacing_low or current_spacing > ref_spacing_high:
                            continue
                        else:
                            current_array_key.append(key)
                            current_array_state.append(state)
                    # remove all the nucleosome in array
                    #     if current_spacing < ref_spacing_low or current_spacing > ref_spacing_high:
                    #         current_array_key.append(0)
                    #         current_array_state.append(0)
                    #         continue
                    #     else:
                    #         current_array_key.append(key)
                    #         current_array_state.append(line_state)
                    #         continue
                    else:
                        if 0 not in current_array_key:
                            for idx,state_key in enumerate(current_array_state):
                                nuc_cell_filtered_dict[celltype][state_key].append(current_array_key[idx])
                        current_array_state = [line_state]
                        current_array_key = [key]
                        line_last_dyad = line_dyad
                        line_last_chr = line_chr
                        array_count = 1

    # step4. re-calculate the regualrity and spacing. nuc_filtered_dict with key:state value:a list of chrom_coor1_coor2
    # step4.1 re-create coverage file
    nuc_filt_count = 0
    for celltype in celltypelist:
        nuc_filt_tmp_count = write_multi_statesfile_from_dict(nuc_cell_filtered_dict[celltype],out_suffix,celltype)
        nuc_filt_count += nuc_filt_tmp_count

    create_file_for_regs_spacing('.',out_suffix,'_gl_an_resp_region.bed','_bin300_resp_total.txt','_filt',
                                 celltypelist,like_wig_files_list,numcpu,states_list,rmtmpfile)

    filt_ave_spacing_dict, filt_detail_spacing_dict, filt_regularity_score_dict \
        = cal_regularity_spacing(states_list,celltypelist, array_up,array_down,hz_up,hz_down,plotmark,regularity_method,
                                 '_bin300_resp_total.txt','filt')

    info_dict = defaultdict(list)
    if writeinfo:
        # write regularity score and average spacing
        info_dict_name = 'Pre_Post_Regularity_Spacing_' + get_time() +'.txt'
        info_dict_cols = ['Pre.Ave. Spacing','Post.Ave. Spacing','Pre.Regularity Score','Post.Regularity Score']
        for state in states_list:
            info_dict['Pre.Ave. Spacing'].append(ave_spacing_dict['state_'+state])
            info_dict['Post.Ave. Spacing'].append(filt_ave_spacing_dict['state_'+state])
            info_dict['Pre.Regularity Score'].append(regularity_score_dict['state_'+state])
            info_dict['Post.Regularity Score'].append(filt_regularity_score_dict['state_'+state])
        df_info = pd.DataFrame(info_dict,index=['S'+state for state in states_list])
        df_info = df_info[info_dict_cols].round(2)
        df_info.to_csv(info_dict_name,sep='\t')

        # write detailed spacing information, might be deprecate in the formal version
        detail_spacing_name = 'Pre_post_Detail_Spacing_' + get_time() + '.txt'
        f_spacing = open(detail_spacing_name,'w')
        f_spacing.write('Pre-filter detail spacing information:\n')
        for state in states_list:
            output_number = [str(value) for value in detail_spacing_dict['state_'+state]]
            f_spacing.write('S'+state+'\t'+'\t'.join(output_number)+'\n')
        f_spacing.write('Post-filter detail spacing information:\n')
        for state in states_list:
            output_filt_number = [str(value) for value in filt_detail_spacing_dict['state_'+state]]
            f_spacing.write('S'+state+'\t'+'\t'.join(output_filt_number)+'\n')
        f_spacing.close()

    return nuc_count,nuc_filt_count,info_dict


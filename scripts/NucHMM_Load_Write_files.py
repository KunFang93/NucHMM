import sys
import pandas as pd
import numpy as np
from NucHMM_common_data import hg19
from NucHMM_utilities import ismember,sort_chrom_coor1_coor2,ShowProcess
from collections import defaultdict

def load_rawhmm(rawhmmfile):
    '''
    Load rawhmm file.
    rawhmm file is compose of four parts: the first line is numstate and numoutput (numobserve);
    the second part is transition matrix; the third part is emission matrix and the last line is the initial parameters.
    :param rawhmmfile: .rawhmm file path and name.
    :return: trans_matrix (python list format,[[]*numstates]);
    emit_matrix,numstate,numoutput
    '''
    inputfile = rawhmmfile
    trans_matrix = []
    emit_matrix = []
    numstate = 0
    numoutput = 0
    with open(inputfile,'r') as input_file:
        count = 0
        for line in input_file:
            L = line.strip().split()
            '''number of state and number of output'''
            if count == 0:
                numstate = int(L[0])
                numoutput = int(L[1])
                count += 1
                continue
            if count >=1 and count <= numstate:
                tmp = []
                for value in L:
                    tmp.append(round(float(value),3))
                trans_matrix.append(tmp)
                count += 1
                continue
            if count >= numstate + 1 and count <= 2 * numstate:
                tmp = []
                for value in L:
                    tmp.append(round(float(value),3))
                emit_matrix.append(tmp)
                count += 1
                continue

    print("Rawhmm File Loaded!")
    input_file.close()
    return trans_matrix,emit_matrix,numstate,numoutput

def load_histonefile(histonelistfile):
    '''
    Load histone mark file
    :param histonelistfile: the path and name of input histone marks file
    :return: histone_list (python list)
    '''
    histonefile = histonelistfile
    histone_list = []
    with open(histonefile,'r') as histone_file:
        '''The order of histones in histone file must be same with the inputlistfile for tfassign'''
        for line in histone_file:
            L = line.strip()
            histone_list.append(L)

    histone_file.close()
    return histone_list

def load_fq_list(fqlistfile):
    '''
    Load fastq list file
    :param fqlistfile: the path and name of input fastq file, Input format: for each line \'FILE READ-LEN SEQ-TYPE(chip/mnase) PEAK-TYPE(narrow?broad/none)\'
    :return: fqfile_list (python list), if pe_mark is True, list format is like [[F1_R1.fq,F1_R2.fq],[F2_R1,F2_R2.fq]]

    '''
    result_dict = defaultdict(list)
    with open(fqlistfile,'r') as fastq_file:
        for line in fastq_file:
            L = line.strip().split()
            if len(L) == 5:
                result_dict['fq_list'].append(L[:2])
                result_dict['read_len_list'].append(int(L[2]))
                result_dict['seq_type_list'].append(L[3])
                result_dict['peak_type_list'].append(L[4])
                result_dict['pe_mark_list'].append('PE')
            elif len(L) == 4:
                result_dict['fq_list'].append([L[0]])
                result_dict['read_len_list'].append(int(L[1]))
                result_dict['seq_type_list'].append(L[2])
                result_dict['peak_type_list'].append(L[3])
                result_dict['pe_mark_list'].append('SE')
            else:
                print("Wrong input format, Please check the input fastq file list format.")
                exit(1)
    fastq_file.close()
    return result_dict

def load_bam_list(bamfilelist):
    result_dict = defaultdict(list)
    with open(bamfilelist,'r') as bam_file:
        for line in bam_file:
            L = line.strip().split()
            result_dict['aligned_file_list'].append(L[0])
            result_dict['pe_mark_list'].append(L[1])
            result_dict['seq_type_list'].append(L[2])
            result_dict['peak_type_list'].append(L[3])
    bam_file.close()
    return result_dict

def load_states_file(statefile,background_state):
    with open(statefile,'r') as inputfile:
        diction = defaultdict(list)
        line_count = 1
        print("Read input file:")
        for line in inputfile:
            line_count += 1
            sys.stdout.write('\rLine: ' + str(line_count))
            L = line.strip().split()
            L_start = L[1]
            L_end = L[2]
            L_state = L[4]
            L_chr = L[0]
            str1 = [L_chr,L_start,L_end]
            key = "_".join(str1)
            if key in diction:
                if L_state in diction[key]:
                    continue
                else:
                    # skip background state if there is any functional states
                    if L_state in background_state:
                        continue
                    else:
                        diction[key].append(L_state)
            else:
                diction[key].append(L_state)
        print('\n')
    inputfile.close()
    return diction

def load_output_file(outputfile):
    '''
    Load the output file from hmm_calling subcommand
    :param outputfile: the path and name of the output file
    :return: a dictionary: key format 'chr1_1_10', value is its output.
    '''
    diction_output = {}
    count = 0
    print('\n')
    with open(outputfile) as out_file:
        for line in out_file:
            count += 1
            sys.stdout.write('\rLoading output file:'+str(count))
            line_info = line.strip().split()
            line_chr = line_info[0]
            line_start = line_info[1]
            line_end = line_info[2]
            line_output = int(line_info[4])
            key = '_'.join([line_chr,line_start,line_end])
            diction_output[key] = line_output
    out_file.close()
    print('\n')
    return diction_output

def load_genome_size_file(genomesizefile):
    '''
    Load genome size file
    :param genomesizefile: the path and name of the genome size file
    :return: the dictionary with key: chromosome (e.g 'chr1') and value: the size of chromosomes
    '''
    ref_dict = {}
    with open(genomesizefile, 'r') as input_file:
        for line in input_file:
            L = line.strip().split()
            ref_dict[L[0]] = int(L[1])
    input_file.close()
    return ref_dict

def load_nuc_position_file(LSEG, POS):
    '''return: pos_diction include index of certain position at begin, and end, specifical order is start_h, start_m
    end_h,end_m, e.g 69090:[start_h, start_m, end_h, end_m]'''
    pos_diction = defaultdict(list)
    MIX = {}
    Pos_Mix = {}
    Start_Pos = {}
    End_Pos = {}
    number_of_cell = len(POS)
    for keys in LSEG:
        tmp_seg = sorted(list(set(LSEG[keys])))
        #tmp_seg.remove(0)
        for pos_index in range(number_of_cell):
            #print('POS:'+str(len(POS))+'\t'+str(len(POS[pos_index]))+'\t'+str(index))
            try:
                MIX[pos_index] = sorted(tmp_seg + POS[pos_index][keys])
                Pos_Mix[pos_index] = ismember(tmp_seg, MIX[pos_index])
                sub_index = ismember(tmp_seg,tmp_seg)
                Start_Pos[pos_index] = list(map(lambda x, y: x - y, Pos_Mix[pos_index], sub_index))
                End_Pos[pos_index] = list(map(lambda x, y: x - y - 1, Pos_Mix[pos_index], sub_index))
            except IndexError:
                print('position diction IndexError'+'\n')
                continue

        for index2, i in enumerate(tmp_seg):
            for cellnum in range(number_of_cell):
                if str(i)+'_'+keys not in pos_diction:
                    pos_diction[str(i)+'_'+keys] = [0]*(2*number_of_cell)
                    '''keep '0_0':[0,0,0,0]'''
                if str(i)+'_'+keys != str(0)+'_'+keys:
                    pos_diction[str(i)+'_'+keys][cellnum] = Start_Pos[cellnum][index2]
                    pos_diction[str(i) + '_' + keys][cellnum+number_of_cell] = End_Pos[cellnum][index2]
    return pos_diction

def load_assignedtf_file(assignedfile):
    POS = defaultdict(list)
    READ = defaultdict(list)
    Start_End = defaultdict(list)
    count_input = 0
    for line in assignedfile:
        try:
            count_input += 1
            sys.stdout.write('\rFile line: ' + str(count_input))
            L = line.strip().split()
            L_start = int(L[1])
            L_end = int(L[2])
            '''minus one to match the later index'''
            L_mid = (L_start + L_end) // 2
            L_chr = L[0]
            L_read = L[3]
            # d_mid_count_hct[str(L_mid)+'_'+L_chr] = L_read
            Start_End[L_chr + '_' + str(L_mid)].append(L_start)
            Start_End[L_chr + '_' + str(L_mid)].append(L_end)
            POS[L_chr].append(L_mid)
            READ[L_chr].append(L_read)
        except IndexError:
            # skip the blank line in cis_element process that produced by paste..
            continue
    assignedfile.close()
    print('\n')
    return POS, READ, Start_End


def load_nuc_detail_file(nuc_detail):
    '''
    Load the nucleosome detail file, which is the result file of iNPS
    :param nuc_detail:
    :return: a dictionary with key: 'chr1_12_13' and value a list of it's detailed information
    '''

    print('Loading the nucleosome detailed file..')
    count_line = 1
    nuc_detail_dict = {}
    with open(nuc_detail,'r') as nuc_file:
        for line in nuc_file:
            sys.stdout.write('\rReading Line:' + str(count_line))
            count_line += 1
            line_info = line.strip().split()
            # skip the header
            try:
                int(line_info[4])
            except IndexError and ValueError:
                continue
            line_chr = line_info[0]
            line_start = line_info[1]
            line_end = line_info[2]
            line_detail = line_info[3:]
            key = line_chr + '_' + line_start + '_' + line_end
            nuc_detail_dict[key] = line_detail
    print('\nFinish Loading!')
    return nuc_detail_dict

def write_binnum_file(precomp_file,output):
    print("Writing binnum file..")
    with open(precomp_file,'r') as inputfile1, open(output, 'w+') as outputfile1:
        d = {}
        for line in inputfile1:
            L = line.strip().split()
            try:
                L_chr = int(L[0]) + 1
            except ValueError:
                print("Abnormal Chr Idx, Skipped")
                continue
            if L_chr in d:
                d[L_chr] += 1
            else:
                d[L_chr] = 1

        x = sorted([key for key in d])
        for j in x:
            if j == 23:
                outputfile1.write('chrX' + '\t' + str(d[j]) + '\n')
            elif j == 24:
                outputfile1.write('chrY' + '\t' + str(d[j]) + '\n')
            else:
                outputfile1.write('chr' + str(j) + '\t' + str(d[j]) + '\n')
    inputfile1.close()
    outputfile1.close()


def write_precomp_file(direction, dict_seg, seg_name, count_out, L4_mark, cell_num,
                   number_of_cell, chr_index, num_seg ,tmp_read,tmp_pos,outputfile):
    '''dict_seg[seg_name][2+2*number_of_cell+cell_num][index][0] is the start of this nuc and
    dict_seg[seg_name][2+2*number_of_cell+cell_num][index][1] is the end of this nuc'''

    try:
        if direction == 'for':
            diff = dict_seg[seg_name][2+cell_num] - max(dict_seg[seg_name][2:2+number_of_cell])
            if diff >= 0:
                for index, read in enumerate(dict_seg[seg_name][2+number_of_cell+cell_num]):
                    count_out += 1
                    outputfile.write(str(chr_index) + '\t' + str(count_out) + '\t' + str(read) + '\t'
                                      + str(dict_seg[seg_name][2+2*number_of_cell+cell_num][index][0]) + '\t' +
                                      str(dict_seg[seg_name][2+2*number_of_cell+cell_num][index][1]) + '\n')
            else:
                if L4_mark == '+':
                    for index, read in enumerate(dict_seg[seg_name][2+number_of_cell+cell_num]):
                        count_out += 1
                        outputfile.write(str(chr_index) + '\t' + str(count_out) + '\t' + str(read) + '\t'
                                          + str(dict_seg[seg_name][2+2*number_of_cell+cell_num][index][0]) + '\t' +
                                          str(dict_seg[seg_name][2+2*number_of_cell+cell_num][index][1]) + '\n')
                    for i in range(-diff):
                        count_out += 1
                        outputfile.write(str(chr_index) + '\t' + str(count_out) + '\t' + str(0) + '\n')
                else:
                    for i in range(-diff):
                        count_out += 1
                        outputfile.write(str(chr_index) + '\t' + str(count_out) + '\t' + str(0) + '\n')
                    for index, read in enumerate(dict_seg[seg_name][2+number_of_cell+cell_num]):
                        count_out += 1
                        outputfile.write(str(chr_index) + '\t' + str(count_out) + '\t' + str(read) + '\t'
                                          + str(dict_seg[seg_name][2+2*number_of_cell+cell_num][index][0]) + '\t' +
                                          str(dict_seg[seg_name][2+2*number_of_cell+cell_num][index][1]) + '\n')
        else:
            diff = dict_seg[seg_name][2 + cell_num] - max(dict_seg[seg_name][2:2 + number_of_cell])
            if num_seg == 1 or num_seg == 2:
                if diff >= 0:
                    dict_seg[seg_name][2+number_of_cell+cell_num].reverse()
                    dict_seg[seg_name][2+2*number_of_cell+cell_num].reverse()
                    tmp_read.append(dict_seg[seg_name][2+number_of_cell+cell_num])
                    tmp_pos.append(dict_seg[seg_name][2+2*number_of_cell+cell_num])
                else:
                    fill = []
                    fill2 = []
                    dict_seg[seg_name][2+number_of_cell+cell_num].reverse()
                    dict_seg[seg_name][2+2*number_of_cell+cell_num].reverse()
                    for k in range(-diff):
                        fill.append(0)
                        fill2.append('')
                    tmp_read.append(dict_seg[seg_name][2+number_of_cell+cell_num] + fill)
                    tmp_pos.append(dict_seg[seg_name][2+2*number_of_cell+cell_num] + fill2)

            if num_seg == 3:
                if diff >= 0:
                    dict_seg[seg_name][2+number_of_cell+cell_num].reverse()
                    dict_seg[seg_name][2+2*number_of_cell+cell_num].reverse()
                    for index, read in enumerate(dict_seg[seg_name][2+number_of_cell+cell_num]):
                        count_out += 1
                        outputfile.write(str(chr_index) + '\t' + str(count_out) + '\t' + str(read) + '\t'
                                          + str(dict_seg[seg_name][2+2*number_of_cell+cell_num][index][0]) + '\t' +
                                          str(dict_seg[seg_name][2+2*number_of_cell+cell_num][index][1]) + '\n')
                else:
                    for i in range(-diff):
                        # print(count_out)
                        count_out += 1
                        outputfile.write(str(chr_index) + '\t' + str(count_out) + '\t' + str(0) + '\n')
                    dict_seg[seg_name][2+number_of_cell+cell_num].reverse()
                    dict_seg[seg_name][2+2*number_of_cell+cell_num].reverse()
                    for index, read in enumerate(dict_seg[seg_name][2+number_of_cell+cell_num]):
                        try:
                            count_out += 1
                            outputfile.write(str(chr_index) + '\t' + str(count_out) + '\t' + str(read) + '\t'
                                              + str(dict_seg[seg_name][2+2*number_of_cell+cell_num][index][0]) + '\t' +
                                              str(dict_seg[seg_name][2+2*number_of_cell+cell_num][index][1]) + '\n')
                        except IndexError:
                            print(str(seg_name)+'\t'+str(dict_seg[seg_name])+'\t'+str(index)+'\t'+str(cell_num)+
                                  '\t'+str(number_of_cell)+'\n')
                            exit(1)

                if len(tmp_read) == 0:
                    print('\n'+seg_name+'\n')
                    print(tmp_read)

                for index, j in enumerate(tmp_read[1]):
                    if tmp_pos[1][index] == '':
                        count_out += 1
                        outputfile.write(str(chr_index) + '\t' + str(count_out) + '\t' + str(0) + '\n')
                    else:
                        count_out += 1
                        outputfile.write(str(chr_index) + '\t' + str(count_out) + '\t' + str(j) + '\t'
                                         + str(tmp_pos[1][index][0]) + '\t' +
                                          str(tmp_pos[1][index][1]) + '\n')

                for index, j in enumerate(tmp_read[0]):
                    if tmp_pos[0][index] == '':
                        count_out += 1
                        outputfile.write(str(chr_index) + '\t' + str(count_out) + '\t' + str(0) + '\n')
                    else:
                        count_out += 1
                        outputfile.write(str(chr_index) + '\t' + str(count_out) + '\t' + str(j) + '\t'
                                          + str(tmp_pos[0][index][0]) + '\t' +
                                          str(tmp_pos[0][index][1]) + '\n')

                tmp_read = []
                tmp_pos = []
    except IndexError:
        outputfile.write('indexerror'+'\n')

    return tmp_read,tmp_pos,count_out


def write_segment_file(sorted_key,diction,outputfile):
    '''
    Write the segment file (used in cal_state_coverage)
    :param sorted_key:
    :param diction:
    :return:
    '''
    with open(outputfile, 'w') as out_file:
        for key in sorted_key:
            key_out = key.split('_')
            out = key_out + diction[key]
            out_file.write('\t'.join(out)+'\n')

    out_file.close()


def write_multi_statesfile_from_dict(diction,outputfile_suffix,cell_type):
    '''notice the diction should have key:state and value: a list of chrom_coor1_coor2'''
    nuc_count = 0
    for state in diction:
        # sort the key
        try:
            sorted_key = sort_chrom_coor1_coor2(diction[state])
        except ValueError:
            print(diction[state])
        filename = cell_type + '_state_' + state + outputfile_suffix
        f = open(filename,'w')
        for key in sorted_key:
            nuc_count += 1
            line_info = key.split('_')
            outline = line_info[:3] + [state] + [line_info[-1]]
            f.write('\t'.join(outline)+'\n')
        f.close()
    return nuc_count

def write_prep_resultlist(fq_info_dict,outfilename):
    with open(outfilename,'w') as out_file:
        for file in fq_info_dict['final_file_list']:
            out_file.write(file+'\n')
    out_file.close()

# def finalize_output(gl_an_resp_filt_file,cutoffdist):
#     # merge individual nucleosoem state to nucleosome state region in gl_an_resp_filt_file
#     print("Loading..")
#     input_df = pd.read_csv(gl_an_resp_filt_file,sep='\t',header=None)
#     print("Processing..")
#     input_df.loc[:,'dyad'] =(input_df.loc[:,1] + input_df.loc[:,2])/2
#     input_df.loc[:,'spacing'] = input_df.loc[:,'dyad'].diff()
#     # The spacing of the 1st nucleosome in a chromosome should be same as the spacing value of the 2nd nucleosome
#     input_df = input_df.fillna(0)
#     input_df.loc[input_df.loc[:,'spacing'] < 0,'spacing'] = \
#         list(input_df.loc[input_df.loc[input_df.loc[:,'spacing'] < 0,'spacing'].index.values + 1,'spacing'])
#     input_df.loc[0,'spacing']= input_df.loc[1,'spacing']
#     split_idx = input_df.loc[input_df.loc[:,'spacing'] > cutoffdist,:].index.values
#     input_df_split = np.split(input_df,split_idx)
#     final_dict = {'chrom':[],'start':[],'end':[],'state':[],'Ave.Local.pos':[],'Ave.Local.spa':[]}
#     for idx, split_df in enumerate(input_df_split):
#         sys.stdout.write("\r{}/{}".format(idx,len(input_df_split)))
#         if idx == 0:
#             # always empty df
#             continue
#         else:
#             gb = split_df.groupby(3)
#             states_array = [gb.get_group(x) for x in gb.groups]
#             for idy,state_array in enumerate(states_array):
#                 current_chrom = state_array.loc[state_array.index[0],0]
#                 current_start = state_array.loc[state_array.index[0],1]
#                 current_end = state_array.loc[state_array.index[-1],2]
#                 current_state = state_array.loc[state_array.index[0],3]
#                 current_pos = state_array[4].mean()
#                 if idy == 0:
#                     current_spa = state_array['spacing'][1:].mean()
#                 else:
#                     current_spa = state_array['spacing'].mean()
#                 final_dict['chrom'].append(current_chrom)
#                 final_dict['start'].append(current_start)
#                 final_dict['end'].append(current_end)
#                 final_dict['state'].append(current_state)
#                 final_dict['Ave.Local.pos'].append(current_pos)
#                 final_dict['Ave.Local.spa'].append(current_spa)
#     final_df = pd.DataFrame.from_dict(final_dict)
#     final_df.loc[:,'Ave.Local.pos'] = np.round(final_df.loc[:,'Ave.Local.pos'])
#     final_df.loc[:,'Ave.Local.spa'] = np.round(final_df.loc[:,'Ave.Local.spa'])
#     return final_df

def write_final_array(file,outname):
    '''
    merged file cols: chr start end state Ave.Local.pos Ave.Local.spa
    :param file: 
    :param outname: 
    :return: 
    '''
    with open(file,'r') as input_file, open(outname, 'w') as out_file:
        for line in input_file:
            line_info = line.strip().split()
            line_chr = line_info[0]
            line_start = line_info[1]
            line_end = line_info[2]
            line_state = line_info[3]
            line_pos = float(line_info[4])
            line_spa = float(line_info[5])
            line_array_mark = line_info[6]
            if line_array_mark == 'single':
                # single type don't have spacing information
                outline = [line_chr,line_start,line_end,line_state,str(line_pos),-1]
                out_file.write('\t'.join(outline) +'\n')
            elif line_array_mark == 'array-start':
                array_pos = [line_pos]
                array_spa = [line_spa]
                array_chr = line_chr
                array_start = line_start
                array_state = line_state
            elif line_array_mark.split('-')[-1] != 'end':
                array_pos.append(line_pos)
                array_spa.append(line_spa)
            elif line_array_mark.split('-')[-1] == 'end':
                array_end = line_end
                array_pos.append(line_pos)
                array_spa.append(line_spa)
                pos_ave = np.round(np.mean(array_pos),3)
                spa_ave = np.round(np.mean(array_spa),3)
                outline = [array_chr,array_start,array_end,array_state,str(pos_ave),str(spa_ave)]
                out_file.write('\t'.join(outline) +'\n')
            else:
                print("Cannot recognize array mark")
                print('\t'.join(outline) +'\n')
    input_file.close()
    out_file.close()

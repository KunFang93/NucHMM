#!/usr/bin/env python
#--coding:utf-8 --

import sys
import subprocess
import pandas as pd
import numpy as np
from collections import defaultdict
from NucHMM_Load_Write_files import load_genome_size_file, write_precomp_file,load_assignedtf_file

def hmm(config, refgenome, rawhmmfile, b, i, numstates, numoutputs, precompfiles, bicfile, binoutfile):
    """HMM model for training nucleosome states, typically for the first run"""
    print("Running HMM model")
    software_name = 'NucHMM-learn'
    if b:
        BIC = ' -b '
    else:
        BIC= ''

    subprocess.call(config.hmm_directory + '/' + software_name + BIC +' -i ' + str(i) + ' -h '+ rawhmmfile +
                ' -v ' + binoutfile + ' -g ' + refgenome
                +' '+ str(numstates) + ' ' + str(numoutputs) + ' ' + (' ').join(precompfiles) + ' > '+ bicfile
                ,shell=True)

def hmm_calling(config, refgenome, state_outfile, output_outfile, precompfile, rawhmminput, binoutfile):
    """hmm_calling use Viterbi algorithm to calls the states
    for each nucleosomes within the effective region in each cell type (in state file), and outputs file contains
    raw mark information.
    Usage: NucHMM --hmm-directory=path hmm-calling ...]"""

    software_name = 'NucHMM-output_results'
    subprocess.call(config.hmm_directory + '/' + software_name + ' ' + refgenome + ' ' + binoutfile + ' ' +
                    state_outfile + ' ' + output_outfile + ' ' + precompfile + ' ' + rawhmminput, shell=True)


def modifyhmm(statefiles,rawhmm,outrawhmm,num,nucnum,percentage,outstats):
    '''remove redundant state and equally distribute its probability'''

    # read rawhmm file
    count_raw = 0
    rawhmm_dict = {}
    print('Loading .rawhmm file...')
    with open(rawhmm,'r') as rawhmm_file:
        for line in rawhmm_file:
            line_info = line.strip().split()
            count_raw += 1
            if count_raw == 1:
                numstate = int(line_info[0])
                numobserve = int(line_info[1])
                continue
            rawhmm_dict[count_raw] = line_info
    trans_key = [count for count in range(2,2+numstate)]
    emit_key = [count for count in range(2+numstate,2+2*numstate)]
    init_prob_key = 2 + 2*numstate
    print('Loading Completed!')
    rawhmm_file.close()

    states_list = ['S'+str(state) for state in range(1,numstate+1)]
    # read state bed for each cell type
    states_dict = {}
    for statefile in statefiles:
        count_line = 1
        celltype = statefile.split('/')[-1].split('_')[0]
        states_dict[celltype] = [0] * numstate
        with open(statefile,'r') as state_file:
            for line in state_file:
                sys.stdout.write('\rReading States file:' + str(count_line))
                count_line += 1
                line_info = line.strip().split()
                line_state = int(line_info[4])
                states_dict[celltype][line_state - 1] += 1
        state_file.close()
        print('\n')

    df_out = pd.DataFrame(states_dict,index=states_list)
    if outstats is not None:
        df_out.to_csv(outstats,sep='\t')

    # remove redundant states
    total_state_num = df_out.sum(axis=1).tolist()
    rm_states_list = []
    rm_info = []
    if num:
        print('Use number method!')
        for idx,number in enumerate(total_state_num):
            if number <= nucnum:
                # e.g if idx=0, value<350, we need to remove the state1.
                rm_states_list.append(idx)
                rm_info.append(number)
    else:
        print('Use percentage method!')
        total_nuc = sum(total_state_num)
        print(total_nuc)
        for idx,number in enumerate(total_state_num):
            number = float(number)
            if number/total_nuc <= percentage:
                # if idx=0, value<350, we need to remove the state1.q
                rm_states_list.append(idx)
                rm_info.append(round(number/total_nuc,5))

    print('Remove States:')
    for idx,state in enumerate(rm_states_list):
        print('S'+str(state+1)+':'+str(rm_info[idx]))

    # modify the rawhmm file
    # modify the transition matrix by deleting the redundant state and distrbute its probability evenly to other states
    if outrawhmm is None:
        dt_string = get_time()
        out_rawhmm = 'modihmm_' + dt_string + '.rawhmm'
    else:
        out_rawhmm = outrawhmm

    f = open(out_rawhmm,'w')
    f.write(' '.join([str(numstate-len(rm_states_list)),str(numobserve)])+'\n')
    new_rawhmm_dict = defaultdict(list)
    rm_trans_key = np.array(rm_states_list) + 2
    for key in trans_key:
        if key in rm_trans_key:
            continue
        redun_prob = 0.0
        for idx,value in enumerate(rawhmm_dict[key]):
            value = float(value)
            if idx in rm_states_list:
                redun_prob += value
            else:
                new_rawhmm_dict[key].append(value)
        redist_prob = redun_prob/(numstate-len(rm_states_list))
        new_rawhmm_dict[key] = (np.array(new_rawhmm_dict[key])+redist_prob).tolist()
        f.write(' '.join(str(prob_new) for prob_new in new_rawhmm_dict[key])+'\n')
    # remove redundant state in emit_key
    rm_emit_key = np.array(rm_states_list) + 2 + numstate
    for key in emit_key:
        if key in rm_emit_key:
            continue
        else:
            new_rawhmm_dict[key] = rawhmm_dict[key]
        f.write(' '.join(str(prob_new) for prob_new in new_rawhmm_dict[key])+'\n')
    # remove redundant state in init_prob
    redun_prob_init = 0
    for idx,value in enumerate(rawhmm_dict[init_prob_key]):
        value = float(value)
        if idx in rm_states_list:
            redun_prob_init += value
        else:
            new_rawhmm_dict[init_prob_key].append(value)
    redist_prob_init = redun_prob_init/(numstate-len(rm_states_list))
    new_rawhmm_dict[init_prob_key] = (np.array(new_rawhmm_dict[init_prob_key])+redist_prob_init).tolist()
    f.write(' '.join(str(prob_new) for prob_new in new_rawhmm_dict[init_prob_key])+'\n')
    f.close()

def hmm_second(config, refgenome, out_rawhmm, b, i, numstates, numoutputs,
             precompfiles, bic_file, rawhmminput,binoutfile):
    """HMM model for training nucleosome states, typically for the second run"""
    software_name = 'NucHMM-learn'
    if b:
        BIC = ' -b '
    else:
        BIC= ''

    # read numstate and numoutput if there is no input
    if numstates is None or numoutputs is None:
        with open(rawhmminput,'r') as rawhmm_file:
            for line in rawhmm_file:
                line_info = line.strip().split()
                numstate = int(line_info[0])
                numobserve = int(line_info[1])
                break
        rawhmm_file.close()
    else:
        numstate = numstates
        numobserve = numoutputs

    subprocess.call(config.hmm_directory + '/' +software_name + BIC +' -i ' + str(i) + ' -h '+ out_rawhmm +  ' -H ' + rawhmminput +
                ' -v ' + binoutfile + ' -g ' + refgenome
                +' '+ str(numstate) + ' ' + str(numobserve) + ' ' + (' ').join(precompfiles) + ' > '+ bic_file
                ,shell=True)




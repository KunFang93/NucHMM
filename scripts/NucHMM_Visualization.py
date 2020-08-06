#!/usr/bin/env python
#--coding:utf-8 --

import os.path
import numpy as np
import pandas as pd
import matplotlib.ticker as mtick
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors
from datetime import datetime
from NucHMM_states_coverage import file2dict
from NucHMM_Load_Write_files import load_rawhmm,load_histonefile
from NucHMM_utilities import file_check, query_yes_no

def plot_state_coverage(genenumber,celltype,mapped_upbed,mapped_genebed,mapped_downbed,interval_other_area,up_boundary,
                        down_boundary,samplepoints,mapped_pic,spe_colors,states,background_state):
    '''
    Plot state coverage from three files
    :param genenumber:
    :param celltype:
    :param mapped_upbed:
    :param mapped_genebed:
    :param mapped_downbed:
    :param interval_other_area:
    :param up_boundary:
    :param down_boundary:
    :param samplepoints:
    :param mapped_pic:
    :param spe_colors:
    :param states:
    :param background_state:
    :return:
    '''
    print('Plotting...'+'\n')
    if not os.path.exists(mapped_upbed):
        print("Missing TSS_up_states.bed")
        exit(1)
    if not os.path.exists(mapped_genebed):
        print("Missing genebody_states.bed")
        exit(1)
    if not os.path.exists(mapped_downbed):
        print("Missing TTS_down_states.bed")
        exit(1)
    drawfiles = [mapped_upbed,mapped_genebed,mapped_downbed]
    '''draw distribution'''
    up_sample = up_boundary // interval_other_area
    down_sample = down_boundary // interval_other_area
    gb_sample = samplepoints
    total_bin = up_sample + down_sample + gb_sample
    States_x = file2dict(drawfiles,interval_other_area,samplepoints,up_boundary,down_boundary,
                         background_state,states,genenumber)

    for state in States_x:
        States_x[state] = States_x[state].tolist()

    plt.figure(figsize=(20,15))
    for state in States_x:
    #    x = smooth(x,35)
        pic = plt.plot(States_x[state],label='state'+str(state),color=spe_colors[str(state)])
        plt.title(celltype)
        plt.ylabel('Frequency')
        plt.xticks([0,up_sample,up_sample+gb_sample,total_bin],
                   ['-'+str(up_boundary//1000)+'kb','TSS','TTS',str(down_boundary//1000)+'kb'])

        fmt='%.2f%%'
        yticks = mtick.FormatStrFormatter(fmt)
        plt.gca().yaxis.set_major_formatter(yticks)
        plt.legend()
        plt.savefig(mapped_pic,bbox_inches='tight',dpi=300)
    plt.close()
    print('\n')
    print('Finished!')

def Marknum(mark):
    '''very strong assumption that probability is equally contributed by each mark,
    Detailly, we divide probability by the number or mark.
    e.g. prob of 1010(10): 0.5, for 8 we count 0.25, for 2 we count 0.25'''
    sum = 0
    for i in bin(mark)[2:]:
        sum += int(i)
    return sum


def emit2mark(emitmatrix,numstates,numoutputs,Mark):
    '''
    transfer emit matrix to mark state matrix
    :param emitmatrix: list with emission probability in it. e.g [[0.1, 0.2],[0.3, 0.4]]
    :param numstates: number of states
    :param numoutputs: number of observations
    :param Mark: see definition in Marknum
    :return: out, np.array matrix with dimension numstates * number of histone mark (log2(numoutputs))
    '''
    emit_matrix = np.array(emitmatrix)
    g = int(np.log2(numoutputs))
    out = np.zeros((numstates,g))
    for m in range(int(np.log2(numoutputs))):
        for o in range(numoutputs):
            if 2**m & o > 0:
                #print(m,o,2**m & o)
                if Mark:
                    marknum = Marknum(o)
                    out[:,m] = out[:,m] + emit_matrix[:,o]/marknum
                else:
                    out[:,m] = out[:,m] + emit_matrix[:,o]
    return out


def rawhmm2matrix(rawhmmfile,histonelistfile,transmat,markstatemat,mark,output_mark):
    '''
    rawhmm file transfer to mark-state matrix
    :param rawhmmfile: input .rawhmm file
    :param histonelistfile: the file have histone marks, one for each row.
    Importantly, The order of histone marks in histone file must be same with the inputlistfile for tfassign
    :param transmat: The path and name for output transition_matrix.txt
    :param markstatemat: The path and name for output markstate_matrix.txt
    :param mark: True/False, see Marknum function
    :param output_mark: True/False, to determine whether should output matrix.txt files
    :return: numstate (number of states); out (see emit2mark function); trans_matrix (transition matrix in list format);
    histone_list (The list of the histone marks)
    '''
    trans_matrix,emit_matrix,numstate,numoutput = load_rawhmm(rawhmmfile)
    histone_list = load_histonefile(histonelistfile)

    # mark_state matrix and transition matrix for writing
    out = emit2mark(emit_matrix,numstate,numoutput,mark)
    mark_state_matrix = pd.DataFrame(out,columns=histone_list,index=np.arange(1,numstate+1,1))
    trans_matrix_final = pd.DataFrame(np.array(trans_matrix),columns=np.arange(1,numstate+1,1),
                                      index=np.arange(1,numstate+1,1))
    if output_mark:
        if markstatemat is None:
            mark_state_matrix.to_csv('./mark_state_matrix.txt',sep='\t')
        else:
            mark_state_matrix.to_csv(markstatemat+'.txt',sep='\t')
        if transmat is None:
            trans_matrix_final.to_csv('./trans_matrix.txt',sep='\t')
        else:
            trans_matrix_final.to_csv(transmat+'.txt',sep='\t')

    return numstate,out,trans_matrix,histone_list


def HMM_matrix_visualization(rawhmmfile,histonelistfile,transmat,markstatemat,mark,outputmark,color_pick):
    '''Visualization of Transition matrix and mark-state matrix'''

    # matrix colors
    if color_pick == 0:
        # red and white
        matrix_color =  matplotlib.colors.LinearSegmentedColormap.from_list("", ["#ffffff","#fdc7c7","#ff8686","#fb4444","#ff0000"])
    elif color_pick == 1:
        # red and yellow
        matrix_color = "YlOrRd"
    elif color_pick == 2:
        # red and blue
        matrix_color = "coolwarm"
    else:
        print("Use default color")
        matrix_color = "coolwarm"

    numstate,out,trans_matrix,histone_list = rawhmm2matrix(rawhmmfile,histonelistfile,transmat,markstatemat,mark,outputmark)
    bar_max = out.max()
    fig, ax = plt.subplots(figsize=(20,15))
    show_annot_array = out >=bar_max/4
    y_label = []
    for i in range(numstate):
        y_label.append('S'+ str(i+1))
    emitpic = sns.heatmap(ax=ax, data=out,linewidths=.5,linecolor="black",yticklabels=y_label,annot=out,
                          xticklabels=histone_list, cmap=matrix_color,vmin=0,vmax=bar_max,annot_kws={"fontsize":24,"color":"black"})
    cbar_emit = emitpic.collections[0].colorbar
    cbar_emit.ax.tick_params(labelsize=18)
    for text, show_annot in zip(ax.texts, (element for row in show_annot_array for element in row)):
        text.set_visible(show_annot)
    plt.xticks(fontsize=20,rotation=0)
    plt.yticks(fontsize=20,rotation=0)
    # get the time
    now = datetime.now()
    dt_string = now.strftime("%d:%m:%Y_%H:%M:%S")
    if markstatemat is None:
        plt.savefig('./mark_state_'+dt_string+'.png')
    else:
        try:
            plt.savefig(markstatemat,bbox_inches='tight',dpi=300)
        except ValueError:
            plt.savefig(markstatemat+'.png',bbox_inches='tight',dpi=300)


    '''plot trans prob matrix'''
    transout = np.array(trans_matrix)
    bar_max = transout.max()
    fig, ax = plt.subplots(figsize=(20,15))
    trans_annot_array = transout >=bar_max/4
    tick_labels = []
    for i in range(numstate):
        tick_labels.append('S'+ str(i+1))

    transpic = sns.heatmap(ax=ax, data=transout,linewidths=.5,linecolor="black",yticklabels=tick_labels,annot=transout,
                          xticklabels=tick_labels, cmap=red_yellow,vmin=0,vmax=bar_max,annot_kws={"fontsize":24,"color":"black"})
    cbar_trans = transpic.collections[0].colorbar
    cbar_trans.ax.tick_params(labelsize=18)
    for text, show_annot in zip(ax.texts, (element for row in trans_annot_array for element in row)):
        text.set_visible(show_annot)
    plt.xticks(fontsize=20,rotation=0)
    plt.yticks(fontsize=20,rotation=0)
    if transmat is None:
        plt.savefig('./trans_'+dt_string+'.png')
    else:
        try:
            plt.savefig(transmat,bbox_inches='tight',dpi=300)
        except ValueError:
            plt.savefig(transmat+'.png',bbox_inches='tight',dpi=300)


def plot_state_coverage_from_dict(celltype,State_dict,interval_other_area,up_boundary,down_boundary,samplepoints,
                                  mapped_pic,spe_colors):
    print('Plotting...')
    '''draw distribution'''
    up_sample = up_boundary // interval_other_area
    down_sample = down_boundary // interval_other_area
    gb_sample = samplepoints
    total_bin = up_sample + down_sample + gb_sample
    States_x = State_dict
    # normalize

    for state in States_x:
        States_x[state] = States_x[state].tolist()

    plt.figure(figsize=(20,15))
    ax = plt.subplot(111)
    for state in States_x:
    #    x = smooth(x,35)
        pic = ax.plot(States_x[state],label='S'+str(state),color=spe_colors[str(state)])
    plt.title(celltype)
    plt.ylabel('Frequency')
    plt.xticks([0,up_sample,up_sample+gb_sample,total_bin],
               ['-'+str(up_boundary//1000)+'kb','TSS','TTS',str(down_boundary//1000)+'kb'])

    fmt='%.2f%%'
    yticks = mtick.FormatStrFormatter(fmt)
    plt.gca().yaxis.set_major_formatter(yticks)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])
    # Put a legend below current axis
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1),
              fancybox=True, shadow=True, ncol=6)
    plt.savefig(mapped_pic,bbox_inches='tight',dpi=300)
    plt.close()

    # plot 1): upstream TSS and 2): rescaled genebody and downstream
    plt.figure()
    ax_up = plt.subplot(111)
    for state in States_x:
        ax_up.plot(States_x[state][:up_sample],label='S'+str(state),color=spe_colors[str(state)])
    plt.title(celltype)
    plt.ylabel('Frequency')
    plt.xticks([0,up_sample//4,up_sample//2,(up_sample*3)//4,up_sample],
               ['-'+str(up_boundary//1000)+'kb','-'+str(3*up_boundary//4000)+'kb',
                '-'+str(up_boundary//2000)+'kb','-'+str(up_boundary//4000)+'kb','TSS'])
    fmt='%.1f%%'
    yticks = mtick.FormatStrFormatter(fmt)
    plt.gca().yaxis.set_major_formatter(yticks)
    box = ax_up.get_position()
    ax_up.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])
    # Put a legend below current axis
    ax_up.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1),
                 fancybox=True, shadow=True, ncol=6)
    plt.savefig(celltype+'_distribution_up.png',bbox_inches='tight',dpi=300)
    plt.close()

    plt.figure()
    ax_gb = plt.subplot(111)
    for state in States_x:
        ax_gb.plot(States_x[state][up_sample:],label='S'+str(state),color=spe_colors[str(state)])
    plt.title(celltype)
    plt.ylabel('Frequency')
    plt.xticks([0,gb_sample,gb_sample+down_sample],
               ['TSS','TTS',str(down_boundary//1000)+'kb'])
    fmt='%.1f%%'
    yticks = mtick.FormatStrFormatter(fmt)
    plt.gca().yaxis.set_major_formatter(yticks)
    box = ax_gb.get_position()
    ax_gb.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])
    # Put a legend below current axis
    ax_gb.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1),
                 fancybox=True, shadow=True, ncol=6)
    plt.savefig(celltype+'_distribution_gbdown.png',bbox_inches='tight',dpi=300)
    plt.close()
    print('\n')
    print('Finished!')


def plot_distribution(celltype,States_x,interval_other_area,up_boundary,down_boundary,samplepoints,spe_colors):
    mapped_pic = celltype + '_distribution.png'
    if file_check(mapped_pic):
        query_mark = query_yes_no('Would you like to overwrite the '+celltype+' distribution plot?')
        if query_mark:
            print('plot '+ celltype + ' distribution...')
            plot_state_coverage_from_dict(celltype,States_x,interval_other_area,up_boundary,down_boundary,samplepoints,
                                  mapped_pic,spe_colors)
        else:
            query_mark2 = query_yes_no('Would you like to write a new '+celltype+' distribution plot?')
            if query_mark2:
                mapped_pic = celltype + '_distribution_' + get_time() +'.png'
                plot_state_coverage_from_dict(celltype,States_x,interval_other_area,up_boundary,down_boundary,samplepoints,
                                  mapped_pic,spe_colors)
    else:
        print('plot '+ celltype + ' distribution...')
        plot_state_coverage_from_dict(celltype,States_x,interval_other_area,up_boundary,down_boundary,samplepoints,
                                  mapped_pic,spe_colors)


def plot_array_distribution(all_info,selected_length,States,region_up,region_down,spe_colors):

    selected_region_sum = {state:np.zeros(selected_length) for state in States}
    for celltype, cover_info in all_info.items():
        for state in States:
            selected_region_sum[state] += cover_info[state][region_up:region_down]

    ave_spacing = {}
    for state in States:
        plt.figure()
        peaks,current_spacing = cal_spacing(State_raw[state])
        ave_spacing[state] = current_spacing * 10
        plt.plot(State_raw[state],color=spe_colors[state.split('_')[1]])
        plt.scatter(peaks,State_raw[state][peaks],marker="^")
        plt.xticks([98,198,298],['0','1000','2000'])
        plt.title(state)
    print(ave_spacing)

def plot_violin_nuc_pos(states_pos_dict,states_list,spe_colors,figname):
    # transfer to dataframe that seaborn recognize
    print("Transfering format..")
    trans_dict = {'Pos. Score':[],'States':[]}
    for state in states_pos_dict:
        trans_dict['Pos. Score'] += states_pos_dict[state]
        trans_dict['States'] += [state] * len(states_pos_dict[state])
    colors_p = [spe_colors[state] for state in states_list]
    df_pos_state = pd.DataFrame(trans_dict)
    print("Plotting..")
    ax = sns.violinplot(x="States", y="Pos. Score", order=states_list,data=df_pos_state, palette=colors_p)
    plt.savefig(figname,dpi=300)
    print("Plotting Finish!")

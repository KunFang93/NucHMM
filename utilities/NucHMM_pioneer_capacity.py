import os
import sys
import bisect
import pandas as pd
import numpy as np
import mygene
import csv
import time
import itertools as it
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from collections import defaultdict, Counter
from statsmodels.stats.proportion import proportion_confint
from matplotlib.patches import Polygon
from scipy import stats
from pathlib import Path

try:
    nucfilelist = sys.argv[1]
    genefile = sys.argv[2]
    FPKMfilelist = sys.argv[3]
    nucTFgenefilelist = sys.argv[4]
    nucTFgenefiltlist = sys.argv[5]
    tffilter = sys.argv[6]
except IndexError:
    print("Usage: python NucHMM_pioneer_capacity.py "
          "list_gl_an_resp_pos_filt.txt gene_expand.txt list_FPKM.txt nuc_TF_gene_intersect.unmofilt.txt "
          "nuc_TF_gene_intersect.mofilt.txt TF_filter.txt")
    exit(1)


print("Current version only suitable for three cell types. Remeber to change Hyperparameters for any new application.")
background_state = [8,12]
distal_states = ['S4','S5','S10']
totalstatenum = 13
alpha_value = 0.05
zFPKM_cutoff = -3
# boxplot upper limit and down limit
top = 200
bottom = -5
zero_mark = False
wrmidfile = True
num_distal = len(distal_states)
spe_colors = {'1':'gold','2':'blue', '3':'orange', '4':'green', '5':'red', '6':'purple', '7':'brown', '8':'pink',
              '9':'gray', '10':'olive','11':'cyan','12':'lightblue','13':'peru', '14':'navy','15':'limegreen',
              '16':'darkcyan','17':'orchid','18':'saddlebrown','19':'sliver','20':'black'}

class zFPKM(object):
    def __init__(self,data,n=512):
        self.__data = np.asarray(data)
        self.__n = n
    def get_data(self):
        return self.__data
    def set_data(self,data):
        self.__data = np.asarray(data)
    def zfpkm_cal(self):
        # remove zero to avoid -inf
        data_no0 = self.__data[self.__data != 0]
        # log2 transform for FPKM values
        data_no0_log2 = np.log2(data_no0)
        # Compute Gaussian kernel density estimation
        current_kernel = stats.gaussian_kde(data_no0_log2)
        # Set the maximum point in the density as the mean for the fitted Gaussian
        x_raw = np.linspace(np.min(data_no0_log2),np.max(data_no0_log2),self.__n)
        y_kernel = current_kernel(x_raw)
        y_kernel_max = max(y_kernel)
        # pick the first point that has the maximum value
        max_idx = [i for i, j in enumerate(y_kernel) if j == y_kernel_max][0]
        mu = x_raw[max_idx]
        # Determine the standard deviation
        U = np.mean(data_no0_log2[data_no0_log2 > mu])
        stdev = (U - mu) * np.sqrt(np.pi/2)
        # Compute zFPKM transform
        zfpkm_result = (data_no0_log2 - mu)/stdev
        return x_raw,y_kernel,mu,stdev,zfpkm_result
    def zFPKM_result(self):
        x_raw, y_kernel,mu, stdev, result = self.zfpkm_cal()
        return result
    def zFPKM_plot(self,outfilename):
        x_raw, y_kernel,mu, stdev, result = self.zfpkm_cal()
        fitted_y = stats.norm.pdf(x_raw,mu,stdev)
        max_y_kernel = max(y_kernel)
        max_fitted_y = max(fitted_y)
        scaleFitted_y = fitted_y * (max_y_kernel / max_fitted_y)
        fig,ax = plt.subplots()
        ax.plot(x_raw,scaleFitted_y,linewidth=4)
        ax.plot(x_raw,y_kernel,linewidth=4)
        ax.plot([-3,-3],[0,max(scaleFitted_y)],linestyle='dashed',color='red',linewidth=3)
        ax.set_ylabel('Density')
        ax.set_xlabel('log2FPKM')
        plt.savefig(outfilename)
        plt.close()

def insert(list, n, pos):
    # pos can be 'start' or 'end'
    while n in list:
        if pos == 'start':
            n -= 1
        elif pos == 'end':
            n += 1
        else:
            print('Wrong position mark!')
            exit(1)
    bisect.insort(list, n)
    return list,n

def W_B_CI(df_col,alpha_value):
    # The reference of Wilson/Brown method: Brown et.al Interval Estimation for a Binomial Proportion
    # the column of the df is a gene and the row of the df is a state
    # build low counts table for Wilson modification first
    low_count_wm_table = {0.90:[0.105,0.532,1.102],0.95:[0.051,0.355,0.818],0.99:[0.010,0.149,0.436]}
    counts = np.asarray(df_col)
    low_bound, up_bound = proportion_confint(count=counts,nobs=sum(counts),alpha=alpha_value,method='wilson')
    for idx, count in enumerate(counts):
        if count == 1 or count == 2 or (count == 3 and sum(counts) > 50):
            low_bound[idx] = low_count_wm_table[1-alpha_value][count-1]/sum(counts)
        else:
            continue
    return low_bound,up_bound

def read_gene_expand(genefile):
    # gene_ordered_name store the gene name in the same order of gene_expand file
    gene_ordered_name = []
    # gene_region_dict: key: gene name; value:[chr, start, end]
    gene_region_dict = defaultdict(list)
    with open(genefile,'r') as gene_file:
        for line in gene_file:
            line_info = line.strip().split()
            line_chr = line_info[0]
            line_start = int(line_info[1])
            line_end = int(line_info[2])
            line_gene = line_info[4]
            gene_ordered_name.append(line_gene)
            gene_region_dict[line_gene].append(line_chr)
            gene_region_dict[line_gene].append(line_start)
            gene_region_dict[line_gene].append(line_end)
    gene_file.close()
    return gene_ordered_name,gene_region_dict

def read_nuchmm_nuc(nucfile):
    # nuc_file should be sorted (default sorted if using nuchmm process)
    # nuc_chr_dict: key:chr; value: nucs position:[0,150,400,550]
    nuc_chr_dict = defaultdict(list)
    # nuc_state_info: key:chr; value: nucs state:[1,2,3,10], for the same chromosome,
    # the length of nuc_state_info[chr] should be the half of the nuc_chr_dict[chr]
    nuc_state_info = defaultdict(list)
    with open(nucfile,'r') as nuc_file:
        # the gl_an_resp_pos_filt file has already been sorted
        for line in nuc_file:
            line_info = line.strip().split()
            line_chr = line_info[0]
            line_start = int(line_info[1])
            line_end = int(line_info[2])
            line_state = int(line_info[3])
            nuc_chr_dict[line_chr].append(line_start)
            nuc_chr_dict[line_chr].append(line_end)
            nuc_state_info[line_chr].append(line_state)
    nuc_file.close()
    return nuc_chr_dict,nuc_state_info

def read_list_file(listfile):
    files_list = []
    with open(listfile,'r') as list_file:
        for line in list_file:
            line_info = line.strip()
            files_list.append(line_info)
    list_file.close()
    return files_list

def state_genes_list(gene_ordered_name,gene_region_dict,nuc_chr_dict,nuc_state_info,rowindex,alpha_value):
    state2index = {int(state[1:]):idx for idx,state in enumerate(rowindex)}
    idx2state = {idx:state for idx,state in enumerate(rowindex)}
    state_gene_counts = {}
    for gene in gene_ordered_name:
        current_chr = gene_region_dict[gene][0]
        current_start = gene_region_dict[gene][1]
        current_end = gene_region_dict[gene][2]
        # step1. insert region start and end to the nuc pos list, make sure not change the original list
        current_nuc_pos_list = nuc_chr_dict[current_chr][:]
        current_nuc_pos_list,revised_start = insert(current_nuc_pos_list,current_start,'start')
        current_nuc_pos_list,revised_end = insert(current_nuc_pos_list,current_end,'end')
        # step2. find index of revised_start and revised_end
        start_index = current_nuc_pos_list.index(revised_start)
        end_index = current_nuc_pos_list.index(revised_end)
        # if index is odd, there is a nucleosome across the boundary
        if start_index % 2 == 0:
            origin_s_index = start_index
        else:
            origin_s_index = start_index - 1
        # while slice(2,6) return 2,3,4,5
        if end_index % 2 == 0:
            origin_e_index = end_index
        else:
            origin_e_index = end_index + 1
        # transfer from nuc pos index to nuc state index
        state_s_index = origin_s_index//2
        state_e_index = origin_e_index//2
        sliced_nuc_pos_list = nuc_chr_dict[current_chr][:][slice(origin_s_index,origin_e_index)]
        sliced_nuc_state_list = nuc_state_info[current_chr][:][slice(state_s_index,state_e_index)]
        # initial gene-state count
        state_gene_counts[gene] = [0] * len(rowindex)
        current_counts_dict = dict(Counter(sliced_nuc_state_list))
        for state in current_counts_dict:
            current_idx = state2index[state]
            state_gene_counts[gene][current_idx] = current_counts_dict[state]
        # double check
        if len(sliced_nuc_pos_list) % 2 == 0:
            continue
        else:
            print(len(sliced_nuc_pos_list))
            print(sliced_nuc_pos_list)
            print(origin_s_index,origin_e_index,current_chr,current_start,current_end)
            print('Find odd number of nucleosome coordinates, please check the code.')
            exit(1)
    # df_counts row represents state and column represents genes
    df_counts = pd.DataFrame(state_gene_counts,index=rowindex)
    CI_info = df_counts.apply(lambda x:W_B_CI(x,alpha_value))
    # state distributed according to its ratio to the total nucleosomes' number
    row_sum = df_counts.sum(axis=1)
    total_nuc_number = sum(row_sum)
    CI_cutoff = []
    for rowsum in row_sum:
        CI_cutoff.append(rowsum/total_nuc_number)
    # initial state_gene_list, key is the states like 'S1" and value is a list contains all genes that
    # associated with this state
    states_gene_list = defaultdict(list)
    gene_list = CI_info.index.values
    for idx, gene_CI in enumerate(CI_info):
        current_gene = gene_list[idx]
        current_low_bound = gene_CI[0]
        for idy, state_lb_prob in enumerate(current_low_bound):
            if state_lb_prob >= CI_cutoff[idy]:
                current_state = idx2state[idy]
                states_gene_list[current_state].append(current_gene)
    return states_gene_list

def write2csv(cts_list,out_name,states_list):
    with open(out_name,"w+") as f:
        writer = csv.writer(f)
        writer.writerow([state for state in states_list])
        if sys.version_info > (3,0):
            for values in it.zip_longest(*cts_list):
                writer.writerow(values)
        else:
            for values in it.izip_longest(*cts_list):
                writer.writerow(values)
    f.close()

def read_FPKM_file(FPKM_file):
    outname = FPKM_file.split('/')[-1].split('.')[0] + '.symbol.csv'
    # check availability of the converted symbol file
    symbol_FPKM_dict = {}
    if os.path.exists(outname):
        symbol_FPKM_df = pd.read_csv(outname)
        for idx,row in symbol_FPKM_df.iterrows():
            symbol_FPKM_dict[row['symbol']] = row['fpkm']
    else:
        # Load the FPKM value of the genes
        df_FPKM = pd.read_csv(FPKM_file,sep="\t")
        # transfer gene id to gene symbol, remove the missing id
        print("ID trasfering..")
        mg = mygene.MyGeneInfo()
        trans_mg=mg.querymany(df_FPKM["gene_id"],scopes='entrezgene,ensembl.gene,HGNC',
                             fields='symbol',species='human',returnall=True)['out']
        # build id_FPKM_dict
        id_FPKM_dict = {}
        for index, row in df_FPKM.iterrows():
            id_FPKM_dict[row["gene_id"]] = row["FPKM"]
        # build symbol_FPKM_dict, key: symbol, value:FPKM
        symbol_FPKM_dict2 = {'symbol':[],'fpkm':[]}
        for info in trans_mg:
            current_gene_id = info['query']
            try:
                current_symbol = info['symbol']
                symbol_FPKM_dict[current_symbol] = id_FPKM_dict[current_gene_id]
                symbol_FPKM_dict2['symbol'].append(current_symbol)
                symbol_FPKM_dict2['fpkm'].append(id_FPKM_dict[current_gene_id])
            except KeyError:
                continue
        symbol_FPKM_df = pd.DataFrame(symbol_FPKM_dict2)
        symbol_FPKM_df.to_csv(outname,index=False)
    return symbol_FPKM_dict

def FPKM_stratify(symbol_FPKM_dict,zFPKM_cutoff,plot_mark):
    # divide symbo_FPKM_dict[celltype] into 4 groups: unexpressed, low/medium/high expressed
    # eventually return a dict with structure dict[celltype][symbol] = [its_FPKM, its_zFPKM, 'mark(unexpressed/low/medium/high)']
    symbol_FPKM_stratified = {celltype:{} for celltype in symbol_FPKM_dict}
    for celltype in symbol_FPKM_dict:
        tmp_FPKM = []
        tmp_symbol = []
        for symbol in symbol_FPKM_dict[celltype]:
            current_FPKM = symbol_FPKM_dict[celltype][symbol]
            if current_FPKM - 0 < 0.000001:
                symbol_FPKM_stratified[celltype][symbol] = [current_FPKM, 'NAN', 'unexpressed']
            else:
                tmp_FPKM.append(current_FPKM)
                tmp_symbol.append(symbol)
        tmp_df = pd.DataFrame(tmp_FPKM,index=tmp_symbol)
        # zFPKM transfer
        tmp_zdata = zFPKM(tmp_df)
        tmp_zFPKM = tmp_zdata.zFPKM_result()
        if plot_mark:
            tmp_zdata.zFPKM_plot(celltype+'_log2FPKM_fiited.png')
        tmp_zFPKM_large_cutoff = []
        tmp_symbol_large_cutoff = []
        for idx, zFPKM_value in enumerate(tmp_zFPKM):
            symbol = tmp_symbol[idx]
            if zFPKM_value <= zFPKM_cutoff:
                symbol_FPKM_stratified[celltype][symbol] = [tmp_FPKM[idx], zFPKM_value, 'unexpressed']
            else:
                tmp_zFPKM_large_cutoff.append(zFPKM_value)
                tmp_symbol_large_cutoff.append(symbol)
        low_cutoff = np.percentile(tmp_zFPKM_large_cutoff,30)
        high_cutoff = np.percentile(tmp_zFPKM_large_cutoff,70)
        for idx, symbol in enumerate(tmp_symbol_large_cutoff):
            if tmp_zFPKM_large_cutoff[idx] <= low_cutoff:
                symbol_FPKM_stratified[celltype][symbol] = [symbol_FPKM_dict[celltype][symbol],
                                                            tmp_zFPKM_large_cutoff[idx], 'low']
            elif tmp_zFPKM_large_cutoff[idx] > high_cutoff:
                symbol_FPKM_stratified[celltype][symbol] = [symbol_FPKM_dict[celltype][symbol],
                                                            tmp_zFPKM_large_cutoff[idx], 'high']
            else:
                symbol_FPKM_stratified[celltype][symbol] = [symbol_FPKM_dict[celltype][symbol],
                                                            tmp_zFPKM_large_cutoff[idx], 'medium']
    return symbol_FPKM_stratified

def assign_FPKM2states(states_gene_list,symbol_FPKM_dict,zero_mark):
    # assign FPKM to each states
    state_FPKM_list = defaultdict(list)
    gene_not_list = []
    for state in states_gene_list:
        for gene in states_gene_list[state]:
            try:
                state_FPKM_list[state].append(symbol_FPKM_dict[gene])
            except KeyError:
                # print("%s not in RNA-seq result" % gene)
                gene_not_list.append(gene)
                if zero_mark:
                    state_FPKM_list[state].append(0)
                continue
    print("%d genes are not in RNA-seq" % len(set(gene_not_list)))
    return state_FPKM_list

def df_log_transform(out):
    # in case there is any 0 value
    out += 1
    out_log = out.applymap(np.log2)
    return out_log

def plot_TF_state_map(out,rowname,max_bar,outfilename):
    mpl.rcParams.update(mpl.rcParamsDefault)
    num_bar = len(rowname)
    print(outfilename,num_bar/2)
    # Set up the matplotlib figure
    f, ax = plt.subplots(figsize=(num_bar/2, num_bar/2))
    # Generate a custom diverging colormap
    cmap = sns.diverging_palette(240, 10, as_cmap=True)
    # sns.set(font_scale=10)
    # ax = sns.clustermap(out,linewidths=.5,linecolor="black",col_cluster=False,square=True,
    #                     yticklabels=rowname,cmap=cmap,vmin=0.3,cbar_kws={"shrink": 0.5})
    ax = sns.heatmap(out,linewidths=.5, square=True, linecolor='white',
                         yticklabels=rowname,cmap=cmap,vmin=0,vmax=max_bar,cbar_kws={"shrink": 0.5},ax=ax)
    for i in range(int(out.shape[1]/4)):
        ax.axvline(x=i*4, color='k',linewidth=3)
    ax.axvline(x=out.shape[1], color='k',linewidth=3)
    _, xlabels = plt.xticks()
    _, ylabels = plt.yticks()
    ax.set_xticklabels(xlabels, size=20)
    ax.set_yticklabels(ylabels, size=20)
    plt.yticks(rotation=0)
    # ax.set_xticklabels(ax.get_xticks(), size=15)
    # ax.set_yticklabels(ax.get_yticks(), size=15)
    strFile = outfilename
    plt.savefig(strFile)
    plt.close()

def plot_heatmap(df,outfilename,TFs_list):
    df_log = (df+1).applymap(np.log2)
    max_bar = df_log.max().max()
    # Set up the matplotlib figure
    f, ax = plt.subplots(figsize=(len(TFs_list)/2, len(TFs_list)/2))
    # Generate a custom diverging colormap
    cmap = sns.diverging_palette(240, 10, as_cmap=True)
    ax = sns.heatmap(df_log,linewidths=.5, square=True, linecolor='white',
                         yticklabels=TFs_list,cmap=cmap,vmin=0,vmax=max_bar,cbar_kws={"shrink": 0.5},ax=ax)
    _, xlabels = plt.xticks()
    _, ylabels = plt.yticks()
    ax.set_xticklabels(xlabels, size=20)
    ax.set_yticklabels(ylabels, size=20)
    plt.yticks(rotation=0)
    # ax.set_xticklabels(ax.get_xticks(), size=15)
    # ax.set_yticklabels(ax.get_yticks(), size=15)
    plt.savefig(outfilename)
    plt.close()

def autolabel(rects,ax):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        test = ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')
        test.set_fontsize(6)

def ratio_list(data_list_dict):
    total_num = 0
    for data_list in data_list_dict:
        total_num += np.array(data_list_dict[data_list])
    data_list_dict_ratio = {}
    for data_list in data_list_dict:
        data_list_dict_ratio[data_list] = np.divide(np.array(data_list_dict[data_list]),total_num) * 100
    for data_list in data_list_dict_ratio:
        data_list_dict_ratio[data_list] = [round(i,1) for i in data_list_dict_ratio[data_list]]
    return data_list_dict_ratio

def plot_bar(data_df,outfilename,y_label):
    # labels is the label for value. e.g In data_dict['unexpressed'] = [1,2,3,4,5,6], labels should be the label of
    # [1,2,3,4,5,6]
    mpl.rcParams.update(mpl.rcParamsDefault)
    fig, ax=plt.subplots()
    sns.barplot(data=data_df,x='States',y='Counts',hue='Exp_Level',palette="rocket")
    ax.set_ylabel(y_label)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    # for rect_key in rect_dict:
    #     autolabel(rect_dict[rect_key],ax)
    fig.tight_layout()
    plt.savefig(outfilename)
    plt.close()

def exp_gene_bar_plot(gene_state_dict,cell_types,expression_level,states_list,out_word,y_label):
    for celltype in cell_types:
        data_4plot_cts = {'Counts':[],'States':[],'Exp_Level':[]}
        for exp_level in expression_level:
            for idx,state in enumerate(states_list):
                data_4plot_cts['Counts'].append(len(gene_state_dict[celltype][exp_level][state]))
                data_4plot_cts['States'].append(state)
                data_4plot_cts['Exp_Level'].append(exp_level)
        data_4plot_cts_df = pd.DataFrame(data_4plot_cts)
        data_4plot_cts_df.to_csv(celltype+'_'+out_word+'_barplot.csv',index=False)
        plot_bar(data_4plot_cts_df,'Bar_plot/'+celltype+'_'+out_word+'_barplot.png',y_label)

def TF_TF_cutoff(TF_df):
    # exclude TF itself
    current_df = TF_df.copy()
    tmp_df = np.fill_diagonal(current_df.values, 0)
    # first cutoff
    # majority_value = np.quantile(TF_df.values.flatten(),0.95)
    # print('majority value',majority_value)
    tmp_df_list = current_df.values.flatten()
    tmp_df_list_filt = np.array([count for count in tmp_df_list if count > 0])
    # minor_value = np.quantile(TF_df.values.flatten(),0.98)
    # tmp_df_list_filt2 = np.array([count for count in tmp_df_list_filt if count < minor_value])
    relationship_cutoff = np.quantile(tmp_df_list_filt,0.75)
    return relationship_cutoff

def plot_bar3(df_normal,category,TFs_list,top_tf,outfilename):
    top_TF = top_tf
    plt.subplots(figsize=(len(TFs_list)//2, len(TFs_list)//2))
    ax = plt.bar(df_normal.index.values.tolist()[-top_TF:],df_normal[category].values.tolist()[-top_TF:],
                 color='grey',edgecolor='black')
    plt.xticks(fontsize=len(TFs_list)//2,rotation=80)
    plt.yticks(fontsize=len(TFs_list)//2)
    plt.savefig(outfilename.split('.')[0]+'.png')
    plt.close()

def build_nuc_tf_state_cts_dict(nuc_TF_gene_files,cts_gene_state_dict,cell_types,expression_level,rowindex):
    nuc_TF_state_cts_dict = {celltype:{exp_level:{state:{} for state in rowindex} for exp_level in expression_level} for celltype in cell_types}
    TF_list = {celltype:[] for celltype in cell_types}
    for file in nuc_TF_gene_files:
        count = 0
        celltype = file.split('/')[-1].split('_')[0]
        with open(file, 'r') as TF_file:
            for line in TF_file:
                count += 1
                sys.stdout.write('\rLine:' + str(count))
                line_info = line.strip().split()
                line_state = line_info[3]
                current_state = 'S'+str(line_state)
                line_gene_name = line_info[-1]
                line_TF_name = line_info[-6]
                TF_list[celltype].append(line_TF_name)
                for exp_level in expression_level:
                    if line_gene_name not in cts_gene_state_dict[celltype][exp_level][current_state]:
                        continue
                    else:
                        if line_TF_name not in nuc_TF_state_cts_dict[celltype][exp_level][current_state]:
                            nuc_TF_state_cts_dict[celltype][exp_level][current_state][line_TF_name] = 1
                        else:
                            nuc_TF_state_cts_dict[celltype][exp_level][current_state][line_TF_name] += 1
        TF_file.close()
    return nuc_TF_state_cts_dict,TF_list

def build_cts_ss_gene_state_dict(cts_gene_state_dict,target_states,cell_types,expression_level):
    cts_ss_gene_state_dict = {celltype:{exp_level:{state:[] for state in distal_states}
                                     for exp_level in expression_level} for celltype in cell_types}
    for celltype in cell_types:
        for exp_level in expression_level:
            for state in target_states:
                tmp_cts_dss = cts_gene_state_dict[celltype][exp_level][state]
                for state_y in target_states:
                    if state_y != state:
                        tmp_cts_dss = set(tmp_cts_dss) - set(cts_gene_state_dict[celltype][exp_level][state_y])
                cts_ss_gene_state_dict[celltype][exp_level][state] = list(tmp_cts_dss)
    return cts_ss_gene_state_dict

def build_TNS_matrix(nuc_TF_gene_files,cts_ss_gene_state_dict,TF2index_dict,target_states,cell_types,expression_level):
    # nuc_TF_gene_list is for TTA
    nuc_TF_gene_list_dict = {celltype:{} for celltype in cell_types}
    nuc_TF_state_cts_ss_dict = {celltype:{exp_level:{state:{} for state in target_states} for exp_level in expression_level}
                                 for celltype in cell_types}
    for file in nuc_TF_gene_files:
        count = 0
        celltype = file.split('/')[-1].split('_')[0]
        with open(file, 'r') as TF_file:
            for line in TF_file:
                count += 1
                sys.stdout.write('\rLine:' + str(count))
                line_info = line.strip().split()
                line_chr = line_info[0]
                line_start = line_info[1]
                line_end = line_info[2]
                line_state = line_info[3]
                current_state = 'S'+str(line_state)
                line_gene_name = line_info[-1]
                line_TF_name = line_info[-6]
                # build dict for the later on TFs association index calculation
                key = line_chr + '_'  + line_start + '_' + line_end
                if key not in nuc_TF_gene_list_dict[celltype]:
                    nuc_TF_gene_list_dict[celltype][key] = {'state':line_state,'TFs':[line_TF_name],'genes':[line_gene_name]}
                else:
                    if nuc_TF_gene_list_dict[celltype][key]['state'] == line_state:
                        if line_TF_name not in nuc_TF_gene_list_dict[celltype][key]['TFs']:
                            nuc_TF_gene_list_dict[celltype][key]['TFs'].append(line_TF_name)
                        if line_gene_name not in nuc_TF_gene_list_dict[celltype][key]['genes']:
                            nuc_TF_gene_list_dict[celltype][key]['genes'].append(line_gene_name)
                    else:
                        print('Not unique state for a nucleosome')
                if current_state in target_states:
                    for exp_level in expression_level:
                        if line_gene_name not in cts_ss_gene_state_dict[celltype][exp_level][current_state]:
                            continue
                        else:
                            if line_TF_name not in nuc_TF_state_cts_ss_dict[celltype][exp_level][current_state]:
                                nuc_TF_state_cts_ss_dict[celltype][exp_level][current_state][line_TF_name] = 1
                            else:
                                nuc_TF_state_cts_ss_dict[celltype][exp_level][current_state][line_TF_name] += 1
                else:
                    continue
        TF_file.close()
    Path("TF_state_expression_cts").mkdir(parents=True, exist_ok=True)
    cts_ss_gene_state_TF_total_df = {}
    cts_ss_gene_state_TF = {celltype:{label:[0]*len(TF_list[celltype]) for label in distal_x_label} for celltype in cell_types}
    for celltype in cell_types:
        for exp_level in expression_level:
            for state in target_states:
                for TF in nuc_TF_state_cts_ss_dict[celltype][exp_level][state]:
                    current_idx = TF2index_dict[celltype][TF]
                    current_x_label = state + '_' + exp_level[0]
                    cts_ss_gene_state_TF[celltype][current_x_label][current_idx] += nuc_TF_state_cts_ss_dict[celltype][exp_level][state][TF]
        cts_ss_gene_state_TF_total_df[celltype] = pd.DataFrame(cts_ss_gene_state_TF[celltype],index=TF_list[celltype])
        if wrmidfile:
            cts_ss_gene_state_TF_total_df[celltype].to_csv('TF_state_expression_cts/'+celltype+'_TF_state_exp_cts_dss.csv')
    return nuc_TF_gene_list_dict, cts_ss_gene_state_TF_total_df

def pioneer_capacity(TF_relationship_df_unexp,TF_relationship_df_high,cts_dss_gene_state_TF_df,TFs_list,tf_filters,distal_states,outfilename):
    cols = ['chromatin opening score','TFs complex unexp_score','TFs complex h_score','TFs complex n_score',
            'TFs complex ratio','TFs complex final','dynamic change score']
    df = pd.DataFrame(0.0,index=TFs_list,columns=cols)
    relationship_cutoff_unexp = TF_TF_cutoff(TF_relationship_df_unexp)
    relationship_cutoff_high = TF_TF_cutoff(TF_relationship_df_high)
    # fill the TFs complex score
    for tf in TFs_list:
        if len(TF_relationship_df_unexp[tf][TF_relationship_df_unexp[tf] >= relationship_cutoff_unexp]) == 0:
            df['TFs complex unexp_score'][tf] = 0
        else:
            df['TFs complex unexp_score'][tf] = float(len(TF_relationship_df_unexp[tf][TF_relationship_df_unexp[tf] >= relationship_cutoff_unexp]))
        # avoid relationship_cutoff_unexp larger than relationship_cutoff_high
        if relationship_cutoff_unexp > relationship_cutoff_high:
            h_cutoff = relationship_cutoff_high - 1
        else:
            h_cutoff = relationship_cutoff_unexp - 1
        if len(TF_relationship_df_high[tf][TF_relationship_df_high[tf] >= h_cutoff]) == 0:
            df['TFs complex h_score'][tf] = 0
        else:
            df['TFs complex h_score'][tf] = float(len(TF_relationship_df_high[tf][TF_relationship_df_high[tf] >= h_cutoff]))
        # calculate the ratio between high and unexp (number of related TFs)
        if len(TF_relationship_df_unexp[tf][TF_relationship_df_unexp[tf] >= relationship_cutoff_unexp]) == 0:
            df['TFs complex ratio'][tf] = 0
        else:
            df['TFs complex ratio'][tf] = np.divide(df['TFs complex h_score'][tf],df['TFs complex unexp_score'][tf])
            print('TFs complex ratio',df['TFs complex ratio'][tf],float(df['TFs complex h_score'][tf]),float(df['TFs complex unexp_score'][tf]),
                  float(df['TFs complex ratio'][tf]))
    # in unexp case, the less related TFs, the strongest of the pioneer capacity
    for tf in TFs_list:
        if df['TFs complex unexp_score'][tf] == 0:
            df['TFs complex n_score'][tf] = 0
        else:
            df['TFs complex n_score'][tf] = np.divide(max(df['TFs complex unexp_score']),df['TFs complex unexp_score'][tf])
        df['TFs complex final'][tf] = df['TFs complex n_score'][tf] * df['TFs complex ratio'][tf]
    # fill the chromatin opening score
    df2 = pd.DataFrame(index=cts_dss_gene_state_TF_df.index.values,columns=distal_states)
    sum_count = {state:0 for state in distal_states}
    for state in distal_states:
        selected_col = [state+'_u',state+'_l',state+'_m',state+'_h']
        current_sum_count = cts_dss_gene_state_TF_df[selected_col].sum().sum()
        sum_count[state] = current_sum_count
        df2[state] = cts_dss_gene_state_TF_df[selected_col].sum(axis=1)/current_sum_count
    # openning_cutoff = df2.sum(axis=1).quantile(0.15)
    df['chromatin opening score'] = (2*np.exp(df2['S10']))/(np.exp(df2['S4'])+np.exp(df2['S5'])) * ((df2['S10']+df2['S4']+df2['S5'])/3)
    # fille the dynamic change score
    # plus 1 avoid inf
    df['dynamic change score'] = 0.5*((cts_dss_gene_state_TF_df['S5_m'] + cts_dss_gene_state_TF_df['S5_h'])/\
                                 (cts_dss_gene_state_TF_df['S5_u'] + cts_dss_gene_state_TF_df['S5_l']+1) * df2['S5']+\
                                 (cts_dss_gene_state_TF_df['S10_u'] + cts_dss_gene_state_TF_df['S10_l'])/\
                                 (cts_dss_gene_state_TF_df['S10_m'] + cts_dss_gene_state_TF_df['S10_h']+1) * df2['S10'])
    df.to_excel(outfilename[:-5]+'_raw.xlsx')
    df_normal = df/df.max()
    # make 0.25 penalty to 0 score of TFs complex final
    df_normal['TFs complex final'] = df_normal['TFs complex final'].replace(0,-0.16)
    # avoid negative value
    df_normal['TFs complex final'] = df_normal['TFs complex final'] + 0.16
    df_normal['pioneer capacity'] = df_normal['chromatin opening score'] + 1.5* df_normal['TFs complex final'] + 0.5*df_normal['dynamic change score']
    df_normal = df_normal.sort_values(by='pioneer capacity')
    df_normal = df_normal.round(3)
    df_normal['pioneer capacity'] = df_normal['pioneer capacity']/df_normal['pioneer capacity'].max()*100
    df_normal['tf'] = df_normal.index.values
    df_normal = df_normal[df_normal['tf'].isin(tf_filters)]
    df_normal.to_excel(outfilename)
    # plot top 25
    if len(tf_filters)>=25:
        plot_bar3(df_normal,'pioneer capacity',TFs_list,25,outfilename)
        plot_bar3(df_normal,'chromatin opening score',TFs_list,25,outfilename.split('_')[0]+'_chromatin_opening_index.xlsx')
        plot_bar3(df_normal,'TFs complex final',TFs_list,25,outfilename.split('_')[0]+'_TF_TF_index.xlsx')
        plot_bar3(df_normal,'dynamic change score',TFs_list,25,outfilename.split('_')[0]+'_dynamic_change_index.xlsx')
    else:
        plot_bar3(df_normal,'pioneer capacity',TFs_list,len(tf_filters),outfilename)
        plot_bar3(df_normal,'chromatin opening score',TFs_list,len(tf_filters),outfilename.split('_')[0]+'_chromatin_opening_index.xlsx')
        plot_bar3(df_normal,'TFs complex final',TFs_list,len(tf_filters),outfilename.split('_')[0]+'_TF_TF_index.xlsx')
        plot_bar3(df_normal,'dynamic change score',TFs_list,len(tf_filters),outfilename.split('_')[0]+'_dynamic_change_index.xlsx')


print("Read gene region file..")
gene_order_name,gene_region_dict = read_gene_expand(genefile)
# nucfilelist should have same order with FPKMfileslist
nucfile_list = read_list_file(nucfilelist)
FPKMfile_list = read_list_file(FPKMfilelist)
TFfilter_list = read_list_file(tffilter)

print("Read filtered nucleosome files..")
cell_types = [filename.split('/')[-1].split('_')[0] for filename in nucfile_list]
rowindex = ['S'+str(i) for i in range(1,totalstatenum+1) if i not in background_state]
nuc_chr_dict_total = {}
nuc_state_info_total = {}
state_gene_list_total = {}
symbol_FPKM_total = {}
symbol_FPKM_stratified = {}
for idx,file in enumerate(nucfile_list):
    current_nuc_chr_dict, current_nuc_state_info = read_nuchmm_nuc(file)
    current_celltype = cell_types[idx]
    nuc_chr_dict_total[current_celltype] = current_nuc_chr_dict
    nuc_state_info_total[current_celltype] = current_nuc_state_info
    print('Building ' + current_celltype + ' state gene counts matrix')
    # current_state_gene_list is a dictionary
    # current_states_gene_list[state] = gene_list
    current_states_gene_list = state_genes_list(gene_order_name,gene_region_dict,current_nuc_chr_dict,
                                                                 current_nuc_state_info,rowindex,alpha_value)
    state_gene_list_total[current_celltype] = current_states_gene_list
    print('Building ' + current_celltype +' state FPKM dictionary')
    current_symbol_FPKM_dict = read_FPKM_file(FPKMfile_list[idx])
    symbol_FPKM_total[current_celltype] = current_symbol_FPKM_dict

# divide FPKM into four groups and plot the ratio of each group in each states
# symbol_FPKM_stratified[celltype][gene] = [FPKM,zFPKM,expression_level(unexp/low/med/high)]
symbol_FPKM_stratified = FPKM_stratify(symbol_FPKM_total,zFPKM_cutoff,True)

# identify cts gene based on expression level
# gene_FPKM[celltype][exp_level] = [genes]
# cts_gene_FPKM[celltype][expression_level] = [cts_genes]
expression_level = ['unexpressed','low','medium','high']
gene_FPKM = {celltype:{exp:[] for exp in expression_level} for celltype in cell_types}
cts_gene_FPKM = {celltype:{exp:[] for exp in expression_level} for celltype in cell_types}
common_gene_FPKM = {celltype:{exp:[] for exp in expression_level} for celltype in cell_types}
for celltype in cell_types:
    for gene in symbol_FPKM_stratified[celltype]:
        current_exp_level = symbol_FPKM_stratified[celltype][gene][2]
        gene_FPKM[celltype][current_exp_level].append(gene)

for celltype_x in cell_types:
    for exp_level in expression_level:
        tmp_gene_list = gene_FPKM[celltype_x][exp_level]
        for celltype_y in cell_types:
            if celltype_y != celltype_x:
                tmp_gene_list = set(tmp_gene_list) - set(gene_FPKM[celltype_y][exp_level])
        cts_gene_FPKM[celltype_x][exp_level] = list(tmp_gene_list)
        common_gene_FPKM[celltype_x][exp_level] = list(set(gene_FPKM[celltype_x][exp_level]) - tmp_gene_list)
        # check if correct
        tmp_length = len(cts_gene_FPKM[celltype_x][exp_level]) + len(common_gene_FPKM[celltype_x][exp_level])
        if tmp_length == len(gene_FPKM[celltype_x][exp_level]):
            continue
        else:
            print('recheck the code.')

# save the basic stats of the cts gene
cts_gene_basic_stats = {'celltype':[],'exp_level':[],'total_genes_num':[],'cts_gene_num':[]}
for celltype in cell_types:
    for exp_level in expression_level:
        print(celltype,exp_level,'total number:',len(gene_FPKM[celltype][exp_level]),
              'cts number:', len(cts_gene_FPKM[celltype][exp_level]))
        cts_gene_basic_stats['celltype'].append(celltype)
        cts_gene_basic_stats['exp_level'].append(exp_level)
        cts_gene_basic_stats['total_genes_num'].append(len(gene_FPKM[celltype][exp_level]))
        cts_gene_basic_stats['cts_gene_num'].append(len(cts_gene_FPKM[celltype][exp_level]))
cts_gene_basic_stats_df = pd.DataFrame(cts_gene_basic_stats)
cts_gene_basic_stats_df.to_csv('cts_gene_basic_stats.txt',index=False,sep='\t')

# build cts_gene_state_dict
cts_gene_state_dict = {celltype:{exp_level:{state:[] for state in rowindex} for exp_level in expression_level} for celltype in cell_types}
common_gene_state_dict = {celltype:{exp_level:{state:[] for state in rowindex} for exp_level in expression_level} for celltype in cell_types}
for celltype in cell_types:
    for state in rowindex:
        current_state_gene_list = state_gene_list_total[celltype][state]
        for exp_level in expression_level:
            tmp_cts_gene_state_list = list(set(cts_gene_FPKM[celltype][exp_level]).intersection(set(current_state_gene_list)))
            tmp_common_gene_state_list = list(set(common_gene_FPKM[celltype][exp_level]).intersection(set(current_state_gene_list)))
            cts_gene_state_dict[celltype][exp_level][state] = tmp_cts_gene_state_list
            common_gene_state_dict[celltype][exp_level][state] = tmp_common_gene_state_list

# write cts_gene_state_dict
Path("cts_genes_list").mkdir(parents=True, exist_ok=True)
for celltype in cell_types:
    for exp_level in expression_level:
        tmp_list = []
        for state in rowindex:
            tmp_list.append(cts_gene_state_dict[celltype][exp_level][state])
        write2csv(tmp_list,'cts_genes_list/'+celltype+'_'+exp_level+'_cts_gene_state.csv',rowindex)

# plot the basic stats of cts_gene_state and common_gene_state (bar plot show gene counts)
Path("Bar_plot").mkdir(parents=True, exist_ok=True)
for celltype in cell_types:
    data_4plot_cts = {'Counts':[],'States':[],'Exp_Level':[]}
    data_4plot_common = {'Counts':[],'States':[],'Exp_Level':[]}
    for exp_level in expression_level:
        for idx,state in enumerate(rowindex):
            data_4plot_cts['Counts'].append(len(cts_gene_state_dict[celltype][exp_level][state]))
            data_4plot_common['Counts'].append(len(common_gene_state_dict[celltype][exp_level][state]))
            data_4plot_cts['States'].append(state)
            data_4plot_common['States'].append(state)
            data_4plot_cts['Exp_Level'].append(exp_level)
            data_4plot_common['Exp_Level'].append(exp_level)
    data_4plot_cts_df = pd.DataFrame(data_4plot_cts)
    data_4plot_common_df = pd.DataFrame(data_4plot_common)
    data_4plot_cts_df.to_csv(celltype+'_cts_gene_state_count.csv',index=False)
    data_4plot_common_df.to_csv(celltype+'_common_gene_state_count.png',index=False)
    plot_bar(data_4plot_cts_df,'Bar_plot/'+celltype+'_cts_gene_state_count.png','Number of genes')
    plot_bar(data_4plot_common_df,'Bar_plot/'+celltype+'_common_gene_state_count.png','Number of genes')

# find TF association with cts genes
# nuc_TF_gene is (nuc & TF) & gene, & means bedtools intersect
# TF heatmap
nuc_TF_gene_files = read_list_file(nucTFgenefilelist)
# motif filtered intersect peak files
nuc_TF_gene_mofilt_files = read_list_file(nucTFgenefiltlist)
nuc_TF_state_cts_dict, TF_list = build_nuc_tf_state_cts_dict(nuc_TF_gene_files,cts_gene_state_dict,cell_types,expression_level,rowindex)
nuc_TF_state_cts_mofilt_dict, TF_list2 = build_nuc_tf_state_cts_dict(nuc_TF_gene_mofilt_files,cts_gene_state_dict,cell_types,expression_level,rowindex)
# build gene state TF counts matrix
print("Processing gene-state-TF information")
TF2index_dict = {celltype:{} for celltype in cell_types}
for celltype in cell_types:
    TF_list[celltype] = list(set(TF_list[celltype]))
    TF_list[celltype].sort()
    print(celltype,'TF number',len(TF_list[celltype]))
    for idx, TF in enumerate(TF_list[celltype]):
        TF2index_dict[celltype][TF] = idx

Path("Heatmap").mkdir(parents=True, exist_ok=True)
distal_x_label = []
for state in distal_states:
    for exp_level in expression_level:
        distal_x_label.append(state+'_'+exp_level[0])

# cts and state specific
# cts_gene_state_dict[celltype][exp_level][state] is only cts
# cts & distal state specific (dss)
cts_dss_gene_state_dict = build_cts_ss_gene_state_dict(cts_gene_state_dict,distal_states,cell_types,expression_level)
for celltype in cell_types:
    for exp_level in expression_level:
        tmp_list = []
        for state in distal_states:
            tmp_list.append(cts_dss_gene_state_dict[celltype][exp_level][state])
        write2csv(tmp_list,'cts_genes_list/'+celltype+'_'+exp_level+'_cts_dss_gene_state.csv',distal_states)

# bar plot show gene counts based on expression level
exp_gene_bar_plot(cts_dss_gene_state_dict,cell_types,expression_level,distal_states,'cts_dss_gene_state_count','Number of genes')

# TNS and TTA matrix
# nuc_TF_gene_list_dict is exactly same as nuc_TF_gene_list_dict2
nuc_TF_gene_list_dict, cts_dss_gene_state_TF_total_nofilt_df \
    = build_TNS_matrix(nuc_TF_gene_files,cts_dss_gene_state_dict,TF2index_dict,distal_states,cell_types,expression_level)
nuc_TF_gene_list_dict2, cts_dss_gene_state_TF_total_mofilt_df = \
    build_TNS_matrix(nuc_TF_gene_mofilt_files,cts_dss_gene_state_dict,TF2index_dict,distal_states,cell_types,expression_level)

for celltype in cell_types:
    out_cts_dss_nofilt = cts_dss_gene_state_TF_total_nofilt_df[celltype]
    out_cts_dss_nofilt_log = df_log_transform(out_cts_dss_nofilt)
    max_value_nofilt = out_cts_dss_nofilt_log.max().max()
    plot_TF_state_map(out_cts_dss_nofilt_log,TF_list[celltype],max_value_nofilt,'Heatmap/'+celltype + '_cts_dss_state_exp_TFs_heatmap_log2.png')
    # motif filt tns heatmap
    out_cts_dss_mofilt = cts_dss_gene_state_TF_total_mofilt_df[celltype]
    out_cts_dss_mofilt_log = df_log_transform(out_cts_dss_mofilt)
    max_value_mofilt = out_cts_dss_mofilt_log.max().max()
    plot_TF_state_map(out_cts_dss_mofilt_log,TF_list[celltype],max_value_mofilt,'Heatmap/'+celltype + '_cts_dss_state_exp_TFs_mofilt_heatmap_log2.png')

print('\n')
print('Building TF-TF association matrix')
TF_relationship_df_distal_dict = {celltype:{state:{exp_level:pd.DataFrame(0,index=TF_list[celltype],columns=TF_list[celltype])
                                            for exp_level in expression_level} for state in distal_states} for celltype in cell_types}
for celltype in cell_types:
    nuc_count = 0
    print('\n',celltype)
    for nuc in nuc_TF_gene_list_dict[celltype]:
        nuc_count += 1
        current_state = 'S' + nuc_TF_gene_list_dict[celltype][nuc]['state']
        if current_state not in distal_states:
            continue
        else:
            current_genes_list = nuc_TF_gene_list_dict[celltype][nuc]['genes']
            assign_exp = ''
            for exp_level in expression_level:
                if len(set(cts_dss_gene_state_dict[celltype][exp_level][current_state]).intersection(set(current_genes_list))) >=1:
                    for tf_x in nuc_TF_gene_list_dict[celltype][nuc]['TFs']:
                        for tf_y in nuc_TF_gene_list_dict[celltype][nuc]['TFs']:
                            TF_relationship_df_distal_dict[celltype][current_state][exp_level][tf_y][tf_x] += 1

Path("TTA_heatmap").mkdir(parents=True, exist_ok=True)
for celltype in cell_types:
    for state in distal_states:
        for exp_level in expression_level:
            TF_relationship_df_distal_dict[celltype][state][exp_level].to_csv('TTA_heatmap/'+celltype+'_'+state+'_'+exp_level+'_TF_heatmap_raw.csv')
            plot_heatmap(TF_relationship_df_distal_dict[celltype][state][exp_level],
                         'TTA_heatmap/'+celltype+'_'+state+'_'+exp_level+'_TF_heatmap',TF_list[celltype])

# for celltype in cell_types:
#     TF_relationship_df_distal_dict[celltype]['S10']['unexpressed'].to_csv(celltype+'_S10_unexp.csv')
TF_filters = {celltype:[] for celltype in cell_types}
for tffilt_file in TFfilter_list:
    celltype = tffilt_file.split('/')[-1].split('_')[0]
    with open(tffilt_file,'r') as input_f:
        for line in input_f:
            line_tf = line.strip()
            TF_filters[celltype].append(line_tf)
    input_f.close()

print('/n')
for celltype in cell_types:
    print(celltype)
    pioneer_capacity(TF_relationship_df_distal_dict[celltype]['S10']['unexpressed'].copy(),
                     TF_relationship_df_distal_dict[celltype]['S10']['high'].copy(),
                     cts_dss_gene_state_TF_total_mofilt_df[celltype].copy(),TF_list[celltype].copy(),TF_filters[celltype],
                     distal_states,celltype+'_pioneer_capacity.xlsx')

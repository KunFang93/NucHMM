import sys
import numpy as np
import pandas as pd
import subprocess
from pybedtools import BedTool
from collections import defaultdict,Counter
import matplotlib.pyplot as plt
from io import StringIO
from matplotlib.patches import Polygon

try:
    nuceventfile = sys.argv[1]
    nucinfo = sys.argv[2]
    srt_uniq = sys.argv[3]
    likewiglist = sys.argv[4]
except IndexError:
    print("Usage: python NucHMM_SE_influence.py nuc_event_intersect.bed nuc_info.bed srt_uniq.bed likewig_list.txt")
    exit(1)


def read_list_file(listfile):
    files_list = []
    with open(listfile,'r') as list_file:
        for line in list_file:
            line_info = line.strip()
            files_list.append(line_info)
    list_file.close()
    return files_list

def outlier_threshold(target_list, k):
    '''
    Use q3 + k * iqr to identify outlier
    :param target_list:
    :param k: coefficient
    :return: up cutoff for the target_list
    '''
    data = target_list
    quartile1 = np.quantile(data, 0.25)
    quartile3 = np.quantile(data, 0.75)
    interval_qr = quartile3 - quartile1
    upb = quartile3 + k * interval_qr
    return upb


def nuc_positioning(inputfile,k):

    # value: list[0]:Nucleosome_index; list[1]:Width_between_inflection; list[2]:Peak_height; list[3]:Area_under_curve
    # list[4]:Physical_property; list[5]:-log10(Pvalue_of_peak); list[6]:-log10(Pvalue_of_valley)
    nuc_pos_result = {}
    nucs_width = []
    nucs_height = []
    nucs_area = []
    nucs_coef = []
    nucs_pvalpeak = []
    nucs_pvalvalley = []
    nucs_coor = []
    # read the post-array-filtered file
    nuc_count = 0
    print('Read input file..')
    with open(inputfile, 'r') as input_file:
        for line in input_file:
            line_info = line.strip().split()
            line_chr = line_info[0]
            line_start = line_info[1]
            line_end = line_info[2]
            nuc_count += 1
            sys.stdout.write('\rRead Line:'+str(nuc_count))
            nucs_coor.append(line_chr+'_'+line_start+'_'+line_end)
            line_width = float(line_info[4])
            line_height = float(line_info[5])
            line_area = float(line_info[6])
            # actually not use it in the equation
            if line_info[7] == "MainPeak":
                line_coef = 1
            else:
                line_coef = 0.5
            line_ppeak = round(float(line_info[8]), 2)
            line_pvalley = round(float(line_info[9]), 2)
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
    print('Positioning score equation is (height*log2(PP*PV+1)*sqrt(area))/(width**3).')

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

    for idx, coors in enumerate(nucs_coor):
        nuc_pos_result[coors] = positioning_score[idx]

    return nuc_pos_result


def filter_heter_event(srt_uniq_file,event_dict,states_list):
    final_event = {}
    final_event_nuc = {state:{} for state in states_list}
    for state in states_list:
        print("processing state {}".format(state))
        event_bed = BedTool.from_dataframe(event_dict[state])

        # filter some events that do not have state homogeneous 11,7,7,7,1,11
        intersect_result = BedTool.intersect(srt_uniq_file,event_bed,wa=True,wb=True)
        # print(len(SE_intersect_result),len(noSE_intersect_result))

        event_keeped_dict = defaultdict(list)
        event_nuc_dict = defaultdict(list)
        for line in intersect_result:
            line_state = line[3]
            line_key = line[4] + '_' + line[5] + '_' + line[6]
            current_nuc = line[0] + '_' + line[1] + '_' +line[2]
            event_keeped_dict[line_key].append(line_state)
            event_nuc_dict[line_key].append(current_nuc)

        tmp = {'chrom':[],'start':[],'end':[]}
        for key in event_keeped_dict:
            num_nuc = len(event_keeped_dict[key])
            # if target state outnumber other state in this event, consider as meaningful event for downstream analysis

            if Counter(event_keeped_dict[key])[str(state)] == 0:
                print("check state type, str or int")
                print(event_keeped_dict[key])

            if Counter(event_keeped_dict[key])[str(state)] <= num_nuc/2:
                # print(key,Counter(event_keeped_dict[key]))
                continue
            else:
                final_event_nuc[state][key] = event_nuc_dict[key]
                tmp['chrom'].append(key.split('_')[0])
                tmp['start'].append(key.split('_')[1])
                tmp['end'].append(key.split('_')[2])
        final_event[state] = pd.DataFrame(tmp,index=np.arange(len(tmp['chrom'])))
    return final_event,final_event_nuc

def colored_box_plot(data_plot,rowindex,title_name,colors,outfilename):
    # data_plot:[[list1 of FPKM value],[list2 of FPKM value]]
    fig, ax1 = plt.subplots(figsize=(len(data_plot)/2, 6))
    fig.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)

    #bp = ax1.boxplot(data_plot, notch=0, sym='+', vert=1, whis=1.5)
    bp = ax1.boxplot(data_plot, notch=0, vert=1, whis=1.5, showfliers=False)
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black')
    # plt.setp(bp['fliers'], color='red', marker='+')

    ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
                   alpha=0.5)

    ax1.set_axisbelow(True)
    ax1.set_title(title_name)
    ax1.set_xlabel('States')
    ax1.set_ylabel('FPKM Value')

    # fill the boxes with desired colors
    box_colors = colors
    num_boxes = len(data_plot)
    medians = np.empty(num_boxes)
    for i in range(num_boxes):
        box = bp['boxes'][i]
        boxX = []
        boxY = []
        for j in range(5):
            boxX.append(box.get_xdata()[j])
            boxY.append(box.get_ydata()[j])
        box_coords = np.column_stack([boxX, boxY])
        ax1.add_patch(Polygon(box_coords, facecolor=box_colors[i]))
        # Now draw the median lines back over what we just filled in
        med = bp['medians'][i]
        medianX = []
        medianY = []
        for j in range(2):
            medianX.append(med.get_xdata()[j])
            medianY.append(med.get_ydata()[j])
            ax1.plot(medianX, medianY, 'k')
        medians[i] = medianY[0]
        # Finally, overplot the sample averages, with horizontal alignment
        # in the center of each box
        ax1.plot(np.average(med.get_xdata()), np.average(data_plot[i]),
                 color='w', marker='*', markeredgecolor='k')

    # Set the axes ranges and axes labels
    ax1.set_xlim(0.5, num_boxes + 0.5)
    # ax1.set_ylim(bottom, top)
    ax1.set_xticklabels(rowindex,rotation=45, fontsize=8)

    # Due to the Y-axis scale being different across samples, it can be
    # hard to compare differences in medians across the samples. Add upper
    # X-axis tick labels with the sample medians to aid in comparison
    # (just use two decimal places of precision)
    pos = np.arange(num_boxes) + 1
    upper_labels = [str(np.round(s, 2)) for s in medians]
    weights = ['bold', 'semibold']
    color_idx = 0
    for tick, label in zip(range(num_boxes), ax1.get_xticklabels()):
        k = tick % 2
        ax1.text(pos[tick], .97, upper_labels[tick],
                 transform=ax1.get_xaxis_transform(),
                 horizontalalignment='center', size='x-small',
                 weight=weights[k], color=box_colors[color_idx])
        color_idx += 1

    # Finally, add a basic legend
    fig.text(0.80, 0.025, '*', color='white', backgroundcolor='silver',
             weight='roman', size='medium')
    fig.text(0.815, 0.023, ' Average Value', color='black', weight='roman',
             size='x-small')

    plt.savefig(outfilename)
    plt.close()


def _c(ca, i, j, p, q):

    if ca[i, j] > -1:
        return ca[i, j]
    elif i == 0 and j == 0:
        ca[i, j] = np.linalg.norm(p[i]-q[j])
    elif i > 0 and j == 0:
        ca[i, j] = max(_c(ca, i-1, 0, p, q), np.linalg.norm(p[i]-q[j]))
    elif i == 0 and j > 0:
        ca[i, j] = max(_c(ca, 0, j-1, p, q), np.linalg.norm(p[i]-q[j]))
    elif i > 0 and j > 0:
        ca[i, j] = max(
            min(
                _c(ca, i-1, j, p, q),
                _c(ca, i-1, j-1, p, q),
                _c(ca, i, j-1, p, q)
            ),
            np.linalg.norm(p[i]-q[j])
            )
    else:
        ca[i, j] = float('inf')

    return ca[i, j]


def frdist(p, q):
    """
    Computes the discrete Fréchet distance between
    two curves. The Fréchet distance between two curves in a
    metric space is a measure of the similarity between the curves.
    The discrete Fréchet distance may be used for approximately computing
    the Fréchet distance between two arbitrary curves,
    as an alternative to using the exact Fréchet distance between a polygonal
    approximation of the curves or an approximation of this value.
    This is a Python 3.* implementation of the algorithm produced
    in Eiter, T. and Mannila, H., 1994. Computing discrete Fréchet distance.
    Tech. Report CD-TR 94/64, Information Systems Department, Technical
    University of Vienna.
    http://www.kr.tuwien.ac.at/staff/eiter/et-archive/cdtr9464.pdf
    Function dF(P, Q): real;
        input: polygonal curves P = (u1, . . . , up) and Q = (v1, . . . , vq).
        return: δdF (P, Q)
        ca : array [1..p, 1..q] of real;
        function c(i, j): real;
            begin
                if ca(i, j) > −1 then return ca(i, j)
                elsif i = 1 and j = 1 then ca(i, j) := d(u1, v1)
                elsif i > 1 and j = 1 then ca(i, j) := max{ c(i − 1, 1), d(ui, v1) }
                elsif i = 1 and j > 1 then ca(i, j) := max{ c(1, j − 1), d(u1, vj) }
                elsif i > 1 and j > 1 then ca(i, j) :=
                max{ min(c(i − 1, j), c(i − 1, j − 1), c(i, j − 1)), d(ui, vj ) }
                else ca(i, j) = ∞
                return ca(i, j);
            end; /* function c */
        begin
            for i = 1 to p do for j = 1 to q do ca(i, j) := −1.0;
            return c(p, q);
        end.
    Parameters
    ----------
    P : Input curve - two dimensional array of points
    Q : Input curve - two dimensional array of points
    Returns
    -------
    dist: float64
        The discrete Fréchet distance between curves `P` and `Q`.
    Examples
    --------
    >>> from frechetdist import frdist
    >>> P=[[1,1], [2,1], [2,2]]
    >>> Q=[[2,2], [0,1], [2,4]]
    >>> frdist(P,Q)
    >>> 2.0
    >>> P=[[1,1], [2,1], [2,2]]
    >>> Q=[[1,1], [2,1], [2,2]]
    >>> frdist(P,Q)
    >>> 0
    """
    p = np.array(p, np.float64)
    q = np.array(q, np.float64)

    len_p = len(p)
    len_q = len(q)

    if len_p == 0 or len_q == 0:
        raise ValueError('Input curves are empty.')

    if len_p != len_q or len(p[0]) != len(q[0]):
        raise ValueError('Input curves do not have the same dimensions.')

    ca = (np.ones((len_p, len_q), dtype=np.float64) * -1)

    dist = _c(ca, len_p-1, len_q-1, p, q)
    return dist


spe_colors = {'1':'gold','2':'blue', '3':'orange', '4':'green', '5':'red', '6':'purple', '7':'brown', '8':'pink',
              '9':'gray', '10':'olive','11':'cyan','12':'lightblue','13':'peru', '14':'navy','15':'limegreen',
              '16':'darkcyan','17':'orchid','18':'saddlebrown','19':'sliver','20':'black'}
states_list = [1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 13]
K79me2_state_list = [1,6,7,11,13]
celltype = nuceventfile.split('/')[-1].split('_')[0]
# read nuc_event and divide to SE and noSE according to mark col
nuc_event = pd.read_table(nuceventfile, sep='\t', header=None)
# expand event range to +-500 bp
ext_size = 500

event_mid = ((nuc_event[6] + nuc_event[7])/2).astype(int)
nuc_event[6] = event_mid - ext_size
nuc_event[7] = event_mid + ext_size
SE_event = nuc_event[nuc_event[15]!=0]
noSE_event = nuc_event[nuc_event[15]==0]

# separate by states
SE_event_dict = {}
noSE_event_dict = {}
for state in states_list:
    SE_event_dict[state] = SE_event[SE_event[3] == state].iloc[:,5:]
    noSE_event_dict[state] = noSE_event[noSE_event[3] == state].iloc[:,5:]

# read like-wig files
likewig_list = read_list_file(likewiglist)
likewig_dict = {}
for file in likewig_list:
    chrom = file.split('/')[-1].split('.')[0].split('_')[1]
    print("\rLoading {} Like-wig file.".format(chrom))
    shell_command = ["tail -n +4 {0}|awk -v var={1} 'BEGIN{{OFS=\"\\t\"}}{{print var,$1,$1+9,$3}}'".format(file,chrom)]
    p = subprocess.Popen(shell_command,stdout=subprocess.PIPE,shell=True)
    df = pd.read_csv(StringIO(p.communicate()[0].decode('utf-8')), sep="\t")
    likewig_dict[chrom] = df

# read srt uniq file
all_nuc = BedTool(srt_uniq)
final_event_SE,final_nuc_event_SE = filter_heter_event(all_nuc,SE_event_dict,states_list)
final_event_noSE, final_nuc_event_noSE = filter_heter_event(all_nuc,noSE_event_dict,states_list)

# get Gaussian smoothed value from likewig file according to the event coordinates
distribution_dict_SE = {state:np.zeros(100) for state in states_list}
distribution_dict_noSE = {state:np.zeros(100) for state in states_list}
for chrom in likewig_dict:
    print("process {}".format(chrom))
    intersect_file1 = BedTool.from_dataframe(likewig_dict[chrom])
    for state in K79me2_state_list:
        intersect_file2 = BedTool.from_dataframe(final_event_SE[state])
        intersect_file3 = BedTool.from_dataframe(final_event_noSE[state])
        G_intersect_result_SE = BedTool.intersect(intersect_file1,intersect_file2,wa=True,wb=True)
        G_intersect_result_noSE = BedTool.intersect(intersect_file1,intersect_file3,wa=True,wb=True)
        count = 0
        last_event = ''

        tmp_dict_SE = defaultdict(list)
        tmp_dict_noSE = defaultdict(list)
        for line in G_intersect_result_SE:
            event_name = line[4] + '_' + line[5] + '_' + line[6]
            tmp_dict_SE[event_name].append(float(line[3]))

        for line in G_intersect_result_noSE:
            event_name = line[4] + '_' + line[5] + '_' + line[6]
            tmp_dict_noSE[event_name].append(float(line[3]))

        for event in tmp_dict_SE:
            try:
                distribution_dict_SE[state] = np.add(distribution_dict_SE[state],np.array(tmp_dict_SE[event])[0:100])
            except ValueError:
                print(event,len(tmp_dict_SE[event]))

        for event in tmp_dict_noSE:
            try:
                distribution_dict_noSE[state] = np.add(distribution_dict_noSE[state],np.array(tmp_dict_noSE[event])[0:100])
            except ValueError:
                print(event,len(tmp_dict_noSE[event]))

# calculate frechetdist
frdist_dict = {}
for state in K79me2_state_list:
    reform_SE = []
    reform_noSE = []
    for idx,value in enumerate(distribution_dict_SE[state]/len(final_event_SE[state])):
        reform_SE.append([idx,value])
    for idx, value in enumerate(distribution_dict_noSE[state]/len(final_event_noSE[state])):
        reform_noSE.append([idx,value])
    frdist_dict[state] = frdist(reform_SE,reform_noSE)

print('frdist',frdist_dict)

plt.figure()
for state in K79me2_state_list:
    print(distribution_dict_SE[state])
    plt.plot(np.arange(100),distribution_dict_SE[state]/len(final_event_SE[state]),
             label=state,color=spe_colors[str(state)])
    plt.plot(np.arange(100),distribution_dict_noSE[state]/len(final_event_noSE[state]),
             label=state,color=spe_colors[str(state)],ls='--')
plt.legend()
plt.savefig(celltype+'_nuc_SE_event_distribution.png')

plt.figure()
for state in K79me2_state_list:
    print(distribution_dict_SE[state])
    plt.plot(np.arange(100),distribution_dict_noSE[state]/len(final_event_noSE[state]),
             label=state,color=spe_colors[str(state)])
plt.legend()
plt.savefig(celltype+'_nuc_noSE_event_distribution.png')

nuc_pos_dict = nuc_positioning(nucinfo,1.5)
nuc_pos_for_plot = []
plot_index = []
plot_color = []
for state in K79me2_state_list:
    tmp_SE = []
    tmp_noSE = []
    plot_index.append(str(state)+'_SE')
    plot_index.append(str(state)+'_noSE')
    plot_color.append(spe_colors[str(state)])
    plot_color.append(spe_colors[str(state)])
    for event in final_nuc_event_SE[state]:
        for nuc in final_nuc_event_SE[state][event]:
            tmp_SE.append(nuc_pos_dict[nuc])
    for event in final_nuc_event_noSE[state]:
        for nuc in final_nuc_event_noSE[state][event]:
            tmp_noSE.append(nuc_pos_dict[nuc])
    nuc_pos_for_plot.append(tmp_SE)
    nuc_pos_for_plot.append(tmp_noSE)

colored_box_plot(nuc_pos_for_plot,plot_index,'Nuc pos box',plot_color,celltype+'_nuc_pos_boxplot.png')


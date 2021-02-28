#!/usr/bin/env python
#--coding:utf-8 --

import re
import sys
import subprocess
import pandas as pd
from itertools import islice
from contextlib import contextmanager
from collections import defaultdict
from NucHMM_common_data import hg19,hg38,val2chr,info2strand
from NucHMM_utilities import ismember, get_time, process_read, ShowProcess
from NucHMM_Load_Write_files import load_genome_size_file,load_assignedtf_file,load_nuc_position_file,write_precomp_file, load_histonefile


def order_type_dict(inputfilelist):
    input_files = []
    diction_order2moditype = {}
    diction_moditype2order = {}
    count_order_output = 0
    if isinstance(inputfilelist,str):
        with open(inputfilelist, 'r') as inputfile_list:
            for line2 in inputfile_list:
                L2 = line2.strip()
                try:
                    modi_type = re.findall("(?i)H\d{1}K\d+(?:me|ac)+\d",line2)[0]
                except IndexError:
                    modi_type = L2.split('/')[-1].split('_')[1]
    
                diction_order2moditype[count_order_output] = modi_type
                diction_moditype2order[modi_type] = count_order_output
                input_files.append(r'''%s''' % (L2))
                count_order_output += 1
    elif isinstance(inputfilelist,list):
        for i in inputfilelist:
            with open(i, 'r') as inputfile_list:
                for line2 in inputfile_list:
                    L2 = line2.strip()
                    try:
                        modi_type = re.findall("(?i)H\d{1}K\d+(?:me|ac)+\d",line2)[0]
                    except IndexError:
                        modi_type = L2.split('/')[-1].split('_')[1]
        
                    diction_order2moditype[count_order_output] = modi_type
                    diction_moditype2order[modi_type] = count_order_output
                    input_files.append(r'''%s''' % (L2))
                    count_order_output += 1

    return diction_moditype2order, diction_order2moditype, input_files

def transcript_factor_assign(intersect_cutoff, gap, inputfilelist, outputfile, nucposition, mark_start):
    """tfassign is the sub-command that assign histone modification
    peaks to nucleosomes"""

    # celltype =
    # if outputfile is None:
    #

    tag = {}
    filename = []
    @contextmanager
    def open_many(files=None, mode='r'):
        if files is None:
            files = []
        try:
            fds = []
            count_t = 0
            for f in files:
                tmp_name = f.split('/')[-1].split("_")[1]
                filename.append(tmp_name)
                tag[tmp_name] = (1 << count_t + mark_start)
                count_t += 1
                fds.append(open(f, mode))
            yield fds
        except ValueError as e:
            print(e)
        finally:
            for fd in fds:
                fd.close()

    diction_moditype2order, diction_order2moditype, input_files = order_type_dict(inputfilelist)
    with open_many(input_files,'r') as bedfiles, \
         open(nucposition,'r') as pos_file, \
         open(outputfile,'w+') as outputfile:

        # BED[0][0][0] represent first bedfile, chromosome 1 and first read's middle position
        BED = []
        for i, f in enumerate(bedfiles):
            count_bed_line = 0
            diction_chr = defaultdict(list)
            print('\n'+filename[i]+'\t')
            process_bar = ShowProcess(input_files[i])
            for line in islice(f,0,None):
                count_bed_line += 1
                process_bar.show_process(count_bed_line)
                # sys.stdout.write('\rLoading bed line: ' + str(count_bed_line))
                L = line.strip().split()
                L_chr = L[0]
                try:
                    L_start = int(L[1])
                    L_end = int(L[2])
                except ValueError:
                    if count_bed_line == 1:
                        # default skip header
                        continue
                    else:
                        print("\nLine {} column 2 or 3 is not genomic coordinates".format(count_bed_line))
                        continue

                # skip chrM and random scaffoldd
                if L_chr == 'chrM' or len(L_chr) > 5:
                    continue

                if L_start in diction_chr[L_chr]:
                    continue
                else:
                    diction_chr[L_chr].append(L_start)
                    diction_chr[L_chr].append(L_end)
            BED.append(diction_chr)
            f.close()

        print('\nRead position file')
        pos_file.seek(0)
        count_pos_line = 0
        diction_pos = defaultdict(list)
        diction_pos_start = defaultdict(list)
        diction_pos_end = defaultdict(list)
        number_of_BED = len(BED)
        process_bar = ShowProcess(nucposition)
        for pos in pos_file:
            count_pos_line += 1
            process_bar.show_process(count_pos_line)
            # sys.stdout.write('\rPosition file line: ' + str(count_pos_line))
            P = pos.strip().split()
            P_chr = P[0]
            P_start = int(P[1])
            P_end = int(P[2])
            diction_pos_start[P_chr+'_'+str(P_start)] = [0] * number_of_BED
            diction_pos_end[P_chr+'_'+str(P_end)] = [0] * number_of_BED
            if P_chr == 'chrM' or len(P_chr) > 5:
                continue
            diction_pos[P_chr].append(P_start)
            diction_pos[P_chr].append(P_end)

        print('\ncalculate the reads')
        for index,his_modification in enumerate(BED):
            for keys in his_modification:
                Mix_seg = sorted(his_modification[keys] + diction_pos[keys])
                peak_index = ismember(his_modification[keys],Mix_seg)
                for peak_i in range(len(peak_index)//2):
                    '''reason +1: a = [2,3,4,5], a[2:4] = [3,4]'''
                    peak_i_start_index = peak_index[peak_i*2]+1
                    peak_i_end_index = peak_index[peak_i*2+1]
                    peak_i_start = Mix_seg[peak_i_start_index-1]
                    peak_i_end = Mix_seg[peak_i_end_index]
                    for position in Mix_seg[slice(peak_i_start_index,peak_i_end_index)]:
                        if keys+'_'+str(position) in diction_pos_start:
                            diction_pos_start[keys+'_'+str(position)][index] = peak_i_end - position
                        elif keys+'_'+str(position) in diction_pos_end:
                            diction_pos_end[keys+'_'+str(position)][index] = position - peak_i_start
                        else:
                            print('Please sort -k1,1V -k2,2n -k3,3n -u your input peak file and make sure there is no overlapped peaks (use bedtools merge)')
                            print("Detected overlapped/duplicated peak {} {} {} {}".format(input_files[index],keys,peak_i_start,peak_i_end))
                            exit(1)

        print('Histone marks: transfered mark')
        for keys in tag:
            print(keys+': '+str(tag[keys]))

        process_bar = ShowProcess(nucposition)
        print('\nWrite output file')
        count_out_line = 0
        pos_file.seek(0)
        L_last_end = 0
        for line in pos_file:
            L = line.strip().split()
            Chr = L[0]
            start = int(L[1])
            end = int(L[2])
            Nuc_length = end - start
            mark = 0
            if count_out_line == 0:
                L_last_end = end
            for i in range(number_of_BED):
                if diction_pos_start[Chr+'_'+str(start)][i] > 0 and diction_pos_end[Chr+'_'+str(end)][i] > 0:
                    mark += pow(2,i+mark_start)
                elif diction_pos_start[Chr+'_'+str(start)][i] >= Nuc_length * intersect_cutoff or \
                    diction_pos_end[Chr+'_'+str(end)][i] >= Nuc_length * intersect_cutoff:
                    mark += pow(2,i+mark_start)
                else:
                    mark += 0

            process_bar.show_process(count_out_line)
            # sys.stdout.write('\routput line: ' + str(count_out_line))
            if gap:
                nucleosome_spacing_mark = 2 ** (len(BED)+mark_start)
                if count_out_line == 0:
                    outputfile.write(Chr + '\t' + str(start) + '\t' + str(end) + '\t' + str(mark) + '\n')
                    continue
                else:
                    outputfile.write(Chr + '\t' + str(L_last_end) + '\t' + str(start) + '\t' + str(
                        nucleosome_spacing_mark) + '\n')
                    outputfile.write(Chr + '\t' + str(start) + '\t' + str(end) + '\t' + str(mark) + '\n')
            else:
                if mark_start == 0:
                    outputfile.write(Chr + '\t' + str(start) + '\t' + str(end) + '\t' + str(mark) + '\n')
                else:
                    outputfile.write(Chr + '\t' + str(start) + '\t' + str(end) + '\t' + str(mark + 2 ** (len(BED)+mark_start)) + '\n')
            count_out_line += 1

    pos_file.close()
    outputfile.close()
    return

def gap_position(nucposition):
    nuc_pos = pd.read_csv(nucposition,sep='\t',header=None)[[0,1,2]]
    gap_pos = pd.concat([nuc_pos[[0]][:-1],nuc_pos[[2]][:-1],nuc_pos[[1]][1:].reset_index()[[1]]],axis=1)
    return gap_pos
    
def ciselements_hm_assgin(intersect_cutoff, inputfilelist, cislist, outputfile, nucposition):
    """assign histone modification peaks to nucleosomes and cis-element peaks to nucleosomes and gaps"""

    # calculate nucleosome marks
    print("Assign marks to nucleosomes..")
    tmpinputfile1 = [inputfilelist, cislist]
    tmpoutfile1 = 'tmp_nuc_assign.bed'
    transcript_factor_assign(intersect_cutoff, False, tmpinputfile1, tmpoutfile1, nucposition, 0)

    # calculate gap marks
    print("\nAssgin marks to gaps")
    tmpinputfile2 = cislist
    tmpoutfile2 = 'tmp_gap_assign.bed'
    gap_pos = gap_position(nucposition)
    tmpgappos = 'tmp_gap_pos.bed'
    gap_pos.to_csv(tmpgappos, header=False, index=False, sep='\t')
    transcript_factor_assign(intersect_cutoff, False, tmpinputfile2, tmpoutfile2, tmpgappos, len(load_histonefile(inputfilelist)))

    print("\nFinalize the process..")
    # combine nucleosome mark file and gap mark file
    subprocess.call("paste -d \\\\n {} {} > {}".format(tmpoutfile1,tmpoutfile2,outputfile),shell=True)
    # remove tmp files
    subprocess.call("rm {} {} {}".format(tmpoutfile1, tmpoutfile2, tmpgappos), shell=True)

def segment_gene(genefile, segoutfile, up, down, refgene):
    """geneseg is sub-command for selecting effect region
    and segment it to three part with directionality label."""

    up_boundary = up
    down_boundary = down
    gene_id = defaultdict(list)
    gene_id2 = defaultdict(list)


    with open(genefile,'r') as inputfile, open(segoutfile,'w') as outputfile:
        d_count ={}
        L_seg =[]
        L_SEG = []
        L_last_chr = 'chr1'
        count = 0
        process_bar = ShowProcess(genefile)
        for line in inputfile:
            count += 1
            process_bar.show_process(count)
            # sys.stdout.write('\rread input: '+str(count) )
            L = line.strip().split()
            L_chr = L[0]
            L_strand = L[3]
            L_name = L[4]
            d_count[L_name] = 0
            if L_strand == '+':
                L_start = int(L[1])
                L_end = int(L[2])
                L_up = L_start - up_boundary
                if L_up < 0:
                    L_up = 0
                L_down = L_end + down_boundary
                if L_down > refgene[L_chr]:
                    L_down = refgene[L_chr]
                str1 = [L_name,L_strand]
                gene_id[str(L_start)+'_'+L_chr].append('_'.join(str1))
                gene_id2[str(L_end)+'_'+L_chr].append('_'.join(str1))
            else:
                L_start = int(L[2])
                L_end = int(L[1])
                L_up = L_start + up_boundary
                if L_up > refgene[L_chr]:
                    L_up = refgene[L_chr]
                L_down = L_end - down_boundary
                if L_down < 0:
                    L_down = 0
                str1 = [L_name,L_strand]
                gene_id[str(L_start)+'_'+L_chr].append('_'.join(str1))
                gene_id2[str(L_end)+'_'+L_chr].append('_'.join(str1))
            if L_chr == L_last_chr:
                if L_strand == '+':
                    L_seg.extend([L_up, L_start, L_end, L_down])
                else:
                    L_seg.extend([L_down,L_end,L_start,L_up])
            else:
                L_seg.append(refgene[L_last_chr])
                L_SEG.append(L_seg)
                L_seg = []
                if L_strand == '+':
                    L_seg.extend([L_up, L_start, L_end, L_down])
                else:
                    L_seg.extend([L_down,L_end,L_start,L_up])
                L_last_chr = L_chr

        L_seg.append(refgene['chrY'])
        L_SEG.append(L_seg)

        N_Lseg = []
        for index,chrom in enumerate(L_SEG):
            # chrom_ss = sorted(chrom)
            if chrom[-1] > refgene[val2chr[index+1]]:
                chrom.remove(chrom[-1])
            if chrom[-1] > refgene[val2chr[index+1]]:
                chrom.remove(chrom[-1])
            N_Lseg.append(chrom)

        count_num = 0
        dessert_count = 0
        for index,chrom in enumerate(N_Lseg):
            for i in range(len(N_Lseg[index])-1):
                # chrom_s = sorted(chrom)
                tmp_start = chrom[slice(i,i+2)][0]
                start = str(tmp_start)+'_'+val2chr[index+1]
                tmp_end = chrom[slice(i,i+2)][1]
                end = str(tmp_end)+'_'+val2chr[index+1]
                if start in gene_id:
                    info = gene_id[start][0].split('_')
                    if info[1] == '-':
                        d_count[info[0]] += 1
                        outputfile.write(val2chr[index+1]+'\t'+str(tmp_start)+'\t'+str(tmp_end)+
                                         '\t'+'+'+'\t'+info[0]+'_'+info2strand[info[1]]+'_'+str(d_count[info[0]])+'\n')
                    else:
                        d_count[info[0]] += 1
                        outputfile.write(val2chr[index+1]+'\t'+str(tmp_start)+'\t'+str(tmp_end)+
                                         '\t'+info[1]+'\t'+info[0]+'_'+info2strand[info[1]]+'_'+str(d_count[info[0]])+'\n')
                elif end in gene_id:
                    info = gene_id[end][0].split('_')
                    d_count[info[0]] += 1
                    outputfile.write(val2chr[index+1]+'\t'+str(tmp_start)+'\t'+str(tmp_end)+
                                     '\t'+'-'+'\t'+info[0]+'_'+info2strand[info[1]]+'_'+str(d_count[info[0]])+'\n')
                elif start in gene_id2:
                    info = gene_id2[start][0].split('_')
                    d_count[info[0]] += 1
                    outputfile.write(val2chr[index+1]+'\t'+str(tmp_start)+'\t'+str(tmp_end)+
                                     '\t'+info[1]+'\t'+info[0]+'_'+info2strand[info[1]]+'_'+str(d_count[info[0]])+'\n')
                elif end in gene_id2:
                    info = gene_id2[end][0].split('_')
                    d_count[info[0]] += 1
                    if info[1] == '+':
                        outputfile.write(val2chr[index+1]+'\t'+str(tmp_start)+'\t'+
                                         str(tmp_end)+'\t'+'-'+'\t'+info[0]+'_'+info2strand[info[1]]+'_'+str(d_count[info[0]])+'\n')
                    else:
                        outputfile.write(val2chr[index+1]+'\t'+str(tmp_start)+
                                         '\t'+str(tmp_end)+'\t'+'-'+'\t'+info[0]+'_'+info2strand[info[1]]+'_'+str(d_count[info[0]])+'\n')
                elif tmp_end - tmp_start > 0:
                    dessert_count += 1
                    outputfile.write(val2chr[index+1]+'\t'+str(tmp_start)+
                                     '\t'+str(tmp_end)+'\t'+' '+'\t'+'dessert'+'_'+str(dessert_count)+'\n')
    inputfile.close()
    outputfile.close()

def precompile_step(assignedfilelist, segfile, outputfilelist, refgene):
    """precompile is the sub-command that prepare the input file
    for the HMM model, this should be done after finishing
    geneseg and cell type tfassign.
    """
    @contextmanager
    def open_many(files=None, mode='r'):
        if files is None:
            files = []
        try:
            fds = []
            for f in files:
                fds.append(open(f, mode))
            yield fds
        except ValueError as e:
            print(e)
        finally:
            for fd in fds:
                fd.close()

    output_file =[]
    diction_order2celltype = {}
    diction_celltype2order ={}
    count_order_output = 0
    with open(outputfilelist, 'r') as outputfile_list:
        for line2 in outputfile_list:
            L2 = line2.strip()
            cell_type = line2.split('/')[-1].split('_')[0]
            diction_order2celltype[count_order_output] = cell_type
            diction_celltype2order[cell_type] = count_order_output
            output_file.append(r'''%s'''%(L2))
            count_order_output += 1

    with open_many(output_file,'w') as outputfiles, open(assignedfilelist, 'r') as inputfile, \
            open(segfile, 'r') as seg_file:
        POS = {}
        READ = {}
        Start_End = {}
        for line in inputfile:
            L = line.strip()
            with open(r'''%s'''%(L), 'r') as raw_read_count:
                cell_type = line.split('/')[-1].split('_')[0]
                '''read the nucleosomes position(mid) and its modification mark'''
                sys.stdout.write(str(cell_type) + '\n')
                POS[diction_celltype2order[cell_type]], READ[diction_celltype2order[cell_type]], \
                Start_End[diction_celltype2order[cell_type]] = load_assignedtf_file(raw_read_count)

        number_of_cell = len(POS)

        '''read the gene segment start and end'''
        d_mark = {}
        dict_seg = defaultdict(list)
        L3_last_chr = 'chr1'
        L_SEG = defaultdict(list)
        count_input3 = 0
        process_bar = ShowProcess(segfile)
        for line3 in seg_file:
            count_input3 += 1
            process_bar.show_process(count_input3)
            # sys.stdout.write('\rSegment file line: ' + str(count_input3))
            L3 = line3.strip().split()
            L3_chr = L3[0]
            L3_start = int(L3[1])
            L3_end = int(L3[2])
            if len(L3) == 5:
                L3_mark = L3[3]
                L3_name = L3[4]
            else:
                L3_name = L3[3]
                L3_mark = 'NAN'

            L_SEG[L3_chr].append(L3_start)
            L_SEG[L3_chr].append(L3_end)
            str1 = [L3_name, L3_mark]
            seg_name = '_'.join(str1)
            d_mark[L3_name] = L3_mark
            dict_seg[seg_name].append(L3_start)
            dict_seg[seg_name].append(L3_end)

        for keys in refgene:
            L_SEG[keys].append(refgene[keys])

        print("\nExtract nucleosomes information from the selected region")
        seg_file.seek(0)
        count_input4 = 0
        pos_dict = load_nuc_position_file(L_SEG, POS)
        process_bar = ShowProcess(segfile)
        for line in seg_file:
            count_input4 += 1
            process_bar.show_process(count_input4)
            # sys.stdout.write('\rReading line: ' + str(count_input4))
            L3 = line.strip().split()
            L3_chr = L3[0]

            if L3_chr == "chrX":
                chr_index = 22
            elif L3_chr == "chrY":
                chr_index = 23
            else:
                chr_index = int(re.findall(r"\d+", L3_chr)[0]) - 1

            if len(L3) == 5:
                L3_mark = L3[3]
                L3_name = L3[4]
            else:
                continue

            str1 = [L3_name, L3_mark]
            seg_name = '_'.join(str1)
            d_mark[L3_name] = L3_mark
            geneseg = dict_seg[seg_name][0:2]
            Num, Read_seq, S_and_E = process_read(geneseg, pos_dict, POS, READ, Start_End, L3_chr)
            for i in range(number_of_cell*3):
                dict_seg[seg_name].append(0)
            for cellnum in range(number_of_cell):
                dict_seg[seg_name][2+cellnum] = Num[cellnum]
                dict_seg[seg_name][2+number_of_cell+cellnum] = Read_seq[cellnum]
                dict_seg[seg_name][2+2*number_of_cell+cellnum] = S_and_E[cellnum]

        '''dict_seg:[segment_start,segment_end,cell1_Num_in_seg,cell2_Num_in_seg,cell3_Num_in_seg,...,cell1_Read_seq(mark),
        cell2_read_seq,...,cell1_Start_and_END_for_the_mark,cell2_Start_and_END_for_the_mark]'''
        for cellnum, f in enumerate(outputfiles):
            seg_file.seek(0)
            count_out = 0
            count_read = 0
            last_chr_index = 0
            count_line = 0
            tmp_read = []
            tmp_pos = []
            process_bar = ShowProcess(segfile)
            print("Writing output file:")
            for line4 in seg_file:
                count_line += 1
                process_bar.show_process(count_line)
                # sys.stdout.write("\rWriting output %: " + str(count_line))
                L4 = line4.strip().split()
                L4_chr = L4[0]
                if L4_chr == "chrX":
                    chr_index = 22
                elif L4_chr == "chrY":
                    chr_index = 23
                else:
                    chr_index = int(re.findall(r"\d+", L4_chr)[0]) - 1

                if len(L4) == 5:
                    L4_mark = L4[3]
                    L4_name = L4[4]
                    direction = L4_name.split('_')[1]
                    gene_name = L4_name.split('_')[0]
                    num_seg = int(L4_name.split('_')[2])
                else:
                    continue

                str1 = [L4_name, L4_mark]
                seg_name = '_'.join(str1)
                if last_chr_index != chr_index:
                    count_read = 0
                    count_out = 0
                    count_out_m = 0
                    count_read_m = 0
                    last_chr_index = chr_index

                temp_read,temp_pos,out_count = write_precomp_file(direction, dict_seg, seg_name, count_out, L4_mark,
                                                                  cellnum, number_of_cell, chr_index, num_seg, tmp_read,
                                                                  tmp_pos,f)
                tmp_read = temp_read
                tmp_pos = temp_pos
                count_out = out_count

    for file in outputfiles:
        file.close()

    seg_file.close()
    inputfile.close()
    outputfile_list.close()

import os
import re
import sys
import subprocess


def trimmed_name(file,pe_mark,R_mark):
    if pe_mark:
        current_name = file.split('/')[-1]
        current_name_mainbody = re.findall('\w+(?=.R1|.R2)',current_name)[0]
        if current_name.split('.')[-1] == 'gz':
            current_trimmed_name = 'fq_trim_result/' + current_name_mainbody \
                                   + '_R' + R_mark + '_val_' + R_mark + '.fq.gz'
        else:
            current_trimmed_name = 'fq_trim_result/' + current_name_mainbody \
                                   + '_R' + R_mark + '_val_' + R_mark + '.fq'
    else:
        current_name = file.split('/')[-1]
        current_name_mainbody = re.findall('\w+(?=.fq|.fastq)',current_name)[0]
        current_trimmed_name = 'fq_trim_result/' + current_name_mainbody \
                               + '_trimmed' + current_name.replace(current_name_mainbody,'')
        if not os.path.exists(current_trimmed_name):
            current_trimmed_name = 'fq_trim_result/' + current_name_mainbody \
                               + '_trimmed' + current_name.replace(current_name_mainbody,'').replace('fastq','fq')
        if not os.path.exists(current_trimmed_name):
            print("Trimmed_name Error, please mannually change trimmed name to {}".format(current_trimmed_name))
            exit(1)
    return current_trimmed_name

def bam_name(file,pe):
    if pe:
        current_name = file[0].split('/')[-1]
        current_name_mainbody =  re.findall('\w+(?=.R1|R2)',current_name)[0]
        final_bam_name = 'bam_result/' + current_name_mainbody + '.bam'
    else:
        current_name = file.split('/')[-1]
        current_name_mainbody = re.findall('\w+(?=.trimmed)',current_name)[0]
        final_bam_name = 'bam_result/' + current_name_mainbody + '.bam'
    return final_bam_name

def QC_step(info_dict,threads):
    '''process raw fastq file'''
    if not os.path.exists('fq_trim_result'):
        os.makedirs('fq_trim_result')
    filt_files_list = []
    print("Starting QC trimming process..")
    inputfiles_list = info_dict['fq_list']
    for files in inputfiles_list:
        if len(files) == 2:
            subprocess.call('trim_galore --trim-n --paired -j '+ str(threads) + ' -o fq_trim_result ' + files[0] + ' ' + files[1], shell=True)
            current_R1_trim = trimmed_name(files[0],True,'1')
            current_R2_trim = trimmed_name(files[1],True,'2')
            filt_files_list.append([current_R1_trim,current_R2_trim])
        else:
            print(files[0])
            subprocess.call('trim_galore --trim-n -j '+ str(threads) + ' -o fq_trim_result ' + files[0], shell=True)
            current_trim = trimmed_name(files[0],False,None)
            filt_files_list.append([current_trim])
    print("QC trimming finished!")
    info_dict['filt_fq_list'] = filt_files_list
    return info_dict

def Mapping_step(fqinfodict,Bindex,B2index,threads,mapq_cutoff):
    '''align reads to reference genome, inputfileslist is trimmed fq list or fq list, Bindex is Bowtie index path,
    B2index is Bowtie2 index'''

    # determine which aligner to use
    aligner = []
    readlenlist = fqinfodict['read_len_list']
    for len_info in readlenlist:
        if len_info >= 50:
            aligner.append('Bowtie2')
        else:
            aligner.append('Bowtie')
    aligner = list(set(aligner))
    Bindexpath = ''
    B2indexpath = ''
    if len(aligner) == 1:
        if aligner[0] == 'Bowtie':
            if Bindex is None:
                Bindexpath = input("Please input the full path of the Bowtie index path:")
            else:
                Bindexpath = Bindex
        else:
            if B2index is None:
                B2indexpath = input("Please input the full path of the Bowtie2 index path:")
            else:
                B2indexpath = B2index
    else:
        if B2index is None:
            B2indexpath = input("Please input the full path of the Bowtie2 index path:")
        else:
            B2indexpath = B2index
        if Bindex is None:
            Bindexpath = input("Please input the full path of the Bowtie index path:")
        else:
            Bindexpath = Bindex

    if not os.path.exists('bam_result'):
        os.makedirs('bam_result')

    print("Starting mapping process..")
    inputfileslist = fqinfodict['filt_fq_list']
    aligned_file_list = []
    for idx,files in enumerate(inputfileslist):
        print('mapping','\t'.join(files))
        if len(files) == 2:
            # paired-end
            out_bam_name = bam_name(files,True)
            print("out bam file:", out_bam_name)
            if readlenlist[idx] >= 50:
                subprocess.call('bowtie2 -p ' + str(threads) + ' -x ' + B2indexpath + ' -1 ' + files[0] + ' -2 ' +
                                files[1] + ' -S |samtools view -F 1804 -f 2 -q ' + str(mapq_cutoff) +
                                ' -b - | samtools sort -O BAM -o ' + out_bam_name + ' - ', shell=True)
            else:
                subprocess.call('bowtie -p ' + str(threads) + ' ' + Bindexpath + ' -1 ' + files[0] + ' -2 ' + files[1] +
                                ' --chunkmbs 200 -S |samtools view -F 1804 -f 2 -q ' + str(mapq_cutoff) +
                                ' -b - | samtools sort -O BAM -o ' + out_bam_name + ' - ', shell=True)
            aligned_file_list.append(out_bam_name)
        else:
            # single-end
            out_bam_name = bam_name(files[0],False)
            if readlenlist[idx] >= 50:
                subprocess.call('bowtie2 -p ' + str(threads) + ' -x ' + B2indexpath + ' -U ' + files[0] +
                                ' -S |samtools view -F 1804 -q ' + str(mapq_cutoff) + ' -b - | samtools sort -O BAM -o ' +
                                out_bam_name + ' - ', shell=True)
            else:
                subprocess.call('bowtie -m 1 -p ' + str(threads) + ' ' + Bindexpath + ' -q ' + files[0] +
                                ' -S |samtools view -F 1804 -q ' + str(mapq_cutoff) + ' -b - | samtools sort -O BAM -o ' +
                                out_bam_name + ' - ', shell=True)
            aligned_file_list.append(out_bam_name)
    fqinfodict['aligned_file_list'] = aligned_file_list
    return fqinfodict

def Nuc_calling(bamfile,pe,inpspath):
    # filter reads within the fragment size range [130,180]
    bamfile_body = re.findall('[^\/]+$',bamfile)[0][:-4]
    if pe=='PE':
        print("Screen out unfitted reads..")
        filt_bam = bamfile_body + '_screened.bam'
        if os.path.exists(bamfile+'.bai'):
            subprocess.call('alignmentSieve -p 4 -b ' + bamfile + ' -o ' + filt_bam +
                            ' --minFragmentLength 135 --maxFragmentLength 175',shell=True)
        else:
            subprocess.call('samtools index -@ 4 '+ bamfile,shell=True)
            subprocess.call('alignmentSieve -p 4 -b ' + bamfile + ' -o ' + filt_bam +
                            ' --minFragmentLength 135 --maxFragmentLength 175',shell=True)
    elif pe=='SE':
        filt_bam = bamfile
    else:
        print("Please indicate paired-end or single-end by PE or SE.")
        exit(1)
    # transfer bamtobed
    print("Transfer bam to bed..")
    bed_file = bamfile_body + '.bed'
    if os.path.exists(bed_file):
        print("Use the existing bed file")
    else:
        subprocess.call('bedtools bamtobed -i ' + filt_bam + ' > ' + bed_file, shell=True)
    # use iNPS calling nucleosome core location
    print("iNPS Calling..")
    if pe=='PE':
        subprocess.call('python ' + inpspath + ' -i ' + bed_file +
                        ' -o nuc_calling_result/' + bamfile_body + ' --s_p p',shell=True)
    else:
        subprocess.call('python ' + inpspath + ' -i ' + bed_file +
                        ' -o nuc_calling_result/' + bamfile_body + ' --s_p s',shell=True)
    # remove header of Gathering.like_bed
    gathering_file = 'nuc_calling_result/' + bamfile_body + '_Gathering.like_bed'
    count = 1
    count2= 1
    with open(gathering_file,'r') as g_file:
        for line in g_file:
            line_info = line.strip().split()
            try:
                line_start = int(line_info[1])
            except ValueError:
                count += 1
                continue
            count2 += 1
            if count2 >100:
                break
    g_file.close()
    count += 1
    final_nuc_file = 'nuc_calling_result/' + bamfile_body + '_nucleosome_location.bed'
    subprocess.call('tail -n +' + str(count) + ' ' + gathering_file + ' > ' + final_nuc_file,shell=True)
    return final_nuc_file

def Peak_Nuc_calling_step(fqinfodict,inpspath):
    '''using macs2 to call narrow peaks and epic2 to call broad peaks'''
    inputbamfile = fqinfodict['aligned_file_list']
    peak_mark_list = fqinfodict['peak_type_list']
    pe_mark_list = fqinfodict['pe_mark_list']
    output_file_list = []
    for idx,bamfile in enumerate(inputbamfile):
        out_peak_name_mainbody = bamfile.split('/')[-1].split('.')[0]
        if peak_mark_list[idx] == 'broad':
            print("Calling " + bamfile)
            if not os.path.exists('peakcalling_result'):
                 os.makedirs('peakcalling_result')
            subprocess.call('epic2 -t '+ bamfile +' -bin 100 -g 2 -fdr 0.05 -fs 200 -o peakcalling_result/' + out_peak_name_mainbody
                            + '-epic2-bin100-g2-fdr005.bed 2> ' + out_peak_name_mainbody+'.log', shell= True)
            output_file_list.append('peakcalling_result/' + out_peak_name_mainbody + '-epic2-bin100-g2-fdr005.bed')
        elif peak_mark_list[idx] == 'narrow':
            print("Calling " + bamfile)
            if not os.path.exists('peakcalling_result'):
                os.makedirs('peakcalling_result')
            # print("MACS2 command: ",'macs2 callpeak -t ' + bamfile + ' -f BAM -q 0.01 -g hs -n ' + out_peak_name_mainbody +
            #                 ' --outdir peakcalling_result/' + out_peak_name_mainbody +'-peaks-q001 2> ' + out_peak_name_mainbody +'.log')
            if pe_mark_list[idx] == 'SE':
                subprocess.call('macs2 callpeak -t ' + bamfile + ' -f BAM -q 0.01 -g hs -n ' + out_peak_name_mainbody +
                                ' --outdir peakcalling_result/' + out_peak_name_mainbody +'-peaks-q001 2> ' + out_peak_name_mainbody +'.log', shell=True)
            elif pe_mark_list[idx] == 'PE':
                subprocess.call('macs2 callpeak -t ' + bamfile + ' -f BAMPE -q 0.01 -g hs -n ' + out_peak_name_mainbody +
                                ' --outdir peakcalling_result/' + out_peak_name_mainbody +'-peaks-q001 2> ' + out_peak_name_mainbody +'.log', shell=True)
            else:
                print("Please indicate paired-end or single-end by PE or SE.")
                exit(1)
            output_file_list.append('peakcalling_result/'+out_peak_name_mainbody+'-peaks-q001/'+out_peak_name_mainbody+'_peaks.narrowPeak')
        elif peak_mark_list[idx] == 'none':
            print("Calling " + bamfile)
            if not os.path.exists('nuc_calling_result'):
                os.makedirs('nuc_calling_result')
            if inpspath == None:
                inps_path = input("Please input iNPS.py full path (includes iNPS.py): ")
            else:
                inps_path = inpspath
            nuc_file = Nuc_calling(bamfile,pe_mark_list[idx],inps_path)
            output_file_list.append('nuc_calling_result/'+nuc_file)
        else:
            print("Unknown peak format, please check the input file.")
            exit(1)
        # sys.std.write('\r'+str(idx+1)+"/",len(inputbamfile))
    fqinfodict['final_file_list'] = output_file_list
    return fqinfodict

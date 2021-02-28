#!/usr/bin/env python
#--coding:utf-8 --

import sys, time
import os.path
import subprocess
import numpy as np
from datetime import datetime


def ismember(a, b):
    '''
    find the index of the element that co-exist in a and b, e.g a =[1,2] ,b =[1,2,3] ismember(a,b) => [0,1]
    :param a: a son list and normally is subset of b
    :param b: a parent list
    :return: a list of index
    '''
    bind = {}
    for i, elt in enumerate(b):
        if elt not in bind:
            bind[elt] = i
    return [bind.get(itm, None) for itm in a]


def get_time():
    # get the current time
    now = datetime.now()
    dt_string = now.strftime("%d.%m.%Y_%H.%M.%S")
    return dt_string

def Chr2Num (Chr):
    '''
    Transfer chromosome to number
    :param Chr: input chromosome, e.g chr1
    :return: corresponding number of chromosomes
    '''
    if Chr:
        New = Chr[3:]
        if New == 'X': New = 23
        elif New == 'Y': New = 24
        elif New == 'M': New = 25
        else: New = int(New)
    else:
        New = 0
    return New

def process_read(geneseg, pos_diction, POS, READ, Start_End, chr):
    '''
    orignial in precompile command
    :param geneseg:
    :param pos_diction:
    :param POS:
    :param READ:
    :param Start_End:
    :param chr:
    :return:
    '''
    Pos_index = {}
    Num = {}
    Mid = {}
    S_and_E = {}
    Read_seq = {}
    number_of_cell = len(POS)
    for cellnum in range(number_of_cell):
        Pos_index[cellnum] = [pos_diction[str(geneseg[0])+'_'+chr][cellnum], pos_diction[str(geneseg[1])+'_'+chr][cellnum+number_of_cell]]
        Num[cellnum] = Pos_index[cellnum][1] - Pos_index[cellnum][0] + 1
        Mid[cellnum] = POS[cellnum][chr][slice(Pos_index[cellnum][0], Pos_index[cellnum][0] + Num[cellnum])]
        S_and_E[cellnum] = [Start_End[cellnum][chr + '_' + str(k)] for k in Mid[cellnum]]
        for k in Mid[cellnum]:
            S_and_E[cellnum].append(Start_End[cellnum][chr + '_' + str(k)])
            if Start_End[cellnum][chr + '_' + str(k)] ==[]:
                print('pay attention: '+str(k)+'\t'+str(chr)+'\t'+str(Mid[cellnum])+'\t'+str(slice(Pos_index[cellnum][0], Pos_index[cellnum][0] + Num[cellnum]))+'\n')
                exit(1)

        Read_seq[cellnum] = READ[cellnum][chr][slice(Pos_index[cellnum][0], Pos_index[cellnum][0] + Num[cellnum])]
    return Num, Read_seq, S_and_E

def sort_chrom_coor1_coor2(input):
    '''
    Sort the key in dictionary according the chromosome -> start -> end (chr1,chr2...chr10...)
    :param input: input diction or list with key or element format 'chr1_1_2'
    :return: sorted_key_list
    '''
    sorted_key = sorted(input, key=lambda x: (Chr2Num(x.partition('_')[0]),int(x.partition('_')[2].partition('_')[0]),
                                                int(x.partition('_')[2].partition('_')[2].partition('_')[0])))
    return sorted_key


def query_yes_no(question, default="yes"):
    # Idea from https://stackoverflow.com/questions/3041986/apt-command-line-interface-like-yes-no-input
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        if sys.version_info <(3,0):
            choice = raw_input().lower()
        else:
            choice = input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")


def file_check(file):
    '''
    check if the file exists
    :param file: The path and name of the file
    :return: True/False
    '''
    return os.path.exists(file)

def input_new_name(question):
    if sys.version_info < (3,0):
        new_name = raw_input(question)
    else:
        new_name = input(question)
    return new_name


def boundary_check(boundary, chrom, ref):
    '''
    Check if the coordinate outside the chromosome range
    :param boundary: input coordinates
    :param chrom: chromosome name
    :param ref: a dictionary contain chromosome size information
    :return: checked boundary
    '''
    if boundary < 0:
        new_boundary = 0
    elif boundary > ref[chrom]:
        new_boundary = ref[chrom]
    else:
        new_boundary = boundary
    return new_boundary


def bedops_map(basefile,tobeassignedfile,mappedfile):
    '''
    Use bedmap assign state to segmented region file
    :param basefile: segmented region file
    :param tobeassignedfile: nucleosome-state file
    :param mappedfile: mapped final result
    '''
    subprocess.call("bedmap --echo --echo-map-id --delim '\\t' " + basefile + " " + tobeassignedfile +
                        " > " + mappedfile, shell=True)

def bedops_diff(basefile,filterfiles,filteredfile):
    '''
    Use bedop -d to filter out some regions in the basefile
    :return:
    '''
    subprocess.call("bedops -d " + basefile + ' ' + ' '.join(filterfiles) + " > " + filteredfile, shell=True)

def count_file_rows(inputfile):
    '''
    count how many lines in the file
    :param inputfile:
    :return:
    '''
    count = int(subprocess.check_output("wc -l "+inputfile, shell=True).split()[0])
    return count

def filter_info_print(feature,count_nuc,count_filt_nuc):
    print("%d nucleosomes before %s filtering " % (count_nuc, feature))
    print("%d nucleosomes after %s filtering " % (count_filt_nuc, feature))
    ratio = round(float(count_filt_nuc)/count_nuc,4)*100
    print("%f percent nucleosomes remaining " % ratio)

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

def which(program):
    """
    Check if program exists on path. Adapted from MISO.
    """
    def is_exe(fpath):
        if not os.path.isfile(fpath):
            return False
        elif not os.access(fpath, os.X_OK):
            # If the file exists but is not executable, warn
            # the user
            print("WARNING: Found %s but it is not executable." %(fpath))
            print("Please ensure %s is executable." %(fpath))
            print("On Unix, use something like: ")
            print("  chmod +x %s" %(fpath))
            time.sleep(10)
            return False
        return True
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None

class ShowProcess:
    """
    The class of handling the process
    """
    # init function, need to know the max steps that need to process
    def __init__(self, file, max_arrow=50, infoDone = ''):
        self.max_steps = int(subprocess.check_output('wc -l ' + file, shell=True).split()[0])
        self.i = 0
        self.infoDone = infoDone
        self.max_arrow = max_arrow

    # the function of showing the process
    # show like [>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>]100.00%
    def show_process(self, i=None):
        if i is not None:
            self.i = i
        else:
            self.i += 1
        num_arrow = int(self.i * self.max_arrow / self.max_steps) #calculate how many '>' to display
        num_line = self.max_arrow - num_arrow # calculate how many '-' to display
        percent = self.i * 100.0 / self.max_steps # calculate the finished part，format xx.xx%
        process_bar = '[' + '>' * num_arrow + '-' * num_line + ']'\
                      + '%.2f' % percent + '%' + '\r' # output
        sys.stdout.write(process_bar)
        sys.stdout.flush()
        if self.i >= self.max_steps:
            self.close()

    def close(self):
        print('')
        print(self.infoDone)
        self.i = 0


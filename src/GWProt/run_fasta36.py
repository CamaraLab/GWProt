
"""
This is for running programs in the fasta36 package
and parsing the results to get the aligned residues in each sequence

guide:
https://fasta.bioch.virginia.edu/fasta_www2/fasta_guide.pdf

Citation:
T. F. Smith and M. S. Waterman, (1981) J. Mol. Biol. 147:195-197;
W.R. Pearson (1991) Genomics 11:635-650

copied 4/1/2024
"""

import numpy as np
import os
import csv
import Bio
import re
import subprocess

def cigar_parser(cigar_string, ref_start = 0, read_start = 0, allow_mismatch = True): #written myself
    ref_indices = []
    ref_current_index = ref_start
    read_indices = []
    read_current_index = read_start
    cigar_pairs = re.findall( r'([0-9]+[MIDNSXHP=])', cigar_string)
    
    for p in cigar_pairs:
        c = p[-1]
        n = int(p[:-1])
        if c == 'M' or c =='=':
            for i in range(n):
                ref_indices.append(ref_current_index)
                read_indices.append(read_current_index)
                ref_current_index +=1
                read_current_index+=1
        elif c == 'I':
            read_current_index += n  #these two might be backwards
        elif c =='D' or c == 'N':
            ref_current_index +=n

        elif c =='X' or c == 'S':
            for i in range(n):
                if allow_mismatch:
                    ref_indices.append(ref_current_index)
                    read_indices.append(read_current_index)
                ref_current_index +=1
                read_current_index+=1
    return ref_indices, read_indices


def alncode_parser(code, ref_start=0, read_start=0, allow_mismatch = True):
    # code of the form =23+9=13-2=10-1=3+1=5
    code_pairs = re.findall( r'([-\+=][0-9]+)', code)

    ref_indices = []
    ref_current_index = ref_start
    read_indices = []
    read_current_index = read_start

    for p in code_pairs:
        type = p[0]
        n = int(p[1:])
        if type =='=':
            for i in range(n):
                ref_indices.append(ref_current_index)
                read_indices.append(read_current_index)
                ref_current_index +=1
                read_current_index+=1
        elif type == 'x':
            for i in range(n):
                if allow_mismatch:
                    ref_indices.append(ref_current_index)
                    read_indices.append(read_current_index)
                ref_current_index +=1
                read_current_index+=1
        elif type == '+':
            ref_current_index += n
        elif type == '-':
            read_current_index +=n
    return ref_indices, read_indices
    
        
def parse_ssearch_output(input, code_type, allow_mismatch = True):

    pattern = r'aln_code\n.+\.pdb\s+\(\s?\d+\)(?:[\s]+[\d\.]+){7}\s+(\d+)\s+\d+\s+(\d+)'
    # needs to be double checked
    matches = re.findall(pattern, input)
    #print(matches)#debugging
    a0 = int(matches[0][0])
    a1 = int(matches[0][1])
    if code_type == 'CIGAR':#9C, 9D
        code_pattern = r'aln_code\n.+\.pdb\s+\(\s?\d+\)(?:[\s]+[\d\.]+){16}\s+([MIDNSXHP=\d]+)'
        aln_code = re.findall(code_pattern, input)[0]
        #print(aln_code)#debugging
        return cigar_parser(cigar_string = aln_code, ref_start= a0 - 1, read_start=a1 - 1, allow_mismatch= allow_mismatch)
        
    if code_type == '9c':
        code_pattern = r'aln_code\n.+\.pdb\s+\(\s?\d+\)(?:[\s]+[\d\.]+){16}\s+([\d\+\=-]+)'
        alncode = re.findall(code_pattern, input)[0]
        #print(alncode) #debugging
        return alncode_parser(code = alncode, ref_start= a0 - 1, read_start=a1 - 1, allow_mismatch= allow_mismatch)




def run_ssearch_cigar(fasta1, #filepath
                      fasta2, #filepath
                      allow_mismatch = True):
    
    #first check that fasta36 is installed and in the path
    a = subprocess.run(['ssearch36'] , stdout = subprocess.PIPE, text=True)
    if "ssearch36 [-options] query_file library_file" not in a.stdout:
        raise FileNotFoundError("ssearch36 was not found in PATH. fasta36 package must be installed and in the PATH.")


    cigar_command =  '-s BP62 -p -T 1 -b 1 -f 0 -g 0 -z -1 -m 9C'
    a = subprocess.run(['ssearch36'] + cigar_command.split(' ') + [fasta1, fasta2], text=True,  stdout = subprocess.PIPE)
    cigar_result = a.stdout

    return parse_ssearch_output(input = cigar_result, code_type = 'CIGAR', allow_mismatch = allow_mismatch)



def parse_ssearch_outputv2(input,  allow_mismatch = True):
    #for parsing the output from the new run_ssearch code, assumes CIGAR format
    pattern = 'aln_code\\n\s+\(\s*\d+\)(?:[\s\\t]+[\d\.]+){7}\s+(\d+)\s+\d+\s+(\d+)'
    # needs to be double checked
    matches = re.findall(pattern, input)
    #print(matches)#debugging
    a0 = int(matches[0][0])
    a1 = int(matches[0][1])
    code_pattern = 'aln_code\\n\s+\(\s*\d+\)(?:[\s\\t]+[\d\.]+){16}[\s\\t]*((?:\d+[MIDNSXHP=])+)\\n'

    aln_code = re.findall(code_pattern, input)[0]
    #print(aln_code)#debugging
    return cigar_parser(cigar_string = aln_code, ref_start= a0 - 1, read_start=a1 - 1, allow_mismatch= allow_mismatch)
    




def run_ssearch_cigar_Ram(fasta1, #string
                      fasta2, #string
                      allow_mismatch = True):

    #first check that fasta36 is installed and in the path
    a = subprocess.run(['ssearch36'] , stdout = subprocess.PIPE, text = True)
    if "ssearch36 [-options] query_file library_file" not in a.stdout:
        raise FileNotFoundError("ssearch36 was not found in PATH. fasta36 package must be installed and in the PATH.")


    command2 = f'''''/bin/bash -c "ssearch36 -s BP62 -p -T 1 -b 1 -f 0 -g 0 -z -1 -m 9C <(echo '{ fasta1 }') <(echo '{fasta2}')" '''
    
    result = subprocess.run(command2, shell=True,stdout = subprocess.PIPE, text=True)
    
    cigar_result = result.stdout

    return parse_ssearch_outputv2(input = cigar_result, allow_mismatch = allow_mismatch)
    
#WARNING - when mismatch == False there sometimes are aligned residues which are different


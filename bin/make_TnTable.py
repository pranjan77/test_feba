#!/usr/bin/env python

import sys
import readSeq
import argparse
import subprocess


class ReadParser():

def get_usr_args():
    desc = "Create a barcode to genome position table"
    parser = argparse.ArgumentParrser(usage="%(prog)s outbase ref fastq1 fastq2 .... fastqN [options]", description=desc)
    
    parser.add_argument("outbase", nargs=1, type=str, help="the basename for output files")
    parser.add_argument("ref", nargs=1, type=str, help="The reference to use to create table")
    parser.add_argument("fastq", nargs='+', type=str, help="TnSeq fastq(s) to create table from")
    parser.add_argument("-m", "--min_count", type=int, default=2, help="minimum number mappings to a location before considering position")
    parser.add_argument("-f", "--max_frac", type=float, default=0.25, help="maximum fraction of second best mapping")
    io_group = parser.add_argument_group(title="Input/Output Arguments")
    io_group.add_argument("-o", "--outdir", default=".", help="the directory to store output [cwd]")
    
        

#end get_usr_args

def index_ref(ref, index_base):
    index_cmd = [ "bowtie2-build", ref, index_base ] 
    index_proc = subprocess.Popen(index_cmd)
    return index_proc.wait()
#end index_ref

def get_reader(fastq):
    file_type = subprocess.check_output(["file", fastq])
#




def map_reads(fastq, index):
    bowtie2_cmd = [ 'bowtie2', '-x', index, '-1', '-', 
    map_proc = subprocess.Popen(

def main():
    
    args = get_usr_args()
    index_base = "%s/%s" % (args.outdir,args.outbase)
    
    index_retcode = index_ref (args.ref,index_base)
    
    if index_retcode != 0:
        sys.stderr.write("Indexing of reference %s returned %d exit status\n" % (args.ref, index_retcode))
        sys.exit(-1)

    
    
    
    

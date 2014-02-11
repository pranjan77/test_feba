#!/usr/bin/env python

import sys
import readSeq
import argparse
import subprocess
import multiprocessing


class ReadParser():
    def __init__(self,read_model):
        # TODO: parse read model from a fasta file and store information in object
        pass
    #end __init__

    def parse_read(sef,hdr,seq,qseq):
        # TODO: parse read sequence based on read model return new 'read'
        pass
    #end parse_read
#end ReadParser

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


def process_reads(fastq, read_parser, mapper_stdin):
    # process reads with readSeq, and print to mapper stdin
    for r1,r2 in readSeq(fastq):
        (bc_hdr,bc_seq,bc_qseq) = read_parser.parse_read(r1.hdr,r1.seq,r1.qseq)
        mapper_stdin.write("@%s\n%s\n+%s\n%s\n" % (bc_hdr,bc_seq,bc_qseq))
# end process_reads

def get_barcode_location(sam_line):
    # TODO: do some stuff to parse the SAM line
    pass
#end get_barcode_and_position

def process_mapping_results(mapper_stdout, counter):
    for line in mapper_stdout:
        if line[0] == '@':
            # TODO: do I need to do anything here to keep track of reference sequences?
            pass
        else:
            barcode, contig, strand, pos = get_barcode_location(line)
            if barcode not in counter:
                counter[barcode] = dict()
            if contig not in counter[barcode]
                counter[barcode][contig] = { '+': collecitions.Counter(), '-': collecitions.Counter() }
            
            # TODO: do some stuff to count the locations
#end process_mapping_results



def map_reads(index):
    bowtie2_cmd = [ 'bowtie2', '-x', index, '-1', '-' ]
    map_proc = subprocess.Popen(bowtie2_cmd, stdin=PIPE, stdout=PIPE)
    return map_proc
#end map_reads

def main():
    
    args = get_usr_args()
    index_base = "%s/%s" % (args.outdir,args.outbase)
    
    index_retcode = index_ref (args.ref,index_base)
    
    if index_retcode != 0:
        sys.stderr.write("Indexing of reference %s returned %d exit status\n" % (args.ref, index_retcode))
        sys.exit(-1)
    
    map_proc = map_reads(index_base)
    
    barcode_counter = None # TODO: still need to figure out what this data structure should be 
    
    read_in_proc = multiprocessing.Process(target=process_reads, args=(args.fastq,read_parser,map_proc.stdin))
    map_out_proc = multiprocessing.Process(target=process_mapping_results, args=(map_proc.stdout, barcode_counter))

    read_in_proc.start()
    map_out_proc.start()

    map_proc.wait()
    map_out_proc.join()
     
    # TODO: process counts 

    
    
    
    

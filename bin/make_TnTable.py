#!/usr/bin/env python

import sys
import readSeq
import argparse
import subprocess
import multiprocessing

NEG_STRAND = 0
POS_STRAND = 1


class ReadParser():
    def __init__(self,read_model):
        # parse read model from a fasta file and store information in object
        model = None
        with open(read_model, 'r') as read_model_file:
            read_model_file.readline()
            model = read_model_file.readline()[:-1]
        
        N_match = re.search('(N+)',model)
        n_match = re.search('(n+)',model)
        if len(N_match.groups()) != 1 or len(n_match.groups()) != 1:
            self.logger.error('Invalid Model %s ' % model)
            sys.exit(2)
        # not sure if I need to store these, but will for now just in case
        self.preseq_start = n_match.end()
        self.barcode_start = N_match.start()
        self.postseq_start = N_match.end()
        self.genome_start = len(model)

        # use these to search each read to parse out the barcode and genomic sequence
        self.preseq = model[n_match.end():N_match.start()]
        self.postseq = model[N_match.end():]
        
        # use this to count the start of key sequences in reads
        self.barcode_start_counts = collections.Counter()
        self.genomic_start_counts = collections.Counter()
    #end __init__

    def parse_read(self,hdr,seq,qseq):
        preseq_start = string.find(seq, self.preseq)
        postseq_start = string.find(seq, self.postseq, preseq_start)

        # get barcode sequence
        barcode_start = preseq_start+len(self.preseq)
        barcode_seq = seq[barcode_start:postseq_start]

        # get genomic sequence from the read and quality sequence 
        genomic_start = postseq_start+len(self.postseq)
        genomic_seq = seq[genomic_start:]
        quality_seq = qseq[genomic_start:]

        # create a new header that contains the original header with the barcode sequence
        new_header = "%s:%s" % (hdr,barcode_seq)
    
        # tally starts of these features for QC purposes
        self.barcode_start_counts[barcode_start] += 1
        self.genomic_start_counts[genomic_start] += 1

        return new_header, genomic_seq, quality_seq 
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
        (genomic_hdr,genomic_seq,genomic_qseq) = read_parser.parse_read(r1.hdr,r1.seq,r1.qseq)
        mapper_stdin.write("@%s\n%s\n+%s\n%s\n" % (genomic_hdr,genomic_seq,genomic_qseq))
# end process_reads

cigar_re = re.compile('([MIDNSHPX=])')
cigar_len_chars = {'M','D','=','X','N'} # Sum of the lengths of these operations gets 
def get_barcode_location(sam_line):
    # do some stuff to parse the SAM line
    sam_ar = sam_line.split('\t')
    hdr = sam_ar[0]
    ref = sam_ar[2]
    pos = int(sam_ar[3]) if sam_ar[3] != '*' else None
    if not pos:
        return None, None, None, None
    strand = (int(sam_ar[1]) >> 4) % 2
    
    # if on negative strand, need to compute mapping length to get 5prime start
    if strand == 1:
        cigar = sam_ar[5]
        cigar_ar = cigar_re.split(cigar)
        for i in xrange(1,len(cigar_ar),2):
            if cigar_ar[i] in cigar_len_chars:
                pos += int(cigar_ar[i-1])
        strand = NEG_STRAND
    else:
        strand = POS_STRAND
    
    barcode = hdr.split(":")[-1]
    
    return barcode, ref, strand, pos
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
                counter[barcode][contig] = { NEG_STRAND: collecitions.Counter(), POS_STRAND: collecitions.Counter() }
            counter[barcode][contig][strand][pos] += 1
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
    
    # a dictionary that gets populated with barcodes mapped to location counts 
    # indexed like so:
    #   barcode_counter[barcode][contig][strand][position]
    barcode_counter = dict() 
    
    read_in_proc = multiprocessing.Process(target=process_reads, args=(args.fastq,read_parser,map_proc.stdin))
    map_out_proc = multiprocessing.Process(target=process_mapping_results, args=(map_proc.stdout, barcode_counter))

    read_in_proc.start()
    map_out_proc.start()

    map_proc.wait()
    map_out_proc.join()
     
    # TODO: process counts 

    
    
    
    

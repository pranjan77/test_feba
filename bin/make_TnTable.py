#!/usr/bin/env python

import sys
import re
import os
import argparse
import subprocess
import threading
import collections

import numpy as np

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
plt.ioff()

import readSeq
import feba


NEG_STRAND = 0
POS_STRAND = 1
strand_to_symbol = ('-','+')

### BEGIN ARGUMENT PARSING

def ref_base(fasta):
    start_idx = fasta.rfind("/")+1
    if fasta.endswith(".fasta") or fasta.endswith(".FASTA"):
        return fasta[start_idx:-6]
    elif fasta.endswith(".fa"):
        return fasta[start_idx:-3]
    elif fasta.endswith(".fna"):
        return fasta[start_idx:-4]
    else:
        return fasta[start_idx:]
#ref_base

# strips off any compression extension and then any FastQ extension of a filename
def get_fastq_basename(fastq):
    ret=""
    start = fastq.rfind("/")+1
    if fastq.endswith(".gz"):
        ret = fastq[start:-3]
    elif fastq.endswith(".bz2"):
        ret = fastq[start:-4]
    else:
        ret = fastq[start:]
    if ret.endswith(".fastq"):
        ret = ret[0:-6]
    elif ret.endswith(".fq"):
        ret = ret[0:-3]
    return ret
#end get_fastq_basename

'''
Function to ensure the read model is valid input
'''
def __read_model__ (string):
    model = None
    if os.path.isfile(string):
        with open(read_model, 'r') as read_model_file:
            read_model_file.readline()
            model = read_model_file.readline()[:-1]
    else:
        model = string
    if len(re.sub('[nNATCG]', "", model)) > 0:
        raise argparse.ArgumentTypeError("%s is not a valid model.\nMake sure it consists of only nNATCG")
    return model
#end __read_model__
        
__DEFAULT_MODEL__ = "nnnnnnCGCCCTGCAGGGATGTCCACGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTCGACGGCCGGCCAGACCGGGGACTTATCAGCCAACCTGT"

def get_usr_args():
    desc = "Create a barcode to genome position table"
    parser = argparse.ArgumentParser(usage="%(prog)s outbase ref fastq1 fastq2 .... fastqN [options]", description=desc)
    
    parser.add_argument("ref", type=str, help="The reference to use to create table")
    parser.add_argument("fastq", nargs='+', type=str, help="TnSeq fastq(s) to create table from")

    parser.add_argument("-c", "--min_count", type=int, default=2, help="minimum number mappings to a location before considering position")
    parser.add_argument("-f", "--max_frac", type=float, default=0.25, help="maximum fraction of second best mapping")
    parser.add_argument("-m", "--model", type=__read_model__, default=__DEFAULT_MODEL__, help="fasta file contaiing the TnSeq read model")

    parser.add_argument("--debug", default=False, action="store_true", help="run in debug mode. Save SAM file and trimmed reads")

    io_group = parser.add_argument_group(title="Input/Output Arguments")
    io_group.add_argument("-o", "--outdir", default=".", help="the directory to store output [cwd]")
    io_group.add_argument("-b", "--outbase", type=str, help="basename for output files [basename(fastq)]")


    # build the output basename if it hasn't been given
    args = parser.parse_args()
    if not args.outbase:
        args.outbase = "%s.v.%s" % ("-".join([ get_fastq_basename(fq) for fq in args.fastq ]), ref_base(args.ref))

    return args
#end get_usr_args

### END ARGUMENT PARSING

def read_to_string(header,seq,qual):
    return "@%s\n%s\n+\n%s\n" % (header, seq, qual)
#end read_to_string

'''
    This command is a key step. It processes reads by parsing out
    genomic DNA from TnSeq reads, and writes them to the STDIN
    of the aligner.
'''
def process_reads(fastq, read_parser, mapper_stdin, clean_fastq_out=None, unpaired=False):
    PROG = "process_reads"
    # process reads with readSeq, and print to mapper stdin
    for fq in fastq:
        sys.stderr.write("[%s] Parsing %s with readSeq\n" % (PROG,fq))
        for header, seq, qual in read_it(fq, unpaired):
            (genomic_hdr,genomic_seq,genomic_qseq) = read_parser.parse_read(header, seq, qual)
            if not genomic_hdr:
                sys.stderr.write("[%s] Unable to parse read:\n%s" % (PROG,read_to_string(header,seq,qual)))
            read_out = read_to_string(genomic_hdr,genomic_seq,genomic_qseq)
            mapper_stdin.write(read_out)
            if clean_fastq_out:
                clean_fastq_out.write(read_out)
    sys.stderr.write("[%s] Closing stdin to mapper\n"  % (PROG))
    mapper_stdin.close()
    sys.stderr.write("[%s] Done parsing reads\n" % (PROG))
# end process_reads

'''
Helper to process_reads for reading the fastq as paired/unpaired
'''
def read_it(fastq,unpaired):
    if unpaired:
        for r1 in readSeq.readSeq(fastq,paired=False):
            yield r1.header, r1.seq, r1.qual
    else:
        for r1,r2 in readSeq.readSeq(fastq,paired=True):
            yield r1.header, r1.seq, r1.qual
#end read_it

'''
parses a SAM line and returns barcode and the location it mapped to (reference, strand, pos)
'''
cigar_re = re.compile('([MIDNSHPX=])')
cigar_len_chars = set(['M','D','=','X','N']) # Sum of the lengths of these operations gets 
def get_barcode_location(sam_line):
    # do some stuff to parse the SAM line
    sam_ar = sam_line.split('\t')
    hdr = sam_ar[0]
    ref = sam_ar[2]
    pos = int(sam_ar[3]) if sam_ar[3] != '*' else None
    barcode = hdr.split(":")[-1]
    if pos:
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
    else:
        ref = None
        strand = None
        barcode = hdr 
    return barcode, ref, strand, pos
#end get_barcode_and_position


'''
    This command is a key step. It processes results of the alignment,
    by counting barcode locations. This will run in parallel with the 
    process that parses TnSeq reads. 
'''
def process_mapping_results(mapper_stdout, counter, bam_out=None):
    PROG = "process_mapping_results"
    mapped = 0
    unmapped = 0
    for line in mapper_stdout:
        if bam_out:
            bam_out.write(line)
        if line[0] == '@':
            # TODO: do I need to do anything here to keep track of reference sequences?
            pass
        else:
            barcode, contig, strand, pos = get_barcode_location(line)
            if not contig:
                sys.stderr.write("[%s] %s did not align\n" % (PROG,barcode))
                unmapped += 1
                continue
            mapped += 1
            if barcode not in counter:
                counter[barcode] = dict()
            if contig not in counter[barcode]:
                counter[barcode][contig] = { NEG_STRAND: collections.Counter(), POS_STRAND: collections.Counter() }
            counter[barcode][contig][strand][pos] += 1
    total = mapped+unmapped
    sys.stderr.write("[%s] Processed %d total reads. %d mapped and %d unmapped\n" % (PROG,total,mapped,unmapped))
#end process_mapping_results

### BEGIN EXTERNAL COMMANDS

'''
Index reference with bowtie2
'''
def index_ref(ref, index_base):
    index_cmd = [ "bowtie2-build", ref, index_base ] 
    ret_code = 1
    with open("%s.bowtie2-build.log" % index_base, "w") as index_log:
        index_proc = subprocess.Popen(index_cmd, stdout=index_log, stderr=subprocess.STDOUT)
        ret_code = index_proc.wait()
    return ret_code
#end index_ref

'''
create a subprocess for writing reads to and reading SAM output from
'''
def map_reads(index, log_file_path):
    PROG = "map_reads"
    bowtie2_cmd = [ 'bowtie2-align', '--reorder', '-x', index, '-U', '-' ]
    sys.stderr.write("[%s] Running bowtie2: %s\n" % (PROG," ".join(bowtie2_cmd)))
    map_proc = subprocess.Popen(bowtie2_cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=open(log_file_path,"w"))
    return map_proc
#end map_reads

'''
Create a subprocess for converting SAM text to BAM
'''
def bam_writer(bam_file_path):
    bam_write_cmd = [ 'samtools', 'view', '-bhS', '-' ]
    bam_write_proc = subprocess.Popen(bam_write_cmd, stdin=subprocess.PIPE, stdout=open(bam_file_path,"w"))
    return bam_write_proc
#end bam_write_proc

'''
create a file handle for writing SAM text
'''
def sam_writer(sam_file_path):
    return open(sam_file_path,"w")
#end bam_write_proc

'''
Create a subprocess for writing fastq to a gzipper
'''
def gzip_fastq_writer(fastq_file_path):
    fastq_write_cmd = [ 'gzip', '-c' ]
    fastq_write_proc = subprocess.Popen(fastq_write_cmd, stdin=subprocess.PIPE, stdout=open(fastq_file_path,"w"))
    return fastq_write_proc
#end gzip_fastq_writer

### END EXTERNAL COMMANDS

def main():
    PROG = os.path.basename(__file__)
    args = get_usr_args()
    
    sys.stderr.write("[%s] Using %s\n" % (PROG,sys.executable))
    sys.stderr.write("[%s] %s\n" % (PROG," ".join(sys.argv)))
    index_base = "%s/%s.ref" % (args.outdir,ref_base(args.ref))
    out_prefix = "%s/%s" % (args.outdir,args.outbase)
 
    # parser for pulling genomic DNA out of TnSeq read
    read_parser = feba.TnSeqReadParser(args.model)

    sys.stderr.write("[%s] Indexing %s with bowtie2. Logging to %s.bowtie2-build.log\n" % (PROG, args.ref, index_base))
    index_retcode = index_ref (args.ref,index_base)
    if index_retcode != 0:
        sys.stderr.write("[%s] Indexing of reference %s returned %d exit status\n" % (PROG, args.ref, index_retcode))
        sys.exit(1)
    
    
    # process for mapping reads. Map reads by writing to process's STDIN
    bowtie2_log = "%s.bowtie2.log" % out_prefix
    sys.stderr.write("[%s] Running alignment to %s with bowtie2\n" % (PROG, args.ref))
    sys.stderr.write("[%s] Index stored in %s. Logging to %s\n" % (PROG, index_base, bowtie2_log))
    map_proc = map_reads(index_base, bowtie2_log)

    # a dictionary that gets populated with barcodes mapped to location counts 
    # indexed like so:
    #   barcode_counter[barcode][contig][strand][position]
    barcode_counter = dict() 

    
    # check if we're in debug mode. If so, set up necessary processes for writing 
    # the outputs we need to save.
    process_reads_args = None
    process_mapping_results_args = None
    cleaned_reads_write_proc = None
    aligned_reads_write_proc = None 
    
    mapper_stdin = map_proc.stdin
    mapper_stdout = map_proc.stdout
    
    # if debugging, save cleaned reads and alignment output in a BAM file
    if args.debug:
        # Save cleaned reads
        cleaned_reads_fastq = "%s.cleaned_reads.fastq.gz" % out_prefix
        sys.stderr.write("[%s] Saving cleaned reads to %s\n" % (PROG,cleaned_reads_fastq))
        cleaned_reads_write_proc = gzip_fastq_writer(cleaned_reads_fastq)
        process_reads_args = (args.fastq,read_parser,mapper_stdin,cleaned_reads_write_proc.stdin)
        
        # Save aligned clean reads
        aligned_reads_bam = "%s.sam" % out_prefix 
        sys.stderr.write("[%s] Saving aligned clean reads to %s\n" % (PROG,aligned_reads_bam))
        aligned_reads_write_proc = sam_writer(aligned_reads_bam)
        process_mapping_results_args = (mapper_stdout, barcode_counter, aligned_reads_write_proc)
    else:
        process_reads_args = (args.fastq,read_parser,mapper_stdin)
        process_mapping_results_args = (mapper_stdout, barcode_counter)

    # create a process to read fastq and parse reads for piping to bowtie2    
    read_in_proc = threading.Thread(target=process_reads, args=process_reads_args)
    # create a process for reading SAM output from bowtie2 and counting barcodes
    map_out_proc = threading.Thread(target=process_mapping_results, args=process_mapping_results_args)

    # start processes
    sys.stderr.write("[%s] Starting cleaning, mapping, and counting\n" % (PROG))
    read_in_proc.start()
    map_out_proc.start()

    #map_proc.wait()
    map_out_proc.join()
    
    # TODO: process counts 
    sys.stderr.write("[%s] All barcodes counted\n" % (PROG))
    
    
    sys.stderr.write("[%s] Making histogram of counts\n" % (PROG))

        
    all_counts = list()
    filtered_counts = list()

    for barcode in barcode_counter:
        for contig in barcode_counter[barcode]:
            for strand in barcode_counter[barcode][contig]:
                for position in barcode_counter[barcode][contig][strand]:
                    count = barcode_counter[barcode][contig][strand][position]
                    s = strand_to_symbol[strand]
                    sys.stdout.write("%s\t%s\t%s\t%d\t%d\n" % (barcode, contig, s, position, count))
                    all_counts.append(count)
                    if count >= args.min_count:
                        filtered_counts.append(count)
    pdf = PdfPages("%s.plots.pdf" % out_prefix)

    # Add a plot for Barcode counts
    all_counts = np.array(all_counts,dtype=np.int)
    count_hist_fig, count_hist_plot = plt.subplots(nrows=1,ncols=1)
    count_hist_plot.hist(all_counts)
    count_hist_plot.set_title("Raw Counts")
    count_hist_plot.set_xlabel("Count")
    count_hist_plot.set_ylabel("Number of Positions")
    pdf.savefig(count_hist_fig)
    
    # Add a plot for Filtered Barcode counts
    filtered_counts = np.fromiter(filtered_counts,dtype=np.int)
    count_hist_fig, count_hist_plot = plt.subplots(nrows=1,ncols=1)
    count_hist_plot.hist(filtered_counts)
    count_hist_plot.set_title("Filtered Counts (min = %d)" % args.min_count)
    count_hist_plot.set_xlabel("Count")
    count_hist_plot.set_ylabel("Number of Positions")
    pdf.savefig(count_hist_fig)
    
    # Add a barplot of the starting positions for barcodes
    starts_fig, starts_barplot = plt.subplots(nrows=1,ncols=1)
    barcode_starts = np.array(read_parser.barcode_start_counts.keys())
    barcode_starts_counts = np.fromiter((read_parser.barcode_start_counts[x] for x in barcode_starts), dtype=np.int)
    starts_barplot.bar(barcode_starts,barcode_starts_counts, color='r')
    starts_barplot.set_title("Starting positions for barcodes")
    starts_barplot.set_xlabel("Read Position")
    starts_barplot.set_ylabel("Number of Reads")
    pdf.savefig(starts_fig)


    starts_fig, starts_barplot = plt.subplots(nrows=1,ncols=1)
    genomic_starts = np.array(read_parser.genomic_start_counts.keys())
    genomic_starts_counts = np.fromiter((read_parser.genomic_start_counts[x] for x in genomic_starts), dtype=np.int)
    starts_barplot.bar(genomic_starts,genomic_starts_counts, color='b')
    starts_barplot.set_title("Starting positions for barcodes")
    starts_barplot.set_xlabel("Read Position")
    starts_barplot.set_ylabel("Number of Reads")
    pdf.savefig(starts_fig)
    
    pdf.close()

    # Add a barplot of the starting positions for genomic sequence
    
    #count_hist, count_edges = np.histogram(all_counts, density=False, bins=50)
    #diff = np.diff(count_edges)
    #midpoints = np.fromiter((diff[i]+count_edges[i] for i in xrange(len(diff))),np.float)
    #count_hist_plot.plot
    #count_hist_plot.plot(midpoints,count_hist,color="red",linewidth=1.0,linestyle="-",label=")
    
    
    
    
if __name__ == "__main__":
    main()    

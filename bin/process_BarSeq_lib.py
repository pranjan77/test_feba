#!/usr/bin/env python

import sys
import argparse
import multiprocessing
import os
import collections
import copy
import Queue
import string
import datetime

##updating the library look up path
sys.path.append(os.path.join(os.path.dirname(__file__),'../lib'));

import BarSeq
import readSeq

STOP = "NO MORE READS. THIS SHOULD PROBABLY BE SOMETHING ELSE, BUT I'M MAKING IT A STRING AS A PLACE HOLDER"
QUEUE_SIZE = 100000

# run this to read the input fastqs 
def read_reads(read_files, queue, buffer_size=1000, verbose=False):
    proc_name = multiprocessing.current_process().name
    count = 0
    start = datetime.datetime.now()
    buf = [None]*buffer_size
    added = False
    for read_file in read_files:
        for read1, read2 in readSeq.readSeq(read_file, paired=True):
            buf[count] = read1
            count += 1
            if count == buffer_size:
                added = False
                while not added:
                    if queue.qsize() <= QUEUE_SIZE - count:
                        for item in buf:
                            queue.put(item)
                        #map(queue.put,buf)
                        end = datetime.datetime.now()
                        if verbose:
                            sys.stderr.write("[%s] added %d reads to the input queue. Took %s to complete\n" % (proc_name,count, str(end-start)))
                        added = True
                count = 0   
                start = datetime.datetime.now()
    added = False
    while not added:
        if queue.qsize() < QUEUE_SIZE - count:
            for i in xrange(count):
                queue.put(buf[i])
            end = datetime.datetime.now()
            if verbose:
                sys.stderr.write("[%s] added remaining %d reads to the input queue. Took %s to complete\n" % (proc_name,count, str(end-start)))
            added = True
        if verbose:
            sys.stderr.write("[%s] sending STOP signal\n" % (proc_name))
        queue.put(id(STOP))
    if verbose:
        sys.stderr.write("[%s] done reading reads\n" % (proc_name))
    return 0
#end fill_input_read_queue

# many of these will run in parallel to parse barcodes from each read
def process_reads(parser, in_queue, out_queue, buffer_size=1000, verbose=False):
    proc_name = multiprocessing.current_process().name
    barcode_counts = get_counter(parser.get_multiplex_barcodes())
    reads = [None]*buffer_size
    barcodes = [None]*buffer_size
    more_reads = True
    bc_count = 0
    start = datetime.datetime.now()
    mplex_bc, barseq_bc = None, None
    retrieved = False
    while more_reads:
        # fill the buffer if we have reads to read yet.
        retrieved = False
        while not retrieved:
            if not in_queue.empty():
                i = 0
                for i in xrange(buffer_size):
                    reads[i] = in_queue.get()
                    if reads[i] == id(STOP):
                        if verbose:
                            sys.stderr.write("[%s] found STOP signal\n" % (proc_name))
                        in_queue.put(id(STOP))
                        i -= 1
                        break
                if verbose:
                    sys.stderr.write("[%s] retrieved %d reads from the input queue\n" % (proc_name,i+1))
                retrieved = True
        # process each read in the buffer we just filled
        for read in reads:
            if read == id(STOP):
                more_reads = False
                mplex_bc, barseq_bc = None, None
            else:
                mplex_bc, barseq_bc = parser.extract_barcodes(read.seq,read.qual)
            if mplex_bc is not None:
                barcodes[bc_count] = (mplex_bc, barseq_bc)
                bc_count += 1
            if bc_count > 0 and (bc_count == buffer_size or read == id(STOP)):
                for i in xrange(bc_count):
                    barcode_counts[barcodes[i][0]][barcodes[i][1]] += 1
                end = datetime.datetime.now()
                if verbose:
                    sys.stderr.write("[%s] counted %d barcodes. Took %s to complete\n" % (proc_name,bc_count, str(end-start)))
                bc_count = 0
                start = datetime.datetime.now()
            if not more_reads:
                break
    out_queue.put(barcode_counts)
    if verbose:
        #sys.stderr.write("[%s] added counts to out queue\n" % (proc_name))
        sys.stderr.write("[%s] added counts to out queue. \n" % (proc_name, out_queue.qsize(), id(out_queue)))

    return 0
#end process_reads
            
def get_counter(mplex_barcodes):
    counts = dict()
    for bc in mplex_barcodes.keys():
         counts[bc] = collections.Counter()
    return counts
#end get_counter

def get_usr_args():
    description = "Script for processing BarSeq library"
    
    parser = argparse.ArgumentParser(description=description ,add_help=True)
    
    
    parser.add_argument("mplex_barcodes", type=str, help="The barcodes used in this experiment for multiplexing")
    parser.add_argument("pool_file", type=str, help="The pool file (TnSeq output) to use")
    parser.add_argument("fastq", nargs='+', type=str, help="One or more fastqs to count barcodes for")
    parser.add_argument("-eN","--expName", type=str, dest="exp_name", default=None, metavar="", help="Name of the experiment")
    parser.add_argument("-p", "--procs", type=int, dest="procs", default=multiprocessing.cpu_count(), help="the number of processes to use")
    parser.add_argument("-b", "--buffer_size", type=int, dest="buffer_size", default=10000, help="the size of the buffer for processing reads")
    parser.add_argument("-v", "--verbose", action="store_true", default=False, help="be verbose")
    
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    
    return parser.parse_args()
#end get_usr_args


def add_counts(main_counter, counts_to_add):
    for mplex_bc in counts_to_add.keys():
        for barseq_bc in counts_to_add[mplex_bc]:
            main_counter[mplex_bc][barseq_bc] += counts_to_add[mplex_bc][barseq_bc] 
#end add_counts

def read_pool_file(pool_file):
    barcodes = dict()
    with open(pool_file,"r") as fh:
        fh.readline()
        for line in fh:
            ar = line[:-1].split()
            if ar[-1] == "unique":
                barcodes[ar[0]] = (ar[1],ar[2],ar[3])
    return barcodes
#end read_pool_file

def read_morgan_pool_file(pool_file):
    barcodes = dict()
    with open(pool_file,"r") as fh:
        fh.readline()
        for line in fh:
            ar = line[:-1].split()
            if ar[4] == "pastEnd":
                barcodes[ar[0]] = (ar[4],"",ar[5])
            else:
                barcodes[ar[0]] = (ar[4],ar[5],ar[6])
    return barcodes
#end read_morgan_pool_file

complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
def rev_comp(barcode):
    return barcode.translate(complements)[::-1]
#end rev_comp

def main():
   
    #sys.stderr.write("Hello World! 1\n") 
    args = get_usr_args()
    proc_name = multiprocessing.current_process().name

    read_parser = BarSeq.BarSeqBarcodeParser(args.mplex_barcodes)
    mplex_barcodes = read_parser.get_multiplex_barcodes()
    
    procs = [None]*max(2,args.procs)
    #sys.stderr.write("Hello World! 2\n") 

    # a Queue for putting input reads into
    in_queue = multiprocessing.Queue(QUEUE_SIZE)
    # a Queue for getting passing counts back to the master process
    out_queue = multiprocessing.Queue()
    
    reader_buffer_size = args.buffer_size 
    process_buffer_size = reader_buffer_size

    procs[0] = multiprocessing.Process(target=read_reads, args=(args.fastq,in_queue,reader_buffer_size,args.verbose)) 
    procs[1] = multiprocessing.Process(target=process_reads, args=(copy.deepcopy(read_parser),in_queue,out_queue,process_buffer_size,args.verbose))
    for i in xrange(2,len(procs)):
        procs[i] = multiprocessing.Process(target=process_reads, args=(copy.deepcopy(read_parser),in_queue,out_queue,process_buffer_size,args.verbose))



    # start and wait for our processes to finish
    #for p in procs:
    #for i in xrange(0,len(procs),-1):
   
    start = datetime.datetime.now()
    sys.stderr.write("[%s] %s\n" % (proc_name,start.strftime("DATE: %02m/%02d/%04Y %02H:%0M:%02S")))
    for i in xrange(len(procs)):
        sys.stderr.write("[%s] Attempting to start process %s\n" % (proc_name,procs[i].name))
        procs[i].start()    
        sys.stderr.write("[%s] %s started\n" % (proc_name,procs[i].name))

    # combine counts
    counts = get_counter(mplex_barcodes)
    recvd = 0
    while recvd < len(procs) - 1:
        add_counts(counts,out_queue.get())
        recvd += 1

    sys.stderr.write("[%s] Read all counts\n" % (proc_name))
    
    sys.stderr.write("[%s] Loading pool file\n" % (proc_name))
    #barseq_barcodes = read_morgan_pool_file(args.pool_file)
    barseq_barcodes = read_pool_file(args.pool_file)
    sys.stderr.write("[%s] Found %d barcodes in pool file\n" % (proc_name,len(barseq_barcodes)))
    #barseq_barcodes = read_pool_file(args.pool_file)
    out = sys.stdout
    sys.stderr.write("[%s] Writing poolcounts table.\n" % (proc_name))

    out.write("barcode\trcbarcode\tscaffold\tstrand\tpos")
    for mplex_bc in mplex_barcodes.keys():
        out.write("\t%s" % mplex_barcodes[mplex_bc].name)
    out.write("\n")
    for barseq_bc in barseq_barcodes.keys():
        tn_seq, tn_strand, tn_pos = barseq_barcodes[barseq_bc] 
        rc_bc = rev_comp(barseq_bc)
        out.write("%s\t%s\t%s\t%s\t%s" % (barseq_bc,rc_bc,tn_seq,tn_strand,tn_pos))
        for mplex_bc in mplex_barcodes.keys():
            #if barseq_bc in counts[mplex_bc]:
            if rc_bc in counts[mplex_bc]:
                sys.stderr.write("Found %s in table\n" % rc_bc)
                out.write("\t%d" % counts[mplex_bc][rc_bc])
            else:
                out.write("\t0")
        out.write("\n")
    
    end = datetime.datetime.now()
    sys.stderr.write("[%s] %s\n" % (proc_name,end.strftime("DATE: %02m/%02d/%04Y %02H:%0M:%02S")))
    
    sys.stderr.write("[%s] Took %s to complete\n" % (proc_name,str(end-start)))

    sys.exit(0)
    # export counts to stdout 
    for mplex_bc in counts.keys():
        print mplex_barcodes[mplex_bc].name
        for barseq_bc in counts[mplex_bc]:
            print "%s\t%s\t%d" % (mplex_bc,barseq_bc,counts[mplex_bc][barseq_bc]) 

    
    
#end main


main()


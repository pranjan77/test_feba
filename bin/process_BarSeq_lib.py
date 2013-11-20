#!/usr/bin/env python

import sys
import argparse
import BarSeq


STOP = "NO MORE READS. THIS SHOULD PROBABLY BE SOMETHING ELSE, BUT I'M MAKING IT A STRING AS A PLACE HOLDER"

# run this ti read the input fastqs 
def fill_input_read_queue(read_files, queue, lock, buffer_size=1000):
    count = 0
    for read_file in read_files:
        buf = [None]*buffer_size
        for read1, read2 in readSeq.readSeq(read_file, paired=True):
            buf[count] = read1
            count += 1
            if count == buffer_size:
                count = 0   
                lock.acquire()
                map(queue.add,buf)
                lock.release()
    lock.acquire()
    queue.add(STOP)
    lock.release()
#end fill_input_read_queue

# run this to pull barcodes from a queue and tally occurences
def read_barcodes(barcode_counts, queue, lock, buffer_size=1000):
    barcodes = [None]*buffer_size
    total_barcodes_processed = 0
    while True:
        lock.acquire()
        for i in xrange(buffer_size):
            barcodes[i] = queue.get()
            if barcodes[i] is STOP:
                break
        lock.release()
        for bc in barcodes:
            if bc is STOP:
               return total_barcodes_processed
            mplex_bc, barseq_bc = bc
            barcode_counts[mplex_bc][barseq_bc] += 1
            total_barcodes_processed += 1        
#end read_barcodes


# many of these will run in parallel to parse barcodes from each read
def process_reads(parser, in_queue, in_lock, out_queue, out_lock, buffer_size=1000):
    reads = [None]*buffer_size
    barcodes = [None]*buffer_size
    more_reads = True
    bc_count = 0
    while more_reads:
        if no_more_reads:
            break
        # fill the buffer if we have reads to read yet.
        in_lock.acquire()
        for i in xrange(buffer_size):
            reads[i] = in_queue.get()
            if reads[i] is STOP:
                in_queue.put(STOP)
                break
        in_lock.release()
        # process each read in the buffer we just filled
        for read in reads:
            if read is STOP:
                out_lock.acquire()
                out_queue.add(STOP)
                out_lock.release()
                more_reads = False
            (mplex_bc, barseq_bc) = parser.get_barcodes(read.seq,read.qual)
            if mplex_bc is not None:
                barcodes[bc_count] = mplex_bc, barseq_bc)
                bc_count += 1
                if bc_count == buffer_size:
                    bc_count = 0
                    out_lock.acquire()
                    map(out_queue.add,barcodes)
                    out_lock.release()
#end process_reads
            

def get_usr_args():
    description = "Script for processing BarSeq library"
    
    parser = argparse.ArgumentParser(description=desc,add_help=True)
    
    
    parser.add_argument("mplex_barcodes", nargs=1, type=BarSeq.BarSeqBarcodeParser, metavar="", help="The barcodes used in this experiment for multiplexing")
    parser.add_argument("pool_file", nargs=1, type=str, metavar="", help="The pool file (TnSeq output) to use")
    parser.add_argument("fastq", nargs='+', type=str, help="One or more fastqs to count barcodes for")
    parser.add_argument("-eN","--expName", type=str, dest="exp_name", default=None, metavar="", help="Name of the experiment")
    
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    
    return parser.parse_args()
#end get_usr_args




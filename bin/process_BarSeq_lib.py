#!/usr/bin/env python

import sys
import argparse
import multiprocessing
import os
import collections

##updating the library look up path
sys.path.append(os.path.join(os.path.dirname(__file__),'../lib'));

import BarSeq
import readSeq

STOP = "NO MORE READS. THIS SHOULD PROBABLY BE SOMETHING ELSE, BUT I'M MAKING IT A STRING AS A PLACE HOLDER"
QUEUE_SIZE = 100000

# run this to read the input fastqs 
def read_reads(read_files, queue, buffer_size=1000):
    proc_name = multiprocessing.current_process().name
    count = 0
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
                        sys.stderr.write("[%s] added %d reads to the input queue\n" % (proc_name,count))
                        added = True
                count = 0   
    added = False
    while not added:
        if queue.qsize() < QUEUE_SIZE - count:
            for i in xrange(count):
                queue.put(buf[i])
            sys.stderr.write("[%s] added remaining %d reads to the input queue\n" % (proc_name,count))
            added = True
        sys.stderr.write("[%s] sending STOP signal\n" % (proc_name))
        queue.put(id(STOP))
    sys.stderr.write("[%s] done reading reads\n" % (proc_name))
#end fill_input_read_queue

# many of these will run in parallel to parse barcodes from each read
def process_reads(parser, in_queue, out_queue, buffer_size=1000):
    proc_name = multiprocessing.current_process().name
    barcode_counts = get_counter(parser.get_multiplex_barcodes())
    reads = [None]*buffer_size
    barcodes = [None]*buffer_size
    more_reads = True
    bc_count = 0
    mplex_bc, barseq_bc = None, None
    retrieved = False
    while more_reads:
        # fill the buffer if we have reads to read yet.
        retrieved = False
        while not retrieved:
            if not in_queue.empty():
                i = 0
                for i in xrange(buffer_size):
                    reads[i] = in_queue.get(timeout=3)
                    if reads[i] == id(STOP):
                        sys.stderr.write("[%s] found STOP signal\n" % (proc_name))
                        in_queue.put(id(STOP))
                        i -= 1
                        break
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
                sys.stderr.write("[%s] counted %d barcodes\n" % (proc_name,bc_count))
                bc_count = 0
            if not more_reads:
                break
    out_queue.put(barcode_counts,timeout=5)
    sys.stderr.write("[%s] added counts to out queue\n" % proc_name)
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
    
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    
    return parser.parse_args()
#end get_usr_args


def add_counts(main_counter, counts_to_add):
    for mplex_bc in counts_to_add.keys():
        for barseq_bc in counts_to_add[mplex_bc]:
            main_counter[mplex_bc][barsqe_bc] += counts_to_add[barseq_bc] 
#end add_counts
    

def main():
    
    args = get_usr_args()
    proc_name = multiprocessing.current_process().name

    read_parser = BarSeq.BarSeqBarcodeParser(args.mplex_barcodes)
    mplex_barcodes = read_parser.get_multiplex_barcodes()
    
    procs = [None]*max(2,args.procs)

    # a Queue for putting input reads into
    in_queue = multiprocessing.Queue(QUEUE_SIZE)
    # a Queue for getting passing counts back to the master process
    out_queue = multiprocessing.Queue()
    
    procs[0] = multiprocessing.Process(target=read_reads, args=(args.fastq,in_queue)) 
    procs[1] = multiprocessing.Process(target=process_reads, args=(read_parser,in_queue,out_queue,1000))
    for i in xrange(2,len(procs)):
        procs[i] = multiprocessing.Process(target=process_reads, args=(read_parser,in_queue,out_queue,1000))


    # start and wait for our processes to finish
    for p in procs:
        p.start()    


    if not procs[-1].is_alive():
        sys.stderr.write("[%s] %s is dead.\n" % (proc_name, procs[-1].name))
    else:
        sys.stderr.write("[%s] %s is alive.\n" % (proc_name, procs[-1].name))
    
    procs[-1].join()
    sys.stderr.write("[%s] Made it past joining the last process\n" % proc_name)
    for p in procs:
        p.join()
        sys.stderr.write("[%s] Process complete!\n" % p.name)

    # combine counts
    counts = get_counter(mplex_barcodes)
    while not out_queue.empty():
        add_counts(counts,out_queue.get())

    # export counts to stdout 
    for mplex_bc in counts.keys():
        print mplex_barcodes[mplex_bc].name
        for barseq_bc in counts[mplex_bc]:
            print "\t%s: %d" % (barseq_bc,counts[barseq_bc]) 
    
#end main


main()











## run this to pull barcodes from a queue and tally occurences
#def read_barcodes(barcode_counts, queue, lock, buffer_size=1000):
#    barcodes = [None]*buffer_size
#    total_barcodes_processed = 0
#    while True:
#        #lock.acquire()
#        sys.stderr.write("getting %d processed barcodes from the output queue\n" % buffer_size)
#        for i in xrange(buffer_size):
#            barcodes[i] = queue.get()
#            if barcodes[i] is STOP:
#                break
#        #lock.release()
#        for bc in barcodes:
#            if bc is STOP:
#               return total_barcodes_processed
#            mplex_bc, barseq_bc = bc
#            barcode_counts[mplex_bc][barseq_bc] += 1
#            total_barcodes_processed += 1        
##end read_barcodes



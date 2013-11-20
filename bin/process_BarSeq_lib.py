#!/usr/bin/env python

import sys
import argparse
import BarSeq

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



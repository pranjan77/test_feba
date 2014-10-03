#!/usr/bin/env python

import sys
import argparse

import feba
import readSeq

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

def get_args():
    desc = "Process TnSeq reads by removing trimming pre and post sequences and barcode. Reads printed to stdout"
    parser = argparse.ArgumentParser(usage="%(prog)s fastq [options]", description=desc)
    parser.add_argument("fastq", type=str, help="the fastq to process")
    parser.add_argument("-m", "--model", default=__DEFAULT_MODEL__, type=str, help="the model with which to process reads")
    parser.add_argument("-u","--unpaired", action='store_true', default=False, help="reads are single-end")
    
    return parser.parse_args()
#end get_args

args = get_args()

parser = feba.TnSeqReadParser(args.model)
reader = readSeq.readSeq(args.fastq,paired=not args.unpaired)

out = sys.stdout 

for r1,r2 in reader:
    new_header, new_seq, new_qual = parser.parse_read(r1.header, r1.seq, r1.qual)
    if new_header:
        out.write("@%s\n%s\n+\n%s\n" % (new_header, new_seq, new_qual))
    else:
        sys.stderr.write("Unable to parse read:\n@%s\n%s\n+\n%s\n" % (r1.header, r1.seq, r1.qual))


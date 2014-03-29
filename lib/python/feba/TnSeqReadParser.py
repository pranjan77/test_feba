#!/usr/bin/env python

import re
import string
import collections

class TnSeqReadParser():
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
        if preseq_start == -1 or postseq_start == -1:
            return None, None, None

        # get barcode sequence
        barcode_start = preseq_start+len(self.preseq)
        barcode_seq = seq[barcode_start:postseq_start]

        # get genomic sequence from the read and quality sequence 
        genomic_start = postseq_start+len(self.postseq)
        genomic_seq = seq[genomic_start:]
        quality_seq = qseq[genomic_start:]

        # create a new header that contains the original header with the barcode sequence
        new_header = None
        header_index_pos = hdr.find(" ")
        if header_index_pos >= 0:
            new_header = "%s:%s" % (hdr[0:header_index_pos],barcode_seq)
        else:
            new_header = "%s:%s" % (hdr,barcode_seq)
            
    
        # tally starts of these features for QC purposes
        self.barcode_start_counts[barcode_start] += 1
        self.genomic_start_counts[genomic_start] += 1

        return new_header, genomic_seq, quality_seq 
    #end parse_read
    
    def __str__(self):
        s = ""
        s += "preseq_start  = %d\n" % self.preseq_start 
        s += "barcode_start = %d\n" % self.barcode_start 
        s += "postseq_start = %d\n" % self.postseq_start 
        s += "genome_start  = %d\n" % self.genome_start 
        s += "preseq  = %s\n" % self.preseq 
        s += "postseq = %s\n" % self.postseq 
        return s

#end ReadParser

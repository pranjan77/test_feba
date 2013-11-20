#!/usr/bin/env python


# the length of BarSeq barcodes
BARSEQ_BARCODE_LENGTH = 20

class MultiplexBarcode: 

    def __init__(self,name,seq,leading_n):
            self.name = name
            self.seq = seq
            self.leading_n = leading_n
    #end __init__

#end MultiplexBarcode

class BarSeqBarcodeParser:
    
    def __init__(self,mplex_bc_file_path, preseq="CAGCGTACG", postseq="AGAGACCTC", n_pre_expected=9):
        self.mplex_barcodes = dict()
        self.mplex_bc_lengths = set()
        self.preseq = preseq
        self.postseq = postseq
        with open(mplex_bc_file_path,"r") as bc_file:
            for line in bc_file:
                if line[0:4] == "name":
                    continue
                (name,bc_seq) = line.split("\t")
                leading_Ns = 0
                if bc_seq[0] == 'N' or bc_seq[0] == 'n':
                    leading_Ns = sentence.count('n') + sentence.count('N')
                    bc_seq = bc_seq[leading_Ns:]
                self.mplex_barcodes[bc_seq] = BarSeqBarcode(name,bc_seq,leading_Ns)
                self.mplex_bc_lengths.add((leading_Ns,len(bc_seq))
    #end __init__
    
    def get_barcodes(self,seq,qual):
        barseq_read_seq = None
        barseq_read_qual = None
        mplex_bc_seq = None
        # find the multiplexing barcode first and remove it to get the barsq read
        for leading_Ns, bc_len in self.mplex_bc_lengths:
            mplex_bc_seq = seq[leading_Ns:leading_Ns+bc_len]
            if mplex_bc_seq in self.mplex_barcodes:
                barseq_read_seq = seq[leading_Ns+bc_len:]
                barseq_read_qual = qual[leading_Ns+bc_len:]
                break

        # no multiplexing barcode, don't move on
        if barseq_read_seq == None:
            return None

        
        preseq_pos = barseq_read_seq.find(self.preseq)
        # make sure preseq is where we expected, and allow for up to 2 insertion/deletions in the reads
        if preseq_pos < 0 or preseq_pos < self.n_pre_expected-2 or preseq_pos > self.n_pre_expected+2:
            raise InvalidBarSeqReadError("%s has invalid preseq (%s) index: %d" % (seq,self.preseq,seq.find(self.preseq)))
        
    
        barseq_bc_index = preseq_pos+len(self.preseq)
        barseq_bc_seq = barseq_read_seq[barseq_bc_index:barseq_bc_index+BARSEQ_BARCODE_LENGTH]
        barseq_bc_qual = barseq_read_qual[barseq_bc_index:barseq_bc_index+BARSEQ_BARCODE_LENGTH]
        if len(barcode_seq) != BARSEQ_BARCODE_LENGTH or !barcode_seq.isupper():
            raise InvalidBarSeqReadError("%s has invalid barcode sequence: %s" %(seq,barcode_seq)) 
        
        # trim the postseq to be expected if the read is not long enough
        read_postseq = self.postseq[:(len(barseq_read_seq)-preseq_pos)-38] if len(barseq_read_seq) - preseq_pos < 38 else self.postseq
        
            
        
        return (mplex_bc_seq, barseq_bc_seq)
    #end get_barseq_barcode

#end BarSeqBarcodeParser


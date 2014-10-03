#!/usr/bin/env python


import sys
import os
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC
import pprint
import string
import re



def check():
    print "Successfully used the module"
    

def dna_Compliment(seq):
    '''given a dna sequence just returns the compliment on the other strand
    '''
    rev_table = string.maketrans('ACGTacgt','TGCAtgca')
    return seq.translate(rev_table)



def dna_RevCompliment(seq):
    '''given a dna sequence just returns the reverser compliment on the same strand
    mainly done to find the sequence of the other strand in the 5' -> 3' direction
    '''
    rev_seq = seq[::-1]
    rev_table = string.maketrans('ACGTacgt','TGCAtgca')
    return rev_seq.translate(rev_table)


def extract_sequence_from_genome(bedFile,genomeFasta):
    
    """
    fasta = pybedtools.BedTool(genomeFasta)
    bedFile
    """
    pass



def get_fasta_sequence_length(fasta,outFile=None):
    SUB='get_fasta_sequence_length'
    base_path = os.path.dirname(fasta) or '.'
    fasta_name = os.path.basename(fasta)
    if outFile is None:
        outFile = base_path + '/' + fasta_name + '.genomesize.txt'
    fh=open(outFile,'w')
    args = ['/global/homes/a/apratap/softwares/seqtk/seqtk','comp',fasta]
    return_code = subprocess.check_call(args,stdout=fh) 
    if return_code == 0:
        fh.close()
        print '[%s]: Created genome size file %s' % (SUB,os.path.basename(outFile))
        return outFile  
        
    
def get_total_bases_in_fasta(fasta): 
    """
        return the total number of bases in fasta
    """
    return sum([ len(record.seq) for record in SeqIO.parse(fasta, "fasta") ])

        
        

    
    
    


def extract_fasta_sequences(fasta,headers,outFile=None):
    """
    given a fasta file extract the sequences by header into a new file
    
    """
    
    if outFile is None:
        outFile = fasta + '_subsequence.fasta'
    
    
    in_handle = open(fasta,"rU")
    out_handle = open(outFile,"w")
    
    headers_list = headers.split(',')
    
    print 'headers %s' % headers_list
    for record in SeqIO.parse(in_handle, "fasta"):
        
        if record.id.strip() in headers_list:
            print 'Found record..%s: extracting.....' % record.id
            SeqIO.write(record,out_handle,"fasta")
        
    
    
    
    
    
        

def get_kmers_dict(sequence,kmer_len=None,step_size=1):
    """
    given a input sequence generate a dictionary of kmers with length = kmer_len
    
    """
    SUB='get_kmers_dict'
    
    if kmer_len is None:
        print '[%s]: kmer_len not passed \n' % (SUB)
        return False
    
    kmers_dict = {}
    sequence = sequence.strip().upper()
    seq_revcomp = dna_RevCompliment(sequence)
    #print 'Seq:%s\nRev:%s \n' % (sequence,seq_revcomp)
    
    #generate dict of k-mers for sequqence
    [ kmers_dict.update( {sequence[x:x+kmer_len]:1} )  for x in xrange(0, (len(sequence)-kmer_len+1),step_size)]
     
    
    #generate dict of k-mers for rev comp sequqence
    [ kmers_dict.update( {seq_revcomp[x:x+kmer_len]:1} ) for x in xrange(0, (len(sequence)-kmer_len+1),step_size)]
    
    return  kmers_dict



def mutate_kmers(sequence,num_errors):
    """
    given a kmer generate a list of kmers 
    with upto #x errors (AGCTN), controlled by num_errors
    """
    
    for error in xrange(num_errors+1,0,-1):
        
        print 'At error %d' % (error)
        
#        for x in xrange(len(sequence)):
#            #print x
#            prefix = sequence[:x]
#            suffix = sequence[x+1:]
#            [ mutate_kmers(prefix+base+suffix,error) for base in 'ACGTN'] 
#        
#        
#        #base case
#        if error == 0:
#            print prefix+base+suffix
        
        #print 'At %d \n \n %s' %(x,temp_dict)
    

def get_kmers_dict_with_one_error(sequence,kmer_len,step_size=1):
    """
    given a input sequence generate a dictionary of kmers with length = kmer_len and with one error allowed
    
    """
    
    SUB = 'get_kmers_dict_with_one_error'
    
    if kmer_len is None:
        print '[%s]: kmer_len not passed \n' % (SUB)
        return False
    
    kmers_dict = {}
    sequence = sequence.strip().upper()
    seq_revcomp = dna_RevCompliment(sequence)
    #print 'Seq:%s\nRev:%s \n' % (sequence,seq_revcomp)
    
    #generate dict of k-mers for sequqence with one error
    [ kmers_dict.update( mutate_kmer_with_one_error( sequence[x:x+kmer_len] ) )  for x in xrange(0, (len(sequence)-kmer_len+1),step_size) ]
        
    #generate dict of k-mers for rev comp sequqence with one error
    [ kmers_dict.update( mutate_kmer_with_one_error( seq_revcomp[x:x+kmer_len] ) ) for x in xrange(0, (len(sequence)-kmer_len+1),step_size)]
    return  kmers_dict




def mutate_kmer_with_one_error(kmer):
    """
    given a sequence(kmer) generate all possible combintations which will allow matching with one error
    """
    
    SUB='gen_kmers_with_one_error'
    
    temp_dict = {}
    for x in xrange(len(kmer)):
        #print x
        prefix = kmer[:x]
        suffix = kmer[x+1:]
        temp_dict.update([ (prefix+base+suffix,1) for base in 'ACGTN'])
    #print temp_dict
    #print (len(temp_dict))
    return temp_dict


def do_kmer_dict_match(query,kmers_dict,kmer_len,report_all_matches=False):
    """
    given a query string and a kmer_dict and a kmer_len
    the function will search each subsequence of kmer_len in the query to kmer dictionary
    and return the starting position of the subsequence which matches the kmer dict
    or None
    """ 
    
    query=query.upper()
    matches =  []
    position = 0
    max_match_upto = ( len(query)-kmer_len+1 )
    while( position <= max_match_upto):
        if kmers_dict.get(query[position:position+kmer_len],None) is not None:
#            print position,query[position:position+kmer_len]
            match_position = position  #previously was adding one to the match position
            if report_all_matches is False:
                 return (match_position)
            else:
                matches.append(match_position)
            # if a kmer of kmer_len matches at position x then there is no use checking next kmer_len bases as they have already been matched
            # the caveat : if two kmer overlaps  then there is a possiblity that second k-mer will be missed
            position = position + kmer_len
        else:
            position += 1
    if matches:
        return matches
    else:
        return False

def plotLen(lengths):
    pass




        
        
    














    
def predict_orfs(sequence,use_ATG_as_start_codon=True,get_best_ORF=True,min_orf_length=300,both_strands=False):
    '''
    given an input DNA sequence returns the best or all ORFs
    in a list
    '''

    sequence = Seq(sequence.replace(' ','').upper().strip())
    min_protein_length = min_orf_length/3
    
    #exit if minium threshold not met
    if( len(sequence) < min_protein_length*3):
        return (0,0,0,0,'NA')
    
    ORFs = []
    
    seq_list = [(+1,sequence)]
    if both_strands:
        seq_list.append( (-1, sequence.reverse_complement() ) )
    
    for strand,sequence in seq_list:
        
        for frame in range(3):
            orf_start_pos = 0
            orf_end_pos = 0
            current_frame_seq = sequence[frame:]  #get the sequence for the current frame
            
            #translate the DNA to protein in a particular frame
            protein=current_frame_seq.translate(table=1)
            current_start_pos = 0
            current_end_pos = 0
            
            #find the position/s where protein terminates due to stop codon denoted by '*'
            for m in re.finditer(r'\*',str(protein)):
                current_end_pos = m.end()    
                
                #if protein is bigger than min threshold
                if current_end_pos - current_start_pos >= min_protein_length:
                    
                    #DEBUG
                    #print 'BEGIN: Frame: %d Start :%d, End: %d, Length %d' % ( (frame+1)*strand,current_start_pos*3,current_end_pos*3,current_end_pos*3-current_start_pos*3)
                    
                    possible_protein = str(protein[current_start_pos:current_end_pos])
                    
                    if use_ATG_as_start_codon:
                        #the protein should start with M (ATG) only
                        location_of_start_codon = possible_protein.find('M')
                        if location_of_start_codon != -1:
                            orf_start_pos = (current_start_pos + location_of_start_codon) * 3 + frame
                            orf_end_pos = ( current_end_pos * 3 ) + frame
                            
                            #check if the length of trimmed protein due to the requirement to start with M is still > min
                            #else skip
                            if orf_end_pos - orf_start_pos < min_orf_length:
                                current_start_pos = current_end_pos
                                continue
                        
                        #if protein doesnt start with M(ATG) update start,end position and skip 
                        else:
                            current_start_pos = current_end_pos
                            continue
                    else:
                        #M not necessary to be in the beginning
                        orf_start_pos = current_start_pos*3 + frame
                        orf_end_pos = current_end_pos*3  + frame
                    
                    
                    if strand == -1:
                    #IMP
                    # if strand is -1 we are on the opp strand wrt to sequence provided to the orf predictor
                    # we need to alter the coordinates with respect to the strand of the sequence provided
                    #    -------------------------------------> (length 100 bp) (1 strand)
                    
                    
                    #   <-----------------y-------x------------- (-1 strand)  
                    
                    #position of y with respect to strand 1 is 100-y(100 is assumed length of strand 1)
                    # and similarly position of x with respect to strand 1 is 100-x
                    #slightly convoluted but interested logic
                        temp_start = len(sequence) - orf_end_pos  
                        temp_end = len(sequence) - orf_start_pos
                        
                        orf_start_pos = temp_start
                        orf_end_pos = temp_end
                    output_frame = ( frame + 1) * strand
                    orf_length=orf_end_pos-orf_start_pos
                    ORFs.append([orf_length,orf_start_pos,orf_end_pos,output_frame])
                    
                    #DEBUG
                    #print '###ORF: Length %d  Start :%d, End: %d Frame: %d ' % (orf_length,orf_start_pos,orf_end_pos,output_frame)
    
                #one iteration over, the end of this becomes the start of another
                current_start_pos = current_end_pos
    
    #sort the list to get the longest ORF on the top
    ORFs.sort(reverse=True)
    if get_best_ORF:
        #return type(ORFs[0])
        return ORFs[0]
    else:
        #all ORFs
        return ORFs
                
                
                
                
            
            
        



if __name__ == "__main__":
    
    
    #a="TGAGAATAGGTGAAGGATGTTTCATCCCCCGAACCTCTAATCATTCGCTTTACCCGATAAAACTGATCAAGCTCCAGCTATCCTGGGGGAAACTTCGGAGGGAACCAGCTACTAGATGGTTCGATTAGTCTTTCGCCCCTATACCCAAGTCTGAAAAGCGATTTGCACGTCAGCACATCTACGAGCCTACGAGGCATTCTTGTGACAATCTCGTGCGGCTGCTGGCCCTCTGGAATGCCT"
    a="GTTAGCGCGTAGTTTGGCACCTTAACTCGACTATCGGTTCATCCCGCATCGCCAGTTCTGCTTACCAAAAATGGCCCACTTGGAGCTCACATTGAATGTGCCGGTTCAATTAAGCAACCGACACGTCTTACCTATTTAAAGTTTGAGAATAGGTGAAGGATGTTTCATCCCCCGAACCTCTAATCATTCGCTTTACCCGATAAAACTGATCAAGCTCCAGCTATCCTGGGGGAAACTTCGGAGGGAACCAGCTACTAGATGGTTCGATTAGTCTTTCGCCCCTATACCCAAGTCTGAAAAGCGATTTGCACGTCAGCACATCTACGAGCCTACGAGGCATTCTTGTGACAATCTCGTGCGGCTGCTGGCCCTCTGGAATGCCTTTGGAAATTCGGGCATGATACTCATTCAAGGTGGTTATCATGGAGTCAGTGCGGAAGACGGTCAGTCCTGTCAGCGACATCCTCCCAGCCTGCCCTGCACCAACGTCCTGCATCCGCTCGCTACGTTGCACTACCCGCCCACCTCTGCCCTCTCCTTCGCCCGCTTTCCCAGCCAGCCTCCACTCTCCCCATCGCCCTCCCTCCGCGGTCCACTCCTAGAACCTTCACTGCCATCGTCGCCGCCGCCCTCTCTGCCCAGCCCGTCGCGCCAGCCGCCCTCTCTGCCCAGCCCGTCGCGCCCGCCGCCCTCCCGCTCGTGCCTCCGGTTCTGTTCCCTCGCACTGCCGTCCAACCGCCGCCGCCTCCTTTTCCGTCCCCTCCTTCGCGCGCTGTCAGCCTCCAGCCTCCCCACTCGCTTAAGCTCCACTCGACCCAAGCCCCACCTCCTGCGCGCTTACCCAACCAGCCTCCGCTTTGACCACCGCGCTCCCTGGCCAGCCCACATCCACTTCCGCCTCCTCCCACGCATCCCCCGTCCAGCCCACCTCAGCTTTTGCCCTTTCGTCCCGTCATGCTCCTCGCCCAGCCAGCCTCAGCCGCCCTCTCCTTCCGCGCGCTCCATGTGCTCCCCGTCCCAGCCTTCCTATTAACTGCAGAAGCATCGTGCGCGCGCCCCTTCCTGTGGCTGGTGTCAAGGCCCCACAGGTGACGGCGTGCCGGCACGACCGGTGCGCCCCGCCTACGCCTCAAGAACCACCGTGCACGCGCCCTTTCCCGCGGCTGGTGTCAAGGCCGCAGGTGGCAGCATGCCGGCATGTGCGGTGCGCCCCGCCTATGCCTCAAGAACCGCCGTGCACGCGCCCTTTCCCGCGGCTGGTGTCAAGGCCCCGCAGGTGGCGGCGTGCCGGCACGTCAGGTGTGCCCGGCCTACGCCTCAAG"
    b=open("/global/projectb/scratch/apratap/1.txt",'r').read().replace('\n','')
    #b=str(Seq(b).reverse_complement())
    print predict_orfs(b)
    
    
    
    

    
    
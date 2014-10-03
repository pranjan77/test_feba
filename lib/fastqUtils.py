#!/usr/bin/env python




import sys
import os
import copy
import random
import re
import matplotlib
#matplotlib.use('Agg') #background plotting
import matplotlib.pyplot as plt
import numpy as np

##updating the library look up path
dir = os.path.dirname(__file__)
sys.path.append(os.path.join(dir,'../lib'));
import seqUtils
import utils


class fastq(object):
    
    #make it faster and memory efficient using slots
    __slots__ = ['count','header','seq','qual','type']
    
    count = 0
    
    def __init__(self,header=None,seq=None,qual=None,type=None):
        self.header=header
        self.seq=seq
        self.qual=qual
        self.type=type
        
    def get_revComp_read(self):
        '''
        given a fastq read object 
        return a new fastq read object with read seq and qual reverse complemented 
        '''
        rev_comp_read = copy.deepcopy(self)
        rev_comp_read.seq     = seqUtils.dna_RevCompliment(self.seq)
        rev_comp_read.qual  = self.qual[::-1]
        return rev_comp_read
    

    def trim_read(self,trim_length):
        '''
        given a fastq read object
        and trim length
        return the new trimmed_read fastq object
        '''
        trimmed_read = copy.deepcopy(self)
        ##creating the trimmed read
        trimmed_read.seq        =   self.seq[:trim_length]
        trimmed_read.qual    =   self.qual[:trim_length]
        return trimmed_read
       
    
    def __str__(self):
        '''
        given a fastq read object 
        print the read object to stdout 
        '''
        return ('%s\n%s\n%s' % (self.header,self.seq,self.qual))
        
        
    def num_of_Ns(self,pattern=None):
        
        if pattern is None:
            pattern = 'N'
        return len(re.findall(pattern,self.seq))
        
        


def get_pairedfastq_iterator(fastq_files):
    '''
    given a fastq file or files return read1 and read2 as an object of class fastq 
    through a generator
    
    if multiple files are given read them one by one
    '''
    
    SUB='get_pairedfastq_iterator'
    
    
    #open the file handle for each file
    fh = ( get_appropriate_read_fileHandle(fastq_file,mode='r') for fastq_file in fastq_files)
    
    
    #get file handle for each file
    for fh in files:     
        for line in fh:
            if not line.startswith('@'):
                print 'File %s is not a fastq type file ' % fastq_file
                yield(None,None)
            else:
                read1 = fastq()
                read2 = fastq()
                read1.header= line.strip()
                read1.seq    = fh.next().strip()
                fh.next() #skip quality header
                read1.qual   = fh.next().strip()
                read1.type = 1
                
                read2.header= fh.next().strip()
                read2.seq    = fh.next().strip()
                fh.next() #skip quality header
                read2.qual   = fh.next().strip()
                read2.type   = 2                                    
                yield(read1,read2)



def get_genericfastq_iterator(fastq_files):
    '''
    given a fastq return read1 and read2 as a dictionary 
    through a generator
    '''
    
    SUB='get_pairedfastq_iterator'
    
    #open the file handle for each file
    files = ( get_appropriate_read_fileHandle(fastq_file,mode='r') for fastq_file in fastq_files )
    
    #get file handle for each file
    for fh in files:         
        for line in fh:
            if not line.startswith('@'):
                print 'File %s doesnt seem to be paired fastq formatted ' % fastq_file
                yield(None,None)
            else:
                read = fastq()
                read.header= line.strip()
                read.seq    = fh.next().strip()
                fh.next() #skip quality header
                read.qual = fh.next().strip()                                  
                yield(read)



def get_appropriate_read_fileHandle(file,mode='r'):
    '''
    given a file, try to detect file type with suffix 
    and return appropriate file handle, by default in read mode
    '''
    SUB='get_appropriate_read_fileHandle'
    if file.endswith('.gz'):
        import gzip
        fh = gzip.open(file,mode)
        return fh
    
    else:
        fh = open(file,mode)
        return fh
    


def get_guessed_fastq_prefix(fastqFile):
    '''
    given a file, try to detect file type with suffix 
    and return a guessed fastq file name
    '''
    SUB='get_guessed_fastq_prefix'
    fastqFile= str(os.path.basename(fastqFile)).strip()
    
    possible_extensions = ['.fastq.gz','fq.gz','.fastq','.fq' ]    
    for ext in possible_extensions:
        if fastqFile.endswith(ext):
            return fastqFile.replace(ext,'') 
        
    print "Can't detect the file extension as fastq type"
    sys.exit(2)


def writeFastq_read(read,fastq_handle=None,header_suffix=None):
    '''
    1. given a read object of class fastq(), write it to   
    take the read object generated by:
        a. opened filehandle if provided
        b. or STDOUT
    '''
    if header_suffix is not None:
        temp = read.header.split(' ',2)
        read.header = temp[0] + ':' + header_suffix + ' ' + temp[1]

    #composing the full read structure
    full_read = read.header + '\n' + read.seq + '\n+\n' + read.qual + '\n'


    if fastq_handle:
        fastq_handle.write(full_read)
    else:
        sys.stdout.write(full_read)
    



def create_lib_complexity_plot(fastqFile,mer_size=20):
    
    #def calculate_library_complexity(fastqfiles):
    
    report_uniqness_after_N_reads = 100000
    #middle_mer_start 
    
    randm_mer_fastq_complexity = {}
    starting_mer_fastq_complexity = {}
    count_uniq_random_mers =0
    count_uniq_starting_mers = 0
    
    bucket_percent_uniq_random_mers = []
    bucket_percent_uniq_starting_mers = []
    
    
    #output file name
    #fastq_prefix = get_guessed_fastq_prefix(fastqFile)
    fastq_prefix = fastqFile + '_kmer' + str(mer_size)
    baseDir = os.path.dirname(fastqFile) or '.'
    
    lib_complexity_plot_name = baseDir + '/' + fastq_prefix + '_lib_complexity_plot.png'
    lib_complexity_numbers_file = baseDir + '/' + fastq_prefix + '_lib_complexity.data'

    fh = open(lib_complexity_numbers_file,'w')

    for count,read in enumerate(get_genericfastq_iterator(fastqFile)):
        if count % report_uniqness_after_N_reads == 0 and count > 0:
            uniq_random_mer_percent = round( (count_uniq_random_mers/float(count+1))*100, 2)
            uniq_starting_mer_percent = round( (count_uniq_starting_mers/float(count+1))*100, 2)
            bucket_percent_uniq_random_mers.append(uniq_random_mer_percent)
            bucket_percent_uniq_starting_mers.append( uniq_starting_mer_percent  )
            outline = '%s \t %s \t %s \t %s \t %s \n' %  (count, count_uniq_random_mers,
                                                        count_uniq_starting_mers,
                                                        uniq_random_mer_percent,
                                                        uniq_starting_mer_percent)
            fh.write(outline)
            
        starting_mer = read.seq[:mer_size]
        random_mer_loc=random.randint(0,(len(read.seq)-mer_size))
        random_mer = read.seq[random_mer_loc:random_mer_loc+mer_size]
        
        if randm_mer_fastq_complexity.get(random_mer,None) is None:
            count_uniq_random_mers +=1
            randm_mer_fastq_complexity[random_mer] = 1
        
        if starting_mer_fastq_complexity.get(starting_mer,None) is None:
            count_uniq_starting_mers += 1
            starting_mer_fastq_complexity[starting_mer] = 1
        
    fh.close()
    
    #plotting the library complexity
    x_axis = [ x*report_uniqness_after_N_reads for x in xrange(1,len(bucket_percent_uniq_random_mers)+1) ]
    plt.plot(x_axis,bucket_percent_uniq_random_mers,label='rdm_%dmer' % mer_size)
    plt.plot(x_axis,bucket_percent_uniq_starting_mers,label='start_%dmer' % mer_size)
    plt.ylabel('percent unqiue')
    plt.xlabel('Read Counts')
    plt.legend()
    plt.savefig(lib_complexity_plot_name)






def create_fastq_length_histogram(fastqFile,logger=None):
    SUB='create_fastq_length_histogram'
    
    fastq_prefix = get_guessed_fastq_prefix(fastqFile)
    baseDir = os.path.dirname(fastqFile) or '.'
    plot_name = baseDir + '/' + fastq_prefix + '_read_length.png'
    
    
    if logger is None:
        logfile = baseDir + '/' + fastq_prefix + '_read_length.log'
        logger = utils._get_logger(logfile,logger_name=SUB)
    
    
    read_lengths = [ len(read.seq)  for read in (get_genericfastq_iterator(fastqFile))]
    np_read_lengths = np.array(read_lengths,dtype=int)

    #plot the histogram
    fig,ax = plt.subplots()
    bins = ax.hist(np_read_lengths,bins=50,color='darkgreen')
    ax.set_ylabel('read counts')
    ax.set_xlabel('read_length')
    ax.set_title(fastq_prefix)
    plt.savefig(plot_name)
    
    logger.info('created fastq read length histogram for : %s' % (fastqFile))
    logger.info('Read Length stats: %0.2f +/- %0.2f' %(np.mean(np_read_lengths),np.std(np_read_lengths)))
    
    
        
                                     



def filter_reads_with_Ns(fastqFile,max_Ns=1):
    SUB='filter_reads_with_Ns'
    
    fastq_prefix = get_guessed_fastq_prefix(fastqFile)
    baseDir = os.path.dirname(fastqFile) or '.'
    
    fastq_N_filtered_file = fastq_prefix + '_N' + str(max_Ns) + '_filtered.fastq' 
    fh = get_appropriate_read_fileHandle(fastq_N_filtered_file,mode='w')
    
    count_reads_with_Ns = 0 
    
    for count,read in enumerate((get_genericfastq_iterator(fastqFile))):
        
        if count > 0 and count % 100000 == 0:
            print '[%s]: processed %d reads' % (SUB,count) 
        
        if read.num_of_Ns() >= max_Ns:
            count_reads_with_Ns += 1
        else:
            writeFastq_read(read,fh)
    fh.close()
    
    print '[%s]: created filtered fastq %s with %d reads with Ns removed ' % (SUB,fastq_N_filtered_file,count_reads_with_Ns)
         

"""        

----------------------
Read name = HISEQ04:205:C18YHACXX:8:1307:8684:63767
Alignment start = 7 (-)
Cigar = 100M
Mapped = yes
Mapping quality = 37
----------------------
Base = C
Base phred quality = 31
----------------------
X0 = 1
X1 = 0
MD = 1A0T0C96
XG = 0
NM = 3
XM = 3
XO = 0
XT = U
-------------------
Alignment start position = Coliphage:7
TTCTGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTTGATAAAGCAGGAATTACTACTGCTTGTTTACGAATT


Orig read
@HISEQ04:205:C18YHACXX:8:1307:8684:63767 1:N:0:NNNNNN
AATTCGTAAACAAGCAGTAGTAATTCCTGCTTTATCAAGATAATTTTTCGACTCATCAGAAATATCCGAAAGTGTTAACTTCTGCGTCATGGAAGCAGAA
+
@@CFFFBDBFFHFIJJJGHIJIGIIIIJIIJJJJIGGIIHJIEIJJIJIGIJJJJJIGH<FHHGIGBCHI>B?CA?BEFFFF@;C>8>';@@CCDDCCD>
"""
        
    
        
        
    
    
    
    
    
    
     
    

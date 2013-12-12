#!/usr/bin/env python


import pysam
import os
import time
import re
import sys
import math
import matplotlib
import pybedtools
import pandas
#matplotlib.use('Agg') #background plotting
import matplotlib.pyplot as plt
import operator
from Bio import SeqIO


##updating the library look up path
dir = os.path.dirname(__file__)
sys.path.append(os.path.join(dir,'../lib'));

##loading in-house modules
import fastqUtils
import utils
import seqUtils


class TnSeq_Experiment:
    
    #class initiaizer
    def __init__(self,readmodel,expName,refGenome,fastq,
                 flanking_region=7,wobble=2,
                 plasmid = None,
                 refAnnotation = None,
                 residual_vector=None,
                 outPath=None,
                 annotation_model = None,
                 forceRun = False,
                 min_clean_readLen = 15):
        
        if outPath is None:
            outPath = os.path.getcwd()
        self.seq                               = readmodel
        self.expName                           = expName
        self.forceRun                          = forceRun
        self.refGenome                         = refGenome
        self.fastq                             = fastq
        self.refAnnotation                     = refAnnotation
        self.annotation                        = annotation_model
        self.outPath                           = outPath
        self.plasmid                           = plasmid  
        self.len                               = len(readmodel)
        self.min_clean_readLen                 = min_clean_readLen
        self.flanking_region                   = flanking_region
        self.wobble                            = wobble
        self.barcodeStart,self.barcodeEnd      = self._get_barcode_location()
        self.barcodeLen                        = self.barcodeEnd -  self.barcodeStart
        self.preBarcodeSeq,self.postBarcodeSeq = self._get_barcode_flakning_seq()
        self.EndSeq                            = readmodel[len(readmodel)-flanking_region:]
        
        #setting out files
        self.cleanedFastq                      = outPath + '/' + expName +'_cleaned.fastq.gz'
        self.barcodesBed                       = outPath + '/' + expName +'_cleaned_barcodes_best_hit.bed'
        self.barcodesFile                      = outPath + '/' + expName +'_cleaned_barcodes.txt'
        
        
        #counters
        self.counters                          = __TnSeqCounters__()
        
        #logger
        _logfileName  = outPath + '/' + expName + '_TnSeq_processing.log'
        self.logger                            = utils._get_logger(_logfileName)
        
        
        #create plot dir
        plot_dir = outPath + '/plots/'
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)
        self.plots_dir = plot_dir
        
        
        #create re compiled object for matching
        self.preBarcode       = re.compile(self.preBarcodeSeq)
        self.postBarcode      = re.compile(self.postBarcodeSeq)
        self.transposonEnd    = re.compile(self.EndSeq)
    
    #for string representation of the object
    def __str__(self):
        line  = ('expName: %s' % self.expName)
        line1 = ('Model: %s \nBarcode Location Start: %d, End: %d \nlength: %d' % (self.seq,self.barcodeStart,self.barcodeEnd,self.barcodeLen))
        line2 = ('pre barcode seq: %s \npost barcode seq: %s' % (self.preBarcodeSeq, self.postBarcodeSeq))
        line3 = ('Mode end seq %s ' % (self.EndSeq))
        line4 = ('DNA start position %d' % self.len )
        return(line+'\n'+line1+'\n'+line2+'\n'+line3+'\n'+line4)
    
    #get the barcode location
    def _get_barcode_location(self):
        match = re.search('(N+)',self.seq)
        if len(match.groups()) != 1:
            self.logger.error('Invalid Model %s ' % self.seq)
            sys.exit(2)
        else:
            return(match.start(),match.end())
    
    def _get_barcode_flakning_seq(self):
        preBarcodeSeq = self.seq[self.barcodeStart-self.flanking_region:self.barcodeStart]
        postBarcodeSeq = self.seq[self.barcodeEnd:self.barcodeEnd+self.flanking_region]
        return(preBarcodeSeq,postBarcodeSeq)
    
    

class __TnSeqCounters__:
    def __init__(self):
        
        self.read = 0
        
        self.proper_barcode = 0
        self.improper_barcode = 0
        self.no_barcode = 0
        
        self.improper_transposon_end = 0
        self.transposon_end_not_found = 0
        self.proper_transposon_end = 0
        
        self.clean_read = 0
        self.clean_but_short_read = 0
        
        self.uniq_mapped_reads = 0
        self.multi_mapped_reads = 0
        self.unmapped_reads = 0
        self.plasmid_mapped_reads = 0
        
        #count of all barcodes
        self.barcodes       = 0
        
        #barcode maps to just one read
        self.barcodes_singletons    = 0
#END __TnSeqCounters__ class





def clean_TnSeq_fastq(TnSeqExp):
    """
    given the TnSeq_Experiment object
    1. take the fastqs and for each read check if it fits the expected read model

    """
    
    SUB = 'clean_TnSeq_fastq'
    
    #create temp variable for short access
    counters = TnSeqExp.counters
    logger = TnSeqExp.logger
    
    tnSeq_fastq_cleaning_indicator = TnSeqExp.outPath + '/.fastq_cleaning.done'
    #check if analysis was done earlier and complete
    if os.path.exists(tnSeq_fastq_cleaning_indicator) and os.path.exists(TnSeqExp.cleanedFastq) and not TnSeqExp.forceRun:
        logger.info('[%s]: Tn-Seq fastq cleaning already done \t %s' % (SUB,TnSeqExp.cleanedFastq))
        return 
        
    #open the input fastqs
    fastq_iterator = fastqUtils.get_genericfastq_iterator(TnSeqExp.fastq)
    #open a cleaned fastq file for writing filtered/cleaned reads
    fh_cleaned_fastq = fastqUtils.get_appropriate_read_fileHandle(TnSeqExp.cleanedFastq,mode='w')
    #out histogram file name
    out_TnSeq_QC_plot = TnSeqExp.plots_dir + TnSeqExp.expName + '_readClean_hist.png'
    
    ##temp holder for all match locations for plotting histogram
    read_match_locations = []
    
    
    logger.info('[%s]: Starting : cleaning reads in fastq ' % SUB)
    prog_start = time.clock()
    start_loop = prog_start
    for read in fastq_iterator:
        counters.read += 1
        if counters.read % 1000000 == 0:
            end_loop = time.clock()
            time_taken = end_loop - start_loop
            logger.info('[%s]: processed %d reads in %d secs total time: %d' % (SUB,counters.read, time_taken,end_loop-prog_start))
            start_loop = time.clock()
    
        #find the barcode in the read
        barcode = _find_barcode_location_in_read(TnSeqExp,read,read_match_locations)
        if barcode:
        #if barcode find the cleaned read if transposon end is obtained 
            cleaned_read = _get_cleaned_DNA_subread(TnSeqExp,read,read_match_locations)
            if cleaned_read:
                fastqUtils.writeFastq_read(cleaned_read,fh_cleaned_fastq,header_suffix=barcode)
                
    #close the cleaned fastq handle
    fh_cleaned_fastq.close()
    logger.info('[%s]: Finished : cleaning %d reads in fastq ' % (SUB,counters.read))

    
      
    ########
    #plot the histogram
    ########
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(111)
    ax.hist(read_match_locations,bins=50,facecolor='green',alpha=0.75)
    ax.axvline(TnSeqExp.barcodeStart,color='red',linestyle='dashed',linewidth=2,label='barcodeStart')
    ax.axvline(TnSeqExp.barcodeEnd,color='blue',linestyle='dashed',linewidth=2,label='barcodeEnd')
    ax.axvline(TnSeqExp.len,color='black',linestyle='dashed',linewidth=2,label='transposonEnd')
    ax.legend(prop={'size':8})
    plt.savefig(out_TnSeq_QC_plot)
    logger.info('[%s]: created read cleaning histogram %s ' % (SUB,out_TnSeq_QC_plot))
    
    
    #store the read cleaning stats : could be improved
    #########
    #TEMP
    #########
    logger.info('[%s]: read cleaning stats' % SUB)
    tnSeq_read_cleaning_stats = TnSeqExp.outPath + '/' + TnSeqExp.expName + '_readClean_stats.txt'
    temp_fh = open(tnSeq_read_cleaning_stats,'w')
    order_of_vars = ['read','proper_barcode','improper_barcode','no_barcode',
                 'proper_transposon_end','improper_transposon_end','transposon_end_not_found',
                 'clean_read','clean_but_short_read']
    for var in order_of_vars:
        val = counters.__dict__.get(var,'NA')
        temp_fh.write('%s \t %s\n' % (var,val))
        logger.info('%s \t %s' % (var,val))
    temp_fh.close()
    #END temp
    ###########################
    
    

    #indicating the read cleaning is done
    fh=open(tnSeq_fastq_cleaning_indicator,'w')
    fh.close()
    
#END clean_TnSeq_fastq
    



def get_TnSeq_cleaned_barcodes_df(mappedBam,TnSeqExp,**kwargs):
    """
    output
    1. get a cleaned barcode list as pandas data frame
    2. create  bam file of uniquely mapped reads
    """
    SUB = 'get_TnSeq_cleaned_barcodes_df'
    
    #temp vars
    counters = TnSeqExp.counters
    logger = TnSeqExp.logger
    
    logger.info('[%s]: creating barcode data frame from uniquely mapped reads' % SUB)
    
    #open the input bamFile
    in_bam_fh = pysam.Samfile(mappedBam,'rb')
    
    #create a new bamfile for uniquely mapped reads
    uniqMapped_bam_outFile = TnSeqExp.mapDir + '/' + TnSeqExp.expName + '_uniqMapped_reads.bam'
    TnSeqExp.uniqMapped_bam = uniqMapped_bam_outFile
    uniqMapped_bam_outFh = pysam.Samfile(uniqMapped_bam_outFile,'wb',template=in_bam_fh)
    
    
    #barcode counter
    #nested dict
    #delete any existing barcodes dict if present, can create a mess with mutable dicts
    try:
        del barcodes
    except:
        pass
    barcodes = utils.NestedDict()
    
    
    
    #get the headers for the plasmids: read mapping to these will be ignored
    if TnSeqExp.plasmid:
        list_of_contig_names_to_ignore = [ str(record.id).split(' ')[0] for record in SeqIO.parse(open(TnSeqExp.plasmid,'r'), "fasta") ]
        #get tid for each contigs name
        list_of_contig_ids_to_ignore = [ in_bam_fh.gettid(contig) for contig in list_of_contig_names_to_ignore]
    else:
        list_of_contig_ids_to_ignore = []
    
    
    #create a nested strucutre to hold the uniq barcodes
    #and also create uniquely mapped reads
    prog_start = time.clock()
    start_loop = prog_start
    for count,read in enumerate(in_bam_fh):
        
        #for logging
        if count % 1000000 == 0 and count > 0:
            end_loop = time.clock()
            time_taken = end_loop - start_loop
            logger.info('[%s]: processed %d reads in %d seconds, total: %d' % (SUB,count,time_taken,end_loop-prog_start))
            start_loop = time.clock()
        #mainly for debugging
        if kwargs.get('debug',None) and count == 1000000:
            barcodes_df = _process_barcodes(barcodes)
            return barcodes_df
        
        if read.is_unmapped:
            counters.unmapped_reads += 1
            continue
        #ignore reads mapping to plasmid 
        elif read.tid in list_of_contig_ids_to_ignore :
            counters.plasmid_mapped_reads +=1
            continue        
        #uniquely mapped reads
        elif read.mapq > 0:
            
            insertion_point = read.pos
            if read.is_reverse:
                #add in the length of the read when it reverse maps to find the right point of insertion
                insertion_point = read.pos + read.qlen 
            
            strand = '-' if read.is_reverse else '+'
            counters.uniq_mapped_reads += 1
            barcode   = read.qname.split(':')[-1]
            reference = in_bam_fh.getrname(read.rname)
            if insertion_point not in barcodes[barcode][reference][strand]:
                barcodes[barcode][reference][strand][insertion_point] = 1
            else:
                barcodes[barcode][reference][strand][insertion_point] += 1
            uniqMapped_bam_outFh.write(read)
        #if read is multi mapped
        else:
            counters.multi_mapped_reads += 1
    
    logger.info('[%s]: Finished processing reads in total %d seconds' % (SUB,time.clock()-prog_start))
    #process the nested dict to create a list ultimately leading to pandas data frame
    barcodes_df = _process_barcodes(barcodes)
    
    
    #store the bam stats : could be improved
    #########
    #TEMP
    #########
    logger.info('[%s]: read mapping stats' % SUB)
    tnSeq_mappedBam_stats = TnSeqExp.mapDir + '/' + TnSeqExp.expName + '_bamStats.txt'
    temp_fh = open(tnSeq_mappedBam_stats,'w')
    order_of_vars =  ['uniq_mapped_reads','multi_mapped_reads','unmapped_reads','plasmid_mapped_reads']
    for var in order_of_vars:
        val = counters.__dict__.get(var,'NA')
        temp_fh.write('%s \t %s\n' % (var,val))
        logger.info('%s \t %s' % (var,val))
    temp_fh.close()
    #END temp
    ###########################
    
    logger.info('[%s]: Finished creating barcode data frame from uniquely mapped reads' % SUB)

    return barcodes_df            
#END get_TnSeq_cleaned_barcodes_df




def render_basic_TnSeq_QC_plots(TnSeqExp,barcodes_df,plot_suffix=''):
    
    #outpath
    base_path = TnSeqExp.plots_dir + TnSeqExp.expName
    
    ######
    #Plot1
    ######
    #get the total number of reads per barcode
    num_reads_per_barcode = barcodes_df.groupby('barcode',sort=False,as_index=False)['numReads'].sum()
    #histogram : number of reads per barcode
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(111)
    ax.hist(num_reads_per_barcode.numReads.map(lambda x: math.log(x,2)),color='#FF7F00',bins=15,rwidth=.9)
    ax.set_xlabel('#read per barcode')
    ax.set_ylabel('number of such barcodes positions')
    ax.set_title('TnSeq: Number of reads per barcode :%s' % TnSeqExp.expName)
    ax.grid(True)
    ax.set_xlim(0,10)
    ax.set_xticks([0,1,2,3,4,5,6,7,8,9,10])
    ax.set_xticklabels([int(2**tick) for tick in ax.get_xticks()])
    #save plot
    hist_file_num_reads_per_barcode = base_path + '_hist_reads_per_barcode_%s.png' % plot_suffix
    plt.savefig(hist_file_num_reads_per_barcode)
    TnSeqExp.logger.info('Created histogram for number of reads per barcode for %s data' % plot_suffix)


    #######
    #Plot2
    #######
    #get the total number of unique barcodes per genome position
    num_barcodes_per_genome_positions = barcodes_df.groupby(['chr','strand','pos'],sort=False,as_index=False).apply(lambda x: x.barcode.nunique()).reset_index()
    num_barcodes_per_genome_positions['num_barcodes'] = num_barcodes_per_genome_positions[0]
    del num_barcodes_per_genome_positions[0]
    #histogram : number of barcodes per genome position
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(111)
    ax.hist(num_barcodes_per_genome_positions.num_barcodes.map(lambda x: math.log(x,2)),color='#4DAF4A',bins=15,rwidth=.9)
    ax.set_title('TnSeq: Number of barcode per genome location :%s' % TnSeqExp.expName)
    ax.set_xlabel('#barcode hits per genome position')
    ax.set_ylabel('number of genome positions')
    ax.grid(True)
    ax.set_xlim(0,10)
    ax.set_xticks([0,1,2,3,4,5,6,7,8,9,10])
    ax.set_xticklabels([int(2**tick) for tick in ax.get_xticks()])
    #save plot
    hist_file_barcodes_per_genome_location = base_path + '_hist_barcodes_per_genome_location_%s.png' % plot_suffix
    plt.savefig(hist_file_barcodes_per_genome_location)
    
    
    #######
    #Plot3
    #######
    #get the total number of unique insertions per barcode
    num_genome_positions_per_barcode = barcodes_df.groupby(['barcode'],sort=False,as_index=False).apply(lambda x: x.pos.nunique()).reset_index()
    num_genome_positions_per_barcode['num_insertions'] = num_genome_positions_per_barcode[0]
    del num_genome_positions_per_barcode[0]
    #histogram : number of genome positions where barcode has inserted 
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(111)
    ax.hist(num_genome_positions_per_barcode.num_insertions.map(lambda x: math.log(x,2)),color='#984EA3',bins=15,rwidth=.9)
    ax.set_title('TnSeq: Number of uniq insertions per barcode :%s' % TnSeqExp.expName)
    ax.set_xlabel('#uniq insertions per barcode')
    ax.set_ylabel('number of barcodes')
    ax.grid(True)
    ax.set_xlim(0,10)
    ax.set_xticks([0,1,2,3,4,5,6,7,8,9,10])
    ax.set_xticklabels([int(2**tick) for tick in ax.get_xticks()])
    #save plot
    hist_file_insertions_per_barcode = base_path + '_hist_insertions_per_barcode_%s.png' % plot_suffix
    plt.savefig(hist_file_insertions_per_barcode)

#END render_basic_TnSeq_QC_plots


def create_nucleotideVar_plot_around_insertionSite(TnSeq_Exp,flanking_region):
    
    #creating a bedtool object
    unique_barcodes_bedtool = pybedtools.BedTool(TnSeq_Exp.barcodesBed)
    
    uniq_barcodes_with_flanking_region = unique_barcodes_bedtool.slop(l=flanking_region-1,r=flanking_region,
                                                                      g=seqUtils.get_fasta_sequence_length(TnSeq_Exp.refGenome))

    #get the flanking region sequence
    flanking_seq_insertion_sites = uniq_barcodes_with_flanking_region.sequence(TnSeq_Exp.refGenome,s=True,name=True,tab=True)
    
    #creating a pandas data frame
    _flanking_seq = []
    seqs = ( line.split('\t')[1].replace('\n','') for line in open(flanking_seq_insertion_sites.seqfn))
    _ = [ _flanking_seq.append(list(seq)) for seq in seqs ]
    df = pandas.DataFrame(_flanking_seq)
    
    def _get_nucleotide_percent(col):
        return col.value_counts().astype(float)/col.count()
    #calculate the % of each nucleotide at each location
    nuc_variation_around_insertion_site = df.apply(_get_nucleotide_percent,axis=0)
    
    #creating col names
    col_names = range(-flanking_region,0,1)
    col_names.extend(range(1,flanking_region+1,1))
    nuc_variation_around_insertion_site.columns = col_names
    
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(111)
    nuc_variation_around_insertion_site.T.plot(ax=ax)
    ax.set_title('%s: Nucletotide Variation around transposon insertion site' % TnSeq_Exp.expName)
    ax.set_xlabel('region around insertion site(0)')
    ax.set_ylabel('percentage')
    ax.axvline(0,color='black',lw=3)
    plot_file_name = TnSeq_Exp.plots_dir + TnSeq_Exp.expName + '_Nucleotide_variation_transposon_insertion_site.png'
    plt.savefig(plot_file_name)

#END create_nucleotideVar_plot_around_insertionSite



def get_formatted_barcodes_df(cleaned_barcodes_df,num_top_hits=3):
    #group the barcodes
    grouped = cleaned_barcodes_df.groupby('barcode',sort=False)
    #temp holder for all the barcodes
    data_list = []
    #loop over each barcode
    for count,(name,group) in enumerate(grouped):
        _list = [name]
        _topN = group.sort(columns='numReads',ascending=False)[:num_top_hits]
        _topN = _topN[['chr','strand','pos','numReads']]
        [ _list.extend(row) for row in _topN.itertuples(index=False)] 
        data_list.append(_list)
    
    #create a pandas df
    df = pandas.DataFrame(data_list)
    
    #set the column names    
    column_names = ['barcode',
                    'h1_chr','h1_strand','h1_pos','h1_reads',
                    'h2_chr','h2_strand','h2_pos','h2_reads',
                    'h3_chr','h3_strand','h3_pos','h3_reads'
                    ]
    df.columns = column_names
    
    #eliminate barcodes with just one read in the tophit
    df = df[df.h1_reads > 1]

    return df


    
    
    
#NOT Used    
# now keep only the top 3 hits for each barcode
# this is slow
# as each group is sorted and a new data frame is created for each barcode 
def _get_topN_hits_per_barcode(group,N=3):
    _topN = group.sort(columns='numReads',ascending=False)[:N]
    rows  =_topN.shape[0]
    return pandas.DataFrame(_topN[['chr','strand','pos','numReads']].values.reshape(1,rows*4))
    #[ _list.extend(row) for row in _topN[['chr','strand','pos','numReads']].itertuples(index=False)] 

#END _get_topN_hits_per_barcode





def _process_barcodes(barcodes):
    #temp list containers
    cleaned_barcodes = []

    #four nested loop could do better
    #parsing the following structure
    #barcodes[barcode][reference][strand][read.pos]
    for barcode,map_locations in barcodes.iteritems():
        for ref,per_ref_map_locations in map_locations.iteritems():
            for strand,per_ref_per_strand_map_locations in per_ref_map_locations.iteritems():
                for position,num_reads in per_ref_per_strand_map_locations.iteritems():
                    out_line = [barcode,strand,ref,position,num_reads]
                    cleaned_barcodes.append(out_line)
    
    cols = ['barcode','strand','chr','pos','numReads']
    df = pandas.DataFrame.from_records(cleaned_barcodes,columns=cols)
    return df
#END _process_barcodes


def _find_barcode_location_in_read(model,read,read_match_locations,minQuality=20):
   
    #do the matching
    m_start = model.preBarcode.search(read.seq)
    m_end   = model.postBarcode.search(read.seq)
    
    #create temp variable for short access
    logger = model.logger
    counters = model.counters
    
    if m_start and m_end:
        barcode_start = m_start.end()
        barcode_end = m_end.start()
        read_match_locations.extend((barcode_start,barcode_end))
        #check on multiple conditions
        if ( 
             (model.barcodeStart - model.wobble <= barcode_start <= model.barcodeStart + model.wobble) and  
             (model.barcodeEnd   - model.wobble <= barcode_end   <= model.barcodeEnd + model.wobble)   and 
             (barcode_end-barcode_start == model.barcodeLen)
           ):
            counters.proper_barcode += 1
            #return read.seq[barcode_start:barcode_end+1]
            return read.seq[barcode_start:barcode_end]
        else:
            counters.improper_barcode += 1
            #print 'improper barcode : start %d : end : %d' % (barcode_start,barcode_end)
            #print read.seq
            return None
    #no barcode found
    else:
        counters.no_barcode += 1
        return None





def _get_cleaned_DNA_subread(model,read,read_match_locations):
    
    #create temp variable for short access
    logger = model.logger
    counters = model.counters

    #do the matching
    transposon_end_match = model.transposonEnd.search(read.seq)
    if transposon_end_match:
        start = transposon_end_match.start()
        end = transposon_end_match.end()
        read_match_locations.append(end)
        #check if the end of transposon seq is in the expected region in the read
        if ( 
               #extreme left end for start                                  #extreme right end for start        
             ( model.len - model.flanking_region - model.wobble <= start <= model.len - model.flanking_region + model.wobble) and  
             ( model.len - model.wobble <= end   <= model.len + model.wobble)
           ):
            counters.proper_transposon_end += 1
            #check if minimum read length criteria is satisfied
            if len(read.seq) - end >= model.min_clean_readLen:
                counters.clean_read += 1 
                sub_dna_read = fastqUtils.fastq(header=read.header,seq=read.seq[end+1:],qual=read.qual[end+1:])
                return sub_dna_read
            else:
                counters.clean_but_short_read += 1 
                return None
        else:
            counters.improper_transposon_end += 1
            return None
    else:
        counters.transposon_end_not_found += 1
        return None
    



    

#!/usr/bin/env python

import sys
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import os
import time
import re

import argparse
from EnvironmentModules import *
module(['load','biopython'])
from Bio import SeqIO
import numpy as np

import pybedtools



##updating the library look up path
dir = os.path.dirname(__file__)
sys.path.append(os.path.join(dir,'../lib'));


##loading in-house modules
import fastqUtils
import utils
import seqUtils
import Mappers
import TnSeq




def classify_barcodes(barcode):
    """
    ########
    #Basis
    ########
    #1. uniq insertions 
        > 10 reads at hit 1
        #reads-hit1 / reads-hit2 > 10
    
    #2. multi site insertions
        > 10 reads at hit 1
        #reads-hit1 / reads-hit2  < 2

    #3. ambiguous-low coverage
        < 10 reads
    
    #4. ambiguous - likely PCR chimera
        > 10 reads at hit 1
         2  < #reads-hit1 / reads-hit2  < 10
    """
    
    #hit 1 sanity check
    #potential probs
    if np.isnan(barcode['h1_reads']) or barcode['h1_reads'] == 0:
        print 'Error: number of reads in hit 1 is nan or 0'
        sys.exit(2)
    else:
        h1_reads = barcode['h1_reads']

    #if h2_reads is None/Nan then use one
    if np.isnan(barcode['h2_reads']):
        h2_reads = 1
    else:
        h2_reads = barcode['h2_reads']
    
    #get the ratio
    ratio = h1_reads / h2_reads
    
    if h1_reads >= 10 and ratio >= 10:
        barcode_type = 'unique'
    elif h1_reads > 10 and ratio < 2:
        barcode_type = 'multi_site_insertion'
    elif h1_reads < 10:
        barcode_type = 'low_coverage'
    elif h1_reads > 10 and ratio > 2 and ratio < 10:
        barcode_type = 'ambiguous_PCR_chimera'
    else:
        barcode_type = 'unclassified'
        
    return barcode_type




def get_user_args():

    #user command line args processing
    description = '''
                    Script to process the TnSeq library 
                    ***Input: is Single end read only***
                  '''
    parser = argparse.ArgumentParser(description=description,add_help=True)
    parser.add_argument('--expName'      ,'-eN'   ,dest='exp_name'      ,default=None,  metavar='', help='Name of the experiment')
    parser.add_argument('--refGenome'    ,'-refG' ,dest='ref_genome'    ,default=None,  metavar='', help='reference genome')
    parser.add_argument('--refAnnotation','-refA' ,dest='ref_annotation',default=None,  metavar='', help='reference annotation in bed format')
    parser.add_argument('--readModel'    ,'-rM'   ,dest='read_model'    ,default=None,  metavar='', help='model of the read| fasta file')
    parser.add_argument('--plasmidSeq'   ,'-ps'   ,dest='plasmid_seq'   ,default=None,  metavar='', help='plasmid sequence to filter the reads against')
    parser.add_argument('--wobble'       ,'-w'    ,dest='wobble'        ,default=2,     metavar='' ,type = int, help='number of shifts allowed in determining the location of barcode')
    parser.add_argument('--flankRegion'  ,'-fr'   ,dest='flankRegion'   ,default=7,     metavar='' ,type = int, help='#bases to use match the barcode region')    
    parser.add_argument('--num_cores'    ,'-nc'   ,dest='num_cores'     ,default=2,     metavar='' ,type = int, help='number of cores to use for mapping')
    parser.add_argument('--outPath'      ,'-o'    ,dest='outPath'       ,default=os.getcwd(),   metavar='', help='output path default:cwd')
    parser.add_argument('--fastq'        ,'-fastq',dest='fastq'         ,default=None,  metavar='' ,nargs ='+',help='list of fastqs')
    parser.add_argument('--forceRun'     ,'-f'    ,dest='forceRun'      ,default=False, action='store_true', help='Force run even if previously complete')
    
    
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    
    sys.stderr.write("Using the following bwa:")
    sys.stderr.write(subprocess.check_output(['which', 'bwa']))
    sys.stderr.write("Using the following samtools:")
    sys.stderr.write(subprocess.check_output(['which', 'samtools']))
    
    #parse the args
    user_args = parser.parse_args()
    
    #BUG
    #right now the metavar='' and required=True cant be used in conjunction
    #therefor explicit sanity check is needed
    reqd_options = ['exp_name','ref_genome','read_model','fastq',]
    for opt in reqd_options:
        if eval('user_args.%s' % opt) is None:
            print '[ERROR]: Required option "%s" is not present \n\n' % opt
            parser.print_help()
            sys.exit()

    
    #read the read-model and store the first sequence
    for record in SeqIO.parse(open(user_args.read_model,'r'), "fasta"):
        user_args.read_model = '%s' % record.seq
        break

    #check to make sure output directory exists
    if os.path.exists(user_args.outPath):
        if not os.path.isdir(user_args.outPath):
            print '[ERROR]: %s already exists as a file' % user_args.outPath
            parser.print_help()
            sys.exit()
    else:
        os.mkdir(user_args.outPath,0755)

    return user_args
#end get_user_args
    
    


def main():
    
    #get the user args
    user_args = get_user_args()
    
    
    #create a unified fasta for creating an index if plasmid file is present
    if user_args.plasmid_seq:
        fastas = [user_args.ref_genome,user_args.plasmid_seq]
        final_reference_file = utils.merge_files(fastas,outPath=user_args.outPath,outFile='merged_ref.fasta')
        user_args.ref_genome = final_reference_file

    
    #setup the experiment meta-data
    TnSeq_Exp = TnSeq.TnSeq_Experiment(readmodel        = user_args.read_model,
                                       expName          = user_args.exp_name,
                                       refGenome        = user_args.ref_genome,
                                       refAnnotation    = user_args.ref_annotation,
                                       flanking_region  = user_args.flankRegion,
                                       wobble           = user_args.wobble,
                                       outPath          = user_args.outPath,
                                       fastq            = user_args.fastq,
                                       plasmid          = user_args.plasmid_seq,
                                       forceRun         = user_args.forceRun
                                       )
    
    #get the logger object
    logger = TnSeq_Exp.logger
    
    #log the user arguments
    logger.info('###########################')
    logger.info("Command line options used ")
    [ logger.info('%s=%s' %(arg_name,value)) for arg_name,value in user_args._get_kwargs()]
    logger.info('###########################')
    logger.info('\n')
    logger.info('TnSeq Experiment')
    logger.info(str(TnSeq_Exp))
    

    
    


    ########################
    #TnSeq processing starts
    ########################    
    
    ####
    #1. clean the fastq
    ####
    TnSeq.clean_TnSeq_fastq(TnSeq_Exp)
    
    
    ######    
    #2.map the reads
    ######
    #map the reads with BWA
    bwa_mapper = Mappers.BWA(TnSeq_Exp.refGenome,read1=TnSeq_Exp.cleanedFastq,
                             ncores = user_args.num_cores, logger=logger,
                             forceRun = user_args.forceRun, outPath = user_args.outPath)
    #map the reads
    bwa_mapper.map_reads()
    
    #set the mapping dir in TnSeq object
    TnSeq_Exp.mapDir = bwa_mapper.mapDir
    
          
    #####
    #3. Parse the Mapped reads and Analyze Barcodes Summary
    #####
    barcodes_df = TnSeq.get_TnSeq_cleaned_barcodes_df(bwa_mapper.mappedBam,TnSeq_Exp,debug=False)
    #plot the raw barcode data
    TnSeq.render_basic_TnSeq_QC_plots(TnSeq_Exp,barcodes_df,plot_suffix='raw')
    
    
   
    ####
    #4. barcode filtering
    ####
    #select barcodes that have > 1 read associated with them
    num_reads_per_barcode = barcodes_df.groupby('barcode',sort=False,as_index=False)['numReads'].sum()
    selected_barcodes = num_reads_per_barcode[num_reads_per_barcode.numReads > 1]
    cleaned_barcodes_df = barcodes_df[ barcodes_df['barcode'].isin(selected_barcodes.barcode) ]
    #plot the cleaned barcode data
    TnSeq.render_basic_TnSeq_QC_plots(TnSeq_Exp,cleaned_barcodes_df,plot_suffix='filtered')
    
    
    ####
    #5. get the formatted data frame : computationally expensive step 
    ####
    barcodes_clnd_fmt = TnSeq.get_formatted_barcodes_df(cleaned_barcodes_df)
    
    
    ####
    #6. classify the barcodes
    ####
    barcodes_clnd_fmt['type'] = barcodes_clnd_fmt.apply(classify_barcodes,axis=1)
    #save the cleaned barcodes file
    barcodes_clnd_fmt.to_csv(TnSeq_Exp.barcodesFile,sep='\t',na_rep=0,header=True,index=False,float_format='%d')
    
    
    ####
    #7. MORE PLOTS
    ###
    ####
    #get the fraction of barcode types as a pie chart
    ####
    barcode_type_fraction = barcodes_clnd_fmt.type.value_counts().map(lambda x: float(x)/barcodes_clnd_fmt.type.count())
    fig =plt.figure(figsize=(10,8))
    ax = fig.add_subplot(111)
    _ = ax.pie(barcode_type_fraction.values,labels=barcode_type_fraction.index,
            shadow=True,explode=[.05,.10,0.02,0,0],autopct='%1.1f%%',
            colors=['#984EA3','#4DAF4A','#377EB8','#FF7F00','#E41A1C'])
    ax.set_title("Fraction of Barcodes / type")
    plot_file_name = TnSeq_Exp.plots_dir + '/' + TnSeq_Exp.expName + '_barcodes_classification_pieChart.png'
    plt.savefig(plot_file_name)
    
    
    
    #Plot: number of reads in top hit v/s second best loci
    fig =plt.figure(figsize=(10,8))
    ax = fig.add_subplot(111)
    ax.scatter(barcodes_clnd_fmt.h1_reads,barcodes_clnd_fmt.h2_reads,color='tomato')
    ax.grid(True,linestyle='-',color='0.50')
    ax.set_title('#reads on top hit v/s second best barcode hit')
    ax.set_xlabel('#reads barcode top hit')
    ax.set_ylabel('#reads barcode second hit')
    ax.set_xlim( (0,max(barcodes_clnd_fmt.h1_reads)) )
    ax.set_ylim( (0,max(barcodes_clnd_fmt.h2_reads)) )
    ax.plot(barcodes_clnd_fmt.h2_reads,barcodes_clnd_fmt.h2_reads,color='black',lw=3)
    ax.legend(['y = x'], loc='upper right')
    plot_file_name = TnSeq_Exp.plots_dir + TnSeq_Exp.expName + '_barcode_top_vs_second_hit_scatterplot.png'
    plt.savefig(plot_file_name)
    
    
    
    ####
    #8. create a bed file for top uniquely classified barcodes
    ####
    #create the bed file for uniquely classified barcodes
    _temp_df = barcodes_clnd_fmt[barcodes_clnd_fmt.type == 'unique'][['h1_chr','h1_pos','barcode','h1_reads','h1_strand']]
    _temp_df['start'] = _temp_df.h1_pos - 1
    _temp_df['end']   = _temp_df.h1_pos
    _temp_df[['h1_chr','start','end','barcode','h1_reads','h1_strand']].to_csv(TnSeq_Exp.barcodesBed,sep='\t',header=False,index=False)
    
    
    #nucleotide variation plot
    TnSeq.create_nucleotideVar_plot_around_insertionSite(TnSeq_Exp,flanking_region=5)
    
    
    #intersection with reference gene annotation
    annotation_bedtool = pybedtools.BedTool(TnSeq_Exp.refAnnotation)
    transposon_insertions_per_gene = annotation_bedtool.intersect(b=TnSeq_Exp.barcodesBed,c=True)
    out_file = TnSeq_Exp.outPath + '/' + TnSeq_Exp.expName + '_insertions_per_gene.txt'
    transposon_insertions_per_gene.saveas(out_file)
    
    
    
    
if __name__ == "__main__":
    main()
    

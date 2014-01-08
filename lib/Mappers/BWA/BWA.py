#!/usr/bin/env python

import os
import sys
import shutil
import subprocess

#updating the library look up path
dir = os.path.dirname(__file__)
sys.path.append(os.path.join(dir,'../lib'));
#loading in-house modules
import utils
import seqUtils
import fastqUtils




def temp_map_reads_with_BWA(bwa_obj):				
	"""
	WARNING:
	most of the arguments are hard coded
	needs significant improvement
	current only done to get the mapping done for TnSeq project
	"""
	
	SUB = 'temp_map_reads_with_BWA'
	
	#temp files  
	read_1_sai	  = bwa_obj.mapDir + '/' + bwa_obj.fastq_prefix  + '.sai'
	
	#open file handles
	bwa_mapping_status_indicator = bwa_obj.mapDir + '/.' + bwa_obj.fastq_prefix  + '_mapping.done'
	read1_fh = open(read_1_sai,'wb')
	bwa_map_log_fh =open(bwa_obj.bwa_map_log,'w')

	 #check if analysis was done earlier and complete
	#if os.path.exists(bwa_mapping_status_indicator) and os.path.exists(bwa_obj.mappedBam) and not bwa_obj.forceRun:
	#	bwa_obj.logger.info('[%s]: Tn-Seq cleaned read mapping already done \t %s' % (SUB,bwa_obj.mappedBam))
	#	return 
		
	
	#map
	args = ['bwa','aln', '-t', str(bwa_obj.ncores), bwa_obj.reference,bwa_obj.read1]
	#args = ['bwa','aln', '-I', '-t', str(bwa_obj.ncores), bwa_obj.reference,bwa_obj.read1]
	bwa_obj.logger.info('mapping reads with BWA to reference %s' % (bwa_obj.reference))
	bwa_obj.logger.info('mapping options used: %s' % (' '.join(map(str,args))) )
	bwa_map_log_fh.write("[%s] running %s\n" % (SUB, " ".join(map(str,args))))
	bwa_map_log_fh.flush()
	return_code = subprocess.check_call(args,stdout=read1_fh,stderr=bwa_map_log_fh)
	if return_code == 0:
		bwa_obj.logger.info('mapping completed successfully %s' % (bwa_obj.reference))
	else:
		bwa_obj.logger.info('bwa mapping returned non zero exit code %d ..exiting' % (bwa_obj.reference))
		sys.exit(3)


	#binary sai file to bam
	sai_to_sam = ['bwa', 'samse',bwa_obj.reference,read_1_sai,bwa_obj.read1 ]
	sam_to_bam = ['samtools', 'view', '-bS', '-']
	samtools_view_log_fh = open(bwa_obj.mapDir + '/' + bwa_obj.fastq_prefix + '.view.log','w')
	
	bwa_obj.logger.info('converting *sai file to bam.....')
	bwa_map_log_fh.write("[%s] running %s" % (SUB, " ".join(sai_to_sam)))
	p1 = subprocess.Popen(sai_to_sam,stdout=subprocess.PIPE,stderr=bwa_map_log_fh)
	p2 = subprocess.Popen(sam_to_bam,stdin = p1.stdout,stdout=open(bwa_obj.mappedBam,'w'),stderr=samtools_view_log_fh)
	p2.communicate()
	p1.stdout.close()
	bwa_obj.logger.info('bam file produced %s' % bwa_obj.mappedBam)
	
	#close file handles	
	read1_fh.close()
	bwa_map_log_fh.close()
	samtools_view_log_fh.close()
	os.remove(read_1_sai)
	
	
	#indicating the read cleaning is done
	fh=open(bwa_mapping_status_indicator,'w')
	fh.close()



class BWA:
	"""
	Map the reads with BWA
	"""
	
	#currently only works for read 1
	#this is done in the interest of time
	def __init__(self,reference,logger=None,combined_fastq = None,
				 forceRun = False,
				 read1 = None,
				 read2 = None,
				 ncores = 2,
				 outPath=os.getcwd()):
		
		if read1 is None and combined_fastq is None:
			print '[BWA]: Error: No fastq file provided\n'
			sys.exit(2)
		
	
		self.reference		  =   os.path.abspath(reference)
		self.forceRun		   =   forceRun
		self.ncores			 =   ncores
		self.combined_fastq	 =   combined_fastq
		self.read1			  =   read1
		self.read2			  =   read2
		self.outPath			=   outPath
		
		
		
		#create a new mapping directory
		mapDir = outPath + '/mapping'
		if not os.path.exists(mapDir):
			os.makedirs(mapDir)
		self.mapDir			  =  mapDir
		
		#setup the logger
		if logger is None:
			_logfileName  = outPath + '/BWA_processing.log'
			logger = utils._get_logger(_logfileName)
		self.logger = logger							
		
		#mapping related
		self.fastq_prefix			=   fastqUtils.get_guessed_fastq_prefix(os.path.basename(self.read1))
		self.bwa_map_log			 =   self.mapDir + '/' + self.fastq_prefix  + '.bwa.log'
		self.mappedBam			   =   self.mapDir + '/' + self.fastq_prefix  + '.bam'
		
		#create the reference index if needed
		final_refPath = creat_bwa_ref_index(reference,logger,outPath=mapDir)
		self.reference = final_refPath
		
	
	
	def map_reads(self):
		temp_map_reads_with_BWA(self)
	



	
	
	
	
	

def creat_bwa_ref_index(reference,logger,outPath=os.getcwd()):
	
	SUB = 'creat_bwa_ref_index'
	
	#check if the index already exists
	#in that case dont create one
	expected_bwa_index_file_suffixes = ['amb','ann','bwt','pac','sa']
	
	#flag
	index_creation_needed = False
	for suffix in expected_bwa_index_file_suffixes:
		if not os.path.exists(reference + '.' + suffix):
			if suffix == 'fai':
				#genome size file is missing create one
				genome_size_file = reference + '.fai'
				seqUtils.get_fasta_sequence_length(reference,outFile=genome_size_file)
			else:
				index_creation_needed = True
				break
			
	if index_creation_needed:
		#choose the right indexing method based on size of the fastq
		#greater than 2Gb genome choose choose bwtsw method else 'is' method
		size = seqUtils.get_total_bases_in_fasta(reference)
		print 'size is %s' % size
		if size > 2e9:
			index_algo = 'bwtsw'
		else:
			index_algo = 'is'
		#create a new dir for storing index
		index_dir = outPath + '/bwa_index/'
		if not os.path.exists(index_dir):
			os.makedirs(index_dir)
		new_reference_location = index_dir + os.path.basename(reference) 
		shutil.copy(reference,new_reference_location)
		#create the index
		args = ['bwa','index','-a',index_algo,new_reference_location]
		return_code = subprocess.check_call(args) 
		#create the fai file
		genome_size_file = new_reference_location + '.fai'
		seqUtils.get_fasta_sequence_length(new_reference_location,outFile=genome_size_file)
		logger.info('[%s]: created bwa reference index for %s' % (SUB,new_reference_location))
		return new_reference_location
	else:
		logger.info('[%s]: bwa reference index already exists for %s' % (SUB,reference))
		return reference
	
	
	

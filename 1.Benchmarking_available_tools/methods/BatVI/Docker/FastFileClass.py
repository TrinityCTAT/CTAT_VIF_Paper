#!/usr/bin/env python

import pandas as pd
import pysam 
import os, sys, re
import logging
import sys, time
import itertools
import gzip
logging.basicConfig(format='\n %(levelname)s : %(message)s', level=logging.DEBUG)
# logging.basicConfig(format='%(asctime)s %(message)s') 
'''
_Review terms_ :) 

Files:
    fastq : stores sequence records. stores an ID and a sequence. 
    FASTQ : contains a sequence of quality scores for each nucleotide


Python: 

    self            : self represents the instance of the class, used to access attributes and methods of the class. maps attributes with given arguments 
    init            : 'init' is a constructor, method called when object created allows for teh initialization of the attributes of a class 
    attributes      : ex)     self.hairColor = red
    methods         : ex)     def dyeHair(self): 
                            print(blue)

Notes:
source file using exec(open('processFA.py').read())
fa_file = faFile(input_file = "")
'''

def CHECK_reading(true_ids, read_names):
    if len(true_ids) == len(read_names):
        print(f"\t\t\t\tAll {len(read_names)} reads found!")
    else:
        m = f"\t\t\t\t NOT All reads found! {len(true_ids)} out of {len(read_names)} Found."
        logger.warning(m)

        print(set(read_names) - set(true_ids))
        exit()
    
def readME(infile, read_names):
    '''
    Read in the fastQ and iterate over it 
    '''
    total_reads = len(read_names)
    final_list = []
    while True:
        i = list(itertools.islice(infile, 4))
        if not i:
            break 
        ID = i[0][1:].split(" ")[0]
        if ID in read_names:
            final_list.append(i)

            # print Percent
            sys.stdout.write("\r"); sys.stdout.flush()
            sys.stdout.write(f"{len(final_list)/total_reads * 100}%...")
    return final_list

# Create class object faFile
class fastFile:
    # initialize object
    def __init__(self, input_file): # arguments to class instantiation 
        super(fastFile, self).__init__()
        self.input_file = input_file
        self.buffer_size = os.stat(input_file).st_blksize
        self.file = gzip.open(input_file, "rt") if re.search(".gz$", input_file) else open(input_file, "rt")
        self.sequence_ids = []
        self.sequences = {}
        self.num_seqs = 0
        self.output_name = ""


    def faReader(self):
        contig_id = None
        # Read in the file to memory 
        tmp = self.file.read()
        # Split the file by the > character
        ## this seperates the different sequences
        ## skip the first as it will be blank do to .split() 
        tmp = tmp.split(">")[1:]
        # Create the dictionary to hold the ids and sequences 
        dic = {}
        # Run through each sequence and remove the \n characters 
        ## Then add it to the dictionary 
        for i in tmp:
            idx = i.find("\n")
            dic[i[0:idx]] = i[idx:].replace("\n","")
        # Save the dictionary in the object as .sequences 
        self.sequences = dic
        # return the object
        return self 

    def fqReader(self):
        '''
        Read in the FastQ formated file 
        This function assumes the file is in standard 4 line FASTQ format 

        FastQ file format
            line1: @sequence_identifier 
            line2: Raw_Sequence
            line3: option
            line4: Quality values s
        '''

        contig_id = None
        # Read in the file 
        tmp = self.file.read().rstrip()
        # split all the lines by \n 
        tmp = tmp.split("\n")
        # Create the dictionary to hold the ids and sequences 
        dic = {}
        # Run through each sequence and remove the \n characters 
        ## Then add it to the dictionary 
        n = 4 # there are 4 lines in a single entry
        for i in list(range(0,len(tmp),n)):
            dic[tmp[i][1:]] = "\n".join(tmp[i+1:i+n])+"\n"
        # Save the dictionary in the object as .sequences 
        self.sequences = dic
        
        # return the object
        return self 

    def getSeqLengths(self):

        ## Index the fasta file
        # fai_filename = fasta_filename + ".fai"
        # if not os.path.exists(fai_filename):
        #     cmd = "samtools faidx " + fasta_filename
        #     subprocess.check_call(cmd, shell=True)

        seq_lengths = dict()
        with open(self.input_file) as fh:
            for line in fh:
                line = line.rstrip()
                vals = line.split("\t")
                acc, length = vals[0], int(vals[1])
                seq_lengths[acc] = length

        # get only chromosomes and not mitochondrial
        ref_chromosomes = list(filter(lambda x: re.match("chr", x), seq_lengths.keys()))

        return seq_lengths



    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # Subset  
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    def fqReader_subset(self, read_names):
        '''
        Subsets the given Fastq files to only include the given read_names.
        Returns the new fastq as a sting.

        This function assumes the file is in standard 4 line FASTQ format 
        FastQ file format
            line1: @sequence_identifier 
            line2: Raw_Sequence
            line3: option
            line4: Quality values s
        '''

        # variable to tell if the files are old format or not
        older_format = False

        # If files are gzipped
        if self.input_file.endswith(".gz"):
            with gzip.open(self.input_file, "rt") as infile:
                # check the first line 
                # Need to determine if these are older or newer formatted fastqs
                first_line = infile.readline().rstrip()
                first_id = first_line[1:].split(" ")[0]
                if ( (first_id.endswith("/1")) or (first_id.endswith("/2")) ):
                    logger.info("\t\tIdentified older formated Fastq, editing read names...")
                    older_format = True
            # Now iterate over file 
            with gzip.open(self.input_file, "rt") as infile:
                if older_format == True:
                    final_output = readMeOldFormat(infile,read_names)
                else:
                    final_output = readME(infile, read_names)

        else:
            with open(self.input_file, "r") as infile:
                # check the first line 
                first_line = infile.readline().rstrip()
                first_id = first_line[1:].split(" ")[0]
                if ( (first_id.endswith("/1")) or (first_id.endswith("/2")) ):
                    logger.info("\t\tIdentified older formated Fastq, editing read names...")
                    older_format = True

            with open(self.input_file, "r") as infile:
                if older_format == True:
                    final_output = readMeOldFormat(infile,read_names)
                else:
                    final_output = readME(infile, read_names)

        #~~~~~~~~~~~~~~~~~~~~~~~~~
        # CHECK
        #~~~~~~~~~~~~~~~~~~~~~~~~~
        CHECK_reading(final_output, read_names)
        # N_reads = len(final_output)
        # if N_reads == len(read_names):
        #     m = f"All Reads Found {N_reads}/{len(read_names)}"
        #     logger.info(m)
        # else:
        #     m = f"NOT All Reads Found {N_reads}/{len(read_names)}"
        #     logger.warning(m)
        #     quit()


        # Now join the list of lists together 
        final_output_join = list(itertools.chain.from_iterable(final_output))

        # now join the string in the list together 
        final_output_joined = "".join(final_output_join)
        

        return final_output_joined



    def faWriter(self, output = "Output_fa.fa"):
        '''
        Will create a new FA file for what is in self.sequences 
        '''
        # Set the constant 
        self.output_name = output

        # Create the file to write to
        output_file = open(self.output_name, "w")

        for i in self.sequences:
            tmp = ">{}\n{}\n".format(i, self.sequences[i])
            output_file.write(tmp)
        output_file.close()

    def fqWriter(self, output = "Output_fq.fa"):
        # Set the constant 
        self.output_name = output

        # Create the file to write to
        output_file = open(self.output_name, "w")

        for i in self.sequences:
            tmp = f"@{i}\n{self.sequences[i]}"
            output_file.write(tmp)
        output_file.close()




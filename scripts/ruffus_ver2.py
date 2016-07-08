#!/usr/bin/env python
import sys, os

from ruffus import *
import ruffus.cmdline as cmdline
from subprocess import check_call
parser = cmdline.get_argparse(description="Chela's Pipeline")

#                                                                                 .
#   Very flexible handling of input files                                         .
#                                                                                 .
#      input files can be specified flexibly as:                                  .
#                 --input a.fastq b.fastq                                         .
#                 --input a.fastq --input b.fastq                                 .
#                 --input *.fastq --input other/*.fastq                           .
#                 --input "*.fastq"                                               .
#                                                                                 .
#       The last form is expanded in the script and avoids limitations on command .
#           line lengths                                                          .
#                                                                                 .
parser.add_argument('-i', '--input', nargs='+', metavar="FILE", action="append", help = "Fastq files")

options = parser.parse_args()

#  standard python logger which can be synchronised across concurrent Ruffus tasks
logger, logger_mutex = cmdline.setup_logging ("Chela", options.log_file, options.verbose)

#                                                                                .
#   Useful code to turn input files into a flat list                             .
#                                                                                .
from glob import glob
original_data_files = [fn for grouped in options.input for glob_spec in grouped for fn in glob(glob_spec)] if options.input else []
if not original_data_files:
    original_data_files = [["C1W1_R1.fastq.gz", "C1W1_R2.fastq.gz"]]
    #raise Exception ("No matching files specified with --input.")
#   <<<----  pipelined functions go here

#start drmaa_session
from ruffus.drmaa_wrapper import run_job, error_drmaa_job
import drmaa
drmaa_session = drmaa.Session()
drmaa_session.initialize()
#_________________________________________________________________________________
#                                                                                .
#   Group together file pairs and make directories                               .
#   The first step trims and generates fastqc reports                            .
#_________________________________________________________________________________
@mkdir(original_data_files,
       #match input pattern
        formatter("([^/]+)R[12]_001.fastq.gz$"),
        # make qc directory
        "{path[0]}/qc")
@collate(original_data_files,
            # match file name up to the "R1.fastq.gz"
            formatter("([^/]+)R[12]_001.fastq.gz$"),
            # Create output parameter supplied to next task
            ["{path[0]}/{1[0]}R1_001_val_1.fq.gz",  # paired file 1
             "{path[0]}/{1[0]}R2_001_val_2.fq.gz"], # paired file 2
            # Extra parameters qc folder
            "{1[0]}R1_001.fastq.gz", "{1[0]}R2_001.fastq.gz","{path[0]}/qc",
            "{path[0]}",
            logger, logger_mutex)
def trim_fastq(input_files, output_files, basename1, basename2, qc_folder, output_folder , logger, logger_mutex):
  if len(input_files) != 2:
    print "length of input file: " + str(len(input_files))    
    raise Exception("One of read pairs %s missing" % (input_files,))
  #cmd = (" trim_galore --fastqc --paired {input_files[0]} {input_files[1]} "  )
  #check_call(cmd.format(**locals()))
  
  #print input_files
  #print output_files
  input1 = input_files[0]
  input2 = input_files[1]
  cmd = "cd $TMPDIR ; cp " +  input1 + " . ; cp " + input2 + " . ; trim_galore --fastqc --paired " + basename1 + " "  +  basename2 + " ; mv *_val_*.fq.gz " + output_folder + " ; mv *fastqc* " +  qc_folder + " ; mv *report* " + qc_folder     + " ; rm * "
  job_name = "trim_fastqc"
     
  try:
    stdout_res, stderr_ress = run_job(cmd,
                                      job_name,
                                      job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                      job_other_options    = "-S /bin/bash -V -l h_rt=01:00:00 -l mem=4G -l tmpfs=10G -wd /home/sejjctj/Scratch -j yes",
                                      job_environment      = { 'BASH_ENV' : '/home/sejjctj/.bashrc' } ,
                                      retain_job_scripts   = True,
                                      working_directory    = "/home/sejjctj/Scratch",
                                      drmaa_session        = drmaa_session,
                                      logger = logger )                                      
                                      
  except error_drmaa_job as err:
    raise Exception("\n".join(map(str,
                      ["Failed to run:",
                        cmd,
                        err,
                        stdout_res,
                        stderr_res])))
                                    
  with logger_mutex:
    logger.debug("trim_fastq worked")
#_______________________________________________________________________________________________________
# 
#              take trimmer output and align with hisat2
#_______________________________________________________________________________________________________
reads_list = []
@merge(trim_fastq,reads_list)
def merge_trim_output(input_files,reads_list):
  reads1 = []
  reads2 = []
  for item in input_files:
    reads1.append(item[0])
    reads2.append(item[1])
  reads_list.append(reads1)
  reads_list.append(reads2)
  # simply flattens the read lists to create two lists for hisat

#"([^/]+)R[12]_001.fastq.gz$"
@transform(merge_trim_output,
           formatter("([^/]+)L001_R[12]_001_val_[12].fq.gz$"), "{path[0]}/{basename[0]}.bam", "{path[0]}","{path[0]}/qc")
def hisat2(input_files, output_file, path, qc_folder):
  first_reads_list = []
  second_reads_list = []
  def get_file_name_list(file_list):
    temp_list = []
    temp_file_list = []
    for read_file in file_list:
      temp_list.append(os.path.basename(read_file))
    return temp_file_list
# above are getting the file name only
  first_reads_list = get_file_name_list(input_files[0])
  second_reads_list = get_file_name_list(input_files[1])
  # use get_file_name_list function to get list of filenames
  first_reads = ','.join(first_reads_list)
  second_reads = ','.join(second_reads_list)
  hisat_output = output_file.split('/')
  hisat_output = hisat_output[-1]
  job_script_directory = "/home/sejjctj/Scratch/test_dir"
  cmd = "cd $TMPDIR; mkdir reference; cp " + path + "/*val*fq.gz" + " . ; cp $HOME/Scratch/reference/grch38/* ./reference ; hisat2 -x ./reference/genome --known-splicesite-infile ./reference/grch38.splices.txt --dta-cufflinks -1 " + first_reads + " -2 " + second_reads + "  2> " + qc_folder  + "/" + hisat_output + ".log | samtools view -bS - -o " + hisat_output + " ; cp " + hisat_output + " " + path + " ; rm -r *"         
  
  job_name = "hisat"

  try:
    stdout_res, stderr_res = run_job(cmd,
                                      job_name,
                                      job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                      job_other_options    = "-S /bin/bash -V -l h_rt=04:00:00 -l mem=4G -l tmpfs=60G -wd /home/sejjctj/Scratch -j yes ",
                                      job_environment      = { 'BASH_ENV' : '/home/sejjctj/.bashrc' } ,
                                      retain_job_scripts   = True,  # retain job scripts for debuging, they go in Scratch/test_dir
                                      working_directory    = "/home/sejjctj/Scratch",
                                      drmaa_session        = drmaa_session,
                                      logger = logger )

  except error_drmaa_job as err:
    raise Exception("\n".join(map(str,
                      ["Failed to run:",
                        cmd,
                        err,
                        stdout_res,
                        stderr_res])))

  with logger_mutex:
    logger.debug("hisat worked")
#_______________________________________________________________________________________________________
# 
#              take trimmer output and run kallisto to get abundances
#_______________________________________________________________________________________________________

@transform(merge_trim_output,
           formatter("([^/]+)L001_R[12]_001_val_[12].fq.gz$"), 
           "{path[0]}/kallisto/abundance.tsv","{path[0]}", 
           "{path[0]}/kallisto")
def kallisto(input_files, output_file, path, kallisto_folder):
  # flatten list
  list_of_reads = []
  reads_list = [item for sublist in input_files for item in sublist]
  for filename in reads_list:
    list_of_reads.append(os.path.basename(filename))
  list_of_reads = ' '.join(list_of_reads)
  
  job_script_directory = "/home/sejjctj/Scratch/test_dir"
  
  cmd = "cd $TMPDIR; mkdir reference; cp " + path + "/*val*fq.gz" + " . ; cp $HOME/Scratch/reference/homo_transcripts.idx ./reference ; kallisto quant -b 100 -i ./reference/homo_transcripts.idx " + list_of_reads + " -o " + kallisto_folder +  " ; rm -r * ;"         
  print cmd
  job_name = "kallisto"
  
  try:
    stdout_res, stderr_res = run_job(cmd,
                                      job_name,
                                      job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                      job_other_options    = "-S /bin/bash -V -l h_rt=04:00:00 -l mem=4G -l tmpfs=60G -wd /home/sejjctj/Scratch -j yes ",
                                      job_environment      = { 'BASH_ENV' : '/home/sejjctj/.bashrc' } ,
                                      retain_job_scripts   = True,
                                      #touch_only           = True,
                                      #run_locally          = True,
                                      working_directory    = "/home/sejjctj/Scratch",
                                      drmaa_session        = drmaa_session,
                                      logger = logger )

  except error_drmaa_job as err:
    raise Exception("\n".join(map(str,
                      ["Failed to run:",
                        cmd,
                        err,
                        stdout_res,
                        stderr_res])))

  with logger_mutex:
    logger.debug("kallisto worked")
#_______________________________________________________________________________________________________
# 
#              samtools sort hisat output
#_______________________________________________________________________________________________________
  
@transform(hisat2,formatter("([^/]+)L001_R1_001_val_1.fq.bam$"), "{path[0]}/{basename[0]}.sorted.bam", "{basename[0]}", "{path[0]}")
def sort(input_file, output_file, prefix, path):
  filename = os.path.basename(input_file)
  print filename
  print output_file
  job_script_directory = "/home/sejjctj/Scratch/test_dir"
  
  cmd = "cd $TMPDIR; cp " + input_file + " . ; samtools sort " + filename + " " + prefix + ".sorted ; cp " + prefix + ".sorted.bam " + path + " ; rm -r * "         

  job_name = "sort"
   
  try:
    stdout_res, stderr_res = run_job(cmd,
                                      job_name,
                                      job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                      job_other_options    = "-S /bin/bash -V -l h_rt=04:00:00 -l mem=4G -l tmpfs=60G -wd /home/sejjctj/Scratch -j yes ",
                                      job_environment      = { 'BASH_ENV' : '/home/sejjctj/.bashrc' } ,
                                      retain_job_scripts   = True,
                                      working_directory    = "/home/sejjctj/Scratch",
                                      drmaa_session        = drmaa_session,
                                      logger = logger )

  except error_drmaa_job as err:
    raise Exception("\n".join(map(str,
                      ["Failed to run:",
                        cmd,
                        err,
                        stdout_res,
                        stderr_res])))

  with logger_mutex:
    logger.debug("sort worked")
#_______________________________________________________________________________________________________
# 
#              cufflinks generate gtf files from sorted hisat output
#_______________________________________________________________________________________________________

@transform(sort,formatter("([^/]+)L001_R1_001_val_1.fq.sorted.bam$"), "{path[0]}/cufflinks/transcripts.gtf", "{path[0]}","{path[0]}/qc")
def cufflinks(input_file, output_file, path, qc_path):
  filename = os.path.basename(input_file)
  job_script_directory = "/home/sejjctj/Scratch/test_dir"
  cmd = "cd $TMPDIR; mkdir reference; cp " + input_file + " . ; cp $HOME/Scratch/reference/grch38/Homo_sapiens.GRCh38* ./reference  ; cufflinks -p 4 -g ./reference/Homo_sapiens.GRCh38.76.gtf -b ./reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa " + filename + " -o " + path + "/cufflinks/ 2> " + qc_path + "/cufflinks.log ; rm -r * ;"
  job_name = "cufflinks"

  try:
    stdout_res, stderr_res = run_job(cmd,
                                      job_name,
                                      job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                      job_other_options    = "-S /bin/bash -V -l h_rt=06:00:00 -w n -l mem=4G -l tmpfs=60G -pe smp 4 -wd /home/sejjctj/Scratch -j yes ",
                                      job_environment      = { 'BASH_ENV' : '/home/sejjctj/.bashrc' } ,
                                      retain_job_scripts   = True,
                                      working_directory    = "/home/sejjctj/Scratch",
                                      drmaa_session        = drmaa_session,
                                      logger = logger )

  except error_drmaa_job as err:
    raise Exception("\n".join(map(str,
                      ["Failed to run:",
                        cmd,
                        err,
                        stdout_res,
                        stderr_res])))

  with logger_mutex:
    logger.debug("cufflinks worked")

if __name__ == '__main__':
  cmdline.run (options)
  drmaa_session.exit()

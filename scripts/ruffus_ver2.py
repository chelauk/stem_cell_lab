#!/usr/bin/env python
import glob
import sys, os, fnmatch
import re
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
#
#    Add argument for where assembly text file required for cuffmerge is
#
parser.add_argument('--cuffdiff_file', metavar="FILE", help = "cuffdiff comparison instructions")
parser.add_argument('--basedir', metavar="DIR", help = "base directory")

options = parser.parse_args()
cuffdiff_file = options.cuffdiff_file
basedir=options.basedir
# working directory
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
        formatter("([^/]+)$"),
        # make qc directory
        ["{subpath[0][1]}/qc",
        # make trimmed_directory)
        "{subpath[0][1]}/fastq_trimmed",
        "{subpath[0][1]}/bam",
        "{subpath[0][1]}/kallisto",
        "{subpath[0][1]}/cufflinks",
        "{subpath[0][3]}/cuffmerge_out"])
@collate(original_data_files,
        formatter("([^/]+)R[12]_001.fastq.gz$"),
              ["{subpath[0][1]}/fastq_trimmed/{1[0]}R1_001_val_1.fq.gz",  # paired file 1
              "{subpath[0][1]}/fastq_trimmed/{1[0]}R2_001_val_2.fq.gz"],  # paired file 1
              "{1[0]}R1_001.fastq.gz", "{1[0]}R2_001.fastq.gz","{subpath[0][1]}/qc",
             "{subpath[0][1]}/fastq_trimmed", logger, logger_mutex)
def trim_fastq(input_files, output_files, basename1, basename2, qc_folder, output_folder ,logger, logger_mutex):
  if len(input_files) != 2:
    print "length of input file: " + str(len(input_files))    
    raise Exception("One of read pairs %s missing" % (input_files,))
  input1 = input_files[0]
  input2 = input_files[1]
  cmd = "cd $TMPDIR ; cp " +  input1 + " . ; cp " + input2 + " . ; trim_galore --fastqc --paired " + basename1 + " "  +  basename2 + " 2> " + qc_folder + "/trim_galore.log ; mv *_val_*.fq.gz " + output_folder + " ; mv *fastqc* " +  qc_folder + " ; mv *report* " + qc_folder     + " ; rm * "
  job_name = "trim_fastqc"
  #print output_files

  try:
    
    stdout_res, stderr_res = "",""
    stdout_res, stderr_res = run_job(cmd,
                                      job_name,
                                      job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                      job_other_options    = "-S /bin/bash -V -l h_rt=01:00:00 -l mem=4G -l tmpfs=10G -wd /home/sejjctj/Scratch -j yes",
                                      job_environment      = { 'BASH_ENV' : '/home/sejjctj/.bashrc' } ,
                                      retain_job_scripts   = True,
                                      working_directory    = "/home/sejjctj/Scratch",
                                      drmaa_session        = drmaa_session,
                                      logger = logger )                                      
   #                                   
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
  #print input_files
  #reads1 = []
  #reads2 = []
  #for item in input_files:
  #  reads1.append(item[0])
  #  reads2.append(item[1])
  #reads_list.append(reads1)
  
  #reads_list.append(reads2)
  #print reads1
  # simply flattens the read lists to create two lists for hisat

#"([^/]+)R[12]_001.fastq.gz$"
hisat_input_list = {}
hisat_input_list2 = []



@collate(trim_fastq, formatter("^.+/(?P<FILENAME>.*)L00[123456]_R[12]_001_val_[12].fq.gz$"), "{subpath[0][1]}/bam/{FILENAME[0]}L001_R1_001_val_1.fq.sorted.bam", "{path[0]}","{subpath[0][1]}/bam","{subpath[0][1]}/qc")
def hisat2(input_files, out_file, path, outpath,qc_folder):
  print "hisat2 input_list: " + str(input_files)
  flat_list = [item for sublist in input_files for item in sublist]
  #print flat_list
  print "hisat2 output file: " + out_file
  print "qc folder: " + qc_folder
  print "path to copy fastqs from:" + path
  
  first_reads = []
  second_reads =[]
  for i in flat_list:
    #print os.path.basename(i)
    if re.search('val_1', os.path.basename(i)):
      first_reads.append(os.path.basename(i))
    elif re.search('val_2', os.path.basename(i)):
      second_reads.append(os.path.basename(i))
  first_reads = ','.join(first_reads)
  second_reads = ','.join(second_reads)
  print "first reads: " + first_reads
#  reads_list = [item for sublist in input_files for item in sublist]
#  for read in reads_list:
#    print read
#  print output_file
#  first_reads_list = []
#  second_reads_list = []
#  def get_file_name_list(file_list):
#    temp_list = [[],[]]
#    for read_file in file_list:
#      temp_list[0].append(os.path.basename(read_file[0]))
#      temp_list[1].append(os.path.basename(read_file[1]))
#    return temp_list
## above for getting the file name only
#  reads_list = get_file_name_list(input_files)
#  first_reads_list = reads_list[0]
#  second_reads_list = reads_list[1]
#  # use get_file_name_list function to get list of filenames
#  first_reads = ','.join(first_reads_list)
#  second_reads = ','.join(second_reads_list)
#  #print first_reads
  hisat_output = out_file.split('/')
  hisat_output = hisat_output[-1]
#  #sort_output = hisat_output[:-4] + ".sorted.bam"
#  job_script_directory = "/home/sejjctj/Scratch/test_dir"


  cmd = "cd $TMPDIR; mkdir reference; cp " + path + "/*val*fq.gz" + " . ; cp $HOME/Scratch/reference/grch38_snp_tran/genome* ./reference ; hisat2 -p 4 -x ./reference/genome_snp_tran  --dta-cufflinks -1 " + first_reads + " -2 " + second_reads + "  2>" + qc_folder  + "/hisat.log | samtools view -bS - -o temp.bam ; samtools sort -@ 4 temp.bam  " + hisat_output[:-4] + " ; cp " + hisat_output + " " + outpath + " ; rm -r *"         
  print cmd 
  #print output_file
  #print hisat_output 
  job_name = "hisat"

  try:
    stdout_res, stderr_res = run_job(cmd,
                                     job_name,
                                     job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                      job_other_options    = "-S /bin/bash -V -l h_rt=07:00:00 -w n -l mem=4G -l tmpfs=60G -pe smp 4 -wd /home/sejjctj/Scratch -j yes ",
                                      job_environment      = { 'BASH_ENV' : '/home/sejjctj/.bashrc' } ,
  #                                    touch_only           = True,
  #                                    run_locally          = True,
                                      retain_job_scripts   = True,  # retain job scripts for debuging, they go in Scratch/test_dir
                                      working_directory    = "/home/sejjctj/Scratch",
                                      drmaa_session        = drmaa_session,
                                      logger = logger )
#
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


@collate(trim_fastq,formatter("([^/]+)L00[1234]_R[12]_001_val_1.fq.gz$"),"{subpath[0][1]}/kallisto/abundance.tsv","{path[0]}","{subpath[0][1]}/kallisto")
def kallisto(input_files, output_file, path, kallisto_folder):
  # flatten list
  list_of_reads = []
  reads_list = [item for sublist in input_files for item in sublist]
  for filename in reads_list:
    list_of_reads.append(os.path.basename(filename))
  list_of_reads = ' '.join(list_of_reads)
  print "kallisto list of reads: " + list_of_reads 
  job_script_directory = "/home/sejjctj/Scratch/test_dir"
  
  cmd = "cd $TMPDIR; mkdir reference; cp " + path + "/*val*fq.gz" + " . ; cp $HOME/Scratch/reference/homo_transcripts.idx ./reference ; kallisto quant -b 100 -t 4 -i ./reference/homo_transcripts.idx " + list_of_reads + " -o " + kallisto_folder +  " ; rm -r * ;"         
  #print cmd
  #print output_file
  job_name = "kallisto"
  
  try:
    stdout_res, stderr_res = run_job(cmd,
                                      job_name,
                                      job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                      job_other_options    = "-S /bin/bash -V -l h_rt=04:00:00 -l mem=4G -w n -pe smp 4 -l tmpfs=60G -wd /home/sejjctj/Scratch -j yes ",
                                      job_environment      = { 'BASH_ENV' : '/home/sejjctj/.bashrc' } ,
                                      retain_job_scripts   = True,
  #                                    touch_only           = True,
  #                                    run_locally          = True,
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
#              cufflinks generate gtf files from sorted hisat output
#_______________________________________________________________________________________________________


@transform(hisat2,formatter("^.+/(?P<FILENAME>.*)bam$"), "{subpath[0][1]}/cufflinks/transcripts.gtf","{subpath[0][1]}/cufflinks","{subpath[0][1]}/qc")
def cufflinks(input_file, output_file,  path, qc_path):
  
  print "cufflinks input: " + input_file
  print "cufflinks output: " + output_files
  
  job_script_directory = "/home/sejjctj/Scratch/test_dir"
  cmd = "cd $TMPDIR; mkdir reference; cp " + input_file + " . ; cp $HOME/Scratch/reference/grch38/Homo_sapiens.GRCh38.dna.primary_assembly.fa* ./reference  ; cp $HOME/Scratch/reference/grch38/Homo_sapiens.GRCh38.85.gtf ./reference ; cp $HOME/Scratch/reference/grch38/ribosomal_mito_mask.gtf ./reference ; cufflinks -p 4 -g ./reference/Homo_sapiens.GRCh38.85.gtf -b ./reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa --mask-file ./reference/ribosomal_mito_mask.gtf " + input_file + " -o " + path + " 2>" + qc_path + "/cufflinks.log ; rm -r * ;"
  job_name = "cufflinks"
  print cmd
  
  try:
    stdout_res, stderr_res = "",""
    stdout_res, stderr_res = run_job(cmd,
                                      job_name,
                                      job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                      job_other_options    = "-S /bin/bash -V -l h_rt=07:00:00 -w n -l mem=4G -l tmpfs=60G -pe smp 4 -wd /home/sejjctj/Scratch -j yes ",
                                      job_environment      = { 'BASH_ENV' : '/home/sejjctj/.bashrc' } ,
   #                                  touch_only           = True,
   #                                  run_locally          = True,
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

#_______________________________________________________________________________________________________
# 
#              cuffmerge, create assembly text for cuffmerge
#_______________________________________________________________________________________________________

@merge(cufflinks, basedir + "assembly_list.txt" )
def create_cuffmerge_input(input_files, output):
  #print "input files :" + str(input_files)
  #print output
  with open(output,'w') as f:
    for item in input_files:
      f.write("%s\n" % item)



@transform(create_cuffmerge_input, formatter("([^/]+).$"), "{path[0]}/cuffmerge_out/merged.gtf", "{path[0]}/cuffmerge_out/")
def cuffmerge(input_file, output_file, output_path):
  print "cuffmerge input: " + input_file
  print "cuffmerge output: " + output_file

  cmd = "cd $TMPDIR; cp " + input_file + " . ; mkdir ./reference; cp $HOME/Scratch/reference/grch38/Homo_sapiens.GRCh38.dna.primary_assembly.fa* ./reference ; cp $HOME/Scratch/reference/grch38/Homo_sapiens.GRCh38.85.gtf ./reference; cuffmerge -p 4  -g ./reference/Homo_sapiens.GRCh38.85.gtf -s ./reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa " + input_file + "  2>>" + output_path + "/cuffmerge.log ; mv ./merged_asm/* " + output_path + "; rm -r * ;" 
  print cmdi
  job_name = "cuffmerge"
  try:
    stdout_res, stderr_res = "",""
    stdout_res, stderr_res = run_job(cmd,
                                     job_name,
                                     job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                     job_other_options    = "-S /bin/bash -V -l h_rt=07:00:00 -w n -l mem=4G -l tmpfs=60G -pe smp 4 -wd /home/sejjctj/Scratch -j yes ",
                                     job_environment      = { 'BASH_ENV' : '/home/sejjctj/.bashrc' } ,
     #                                touch_only           = True,
     #                                run_locally          = True,
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
    logger.debug("cuffmerge worked")


#______________________________________________________________
#
#   creating pairs for comparison with cuffdiff
#______________________________________________________________


import glob
@follows(cuffmerge)
@transform(basedir + "cuffdiff_samples.txt", formatter("([^/]+)"), "{path[0]}/cuffdiff/test1/genes.fpkm_tracking", "{path[0]}/cuffdiff/")
def prepare_cuffdiff(input_files, output_file, output_path):
  
  with open(input_files, 'r') as f:
    for line in f:
      line = line.rstrip()
      line = line.split(',')
      #print line
      path = output_path + line[0]
      #print path
      try:
        os.mkdir( path , 0755)
      except OSError, e:
        if e.errno == 17:
          pass
      test1 = glob.glob(basedir + line[1] + "/replicate*/bam/*")
      test2 = glob.glob(basedir + line[2] + "/replicate*/bam/*")
      test1 = ','.join(test1)
      test2 = ','.join(test2)
      
      cmd = "cd $TMPDIR; mkdir reference; cp $HOME/Scratch/reference/grch38/Homo_sapiens.GRCh38.dna.primary_assembly.fa* ./reference; cuffdiff  -p 4 -b ./reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa " + basedir + "cuffmerge_out/merged.gtf -o " + path + " " + test1 + " " + test2 + " 2>>" + basedir + "/cuffdiff.log ; rm -r * ;"
      print cmd

  job_name = "cuffdiff"
  try:
    stdout_res, stderr_res = "",""
    stdout_res, stderr_res = run_job(cmd,
                                     job_name,
                                     job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                     job_other_options    = "-S /bin/bash -V -l h_rt=12:00:00 -w n -l mem=4G -l tmpfs=60G -pe smp 4 -wd /home/sejjctj/Scratch -j yes ",
                                     job_environment      = { 'BASH_ENV' : '/home/sejjctj/.bashrc' } ,
     #                                touch_only           = True,
     #                                run_locally          = True,
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
    logger.debug("cuffdiff worked")


#______________________________________________________________

if __name__ == '__main__':
  cmdline.run (options, multithread = options.jobs)
  drmaa_session.exit()

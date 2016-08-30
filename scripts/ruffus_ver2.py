#!/usr/bin/env python
import glob
import sys, os, fnmatch
import re
from ruffus import *
import ruffus.cmdline as cmdline
from subprocess import check_call
parser = cmdline.get_argparse(description="Chela's Pipeline")

#parser.add_argument('-i', '--input', nargs='+', metavar="FILE", action="append", help = "Fastq files")
parser.add_argument('-i', '--input', metavar="FILE", help = "Fastq files")
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
#print "input " + options.input
#print "basedir " + basedir
#original_data_files = [fn for grouped in options.input for glob_spec in grouped for fn in glob(glob_spec)] if options.input else []
#print options.input + "\n\n"

import glob
input_files=[]
def get_data_files(x):
  with open(x, 'r') as f:
    content = [line.rstrip('\n') for line in f] 
    for line in content:
      line = line.rstrip()
      input_files.append(glob.glob(basedir + line + "/replicate*/fastq_raw/*fastq.gz"))  
    return input_files
 
original_data_files=get_data_files(options.input)
original_data_files = [item for sublist in original_data_files for item in sublist]

if not original_data_files:
    raise Exception ("No matching files specified with --input.")

#start drmaa_session

from ruffus.drmaa_wrapper import run_job, error_drmaa_job
import drmaa
drmaa_session = drmaa.Session()
drmaa_session.initialize()

#   <<<----  pipelined functions go here
#_________________________________________________________________________________
#                                                                                .
#   Group together file pairs and make directories                               .
#   The first step trims pairs of fastqs and generates fastqc reports            .
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
        "{subpath[0][3]}/cuffmerge_out",
        "{subpath[0][3]}/cuffdiff",
        "{subpath[0][1]}/cuffquant"])
@collate(original_data_files,
        # input formatter to provide read pairs
        formatter("([^/]+)R[12](.+)gz$"),
        # create output parameter to be supplied to next task
        ["{subpath[0][1]}/fastq_trimmed/{1[0]}R1_001_val_1.fq.gz",
         "{subpath[0][1]}/fastq_trimmed/{1[0]}R2_001_val_2.fq.gz"],
        ["{1[0]}R1_001.fastq.gz", "{1[0]}R2_001.fastq.gz"],     # basenames for trim_galore
        "{subpath[0][1]}/qc",               # qc folder
        "{subpath[0][1]}/fastq_trimmed",    # trimmed_folder
        logger, logger_mutex)
def trim_fastq(input_files, output_files, basenames, qc_folder, output_folder ,logger, logger_mutex):
  if len(input_files) != 2:
    #print "length of input file: " + str(len(input_files))    
    raise Exception("One of read pairs %s missing" % (input_files,))
  cmd = (" cd $TMPDIR ; "
         " cp {input_files[0]} . ;"
         " cp {input_files[1]} . ;"
         " trim_galore --fastqc --paired {basenames[0]}  {basenames[1]}  2> {qc_folder}/trim_galore.log ; "
         " mv *val_*.fq.gz  {output_folder} ; "
         " mv *fastqc*  {qc_folder} ; "
         " mv *report* {qc_folder}; rm * ; " )
  
  job_name = "trim_fastqc"
  ## formats the cmd input to get the variables in the {}
  cmd = cmd.format(**locals())
  
  try:
    
    stdout_res, stderr_res = "",""
    stdout_res, stderr_res = run_job(cmd,
                                      job_name,
                                      job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                      job_other_options    = "-S /bin/bash -V -l h_rt=05:00:00 -l mem=4G -l tmpfs=60G -wd /home/sejjctj/Scratch -j yes",
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


@collate(trim_fastq, formatter("([^/]+)L00[1234]_R[12]_001_val_[12].fq.gz$"), 
                              "{subpath[0][1]}/bam/{1[0]}R1_001_val_1.fq.sorted.bam", 
                              "{path[0]}","{subpath[0][1]}/bam",
                              "{subpath[0][1]}/qc")
def hisat2(input_files, out_file, path, outpath,qc_folder):
  flat_list = [item for sublist in input_files for item in sublist]  
  first_reads = []
  second_reads =[]
  for i in flat_list:
    if re.search('val_1', os.path.basename(i)):
      first_reads.append(os.path.basename(i))
    elif re.search('val_2', os.path.basename(i)):
      second_reads.append(os.path.basename(i))
  first_reads = ','.join(first_reads)
  second_reads = ','.join(second_reads)
  hisat_output = out_file.split('/')
  hisat_output = hisat_output[-1]

  cmd = ( " cd $TMPDIR ; " 
          " mkdir reference ; "
          " cp  {path}"  + "/*val*fq.gz" + " . ; "
          " cp $HOME/Scratch/reference/grch38_snp_tran/genome* ./reference ; "
          " hisat2 -p 4 -x ./reference/genome_snp_tran  --dta-cufflinks "
          " --novel-splicesite-outfile ./novel_splice.txt "
          " --novel-splicesite-infile ./novel_splice.txt "
          " -1 " + first_reads + " -2 " + second_reads + "  2> {qc_folder}/hisat.log | samtools view -bS - -o temp.bam ; "
          " samtools sort -@ 4 temp.bam  " + hisat_output[:-4] + " ; "
          " cp " + hisat_output + " {outpath} ; "
          " cp novel_splice.txt {outpath} ; "
          " rm -r * ; ")
  
  cmd = cmd.format(**locals())         

  try:
    stdout_res, stderr_res = "",""
    stdout_res, stderr_res = run_job(cmd,
                                     job_name = "hisat",
                                     job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                     job_other_options    = "-S /bin/bash -V -l h_rt=07:00:00 -w n -l mem=4G -l tmpfs=60G -pe smp 4 -wd /home/sejjctj/Scratch -j yes ",
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


@collate(trim_fastq,formatter("([^/]+)L00[1234]_R[12]_001_val_[12].fq.gz$"),"{subpath[0][1]}/kallisto/abundance.tsv","{path[0]}","{subpath[0][1]}/kallisto")
def kallisto(input_files, output_file, path, kallisto_folder):
  list_of_reads = []
  reads_list = [item for sublist in input_files for item in sublist]
  for filename in reads_list:
    list_of_reads.append(os.path.basename(filename))
  list_of_reads = ' '.join(list_of_reads)
  job_script_directory = "/home/sejjctj/Scratch/test_dir"
  
  cmd = (" cd $TMPDIR; mkdir reference; "
         " cp {path}/*val*fq.gz   . ; "
         " cp $HOME/Scratch/reference/hg38_ver84_transcripts.idx ./reference ; "
         " kallisto quant -b 100 -t 4 -i ./reference/hg38_ver84_transcripts.idx " + list_of_reads + " -o {kallisto_folder} ; "
         " rm -r * ;")         
  
  cmd = cmd.format(**locals())         
  
  try:
    stdout_res, stderr_res = "",""
    stdout_res, stderr_res = run_job(cmd,
                                     job_name             = "kallisto",
                                     job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                     job_other_options    = "-S /bin/bash -V -l h_rt=04:00:00 -l mem=4G -w n -pe smp 4 -l tmpfs=60G -wd /home/sejjctj/Scratch -j yes ",
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
    logger.debug("kallisto worked")


#_______________________________________________________________________________________________________
# 
#              cufflinks generate gtf files from sorted hisat output
#_______________________________________________________________________________________________________


@transform(hisat2,formatter("([^/]+)bam$"), "{subpath[0][1]}/cufflinks/transcripts.gtf","{1[0]}bam", "{subpath[0][1]}/cufflinks","{subpath[0][1]}/qc")
def cufflinks(input_file, output_file, cuff_input, path, qc_path):
  
  job_script_directory = "/home/sejjctj/Scratch/test_dir"
  cmd = ( " cd $TMPDIR ; "
          " mkdir reference ; "
          " cp {input_file} . ; "
          " cp $HOME/Scratch/reference/grch38/Homo_sapiens.GRCh38.dna.primary_assembly.fa* ./reference  ; "
          " cp $HOME/Scratch/reference/grch38/Hs.GRCh38.84.exon.gtf ./reference ; "
          " cp $HOME/Scratch/reference/grch38/ribosomal_mito_mask.gtf ./reference ; "
          " cufflinks -q -u --no-update-check -p 4 -g ./reference/Hs.GRCh38.84.exon.gtf "
          " -b ./reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa "
          " --mask-file ./reference/ribosomal_mito_mask.gtf {cuff_input} -o  {path}  2>{qc_path}/cufflinks.log ; "
          " rm -r * ; " )
  cmd = cmd.format(**locals())
  print cmd  
  try:
    stdout_res, stderr_res = "",""
    stdout_res, stderr_res = run_job(cmd,
                                     job_name = "cufflinks",
                                     job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                     job_other_options    = "-S /bin/bash -V -l h_rt=30:00:00 -w n -l mem=4G -l tmpfs=60G -pe smp 4 -wd /home/sejjctj/Scratch -j yes ",
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


#_______________________________________________________________________________________________________
# 
#              cuffmerge, create assembly text for cuffmerge
#_______________________________________________________________________________________________________


# creates a templated for cuffmerge
@merge(cufflinks, "assembly_list.txt" )
def create_cuffmerge_input(input_files, output):
  with open( output, "w" ) as f:
    for item in input_files:
      f.write("%s\n" % item)


@transform(create_cuffmerge_input, formatter("([^/]+)txt"), "{path[0]}/cuffmerge_out/merged.gtf", "{path[0]}/cuffmerge_out/")
def cuffmerge(input_file, output_file, output_path):
  cmd = (" cd $TMPDIR; cp " + basedir + "assembly_list.txt . ; "
         " mkdir ./reference ; "
         " cp $HOME/Scratch/reference/grch38/Homo_sapiens.GRCh38.dna.primary_assembly.fa* ./reference ; "
         " cp $HOME/Scratch/reference/grch38/Hs.GRCh38.84.exon.gtf ./reference ; "
         " cuffmerge -q --no-update-check -p 4  -g ./reference/Hs.GRCh38.84.exon.gtf "
         " -s ./reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa assembly_list.txt 2>>{output_path}/cuffmerge.log ; "
         " mv ./merged_asm/* {output_path} ; "
         " rm -r * ; ") 
  
  cmd = cmd.format(**locals())
  
  try:
    stdout_res, stderr_res = "",""
    stdout_res, stderr_res = run_job(cmd,
                                     job_name = "cuffmerge",
                                     job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                     job_other_options    = "-S /bin/bash -V -l h_rt=12:00:00 -w n -l mem=4G -l tmpfs=60G -pe smp 4 -wd /home/sejjctj/Scratch -j yes ",
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
    logger.debug("cuffmerge worked")


@follows(cuffmerge)
@transform(hisat2,formatter("([^/]+)fq.sorted.bam$"), "{subpath[0][1]}/cuffquant/abundances.cxb", "{subpath[0][1]}/cuffquant")
def cuffquant(input_file, output_file, out_path):
  job_script_directory = "/home/sejjctj/Scratch/test_dir"
  cmd = ( " cd $TMPDIR ; "
          " mkdir reference; "
          " cp $HOME/Scratch/reference/grch38/Homo_sapiens.GRCh38.dna.primary_assembly.fa* ./reference  ; "
          " cp $HOME/Scratch/reference/grch38/ribosomal_mito_mask.gtf ./reference ; "
          " cuffquant -q --no-update-check -p 4 --mask-file "
          " ./reference/ribosomal_mito_mask.gtf -b ./reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa -u " + basedir + "cuffmerge_out/merged.gtf " 
          " {input_file}  -o {out_path}  2>{out_path}/cuffquant.log ; "
          " rm -r * ; ")
  
  cmd = cmd.format(**locals())
   
  try:
    stdout_res, stderr_res = "",""
    stdout_res, stderr_res = run_job(cmd,
                                     job_name             = "cuffquant",
                                     job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                     job_other_options    = "-S /bin/bash -V -l h_rt=12:00:00 -w n -l mem=4G -l tmpfs=60G -pe smp 4 -wd /home/sejjctj/Scratch -j yes ",
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
    logger.debug("cuffquant worked")
#______________________________________________________________
#
#   creating pairs for comparison with cuffdiff
#______________________________________________________________

@follows(cuffquant)
@transform(basedir + "cuffdiff_samples.txt", formatter("([^/]+)"), "{path[0]}/cuffdiff/cuffdiff.done", "{path[0]}/cuffdiff/")
def cuffdiff(input_files, output_file, output_path):
  with open(input_files, 'r') as f:
    for line in f:
      line = line.rstrip()
      line = line.split(',')
      print line
      path = output_path + line[0]
      print "directory path: " + path
      try:
        os.makedirs( path , 0755)
      except OSError, e:
        if e.errno == 17:
          pass
      test1 = glob.glob(basedir + line[1] + "/replicate*/cuffquant/*cxb")
      test2 = glob.glob(basedir + line[2] + "/replicate*/cuffquant/*cxb")
      test1 = ','.join(test1)
      test2 = ','.join(test2)
      print test1
      print test2
      cmd = ( " cd $TMPDIR ; "
              " mkdir reference ; "
              " cp $HOME/Scratch/reference/grch38/Homo_sapiens.GRCh38.dna.primary_assembly.fa* ./reference ; "
              " cuffdiff -q --no-update-check "
              " -p 4 -b ./reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa -L " + line[1] + ","+ line[2] + " " + basedir + "cuffmerge_out/merged.gtf "
              " -o " + path + " " + test1 + " " + test2 + " 2>>" + path + "/cuffdiff.log ; "
              " rm -r * ; " )
    cmd =cmd.format(**locals())  
    try:
      stdout_res, stderr_res = "",""
      stdout_res, stderr_res = run_job(cmd,
                                       job_name = "cuffdiff",
                                       job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                       job_other_options    = "-S /bin/bash -V -l h_rt=12:00:00 -w n -l mem=4G -l tmpfs=60G -pe smp 4 -wd /home/sejjctj/Scratch -j yes ",
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
  
  
    with open(output_file, 'w') as f:
      f.write("cuffdiff.done")
    
    with logger_mutex:
      logger.debug("cuffdiff worked")
  
#______________________________________________________________

if __name__ == '__main__':
  cmdline.run (options, multithread = options.jobs)
  drmaa_session.exit()


#!/usr/bin/env python
from glob import glob
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
parser.add_argument('--basedir', metavar="DIR", help = "base directory")

options = parser.parse_args()
basedir=options.basedir

print "Basedir: " + basedir
# working directory
#  standard python logger which can be synchronised across concurrent Ruffus tasks
logger, logger_mutex = cmdline.setup_logging ("Chela", options.log_file, options.verbose)

input_files=[]
files_list=options.input
print files_list
if not files_list:
  raise Exception ("No matching files specified with --input.")

with open(files_list, 'r') as f:
    content = [line.decode('utf-8').rstrip('\n') for line in f] 
    for line in content:
    #print(line)
        line = line.rstrip()
        #print(basedir + line)
        #print glob.glob(basedir + line + "/replicate*/fastq_raw/*gz")
        input_files.append(glob.glob(basedir + line + "/replicate*/fastq_raw/*gz"))  
input_files = [item for sublist in input_files for item in sublist]  
print input_files
#start drmaa_session
from ruffus.drmaa_wrapper import run_job, error_drmaa_job
import drmaa
drmaa_session = drmaa.Session()
drmaa_session.initialize()

#   <<<----  pipelined functions go here
#_________________________________________________________________________________
#                                                                                .
#   The first step trim fastqs and generates fastqc reports                      .
#_________________________________________________________________________________

@mkdir(input_files,
       #match input pattern
        formatter("([^/]+)$"),
        # make qc directory
        ["{subpath[0][1]}/qc",
        # make trimmed_directory)
        "{subpath[0][1]}/fastq_trimmed",
        "{subpath[0][1]}/bam"])
@transform(input_files,
        # input formatter to provide reads
        formatter("(?P<basedir>[/.].+)/(?P<sample>[a-zA-Z0-9_\-\.]+)/(?P<replicate>replicate_[0-9])/(?P<fastq_raw_dir>fastq_raw)/(?P<prefix>[a-zA-Z0-9_\-\.]+)fastq.gz$"),
        # create output parameter to be supplied to next task
	      "{basedir[0]}/{sample[0]}/{replicate[0]}/fastq_trimmed/{prefix[0]}_trimmed.fq.gz",
        "{basedir[0]}/{sample[0]}/{replicate[0]}/qc",               # qc folder
        "{basedir[0]}/{sample[0]}/{replicate[0]}/fastq_trimmed",    # trimmed_folder
        logger, logger_mutex)
def trim_fastq(input_file, output_files, basenames, qc_folder, output_folder ,logger, logger_mutex):
    raw_fastq=os.path.basename(input_file)
    cmd = (" cd $TMPDIR ; "
         " cp {input_file} . ;"
         " trim_galore --fastqc  {raw_fastq} 2> {qc_folder}/trim_galore.log ; "
         " mv *val_*.fq.gz  {output_folder} ; "
         " mv *fastqc*  {qc_folder} ; "
         " mv *report* {qc_folder}; rm * ; " )
  
    job_name = "trim_fastqc"
  ## formats the cmd input to get the variables in the {}
    cmd = cmd.format(**locals())
    #print(cmd)  
    try:
    
      stdout_res, stderr_res = "",""
      stdout_res, stderr_res = run_job(cmd,
                                      job_name,
                                      job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                      job_other_options    = "-w n -S /bin/bash -V -l h_rt=05:00:00 -l mem=4G -l tmpfs=60G -wd /home/sejjctj/Scratch -j yes",
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

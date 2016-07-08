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
                                      job_other_options    = "-S /bin/bash -V -l h_rt=01:00:00 -l mem=4G -l tmpfs=10G -pe smp 4 -wd /home/sejjctj/Scratch -j yes",
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
    logger.debug(trim_fastq worked")


if __name__ == '__main__':
  cmdline.run (options)
  drmaa_session.exit()

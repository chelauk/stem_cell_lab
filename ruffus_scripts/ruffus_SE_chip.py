
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
        formatter("(?P<basedir>[/.].+)/(?P<sample>[a-zA-Z0-9_\-\.]+)/(?P<replicate>replicate_[0-9])/(?P<fastq_raw_dir>fastq_raw)/(?P<prefix>[a-zA-Z0-9_\-]+).fastq.gz$"),
        # create output parameter to be supplied to next task
	      "{basedir[0]}/{sample[0]}/{replicate[0]}/fastq_trimmed/{prefix[0]}_trimmed.fq.gz",
        "{basedir[0]}/{sample[0]}/{replicate[0]}/qc",               # qc folder
        "{basedir[0]}/{sample[0]}/{replicate[0]}/fastq_trimmed",    # trimmed_folder
        logger, logger_mutex)
def trim_fastq(input_file, output_files, qc_folder, output_folder ,logger, logger_mutex):
    raw_fastq=os.path.basename(input_file)
    cmd = (" cd $TMPDIR ; "
         " cp {input_file} . ;"
         " trim_galore --fastqc  {raw_fastq} 2> {qc_folder}/trim_galore.log ; "
         " mv *.fq.gz  {output_folder} ; "
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


#_______________________________________________________________________________________________________
#
#              take trimmer output and align with bowtie2
#_______________________________________________________________________________________________________
@collate(trim_fastq, formatter("([^/]+)_L00[1234]_R1_001_trimmed.fq.gz$"),
                                "{subpath[0][1]}/bam/{1[0]}.fq.sorted.bam",
                                "{path[0]}","{subpath[0][1]}/bam",
                                "{subpath[0][1]}/qc",logger,logger_mutex)
def bowtie2(input_files, out_file, path, outpath,qc_folder,logger, logger_mutex):
    reads = []
    for i in input_files:
      reads.append(os.path.basename(i))
    reads = ','.join(reads)
    print reads
    bowtie2_output = out_file.split('/')
    bowtie2_output = bowtie2_output[-1]
    cmd = ( "cd $TMPDIR \n"
            "mkdir reference \n"
            "cp  {path}"  + "/*fq.gz" + " . \n "
            "ls -l \n"
            "date \n"
            "cp $HOME/Scratch/reference/grch38/bowtie2/*bt2 ./reference \n "
            " bowtie2 -k 4 -X2000 --mm --local --threads 8 "
            " -x  ./reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.bowtie_index "
            " -U {reads} 2> {qc_folder}/bowtie2.log | samtools view -bS - -o temp.bam \n"
            " samtools sort -@ 8 temp.bam -m 2G " + bowtie2_output[:-4] + " 2>{qc_folder}/samtools.log \n"
            " samtools flagstat {bowtie2_output} > {qc_folder}/{bowtie2_output}.mapstats \n"
            " cp {bowtie2_output} {outpath} \n"
            " rm -r * \n")
    cmd = cmd.format(**locals())
    #print cmd
    try:
        stdout_res, stderr_res = "",""
        stdout_res, stderr_res = run_job(cmd,
                                        job_name = "bowtie2",
                                        job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                        job_other_options    = "-w n -S /bin/bash -V -l h_rt=08:00:00 -w n -l mem=2G -l tmpfs=60G -pe smp 8 -wd /home/sejjctj/Scratch -j yes ",
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
        logger.debug("bowtie2 worked")


#_______________________________________________________________________________________________________
#
#              Post-alignment filtering

# =============================
# Remove  unmapped, mate unmapped
# not primary alignment, reads failing platform
# Remove low MAPQ reads
# ==================  

#_______________________________________________________________________________________________________

@transform(bowtie2,formatter("([^/]+)bam$"),"{subpath[0][1]}/bam/{1[0]}filt.nodup.srt.bam","{subpath[0][1]}/bam",
                            "{subpath[0][1]}/qc/filtering.log",
                             logger, logger_mutex )
def post_alignment_filter(input_file, output_file, out_dir,log_file, logger, logger_mutex):
  raw_bam=os.path.basename(input_file)
  prefix=raw_bam[:-4]
  FILT_BAM_PREFIX=prefix + ".filt.srt"
  FILT_BAM_FILE=FILT_BAM_PREFIX +".bam"
  MAPQ_THRESH=30
  TMP_FILT_BAM_FILE=FILT_BAM_PREFIX + "dupmark.bam"
  DUP_FILE_QC=FILT_BAM_PREFIX + ".dup.qc"
  FINAL_BAM_PREFIX=prefix + ".filt.nodup.srt"
  FINAL_BAM_FILE=FINAL_BAM_PREFIX + ".bam"
  FINAL_BAM_INDEX_FILE=FINAL_BAM_PREFIX + ".bai"
  FINAL_BAM_FILE_MAPSTATS=FINAL_BAM_PREFIX + ".flagstat.qc"
  PBC_FILE_QC=FINAL_BAM_PREFIX + ".pbc.qc"
  picard_loc="/shared/ucl/apps/picard-tools/1.136/picard-tools-1.136/"
  cmd=("cd $TMPDIR \n"
       "cp {input_file} . \n"
       "date \n"
       "ls -l \n"
       "\n"
       "samtools sort -@ 4 -m 8G {raw_bam} temporary \n"
       "mv temporary.bam {raw_bam} \n"
       "samtools view -@ 4 -F 1804 -q {MAPQ_THRESH} -b {raw_bam} > {FILT_BAM_FILE} \n"
       "mv temporary_bam.bam {FILT_BAM_FILE} \n"
       "echo \"first filter done\" \n"
       "ls -lh \n"
       "#=========================\n"
       "# Mark Duplicates \n"
       "#==========================\n"
       "\n"
       "java -Xmx8G -jar {picard_loc}picard.jar MarkDuplicates INPUT={FILT_BAM_FILE} \\\n"
       "OUTPUT={TMP_FILT_BAM_FILE} METRICS_FILE={DUP_FILE_QC} VALIDATION_STRINGENCY=LENIENT \\\n"
       "ASSUME_SORTED=true REMOVE_DUPLICATES=false \n"
       "mv {TMP_FILT_BAM_FILE} {FILT_BAM_FILE} \n"
       "echo \"mark duplicates done\" \n"
       "ls -lh \n"
       "date \n"
       "\n"
       "# ============================ \n"
       "# Remove duplicates\n"
       "# Index final position sorted BAM \n"
       "# ============================ \n"
       "\n"
       "samtools view -@ 4 -F 1804 -b {FILT_BAM_FILE} > {FINAL_BAM_FILE} \n"
       "\n"
       "# Index Final BAM file \n"
       "samtools index {FINAL_BAM_FILE} {FINAL_BAM_INDEX_FILE} \n"
       "samtools flagstat {FINAL_BAM_FILE} > {FINAL_BAM_FILE_MAPSTATS} \n"
       "# Compute library complexity \n"
       "# ============================= \n"
       "# sort by position and strand \n"
       "# Obtain unique count statistics \n"
       "\n"
       "PBC_FILE_QC={FINAL_BAM_PREFIX}.pbc.qc \n"
       "# PBC File output \n"
       "echo -e \"TotalReadPairs\\tDistinctReadPairs\\tOneReadPair\\tTwoReadPairs\\tNRF=Distinct/Total\\tPBC1=OnePair/Distinct\\tPBC2=OnePair/TwoPair\" > header \n"
       "bedtools bamtobed -i {FILT_BAM_FILE} | awk 'BEGIN{{OFS=\"\\t\"}}{{print $1,$2,$3,$6}}' | \\\n"
       "grep -v chrM | sort | uniq -c | awk 'BEGIN{{mt=0;m0=0;m1=0;m2=0}} ($1==1){{m1=m1+1}} ($1==2){{m2=m2+1}} \\\n"
       "{{m0=m0+1}} {{mt=mt+$1}} END{{printf \"%d\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n\",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}}' \\\n"
       "> {PBC_FILE_QC} \n"
       "mv {FINAL_BAM_FILE} {out_dir} \n"
       "cat header {PBC_FILE_QC} > temp_file && mv temp_file {PBC_FILE_QC} \n"
       "mv {PBC_FILE_QC} {out_dir} \n"
       "mv {FINAL_BAM_FILE} {out_dir} \n")
  cmd = cmd.format(**locals())
  try:
    stdout_res, stderr_res = "",""
    stdout_res, stderr_res = run_job(cmd,
                                     job_name = "filter_bam",
                                     job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                     job_other_options    = "-S /bin/bash -V -l h_rt=10:00:00 -w n -l mem=8G -l tmpfs=60G -pe smp 4 -wd /home/sejjctj/Scratch -j yes ",
                                     job_environment      = { 'BASH_ENV' : '/home/sejjctj/.bashrc' } ,
                                     retain_job_scripts   = True,  # retain job scripts for debuging, they go in Scratch/test_dir
                                     working_directory    = "/home/sejjctj/Scratch/test_dir",
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
      logger.debug("post_alignment_filter worked")

@transform(post_alignment_filter,formatter("(?P<basedir>[/.].+)/(?P<sample>[a-zA-Z0-9_\-\.]+)/(?P<replicate>replicate_[0-9])/(?P<bam_dir>bam)/(?P<prefix>[a-zA-Z0-9_\-]+).fq.sorted.filt.nodup.srt.bam$"),
          ["{basedir[0]}/{sample}[0]/{replicate[0]}/{bam_dir[0]}/{prefix[0]}.SE.tagAlign.gz",
           "{basedir[0]}/{sample[0]}/{replicate[0]}/{bam_dir[0]}/{prefix[0]}.sample.SE.tagAlign.gz"],
           "{basedir[0]}/{sample[0]}/{replicate[0]}/{bam_dir[0]}",
           "{prefix[0]}",
           logger, logger_mutex)
def bam_to_tagAlign(input_file, output_file, out_dir,prefix, logger, logger_mutex):
  FINAL_BAM_FILE=os.path.basename(input_file)
  FINAL_BAM_PREFIX=prefix
  cmd = ("# =================== \n"
         "# Create tagAlign file \n"
         "# =================== \n"
         "cd $TMPDIR \n"
         "cp {input_file} . \n"
         "FINAL_TA_FILE={FINAL_BAM_PREFIX}.SE.tagAlign.gz \n"
         "bedtools bamtobed -i {FINAL_BAM_FILE} | awk 'BEGIN{{OFS=\"\\t\"}}{{$4=\"N\";$5=\"1000\";print $0}}' | gzip -nc > \"$FINAL_TA_FILE\" \n"
         "# ================================= \n"
         "# Subsample tagAlign file \n"
         "# ================================ \n"
         "NREADS=15000000 \n"
         "SUBSAMPLED_TA_FILE={FINAL_BAM_PREFIX}.sample.SE.tagAlign.gz\n"
         "zcat \"$FINAL_TA_FILE\" | grep -v chrM | shuf -n \"$NREADS\" --random-source=\"$FINAL_TA_FILE\" | gzip -nc > \"$SUBSAMPLED_TA_FILE\" \n"
         "mv \"$SUBSAMPLED_TA_FILE\" {out_dir} \n"
         "mv \"$FINAL_TA_FILE\" {out_dir}")
  cmd = cmd.format(**locals())

  try:
    stdout_res, stderr_res = "",""
    stdout_res, stderr_res = run_job(cmd,
                                     job_name = "bam2tag",
                                     job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                     job_other_options    = "-S /bin/bash -V -l h_rt=10:00:00 -w n -l mem=24G -l tmpfs=30G -wd /home/sejjctj/Scratch -j yes ",
                                     job_environment      = { 'BASH_ENV' : '/home/sejjctj/.bashrc' } ,
                                     retain_job_scripts   = True,  # retain job scripts for debuging, they go in Scratch/test_dir
                                     working_directory    = "/home/sejjctj/Scratch/test_dir",
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
    logger.debug("bam_to_tagAlign worked")

@transform(bam_to_tagAlign, formatter("(?P<basedir>[/.].+)/(?P<sample>[a-zA-Z0-9_\-\.]+)/(?P<replicate>replicate_[0-9])/(?P<bam_dir>bam)/(?P<sub_sampled>[a-zA-Z0-9_]+.filt.nodup.sample.SE.tagAlign.gz$)"),
           "{basedir[0]}/{sample[0]}/{replicate[0]}/{bam_dir[0]}/{sample[0]}.cc.plot.pdf",
           "{basedir[0]}/{sample[0]}/{replicate[0]}/{bam_dir[0]}",
           "{sample[0]}.cc.plot.pdf",
           "{sample[0]}.cc.qc",
           logger, logger_mutex)
def phantom_peak_quals(input_file, output_file, out_dir, outfile1,outfile2,logger, logger_mutex):
  SUBSAMPLED_TA_FILE=os.path.basename(input_file)
  cmd = ("#########################\n"
         "# run  phantompeakquals #\n"
         "#########################\n"
         "cd $TMPDIR \n"
         "mkdir job_temp \n"
         "cp {input_file} . \n"
         "Rscript ~/applications/phantompeakqualtools/run_spp.R "
         " -c={SUBSAMPLED_TA_FILE} -filtchr=chrM "
         " -savp={outfile1} -out={outfile2} "
         " -tmpdir=./job_temp \n"
         "echo -e \"Filename\\tnumReads\\testFragLen\\tcorr_estFragLen\\tPhantomPeak\\tcorr_phantomPeak\\targmin_corr\\tmin_corr\\tphantomPeakCoef\\trelPhantomPeakCoef\\tQualityTag\" > header \n"
         "sed -r 's/,[^\\t]+//g' {outfile2} > temp \n"
         "cat header temp > temporary && mv temporary temp \n"
         "mv temp {outfile2} \n"
         "mv {outfile2} {out_dir}\n"
         "mv {outfile1} {out_dir}")

  cmd = cmd.format(**locals())

  try:
    stdout_res, stderr_res = "",""
    stdout_res, stderr_res = run_job(cmd,
                                     job_name = "phantom",
                                     job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                     job_other_options    = "-S /bin/bash -V -l h_rt=10:00:00 -w n -l mem=24G -l tmpfs=30G -wd /home/sejjctj/Scratch -j yes ",
                                     job_environment      = { 'BASH_ENV' : '/home/sejjctj/.bashrc' } ,
                                     retain_job_scripts   = True,  # retain job scripts for debuging, they go in Scratch/test_dir
                                     working_directory    = "/home/sejjctj/Scratch/test_dir",
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
    logger.debug("create_pseudoreplicates")
'''
@transform(bam_to_tagAlign, formatter("(?P<basedir>[/.].+)/(?P<sample>[a-zA-Z0-9_\-\.]+)/(?P<replicate>replicate_[0-9])/(?P<bam_dir>bam)/(?P<sub_sampled>[a-zA-Z0-9_]+.filt.nodup.sample.SE.tagAlign.gz$)"),
''' 
if __name__ == '__main__':
  cmdline.run (options, multithread = options.jobs)
  drmaa_session.exit()

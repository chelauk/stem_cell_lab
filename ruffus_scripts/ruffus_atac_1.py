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
        line = line.rstrip()
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
#   The first step trims pairs of fastqs and generates fastqc reports            .
#_________________________________________________________________________________

#@mkdir(input_files,
@mkdir(input_files,
       #match input pattern
        formatter("([^/]+)$"),
        # make qc directory
        ["{subpath[0][1]}/qc",
        # make trimmed_directory)
        "{subpath[0][1]}/fastq_trimmed",
        "{subpath[0][1]}/bam"])
@collate(input_files,
        # input formatter to provide read pairs
        formatter("([^/]+)R[12](.+)gz"),
        # create output parameter to be supplied to next task
	["{subpath[0][1]}/fastq_trimmed/{1[0]}R1_001_val_1.fq.gz",
	 "{subpath[0][1]}/fastq_trimmed/{1[0]}R2_001_val_2.fq.gz"],
        ["{1[0]}R1_001.fastq.gz","{1[0]}R2_001.fastq.gz"],    # basename for trim_galore
        "{subpath[0][1]}/qc",               # qc folder
        "{subpath[0][1]}/fastq_trimmed",    # trimmed_folder
        logger, logger_mutex)
def trim_fastq(input_files, output_files, basenames, qc_folder, output_folder ,logger, logger_mutex):
    
    #print input_files
    if len(input_files) !=2:
        raise Exception("One of the reads pairs %s missing" % (input_files,))
    cmd = (" source ~/.bashrc \n" 
         " cd $TMPDIR \n"
         " cp {input_files[0]} . \n"
         " cp {input_files[1]} . \n"
         " trim_galore --fastqc --paired {basenames[0]} {basenames[1]} 2> {qc_folder}/trim_galore.log \n"
         " mv *val_*.fq.gz  {output_folder} \n"
         " mv *fastqc*  {qc_folder} \n"
         " mv *report* {qc_folder}; rm * \n" )
  
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
                                      #job_environment      = { 'BASH_ENV' : '/home/sejjctj/.bashrc' } ,
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
@collate(trim_fastq, formatter("([^/]+)_L00[1234]_R[12]_001_val_[12].fq.gz$"),
                                "{subpath[0][1]}/bam/{1[0]}_L001_R1_001_val_1.fq.sorted.bam",
                                "{path[0]}","{subpath[0][1]}/bam",
                                "{subpath[0][1]}/qc",logger,logger_mutex)
def bowtie2(input_files, out_file, path, outpath,qc_folder,logger, logger_mutex):
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
    bowtie2_output = out_file.split('/')
    bowtie2_output = bowtie2_output[-1]

    cmd = ( " cd $TMPDIR \n"
            " mkdir reference \n"
            " cp  {path}"  + "/*fq.gz" + " . \n "
            " ls -l \n"
            " date \n"
            " cp $HOME/Scratch/reference/grch38/bowtie2/*bt2 ./reference \n "
            " bowtie2 -k 4 -X2000 --mm --local --threads 8 \\\n"
            " -x  ./reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.bowtie_index \\\n"
            " -1 {first_reads} -2 {second_reads} 2> {qc_folder}/bowtie2.log | samtools view -bS - -o temp.bam \n"
            " samtools sort -n -@ 8 temp.bam -m 2G " + bowtie2_output[:-4] + " 2>{qc_folder}/samtools.log \n"
            " cp " + bowtie2_output + " {outpath} \n"
            " rm -r * \n ")
    cmd = cmd.format(**locals())
    #print cmd
    try:
        stdout_res, stderr_res = "",""
        stdout_res, stderr_res = run_job(cmd,
                                        job_name = "bowtie2",
                                        job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                        job_other_options    = "-w n -S /bin/bash -V -l h_rt=08:00:00 -w n -l mem=2G -l tmpfs=60G -pe smp 8 -wd /home/sejjctj/Scratch -j yes ",
                                        #job_environment      = { 'BASH_ENV' : '/home/sejjctj/.bashrc' } ,
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
# Remove unmapped, mate unmapped
# not primary alignment, reads failing platform
# Only keep properly paired reads
# Obtain name sorted BAM file
# Remove orphan reads (pair was removed)
# and read pairs mapping to different chromosomes
# Obtain position sorted BAM
# ============================
# Remove duplicates
# Index final position sorted BAM
# Create final name sorted BAM
# ============================# ==================
# =============================
# Compute library complexity
# =============================
# Sort by name
# convert to bedPE and obtain fragment coordinates
# sort by position and strand
# Obtain unique count statisticsi
# TotalReadPairs  DistinctReadPairs  OneReadPair TwoReadPairs 

#_______________________________________________________________________________________________________

@collate(bowtie2,formatter("([^/]+)bam$"),"{subpath[0][1]}/bam/{1[0]}filtered.bam","{subpath[0][1]}/bam",
                            "{subpath[0][1]}/qc/filtering.log",
                             logger, logger_mutex )
def post_alignment_filter(input_file, output_file, out_dir,log_file, logger, logger_mutex):
    input_file = input_file[0]
    print "input:"
    print input_file
    print "output:"
    print output_file
    output_bam = os.path.basename(output_file)
    print "output_bam:"
    print output_bam
    RAW_BAM_FILE=os.path.basename(input_file)
    OFPREFIX = output_bam[:-4]
    print OFPREFIX
    FILT_BAM_PREFIX=OFPREFIX + ".filt"
    FILT_BAM_FILE=FILT_BAM_PREFIX + ".bam"
    TMP_FILT_BAM_PREFIX="tmp." + FILT_BAM_PREFIX + ".nmsrt"
    TMP_FILT_BAM_FILE=TMP_FILT_BAM_PREFIX + ".bam"
    TMP_FILT_FIXMATE_BAM_FILE=TMP_FILT_BAM_PREFIX + ".fixmate.bam"
    TMP_DUP_BAM_FILE=FILT_BAM_PREFIX + ".dupmark.bam"
    DUP_FILE_QC=FILT_BAM_PREFIX + ".dup.qc"
    FINAL_BAM_PREFIX=OFPREFIX + ".nodup"
    FINAL_BAM_FILE=FINAL_BAM_PREFIX + ".bam"
    FINAL_BAM_INDEX_FILE=FINAL_BAM_FILE + ".bai"
    FINAL_BAM_FILE_MAPSTATS=FINAL_BAM_PREFIX + ".flagstat.qc"
    PBC_FILE_QC=OFPREFIX + ".pbc.qc"
    picard_loc="/shared/ucl/apps/picard-tools/1.136/picard-tools-1.136/"
    cmd = ( " # =============================  \n"
            " # Remove unmapped, mate unmapped \n"
            " # not primary alignment, reads failing platform \n"
            " # Only keep properly paired reads \n"
            " # Obtain name sorted BAM file \n"
            " # ================== \n"
            " source ~/.bashrc \n"
            " cd $TMPDIR \n"
            " cp {input_file} ./ \n"
            " ls -lh \n"
            " date \n"
            " samtools view -F 524 -f 2 -u {RAW_BAM_FILE} \\\n"
            " | sambamba sort -n -m 16G -t 4 /dev/stdin -o {TMP_FILT_BAM_FILE} \n"
            " samtools view -h {TMP_FILT_BAM_FILE} | assign_multimappers.py -k 4 --paired-end \\\n"
            " | samtools fixmate -r /dev/stdin {TMP_FILT_FIXMATE_BAM_FILE} \n"
            " ls -lh \n"
            " date \n"
            " # Remove orphan reads (pair was removed) \n"
            " # and read pairs mapping to different chromosomes \n"
            " # obtain position sorted BAM \n"
            " samtools view -F 1804 -f 2 -u {TMP_FILT_FIXMATE_BAM_FILE} \\\n"
            " | sambamba sort -m 16G -t 4 /dev/stdin -o {FILT_BAM_FILE} \n"
            " rm {TMP_FILT_FIXMATE_BAM_FILE} \n"
            " rm {TMP_FILT_BAM_FILE} \n"
            " ls -lh \n"
            " date \n"
            " # ============= \n"
            " # Mark duplicates \n"
            " # ============= \n"
            " java -Xmx16G -jar {picard_loc}picard.jar MarkDuplicates INPUT={FILT_BAM_FILE} "
            " OUTPUT={TMP_DUP_BAM_FILE} METRICS_FILE={DUP_FILE_QC} "
            " VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true  REMOVE_DUPLICATES=false \n"
            " mv {TMP_DUP_BAM_FILE}  {FILT_BAM_FILE} \n"
            " # ============================ \n"
            " # Remove duplicates \n"
            " # Index final position sorted BAM \n"
            " # Create final name sorted BAM \n"
            " # ============================ \n"
            " samtools view -F 1804 -f 2 -b {FILT_BAM_FILE} > {output_bam} \n"
            " samtools sort -n -m 16G -@ 4 {output_bam} {OFPREFIX}.final_filt_nmsrt \n"
            " # used later on \n"
            " samtools index {output_bam} \n"
            " samtools flagstat {output_bam} > {output_bam}.mapstats \n"
            " mv {output_bam}.mapstats {out_dir}\n"
            " # ============================= \n"
            " # Compute library complexity    \n"
            " # ============================= \n"
            " # Sort by name \n"
            " # convert to bedPE and obtain fragment coordinates \n"
            " # sort by position and strand \n"
            " # Obtain unique count statistics \n"
            " sambamba sort -n -m 16G -t 4 {FILT_BAM_FILE} -o {OFPREFIX}.srt.tmp.bam  \n"
            " echo -e '# PBC File output\n# TotalReadPairs\tDistinctReadPairs\tOneReadPair\tTwoReadPairs\tNRF=Distinct/Total\tPBC1=OnePair/Distinct\tPBC2=OnePair/TwoPair' > header \n"
            " bedtools bamtobed -bedpe -i {OFPREFIX}.srt.tmp.bam \\\n"
            " | awk 'BEGIN{{OFS=\"\\t\"}}{{print $1,$2,$4,$6,$9,$10}}' | grep -v 'chrM' | sort \\\n"
            " | uniq -c | \\\n"
            " awk 'BEGIN{{mt=0;m0=0;m1=0;m2=0}}($1==1){{m1=m1+1}} ($1==2){{m2=m2+1}} {{m0=m0+1}}{{t=mt+$1}}END{{printf\"%d\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n\",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}}' " 
            "  > {PBC_FILE_QC} \n"
            " rm {FILT_BAM_FILE} \n"
            " mv {output_bam} {out_dir} \n"
            " mv {output_bam}.bai {out_dir} \n"
            " cat header header {PBC_FILE_QC} > temporary && mv temporary > {PBC_FILE_QC} \n"
            " mv {PBC_FILE_QC} {out_dir} \n"
            " mv {OFPREFIX}.final_filt_nmsrt.bam {out_dir} ")

    cmd = cmd.format(**locals())
    try:
      stdout_res, stderr_res = "",""
      stdout_res, stderr_res = run_job(cmd,
                                     job_name = "filter_bam",
                                     job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                     job_other_options    = "-S /bin/bash -V -l h_rt=10:00:00 -w n -l mem=16G -l tmpfs=60G -pe smp 4 -wd /home/sejjctj/Scratch -j yes ",
                                     #job_environment      = { 'BASH_ENV' : '/home/sejjctj/.bashrc' } ,
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


###################################################################################################
# Convert PE BAM to tagAlign (BED 3+3 format): function _bam_to_bedpe() and function _bam_to_tag()#
###################################################################################################


@transform(post_alignment_filter,formatter("([^/]+)filtered.bam$"), 
          "{subpath[0][1]}/bam/{1[0]}filtered.final_filt_nmsrt.bedpe.gz", "{subpath[0][1]}/bam", logger, logger_mutex)
def bam_to_tagAlign(input_file, output_file, out_dir, logger, logger_mutex):
  print "\n\ninput_file: " + str(input_file) + "\n\n"
  print "\n\ninput_file: " + str(output_file) + "\n\n"
  FINAL_BAM_FILE=os.path.basename(input_file)
  FINAL_BAM_PREFIX=FINAL_BAM_FILE[:-4]
  BAM_LOC=os.path.dirname(input_file)
  OFPREFIX=FINAL_BAM_FILE[:-4]
  FINAL_NMSRT_BAM=OFPREFIX + ".final_filt_nmsrt.bam"
  FINAL_NMSRT_BAM_PREFIX = FINAL_NMSRT_BAM[:-4]
  FINAL_BEDPE_FILE=FINAL_NMSRT_BAM_PREFIX + ".bedpe.gz"
  FINAL_TA_FILE=FINAL_BAM_PREFIX +".PE2SE.tagAlign.gz"
  NREADS=25000000
  SUBSAMPLED_TA_FILE=OFPREFIX + ".filt.nodup.sample" + str(25) + ".MATE1.tagAlign.gz"
  cmd = ("# =================== \n"
        "# Create tagAlign file \n"
        "# =================== \n"
        "source ~/.bashrc \n"
        "cd $TMPDIR \n"
        "cp {input_file} . \n"
        "cp {BAM_LOC}/{FINAL_NMSRT_BAM} . \n"
        "# Create virtual SE file containing both read pairs \n"
        "bedtools bamtobed -i {FINAL_BAM_FILE} \\\n"
        " | awk 'BEGIN{{OFS=\"\\t\"}}{{$4=\"N\";$5=\"1000\";print $0}}' | gzip -nc > {FINAL_TA_FILE} \n"
        "# ================ \n"
        "# Create BEDPE file \n"
        "# ================ \n"
        "bedtools bamtobed -bedpe -mate1 -i {FINAL_NMSRT_BAM} | gzip -nc > {FINAL_BEDPE_FILE} \n"
        "# ================================= \n"
        "# Subsample tagAlign file \n"
        "# Restrict to one read end per pair for CC analysis \n"
        "# ================================ \n"
        "zcat {FINAL_BEDPE_FILE} | grep -v \"chrM\" | shuf -n {NREADS} --random-source={FINAL_BEDPE_FILE}  \\\n"
        " | awk 'BEGIN{{OFS=\"\\t\"}}{{print $1,$2,$3,\"N\",\"1000\",$9}}' | gzip -nc > {SUBSAMPLED_TA_FILE} \n"
        "mv {FINAL_TA_FILE} {out_dir} \n"
        "mv {SUBSAMPLED_TA_FILE} {out_dir} \n"
        "mv {FINAL_BEDPE_FILE} {out_dir} ")
  cmd = cmd.format(**locals())

  try:
    stdout_res, stderr_res = "",""
    stdout_res, stderr_res = run_job(cmd,
                                     job_name = "bam2tag",
                                     job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                     job_other_options    = "-S /bin/bash -V -l h_rt=10:00:00 -w n -l mem=24G -l tmpfs=30G -wd /home/sejjctj/Scratch -j yes ",
                                     #job_environment      = { 'BASH_ENV' : '/home/sejjctj/.bashrc' } ,
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

@follows(bam_to_tagAlign)
@transform(bam_to_tagAlign,formatter("([^/]+).final_filt_nmsrt.bedpe.gz$"),
          "{subpath[0][1]}/bam/{1[0]}.PE2SE.pr2.tagAlign.gz", "{subpath[0][1]}/bam",
          logger, logger_mutex)
def create_pseudoreplicates(input_file, output_file, out_dir, logger, logger_mutex):
  FINAL_BEDPE_FILE=os.path.basename(input_file)
  PR_PREFIX=FINAL_BEDPE_FILE[:-26]
  PR1_TA_FILE=PR_PREFIX + ".PE2SE.pr1.tagAlign.gz"
  PR2_TA_FILE=PR_PREFIX + ".PE2SE.pr2.tagAlign.gz"
  cmd = ("# ========================\n"
       "# Create pseudoReplicates\n"
       "# =======================\n"
       "source ~/.bashrc \n"
       "cd $TMPDIR \n"
       "cp {input_file} . \n"
       "# Get total number of read pairs \n"
       "nlines=$( zcat {FINAL_BEDPE_FILE} | wc -l ) \n"
       "nlines=$(( (nlines + 1) / 2 )) \n"
        "# Shuffle and split BEDPE file into 2 equal parts \n"
       "zcat {FINAL_BEDPE_FILE} | shuf --random-source={FINAL_BEDPE_FILE} | split -d -l $nlines - {PR_PREFIX} \n"
       "# Will produce {PR_PREFIX}00 and {PR_PREFIX}01 \n"
       "# Convert read pairs to reads into standard tagAlign file \n"
       "awk 'BEGIN{{OFS=\"\\t\"}}{{printf \"%s\\t%s\\t%s\\tN\\t1000\\t%s\\n%s\\t%s\\t%s\\tN\\t1000\\t%s\\n\",$1,$2,$3,$9,$4,$5,$6,$10}}' {PR_PREFIX}00 | gzip -nc > {PR1_TA_FILE} \n"        
       "rm {PR_PREFIX}00 \n"
       "awk 'BEGIN{{OFS=\"\\t\"}}{{printf \"%s\\t%s\\t%s\\tN\\t1000\\t%s\\n%s\\t%s\\t%s\\tN\\t1000\\t%s\\n\",$1,$2,$3,$9,$4,$5,$6,$10}}' {PR_PREFIX}01 | gzip -nc > {PR2_TA_FILE} \n"
       "rm {PR_PREFIX}01 \n"
       "mv {PR1_TA_FILE} {out_dir} \n"
       "mv {PR2_TA_FILE} {out_dir} "
    )
  cmd = cmd.format(**locals())
  try:
    stdout_res, stderr_res = "",""
    stdout_res, stderr_res = run_job(cmd,
                                     job_name = "create_pseudo",
                                     job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                     job_other_options    = "-S /bin/bash -V -l h_rt=01:00:00 -w n -l mem=8G -l tmpfs=30G -wd /home/sejjctj/Scratch -j yes ",
                                     #job_environment      = { 'BASH_ENV' : '/home/sejjctj/.bashrc' } ,
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




@transform(bam_to_tagAlign, formatter("(?P<basedir>[/.].+)/(?P<sample>[a-zA-Z0-9_\-\.]+)/(?P<replicate>replicate_[0-9])/(?P<bam_dir>bam)/(?P<bed_pe>[a-zA-Z0-9_\-\.]+$)"),
           "{basedir[0]}/{sample[0]}/{replicate[0]}/{bam_dir[0]}/{sample[0]}.cc.plot.pdf",
           "{basedir[0]}/{sample[0]}/{replicate[0]}/{bam_dir[0]}",
           "{sample[0]}.cc.plot.pdf",
           "{sample[0]}.cc.qc",
           logger, logger_mutex)
def phantom_peak_quals(input_file, output_file, out_dir, outfile1,outfile2,logger, logger_mutex):
  
  SUBSAMPLED_TA_FILE=os.path.basename(input_file)
  SUBSAMPLED_TA_FILE=SUBSAMPLED_TA_FILE[:-25] + "filt.nodup.sample25.MATE1.tagAlign.gz"
 
  cmd = ("#########################\n"
         "# run  phantompeakquals #\n"
         "#########################\n"
         "source ~/.bashrc \n"
         "cd $TMPDIR \n"
         "mkdir job_temp \n"
         "mv {out_dir}/{SUBSAMPLED_TA_FILE} . \n"
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
                                     #job_environment      = { 'BASH_ENV' : '/home/sejjctj/.bashrc' } ,
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

@transform(bam_to_tagAlign, formatter("(?P<basedir>[/.].+)/(?P<sample>[a-zA-Z0-9_\-\.]+)/(?P<replicate>replicate_[0-9])/(?P<bam_dir>bam)/(?P<prefix>[a-zA-Z0-9_\-\.]+_S[0-9]+_L00[1234]_R[12]_[0-9]+_val_[0-9]+.fq.sorted.filtered)(.final_filt_nmsrt.bedpe.gz$)"),
           "{basedir[0]}/{sample[0]}/{replicate[0]}/{bam_dir[0]}/{prefix[0]}.PE2SE.pr2.tn5.tagAlign.gz",
           "{basedir[0]}/{sample[0]}/{replicate[0]}/{bam_dir[0]}",logger, logger_mutex)  
def tn5_shift(input_file, output_file, out_dir, logger, logger_mutex):
  cmd = ("#==================================\n"
         "# TN5 shift for atac seq           \n"
         "#==================================\n"
         "source ~/.bashrc \n"
         "cd $TMPDIR \n"
         "cp {out_dir}/*tagAlign.gz . \n"
         "for tag in *tagAlign.gz \n"
         "do zcat ""$tag"" | awk -F $'\t' 'BEGIN {{OFS = FS}}{{ if ($6 == \"+\") {{$2 = $2 + 4}} else if ($6 == \"-\") {{$3 = $3 - 5}} print $0}}' | \\\n"
         "gzip -nc > ""${{tag:0:${{#tag}}-12}}.tn5.tagAlign.gz"" \n"
         "done \n"
         "mv *tn5* {out_dir} \n")
  cmd = cmd.format(**locals())
  try:
    stdout_res, stderr_res = "",""
    stdout_res, stderr_res = run_job(cmd,
                                     job_name = "tn5_shift",
                                     job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                     job_other_options    = "-S /bin/bash -V -l h_rt=04:00:00 -w n -l mem=4G -l tmpfs=10G -wd /home/sejjctj/Scratch -j yes ",
                                     #job_environment      = { 'BASH_ENV' : '/home/sejjctj/.bashrc' } ,
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
    logger.debug("tn5_shift")

@transform(tn5_shift, formatter("(?P<basedir>[/.].+)/(?P<sample>[a-zA-Z0-9_\-\.]+)/(?P<replicate>replicate_[0-9])/(?P<bam_dir>bam)/(?P<prefix>[a-zA-Z0-9_\-\.]+).tagAlign.gz$"),
           "{basedir[0]}/{sample[0]}/{replicate[0]}/{bam_dir[0]}/{prefix[0]}.narrowPeak.gz",
           "{basedir[0]}/{sample[0]}/{replicate[0]}/{bam_dir[0]}",logger, logger_mutex)
def macs2(input_file, output_file,out_dir, logger, logger_mutex):
  cmd = ("#===================================\n"
         "#  run mac2 2 on tn5 shifted files  \n"
         "#===================================\n"
         "source ~/.bashrc \n"
         "cd $TMPDIR \n"
         "cp {out_dir}/*tn5.tagAlign.gz . \n"
         "for tag in *tagAlign.gz \n"
         "do \n"
         "prefix=""${{tag:0:${{#tag}}-12}}""   #remove.tagAlign.gz \n" 
         "peakfile=""${{prefix}}"".narrowPeak.gz \n"
         "pval_thresh=0.01 \n"
         "fc_bedgraph=""${{prefix}}"".fc.signal.bedgraph \n"
         "fc_bedgraph_srt=""${{prefix}}"".fc.signal.srt.bedgraph \n"
         "fc_bigwig=""${{prefix}}""_sig.fc.signal.bigwig \n"
         "pval_bedgraph=""${{prefix}}"".pval.signal.bedgraph \n"
         "pval_bedgraph_srt=""${{prefix}}"".pval.signal.srt.bedgraph \n"
         "pval_bigwig=""${{prefix}}_sig.pval.signal.bigwig \n"
         "chrsz=\"/home/sejjctj/Scratch/reference/grch38/hg38.chrom.sizes\" \n"
         "## see https://github.com/taoliu/MACS/issues/145 for choice of --shift and --extsize \n"
         "macs2 callpeak \\\n"
         "-t ""$tag"" -f BED -n ""$prefix"" -g 2700000000 -p $pval_thresh \\\n"
         "--nomodel --shift -100 --extsize 200 -B --SPMR --keep-dup all --call-summits \n"
         "# Sort by Col8 in descending order and replace long peak names in Column 4 with Peak_<peakRank> \n"
         "sort -k 8gr,8gr \"$prefix\"_peaks.narrowPeak | awk 'BEGIN{{OFS=\"\\t\"}}{{$4=""Peak_""NR ; print $0}}' \\\n"
         " | gzip -nc > ""$peakfile"" \n"
         "rm -f \"$prefix\"_peaks.narrowPeak \n"
         "rm -f \"$prefix\"_peaks.xls \n"
         "rm -f \"$prefix\"_summits.bed \n"
         '''
	 "macs2 bdgcmp -t \"$prefix\"_treat_pileup.bdg -c \"$prefix\"_control_lambda.bdg \\\n"
         "--o-prefix ""$prefix"" -m FE \n"
         "slopBed -i \"$prefix\"_FE.bdg -g ""$chrsz"" -b 0 | bedClip stdin ""$chrsz"" ""$fc_bedgraph"" \n"
         "rm -f ""$prefix""_FE.bdg \n"
         "sort -k1,1 -k2,2n ""$fc_bedgraph"" > ""$fc_bedgraph_srt"" \n"
         "bedGraphToBigWig ""$fc_bedgraph_srt"" ""$chrsz"" ""$fc_bigwig"" \n"
         "rm -f ""$fc_bedgraph"" ""$fc_bedgraph_srt"" \n"
         "# sval counts the number of tags per million in the compressed BED file \n"
         "#sval=$(wc -l <(zcat -f \"$tag\" ) | awk '{{printf \"%f\", $1/1000000}}') \n"
         "sval=$(zcat \"$tag\" | wc -l | awk '{{print $1/1000000}}') \n"
         "macs2 bdgcmp \\\n"
         "-t \"$prefix\"_treat_pileup.bdg -c \"$prefix\"_control_lambda.bdg \\\n"
         "--o-prefix ""$prefix"" -m ppois -S ""${{sval}}"" \n"
         "slopBed -i \"$prefix\"_ppois.bdg -g ""$chrsz"" -b 0 | \\\n"
         "bedClip stdin ""$chrsz"" ""$pval_bedgraph"" \n"
         "rm -f \"$prefix\"_ppois.bdg \n"
         "sort -k1,1 -k2,2n ""$pval_bedgraph"" > ""$pval_bedgraph_srt"" \n"
         "bedGraphToBigWig ""$pval_bedgraph_srt"" ""$chrsz"" ""$pval_bigwig"" \n"
         "rm -f ""$pval_bedgraph"" ""$pval_bedgraph_srt"" \n"
         "rm -f \"$prefix\"_treat_pileup.bdg \"$prefix\"_control_lambda.bdg \n"
         '''
	 "mv ./\"$prefix\"* {out_dir} \n"
         "done \n")
  cmd = cmd.format(**locals())
  try:
    stdout_res, stderr_res = "",""
    stdout_res, stderr_res = run_job(cmd,
                                     job_name = "macs2",
                                     job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                     job_other_options    = "-S /bin/bash -V -l h_rt=08:00:00 -w n -l mem=16G -l tmpfs=60G -wd /home/sejjctj/Scratch -j yes ",
                                     #job_environment      = { 'BASH_ENV' : '/home/sejjctj/.bashrc' } ,
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
    logger.debug("mac2_callpeaks")


@transform(macs2, formatter("(?P<basedir>[/.].+)/(?P<sample>[a-zA-Z0-9_\-\.]+)/(?P<replicate>replicate_[0-9])/(?P<bam_dir>bam)/(?P<prefix>[a-zA-Z0-9_\-\.]+)(.tn5.narrowPeak.gz$)"),
           "{basedir[0]}/{sample[0]}/{replicate[0]}/{bam_dir[0]}/{prefix[0]}.PE2SE.pr2.tn5.narrowPeak.filt.gz",
           "{basedir[0]}/{sample[0]}/{replicate[0]}/{bam_dir[0]}",logger, logger_mutex)
def blacklist(input_file, output_file,out_dir, logger, logger_mutex): 
  cmd = ("#===================================\n"
         "#  run mac2 2 on tn5 shifted files  \n"
         "#===================================\n"
         "source ~/.bashrc \n"
         "cd $TMPDIR \n"
         "cp {out_dir}/*narrowPeak.gz . \n"
         "blacklist=\"/home/sejjctj/Scratch/reference/grch38/chipseq_blacklist/hg38.blacklist.bed.gz\" \n"
         "for peak in *narrowPeak.gz \n"
         "do \n"
         "prefix=""${{tag:0:${{#tag}}-14}}""   #remove .narrowPeak.gz \n" 
         "filtered_peak=""${{prefix}}.narrowPeak.filt.gz \n"
         "bedtools intersect -v a ${{peak}} -b ${{blacklist}} \\\n"
         "| awk 'BEGIN{{OFS=\"\\t\"}}{{if($5>1000) $5=1000; print $0}}' \\\n"
         "| grep -P 'chr[\dXY]+[\\t]' | gzip -nc > ${{filtered_peak}} \n"
         "mv ${{filtered_peak}} {out_dir} \n"
         "done \n")

  cmd = cmd.format(**locals())
  try:
    stdout_res, stderr_res = "",""
    stdout_res, stderr_res = run_job(cmd,
                                     job_name = "blacklist",
                                     job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                     job_other_options    = "-S /bin/bash -V -l h_rt=08:00:00 -w n -l mem=16G -l tmpfs=60G -wd /home/sejjctj/Scratch -j yes ",
                                     #job_environment      = { 'BASH_ENV' : '/home/sejjctj/.bashrc' } ,
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
    logger.debug("blacklist")


if __name__ == '__main__':
  cmdline.run (options, multithread = options.jobs)
  drmaa_session.exit()

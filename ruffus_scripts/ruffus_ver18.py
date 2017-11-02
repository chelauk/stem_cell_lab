#!/usr/bin/env python
from glob import glob
import glob
import sys, os, fnmatch
import re
from ruffus import *
import ruffus.cmdline as cmdline
from subprocess import check_call
from ruffus.drmaa_wrapper import run_job, error_drmaa_job
import drmaa

### parse command line for arguments

parser = cmdline.get_argparse(description="Chela's Pipeline")
parser.add_argument('-i', '--input', metavar="FILE", help = "Fastq files")
parser.add_argument('--cuffdiff_file', metavar="FILE", help = "cuffdiff comparison instructions")
parser.add_argument('--basedir', metavar="DIR", help = "base directory")
parser.add_argument('--aligner', metavar="choice", help = "choice of aligner; enter hisat or star")
parser.add_argument('--kallisto', metavar="choice", help = "use kallisto?")

options = parser.parse_args()
cuffdiff_file = options.cuffdiff_file
basedir=options.basedir
aligner=options.aligner
kallisto=options.kallisto

hisat_check=aligner=="hisat"
star_check=aligner=="star"
kallisto_check=kallisto=="yes"

print hisat_check
print star_check

print "Basedir: " + basedir

#  standard python logger which can be synchronised across concurrent Ruffus tasks
logger, logger_mutex = cmdline.setup_logging ("Chela", options.log_file, options.verbose)

if not aligner: 
    raise Exception ("Aligner not selected with --aligner")
files_list=options.input
if not files_list:
  raise Exception ("No matching samples specified with --input.")

input_files=[]
with open(files_list, 'r') as f:
    content = [line.decode('utf-8').rstrip('\n') for line in f] 
    for line in content:
        #print(line)
        line = line.rstrip()
        #print(basedir + line)
        print "input " + str(glob.glob(basedir + line + "/replicate*/fastq_raw/*gz"))
        input_files.append(glob.glob(basedir + line + "/replicate*/fastq_raw/*gz"))  
#                                                                                .
#   Useful code to turn input files into a flat list                             .
#                                                                                .
input_files = [item for sublist in input_files for item in sublist]  
print input_files
# start drmaa_session
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
        formatter("(?P<basedir>[/.].+)/(?P<sample>[a-zA-Z0-9_\-\.]+)/(?P<replicate>replicate_[0-9])/(?P<fastq_dir>fastq_raw)/(?P<fastq_raw>[a-zA-Z0-9_\-\.]+$)"),
        # make qc directory
        ["{basedir[0]}/{sample[0]}/{replicate[0]}/qc",
        # make trimmed_directory)
        "{basedir[0]}/{sample[0]}/{replicate[0]}/fastq_trimmed",
        "{basedir[0]}/{sample[0]}/{replicate[0]}/bam",
        "{basedir[0]}/{sample[0]}/{replicate[0]}/kallisto",
        "{basedir[0]}/{sample[0]}/{replicate[0]}/cufflinks",
        "{basedir[0]}/cuffmerge_out",
        "{basedir[0]}/cuffdiff",
        "{basedir[0]}/{sample[0]}/{replicate[0]}/cuffquant"])
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
    print "OUTPUT FILES!   " + str(output_files)    
    if len(input_files) !=2:
        raise Exception("One of the reads pairs %s missing" % (input_files,))
    cmd = ( " date \n"
            " echo $HOSTNAME \n"
            " cd $TMPDIR \n"
            " cp {input_files[0]} . \n"
            " cp {input_files[1]} . \n"
            " date \n"
            " ls -l \n"
            " trim_galore --fastqc --paired {basenames[0]} {basenames[1]} &> {qc_folder}/trim_galore.log \n"
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
@active_if(hisat_check)
@collate(trim_fastq, formatter("(?P<basedir>[/.].+)/(?P<sample>[a-zA-Z0-9_\-\.]+)/(?P<replicate>replicate_[0-9])/(?P<fastq_trimmed>fastq_trimmed)/(?P<trimmed_fq>[a-zA-Z0-9_\-\.]+)"),
                               "{basedir[0]}/{sample[0]}/{replicate[0]}/bam/{sample[0]}.{replicate[0]}.sorted.bam",
                               "{basedir[0]}/{sample[0]}/{replicate[0]}/{fastq_trimmed[0]}",
                               "{basedir[0]}/{sample[0]}/{replicate[0]}/bam",
                               "{basedir[0]}/{sample[0]}/{replicate[0]}/qc",logger,logger_mutex)
def hisat2(input_files, out_file, path, outpath,qc_folder,logger, logger_mutex):
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

    cmd = ( "cd $TMPDIR \n"
            "mkdir reference \n"
            "cp  {path}/*fq.gz  . \n"
            "cp $HOME/Scratch/reference/grch38_snp_tran/genome* ./reference \n"
            "hisat2 -p 8 -x ./reference/genome_snp_tran  --dta-cufflinks \\\n"
            "--novel-splicesite-outfile ./novel_splice.txt \\\n"
            "--novel-splicesite-infile ./novel_splice.txt \\\n"
            "-1 {first_reads} \\\n"
            "-2 {second_reads} \\\n"
            "2> {qc_folder}/hisat.log | samtools view -bS - -o temp.bam \n"
            "samtools sort -n -@ 8 temp.bam -m 4G " + hisat_output[:-4] + " 2>{qc_folder}/samtools.log \n"
            "mv {hisat_output} {outpath} \n"
            "mv novel_splice.txt {outpath} \n")
    cmd = cmd.format(**locals())
    try:
        stdout_res, stderr_res = "",""
        stdout_res, stderr_res = run_job(cmd,
                                        job_name = "hisat",
                                        job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                        job_other_options    = "-w n -S /bin/bash -V -l h_rt=08:00:00 -w n -l mem=4G -l tmpfs=60G -pe smp 8 -wd /home/sejjctj/Scratch -j yes ",
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
#              take trimmer output and align with star
#_______________________________________________________________________________________________________
@active_if(star_check)
@collate(trim_fastq, formatter("(?P<basedir>[/.].+)/(?P<sample>[a-zA-Z0-9_\-\.]+)/(?P<replicate>replicate_[0-9])/(?P<trimmed_dir>fastq_trimmed)/(?P<fastq>[a-zA-Z0-9_\-\.]+$)"),
                              "{basedir[0]}/{sample[0]}/{replicate[0]}/bam/{sample[0]}.Aligned.sortedByCoord.out.bam",
                              "{basedir[0]}/{sample[0]}/{replicate[0]}/{trimmed_dir[0]}",
                              "{basedir[0]}/{sample[0]}/{replicate[0]}/bam",
		              "{sample[0]}",
                              "{basedir[0]}/{sample[0]}/{replicate[0]}/qc",
                               logger, logger_mutex)
def star(input_files, out_file, path,outpath,sample,qc_folder,logger, logger_mutex):
  flat_list = [item for sublist in input_files for item in sublist]
  print(flat_list)
  first_reads = []
  second_reads =[]
  for i in flat_list:
    if re.search('val_1', os.path.basename(i)):
      first_reads.append(os.path.basename(i))
    elif re.search('val_2', os.path.basename(i)):
       second_reads.append(os.path.basename(i))
  first_reads = ','.join(first_reads)
  second_reads = ','.join(second_reads)
  star_output = out_file.split('/')
  star_output = star_output[-1]
  #print star_output
  cmd = ( "cd $TMPDIR \n "
          "cp {path}/*fq.gz  . \n "
          "~/applications/STAR-2.5.3a/bin/Linux_x86_64/STAR --runThreadN 4 \\\n"
          "--genomeDir ~/Scratch/reference/star_single_cell/index/ \\\n"
          "--readFilesIn " + first_reads + " " + second_reads + " \\\n"
          "--readFilesCommand zcat \\\n" 
          "--outSAMstrandField intronMotif \\\n"
          "--outFilterIntronMotifs RemoveNoncanonical \\\n" ## added for compatibility with
          "--outFileNamePrefix {sample} \\\n"               ## cufflinks
          "--outSAMtype BAM SortedByCoordinate \\\n"  
          "cp *bam {outpath} \\\n"
          "cp *Log.* {qc_folder} ")
  cmd = cmd.format(**locals())
  print cmd
  try:
    stdout_res, stderr_res = "",""
    stdout_res, stderr_res = run_job(cmd,
                                     job_name = "star",
                                     job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                     job_other_options    = "-w n -S /bin/bash -V -l h_rt=02:00:00 -w n -l mem=24G -l tmpfs=60G -pe smp 4 -wd /home/sejjctj/Scratch -j yes ",
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
    logger.debug("star worked")


#_______________________________________________________________________________________________________
################# USING -G flag for cufflinks, no novel Isoforms ######################################
#_______________________________________________________________________________________________________
#
#              cufflinks generate gtf files from sorted hisat/star output
#_______________________________________________________________________________________________________
################# USING -G flag for cufflinks, no novel Isoforms ######################################
@active_if(star_check or hisat_check)
@transform([hisat2,star],formatter("(?P<basedir>[/.].+)/(?P<sample>[a-zA-Z0-9_\-\.]+)/(?P<replicate>replicate_[0-9])/(?P<bam_dir>bam)/(?P<bam_file>[a-zA-Z0-9_\-\.]+$)"),
                                  "{basedir[0]}/{sample[0]}/{replicate[0]}/cufflinks/transcripts.gtf",
                                  "{basedir[0]}/{sample[0]}/{replicate[0]}/cufflinks",
                                  "{basedir[0]}/{sample[0]}/{replicate[0]}/qc",
                                   logger,logger_mutex)
def cufflinks(input_file, output_file, path, qc_path,logger, logger_mutex):
  bam=os.path.basename(input_file)
  cmd = ( "cd $TMPDIR \n"
          "mkdir reference \n"
          "cp {input_file} . \n"
          "samtools sort -@ 8 -m 2G {bam} temp \n"
          "cp $HOME/Scratch/reference/grch38/Homo_sapiens.GRCh38.dna.primary_assembly.fa* ./reference  \n"
          "cp $HOME/Scratch/reference/grch38/Hs.GRCh38.84.exon.gtf ./reference \n"
          "cp $HOME/Scratch/reference/grch38/ribosomal_mito_mask.gtf ./reference \n"
          "cufflinks -q -u --no-update-check -p 8 -G ./reference/Hs.GRCh38.84.exon.gtf \\\n"
          "-b ./reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa \\\n"
          "--mask-file ./reference/ribosomal_mito_mask.gtf temp.bam  \\\n"
          "-o  {path}  \\\n"
          "2>{qc_path}/cufflinks.log \n" )
  cmd = cmd.format(**locals())
  #print cmd
  try:
    stdout_res, stderr_res = "",""
    stdout_res, stderr_res = run_job(cmd,
                                     job_name = "cufflinks",
                                     job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                     job_other_options    = "-w n -S /bin/bash -V -l h_rt=04:00:00 -w n -l mem=4G -l tmpfs=60G -pe smp 8 -wd /home/sejjctj/Scratch -j yes ",
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
#              run QoRTs
#_______________________________________________________________________________________________________
@active_if(star_check or hisat_check)
@collate([hisat2,star],formatter("(?P<basedir>[/.].+)/(?P<sample>[a-zA-Z0-9_\-\.]+)/(?P<replicate>replicate_[0-9])/(?P<bam_dir>bam)/(?P<bam_file>[a-zA-Z0-9_\-\.]+$)"),
                            "{basedir[0]}/{sample[0]}/{replicate[0]}/qc/QC.summary.txt", 
                            "{basedir[0]}/{sample[0]}/{replicate[0]}/qc/qorts.log",
                             logger, logger_mutex )
def qorts(input_file, output_file, log_file, logger, logger_mutex):
    bam=os.path.basename(input_file[0])
    cmd = (" cd $TMPDIR; mkdir tmp \n"
           " cp {input_file[0]} ./ \n"
           " java -Xmx48G -Djava.io.tmpdir=./tmp \\\n"
           " -jar ~/applications/QoRTs/QoRTs.jar QC \\\n" 
           " --minMAPQ 60 \\\n"
           " --maxReadLength 85 \\\n"
           " {bam} \\\n"
           " ~/Scratch/reference/star_single_cell/Hs.GRCh38.84.exon.ercc.gtf \\\n"
           " {output_file} \\\n"
           " 2>{log_file} " )
    cmd = cmd.format(**locals())
    #print cmd
    try:
      stdout_res, stderr_res = "",""
      stdout_res, stderr_res = run_job(cmd,
                                     job_name = "qorts",
                                     job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                     job_other_options    = "-w n -S /bin/bash -V -l h_rt=08:00:00 -w n -l mem=48G -l tmpfs=30G -wd /home/sejjctj/Scratch -j yes ",
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
      logger.debug("qorts worked")



@active_if(kallisto_check)
@collate(trim_fastq,formatter("(?P<basedir>[/.].+)/(?P<sample>[a-zA-Z0-9_\-\.]+)/(?P<replicate>replicate_[0-9])/(?P<trimmed_dir>fastq_trimmed)/(?P<fastq>[a-zA-Z0-9_\-\.]+$)"),
                              "{basedir[0]}/{sample[0]}/{replicate[0]}/kallisto/abundance.tsv",
                              "{basedir[0]}/{sample[0]}/{replicate[0]}/{trimmed_dir[0]}",
                              "{basedir[0]}/{sample[0]}/{replicate[0]}/kallisto",
                              "{basedir[0]}/{sample[0]}/{replicate[0]}/qc" )
def kallisto(input_files, output_file, path,kallisto_folder,qc_folder):
    input_files = [item for sublist in input_files for item in sublist]
    list_of_reads = []
    for filename in input_files:
        list_of_reads.append(os.path.basename(filename))
    list_of_reads = ' '.join(list_of_reads)

    cmd = ("cd $TMPDIR \n"
           "mkdir reference \n"
           "cp {path}/*fq.gz   . \n"
           "cp $HOME/Scratch/reference/hg38_ver84_transcripts.idx ./reference \n"
           "kallisto quant -b 100 -t 4 -i \\\n"
           "./reference/hg38_ver84_transcripts.idx {list_of_reads} \\\n"
           "-o {kallisto_folder}")
    cmd = cmd.format(**locals())
    try:
      stdout_res, stderr_res = "",""
      stdout_res, stderr_res = run_job(cmd,
                                     job_name             = "kallisto",
                                     job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                     job_other_options    = "-w n -S /bin/bash -V -l h_rt=04:00:00 -l mem=8G -w n -pe smp 4 -l tmpfs=60G -wd /home/sejjctj/Scratch -j yes ",
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




if __name__ == '__main__':
  cmdline.run (options, multithread = options.jobs)
  drmaa_session.exit()

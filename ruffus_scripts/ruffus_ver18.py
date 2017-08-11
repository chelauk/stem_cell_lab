#!/usr/bin/env python

#######################################################################################
# An RNA seq pipeline developed by Chela James using the ruffus package               #
# http://www.ruffus.org.uk/                                                           #
# current version:                                                                    #
# * trim galore --hisat-----bam---QORTS (qc and DESEQ2 counts)                        # 
#                 |-STAR----|    |-cufflinks-fpkms-cuffdiff                           #
#                 |-kallisto-tsv(for sleuth)                                          #
# notes:                                                                              #                                                                    
# currently cuffdiff using the -G option for known transcripts from gtf file          #
#######################################################################################
from glob import glob
import glob
import sys, os, fnmatch
import re
from ruffus import *
import ruffus.cmdline as cmdline
from subprocess import check_call

parser = cmdline.get_argparse(description="Standard RNAseq Pipeline")

#parser.add_argument('-i', '--input', nargs='+', metavar="FILE", action="append", help = "Fastq files")
parser.add_argument('-i', '--input', metavar="FILE", help = "Fastq files")
#
#    Add arguments ruffus needs a more useful help file
#
parser.add_argument('--cuffdiff_file', metavar="FILE", help = "cuffdiff comparison instructions")
parser.add_argument('--basedir', metavar="DIR", help = "base directory")
parser.add_argument('--aligner', metavar="choice", help = "choice of aligner")
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
# working directory
#  standard python logger which can be synchronised across concurrent Ruffus tasks
logger, logger_mutex = cmdline.setup_logging ("Chela", options.log_file, options.verbose)

if not aligner: 
    raise Exception ("Aligner not selected with --aligner")
input_files=[]
files_list=options.input
# print files_list
if not files_list:
  raise Exception ("No matching files specified with --input.")
#  Here I need to get more sophisticated and get the user to input a regex
#  to parse and capture the fastq files and use the pattern later


pattern = re.compile(r"/(?P<base_dir>[a-zA-Z0-9 ]+?)(?P<sample>[a-zA-Z0-9 ]+?)/(?P<replicate>[a-zA-Z0-9 ]+?)/(?P<fastq_raw>[a-zA-Z0-9 ]+?)")
s = "/The Prodigy/The Fat Of The Land/04 - Funky Stuff.flac"
m = pattern.search(s)
print m.group('artist')
print m.group('release')
print m.group('track number')
print m.group('title')

with open(files_list, 'r') as f:
    content = [line.decode('utf-8').rstrip('\n') for line in f] 
    for line in content:
    #print(line)
        line = line.rstrip()
        #print(basedir + line)
        #print glob.glob(basedir + line + "/replicate*/fastq_raw/*gz")
        input_files.append(glob.glob(basedir + line + "/replicate*/fastq_raw/*gz"))  
input_files = [item for sublist in input_files for item in sublist]  
#print input_files
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
        "{subpath[0][1]}/bam",
        "{subpath[0][1]}/kallisto",
        "{subpath[0][1]}/cufflinks",
        "{subpath[0][3]}/cuffmerge_out",
        "{subpath[0][3]}/cuffdiff",
        "{subpath[0][1]}/cuffquant"])
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
    cmd = (" cd $TMPDIR ; "
         " cp {input_files[0]} . ;"
         " cp {input_files[1]} . ;"
         " trim_galore --fastqc --paired {basenames[0]} {basenames[1]} 2> {qc_folder}/trim_galore.log ; "
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

#_______________________________________________________________________________________________________
# 
#              take trimmer output and align with hisat2
#_______________________________________________________________________________________________________
## here we have a problem because this works for an eight-file sample (L001-4) but not for single end or
## two-file sample

@active_if(hisat_check)
@collate(trim_fastq, formatter("([^/]+)_L00[1234]_R[12]_001_val_[12].fq.gz$"),
                                "{subpath[0][1]}/bam/{1[0]}_L001_R1_001_val_1.fq.sorted.bam",
                                "{path[0]}","{subpath[0][1]}/bam",
                                "{subpath[0][1]}/qc",logger,logger_mutex)
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

    cmd = ( " cd $TMPDIR ; "
            " mkdir reference ; "
            " cp  {path}"  + "/*fq.gz" + " . ; "
            " cp $HOME/Scratch/reference/grch38_snp_tran/genome* ./reference ; "
            " hisat2 -p 8 -x ./reference/genome_snp_tran  --dta-cufflinks "
            " --novel-splicesite-outfile ./novel_splice.txt "
            " --novel-splicesite-infile ./novel_splice.txt "
            " -1 {first_reads} -2 {second_reads} 2> {qc_folder}/hisat.log | samtools view -bS - -o temp.bam ; "
            " samtools sort -n -@ 4 temp.bam -m 2G " + hisat_output[:-4] + " 2>{qc_folder}/samtools.log ; "
            " cp " + hisat_output + " {outpath} ; "
            " cp novel_splice.txt {outpath} ; "
            " rm -r * ; ")
    cmd = cmd.format(**locals())
    #print cmd
    try:
        stdout_res, stderr_res = "",""
        stdout_res, stderr_res = run_job(cmd,
                                        job_name = "hisat",
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
        logger.debug("hisat worked")
#_______________________________________________________________________________________________________
# 
#              take trimmer output and align with star
#_______________________________________________________________________________________________________
@active_if(star_check)
@collate(trim_fastq, formatter("([^/]+)_L00[1234]_R[12]_001_val_[12].fq.gz$"),
                              "{subpath[0][1]}/bam/{1[0]}.Aligned.sortedByCoord.out.bam",
                              "{path[0]}","{subpath[0][1]}/bam",
		              "{subpath[0][1]}/qc",
                               logger, logger_mutex)
def star(input_files, out_file, path, outpath,qc_folder,logger, logger_mutex):
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
  star_output = out_file.split('/')
  star_output = star_output[-1]
  #print star_output
  cmd = ( "cd $TMPDIR ; "
          "cp  {path}"  + "/*fq.gz" + " . ; "
          "~/applications/STAR-2.5.3a/bin/Linux_x86_64/STAR --runThreadN 4 "
          "--genomeDir ~/Scratch/reference/star_single_cell/index/ "
          "--readFilesIn " + first_reads + " " + second_reads + " "
          "--readFilesCommand zcat " 
          "--outFileNamePrefix " + star_output[:-35] + " "
          "--outSAMtype BAM SortedByCoordinate ; "  
          "cp *bam {outpath} ; "
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
# 
#              run QoRTs
#_______________________________________________________________________________________________________

'''


    with open(files_list, 'r') as f:
        content = [line.decode('utf-8').rstrip('\n') for line in f] 
        for line in content:
            #print(line)
            line = line.rstrip()
            #print(basedir + line)
            #print glob.glob(basedir + line + "/replicate*/bam/*bam")
            input_files.append(glob.glob(basedir + line + "/replicate*/bam/*bam"))  
    bam_list = [item for sublist in input_files for item in sublist]  
    return bam_list
'''
@collate([hisat2,star],formatter("([^/]+)bam$"),"{subpath[0][1]}/qc/QC.summary.txt", 
                            "{subpath[0][1]}/qc/qorts.log",
                             logger, logger_mutex )
def qorts(input_file, output_file, log_file, logger, logger_mutex):
    if hisat_check:
       gtf="~/Scratch/reference/grch38/Homo_sapiens.GRCh38.84.gtf "
    elif star_check:
       gtf="~/Scratch/reference/star_single_cell/Hs.GRCh38.84.exon.ercc.gtf "
    input_file = input_file[0]
    bam=os.path.basename(input_file)
    cmd = ("cd $TMPDIR \n"
           "mkdir tmp \n"
           "cp {input_file} ./ \n"
           "samtools sort -n {bam} -m 24G temp \n"
           "java -Xmx24G -Djava.io.tmpdir=./tmp \\\n"
           "-jar ~/applications/QoRTs/QoRTs.jar QC \\\n" 
           "--minMAPQ 60 \\\n"
           "--maxReadLength 76 \\\n"
           "--keepMultiMapped \\\n"
           "temp.bam \\\n"
           "{gtf}" + " " + "{output_file}" )
    cmd = cmd.format(**locals())
    #print cmd
    try:
      stdout_res, stderr_res = "",""
      stdout_res, stderr_res = run_job(cmd,
                                     job_name = "qorts",
                                     job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                     job_other_options    = "-w n -S /bin/bash -V -l h_rt=05:00:00 -w n -l mem=24G -l tmpfs=30G -wd /home/sejjctj/Scratch -j yes ",
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
@collate(trim_fastq,formatter("([^/]+)_L00[1234]_R[12]_001_val_[12].fq.gz$"),
                              "{subpath[0][1]}/kallisto/abundance.tsv",
                              "{path[0]}",
                              "{subpath[0][1]}/kallisto",
                              "{subpath[0][1]}/qc" )
def kallisto(input_files, output_file, path,kallisto_folder,qc_folder):
    input_files = [item for sublist in input_files for item in sublist]
    list_of_reads = []
    for filename in input_files:
        list_of_reads.append(os.path.basename(filename))
    list_of_reads = ' '.join(list_of_reads)

    cmd = (" cd $TMPDIR; mkdir reference; "
         " cp {path}/*fq.gz   . ; "
         " cp $HOME/Scratch/reference/hg38_ver84_transcripts.idx ./reference ; "
         " kallisto quant -b 100 -t 4 -i ./reference/hg38_ver84_transcripts.idx " + list_of_reads +
         " -o {kallisto_folder} ; "
         " rm -r * ;")
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



#_______________________________________________________________________________________________________
#
#              cufflinks generate gtf files from sorted hisat output
#_______________________________________________________________________________________________________
################# USING -G flag for cufflinks, no novel Isoforms ######################################
'''
@active_if(hisat_check)
@transform(hisat2,formatter("([^/]+)bam$"), "{subpath[0][1]}/cufflinks/transcripts.gtf","{1[0]}bam", "{subpath[0][1]}/cufflinks","{subpath[0][1]}/qc",logger,logger_mutex)
def cufflinks(input_file, output_file, cuff_input, path, qc_path,logger, logger_mutex):
  print input_file

  cmd = ( " cd $TMPDIR ; "
          " mkdir reference ; "
          " cp {input_file} . ; "
          " samtools sort -@ 8 -m 2G {cuff_input} temp ;"
          " cp $HOME/Scratch/reference/grch38/Homo_sapiens.GRCh38.dna.primary_assembly.fa* ./reference  ; "
          " cp $HOME/Scratch/reference/grch38/Hs.GRCh38.84.exon.gtf ./reference ; "
          " cp $HOME/Scratch/reference/grch38/ribosomal_mito_mask.gtf ./reference ; "
          " cufflinks -q -u --no-update-check -p 8 -G ./reference/Hs.GRCh38.84.exon.gtf "
          " -b ./reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa "
          " --mask-file ./reference/ribosomal_mito_mask.gtf temp.bam  -o  {path}  2>{qc_path}/cufflinks.log ; "
          " rm -r * ; " )
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

'''
#_______________________________________________________________________________________________________
#
#              cufflinks generate gtf files from sorted hisat output
#_______________________________________________________________________________________________________
################# USING -G flag for cufflinks, no novel Isoforms ######################################
@active_if(star_check or hisat_check)
@transform([hisat2,star],formatter("([^/]+)bam$"), "{subpath[0][1]}/cufflinks/transcripts.gtf","{1[0]}bam", "{subpath[0][1]}/cufflinks","{subpath[0][1]}/qc",logger,logger_mutex)
def cufflinks(input_file, output_file, cuff_input, path, qc_path,logger, logger_mutex):
  print input_file

  cmd = ( " cd $TMPDIR ; "
          " mkdir reference ; "
          " cp {input_file} . ; "
          " samtools sort -@ 8 -m 2G {cuff_input} temp ;"
          " cp $HOME/Scratch/reference/grch38/Homo_sapiens.GRCh38.dna.primary_assembly.fa* ./reference  ; "
          " cp $HOME/Scratch/reference/grch38/Hs.GRCh38.84.exon.gtf ./reference ; "
          " cp $HOME/Scratch/reference/grch38/ribosomal_mito_mask.gtf ./reference ; "
          " cufflinks -q -u --no-update-check -p 8 -G ./reference/Hs.GRCh38.84.exon.gtf "
          " -b ./reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa "
          " --mask-file ./reference/ribosomal_mito_mask.gtf temp.bam  -o  {path}  2>{qc_path}/cufflinks.log ; "
          " rm -r * ; " )
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



if __name__ == '__main__':
  cmdline.run (options, multithread = options.jobs)
  drmaa_session.exit()

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
parser.add_argument('--gtf', metavar="choice", help = "choice of gtf; enter all_transcripts, all_coding or ercc")
#parser.add_argument('--flat_gff', metavar="choice", help = "flat gff for DEXseq")
parser.add_argument('--kallisto', metavar="choice", help = "use kallisto?")
parser.add_argument('--species', metavar="choice", help = "species" )
parser.add_argument('--stringtie', metavar="choice", help = "use stringtie?")

options = parser.parse_args()
cuffdiff_file = options.cuffdiff_file
basedir=options.basedir
aligner=options.aligner
kallisto=options.kallisto
stringtie=options.stringtie
hisat_check=aligner=="hisat"
star_check=aligner=="star"
kallisto_check=kallisto=="yes"
stringtie_check=stringtie=="yes"
#flat_gff=options.flat_gff
gtf=options.gtf
species=options.species

## could be done better with a dictionary maybe

if species == "human" and gtf=="exon" and hisat_check and not stringtie_check:
	flat_gff="$HOME/Scratch/reference/grch38/Hs.GRCh38.84.exon.flattened.gff"
        gtf="$HOME/Scratch/reference/grch38/Hs.GRCh38.84.exon.gtf"
	hisat_genome_index="$HOME/Scratch/reference/grch38_snp_tran/"
        genome="$HOME/Scratch/reference/grch38/"
	genome_name="Homo_sapiens.GRCh38.dna.primary_assembly.fa"
	mask="$HOME/Scratch/reference/grch38/ribosomal_mito_mask.gtf"
elif species == "human" and gtf=="all_transcripts" and stringtie_check:
        flat_gff="$HOME/Scratch/reference/grch38/gencode.flattened.gff"
	gtf="$HOME/Scratch/reference/grch38/gencode.v28.annotation.gtf"
	hisat_genome_index="$HOME/Scratch/reference/grch38_snp_tran/"
        genome="$HOME/Scratch/reference/grch38/"
	genome_name="Homo_sapiens.GRCh38.dna.primary_assembly.fa"
	mask="$HOME/Scratch/reference/grch38/ribosomal_mito_mask.gtf"
elif species == "human" and gtf=="all_transcripts" and not stringtie_check:
        flat_gff="$HOME/Scratch/reference/grch38/gencode.flattened.gff"
	gtf="$HOME/Scratch/reference/grch38/gencode.v28.annotation.gtf"
	hisat_genome_index="$HOME/Scratch/reference/grch38_snp_tran/"
        genome="$HOME/Scratch/reference/grch38/"
	genome_name="Homo_sapiens.GRCh38.dna.primary_assembly.fa"
	mask="$HOME/Scratch/reference/grch38/ribosomal_mito_mask.gtf"
elif species == "human" and gtf=="all_coding" and hisat_check:
	gtf="$HOME/Scratch/reference/grch38/Hs.GRCh38.84.protein_coding.gtf"
	hisat_genome_index="$HOME/Scratch/reference/grch38_snp_tran/"
	genome="$HOME/Scratch/reference/grch38/"
	genome_name="Homo_sapiens.GRCh38.dna.primary_assembly.fa"
        mask="$HOME/Scratch/reference/grch38/ribosomal_mito_mask.gtf"
elif species == "human" and gtf=="ercc" and star_check:
	gtf="$HOME/Scratch/reference/star_single_cell/Hs.GRCh38.84.exon.ercc.gtf"
        genome="$HOME/Scratch/reference/grch38/"
	mask="$HOME/Scratch/reference/grch38/ribosomal_mito_mask.gtf"
elif species == "mouse" and hisat_check:
	gtf="$HOME/Scratch/reference/GRcm38/gtf/Mus_musculus.GRCm38.84.gtf"
	hisat_genome_index="$HOME/Scratch/reference/GRcm38/hisat2/grcm38_snp_tran"
        genome="$HOME/Scratch/reference/GRcm38/GRCm38/Mus_musculus.GRCm38"
        genome_name="Mus_musculus.GRCm38.dna.primary_assembly.fa"
	mask="$HOME/Scratch/reference/GRcm38/gtf/Mus_musculus.GRCm38.84.ribo.mito.mask.gtf"

if hisat_check:
   print "hisat check: " + str(hisat_check)
   print(str(hisat_genome_index))
if stringtie_check:
   print "stringtie check: " + str(stringtie_check)
print "STAR check: " + str(star_check)
print gtf
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
        # print "input " + str(glob.glob(basedir + line + "/replicate*/fastq_raw/*gz"))
        input_files.append(glob.glob(basedir + line + "/replicate*/fastq_raw/*gz"))  
#                                                                                .
#   Useful code to turn input files into a flat list                             .
#     

input_files = [item for sublist in input_files for item in sublist]  

# for @mkdir 
my_dirs = []
with open(files_list, 'r') as f:
    for line in f:
        #print(line)
        line = line.rstrip()
        print "Dirs: " + str(basedir + line)
        my_dirs.append(glob.glob(basedir + line + "/replicate*"))  

my_dirs = [item for sublist in my_dirs for item in sublist]

print "Dirs:  " +str(my_dirs)
# start drmaa_session
drmaa_session = drmaa.Session()
drmaa_session.initialize()



#   <<<----  pipelined functions go here
#_________________________________________________________________________________
#                                                                                .
#   The first step trims pairs of fastqs and generates fastqc reports            .
#_________________________________________________________________________________


@mkdir(my_dirs,
       #match input pattern
        formatter("(?P<basedir>[/.].+)/(?P<sample>[a-zA-Z0-9_\-\.]+)/(?P<replicate>replicate_[0-9]+)"),
        # make qc directory
        ["{basedir[0]}/{sample[0]}/{replicate[0]}/qc",
        # make trimmed_directory)
        "{basedir[0]}/{sample[0]}/{replicate[0]}/fastq_trimmed",
        "{basedir[0]}/{sample[0]}/{replicate[0]}/bam",
        "{basedir[0]}/{sample[0]}/{replicate[0]}/kallisto",
        "{basedir[0]}/{sample[0]}/{replicate[0]}/cufflinks",
        "{basedir[0]}/{sample[0]}/{replicate[0]}/stringtie",
        "{basedir[0]}/cuffmerge_out",
        "{basedir[0]}/cuffdiff",
        "{basedir[0]}/{sample[0]}/{replicate[0]}/cuffquant"])

@collate(input_files,
        # input formatter to provide read pairs
        formatter("(?P<basedir>[/.].+)/(?P<sample>[a-zA-Z0-9_\-\.]+)/(?P<replicate>replicate_[0-9]+)/fastq_raw/(?P<pair>.+)R[12](?P<fill>.*).fastq.gz"),
        # create output parameter to be supplied to next task
	#["{basedir[0]}/{sample[0]}/{replicate[0]}/fastq_trimmed/{pair[0]}R1_001_val_1.fq.gz",
	# "{basedir[0]}/{sample[0]}/{replicate[0]}/fastq_trimmed/{pair[0]}R2_001_val_2.fq.gz"],
	["{basedir[0]}/{sample[0]}/{replicate[0]}/fastq_trimmed/{pair[0]}R1{fill[0]}_val_1.fq.gz",
	 "{basedir[0]}/{sample[0]}/{replicate[0]}/fastq_trimmed/{pair[0]}R2{fill[0]}_val_2.fq.gz"],
        #["{pair[0]}R1_001.fastq.gz","{pair[0]}R2_001.fastq.gz"],    # basename for trim_galore
        "{basedir[0]}/{sample[0]}/{replicate[0]}/qc",               # qc folder
        "{basedir[0]}/{sample[0]}/{replicate[0]}/fastq_trimmed",    # trimmed_folder
        logger, logger_mutex)
#def trim_fastq(input_files, output_files, basenames, qc_folder, output_folder ,logger, logger_mutex):
def trim_fastq(input_files, output_files, qc_folder, output_folder ,logger, logger_mutex):
    print "OUTPUT FILES!   " + str(output_files)    
    if len(input_files) !=2:
        raise Exception("One of the reads pairs %s missing" % (input_files,))
    cmd = ( " source ~/.bashrc \n"
            " date \n"
            " echo $HOSTNAME \n"
            " cd $TMPDIR \n"
            " cp {input_files[0]} . \n"
            " cp {input_files[1]} . \n"
            " basename1=$(basename {input_files[0]}) \n"
            " basename2=$(basename {input_files[1]}) \n"
            " date \n"
            " ls -l \n"
            #" trim_galore --fastqc --paired {basenames[0]} {basenames[1]} &> {qc_folder}/trim_galore.log \n"
            " trim_galore --fastqc --paired $basename1 $basename2 &> {qc_folder}/trim_galore.log \n"
            " mv *.fq.gz  {output_folder} \n"
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
                                      job_other_options    = "-w n -S /bin/bash -l h_rt=05:00:00 -l mem=4G -l tmpfs=60G -wd /home/sejjctj/Scratch -j yes",
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
#              take trimmer output and align with hisat2
#_______________________________________________________________________________________________________
@active_if(hisat_check)
@collate(trim_fastq, formatter("(?P<basedir>[/.].+)/(?P<sample>[a-zA-Z0-9_\-\.]+)/(?P<replicate>replicate_[0-9])/(?P<fastq_trimmed>fastq_trimmed)/(?P<trimmed_fq>[a-zA-Z0-9_\-\.]+)"),
                               "{basedir[0]}/{sample[0]}/{replicate[0]}/bam/{sample[0]}.{replicate[0]}.sorted.bam",
                               "{basedir[0]}/{sample[0]}/{replicate[0]}/{fastq_trimmed[0]}",
                               "{basedir[0]}/{sample[0]}/{replicate[0]}/bam",
                               "{basedir[0]}/{sample[0]}/{replicate[0]}/qc",
                               hisat_genome_index,
                                logger,logger_mutex)
def hisat2(input_files, out_file, path, outpath,qc_folder,hisat_genome_index,logger, logger_mutex):
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
    
    cmd = ( "source ~/.bashrc \n"
            "cd $TMPDIR \n"
            "mkdir reference \n"
            "cp  {path}/*fq.gz  . \n"
            "cp {hisat_genome_index}/genome* ./reference \n"
            "hisat2 -p 8 -x ./reference/genome_snp_tran  --dta-cufflinks \\\n"
            "--novel-splicesite-outfile ./novel_splice.txt \\\n"
            "--novel-splicesite-infile ./novel_splice.txt \\\n"
            "-1 {first_reads} \\\n"
            "-2 {second_reads} \\\n"
            "2> {qc_folder}/hisat.log | samtools view -bS - -o temp.bam \n"
            "samtools sort -@ 8 temp.bam -m 4G " + hisat_output[:-4] + " 2>{qc_folder}/samtools.log \n"
            "mv {hisat_output} {outpath} \n"
            "mv novel_splice.txt {outpath} \n")
    cmd = cmd.format(**locals())
    try:
        stdout_res, stderr_res = "",""
        stdout_res, stderr_res = run_job(cmd,
                                        job_name = "hisat",
                                        job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                        job_other_options    = "-w n -S /bin/bash -l h_rt=08:00:00 -w n -l mem=4G -l tmpfs=60G -pe smp 8 -wd /home/sejjctj/Scratch -j yes ",
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
        logger.debug("hisat worked")
#_______________________________________________________________________________________________________
# 
#              take trimmer output and align with star
#_______________________________________________________________________________________________________
@active_if(star_check)
@collate(trim_fastq, formatter("(?P<basedir>[/.].+)/(?P<sample>[a-zA-Z1-9_\-\.]+)/(?P<replicate>replicate_[0-9])/(?P<trimmed_dir>fastq_trimmed)/(?P<fastq>[a-zA-Z0-9_\-\.]+$)"),
                              "{basedir[0]}/{sample[0]}/{replicate[0]}/bam/{sample[0]}.Aligned.sortedByCoord.out.bam", #out_file  ~1
                              "{basedir[0]}/{sample[0]}/{replicate[0]}/{trimmed_dir[0]}",                               #path      ~2
                              "{basedir[0]}/{sample[0]}/{replicate[0]}/bam",                                            #outpath   ~3
		              "{sample[0]}.",                                                                            #sample    ~4
                              "{basedir[0]}/{sample[0]}/{replicate[0]}/qc",                                             #qc_folder ~5 
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
  cmd = ( "source ~/.bashrc \n"
          "cd $TMPDIR \n "
          "cp {path}/*fq.gz  . \n "
          "STAR --runThreadN 4 \\\n"
          "--genomeDir ~/Scratch/reference/star_single_cell/index/ \\\n"
          "--readFilesIn " + first_reads + " " + second_reads + " \\\n"
          "--readFilesCommand zcat \\\n" 
          "--twopassMode Basic \\\n" 
          "--outReadsUnmapped None \\\n" 
          "--chimSegmentMin 12 \\\n" 
          "--chimJunctionOverhangMin 12 \\\n"  
          "--alignSJDBoverhangMin 10 \\\n" 
          "--alignMatesGapMax 100000 \\\n" 
          "--alignIntronMax 100000 \\\n" 
          "--chimSegmentReadGapMax 3 \\\n"                                                                                     
          "--alignSJstitchMismatchNmax 5 -1 5 5 \\\n" 
          "--outSAMstrandField intronMotif \\\n"
          "--outFilterIntronMotifs RemoveNoncanonical \\\n" ## added for compatibility with
          "--outFileNamePrefix {sample} \\\n"               ## cufflinks
          "--outSAMtype BAM SortedByCoordinate\n"
          "cp *junction {outpath} \n"
          "cp *bam {outpath} \n"
          "cp *Log.* {qc_folder} ")
  cmd = cmd.format(**locals())
  print cmd
  try:
    stdout_res, stderr_res = "",""
    stdout_res, stderr_res = run_job(cmd,
                                     job_name = "star",
                                     job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                     job_other_options    = "-w n -S /bin/bash -l h_rt=02:00:00 -w n -l mem=24G -l tmpfs=60G -pe smp 4 -wd /home/sejjctj/Scratch -j yes ",
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
    logger.debug("star worked")

#####################################################################################################
#___________________________________________________________________________________________________#
#                          automatic fusion detection if star is used (I might make that default)   #
#####################################################################################################

@active_if(star_check)
@transform(star,formatter("(?P<basedir>[/.].+)/(?P<sample>[a-zA-Z0-9_\-\.]+)/(?P<replicate>replicate_[0-9])/(?P<bam_dir>bam)/(?P<bam_file>[a-zA-Z0-9_\-\.]+$)"),
                          add_inputs("{basedir[0]}/{sample[0]}/{replicate[0]}/bam/*.Chimeric.out.junction"),
                          "{basedir[0]}/{sample[0]}/{replicate[0]}/bam/star-fusion.fusion_predictions.abridged.tsv",
		          "{sample[0]}.",                                                                            #sample    ~4
                          "{basedir[0]}/{sample[0]}/{replicate[0]}/bam",
                          "{basedir[0]}/{sample[0]}/{replicate[0]}/qc",
                          logger, logger_mutex)
def star_fusion(input_files, out_file,sample, outpath,qc_folder,logger, logger_mutex):
  fusion_input = input_files[1]
  fusion_name = os.path.basename(fusion_input)
  cmd = ( "source ~/.bashrc \n"
          "module unload perl \n"
          "module load perl/5.16.0 \n"
          "export PERL5LIB=$PERL5LIB:/home/sejjctj/perl5/lib/perl5 \n"
          "cd $TMPDIR \n "
          "cp {fusion_input}  . \n "
          "awk 'BEGIN{{OFS=\"\\t\"}}{{$1=\"chr\"$1;$4=\"chr\"$4;print $0}}' {fusion_name} > temp && mv temp {fusion_name} \n"
          "STAR-Fusion \\\n"
          "--genome_lib_dir /home/sejjctj/Scratch/reference/star_single_cell/fusion/GRCh38_v27_CTAT_lib_Feb092018/ctat_genome_lib_build_dir \\\n"
          "-J {fusion_name} \\\n"
          "--output_dir {outpath} \n" )

  cmd = cmd.format(**locals())
  #print cmd
  try:
    stdout_res, stderr_res = "",""
    stdout_res, stderr_res = run_job(cmd,
                                     job_name = "star_fusion",
                                     job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                     job_other_options    = "-w n -S /bin/bash -l h_rt=02:00:00 -w n -l mem=24G -l tmpfs=60G  -wd /home/sejjctj/Scratch -j yes ",
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
    logger.debug("star_fusion worked")


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
                                   gtf,genome,mask,genome_name,
                                   logger,logger_mutex)
def cufflinks(input_file, output_file, path,qc_path,gtf,genome,mask,genome_name,logger, logger_mutex):
  bam=os.path.basename(input_file)
  my_mask=os.path.basename(mask)
  cmd = ( "source ~/.bashrc \n"
          "cd $TMPDIR \n"
          "mkdir reference \n"
          "cp {input_file} . \n"
          "cp {genome}*fa* ./reference  \n"
          "cp {gtf} ./reference/gencode.gtf \n"
          "cp {mask} ./reference \n"
          "cufflinks -q -u --no-update-check -p 8 -G ./reference/gencode.gtf \\\n"
          "-b ./reference/{genome_name} \\\n"
          "--mask-file ./reference/{my_mask} {bam} \\\n"
          "-o  {path}  \\\n"
          "2>{qc_path}/cufflinks.log \n" )
  cmd = cmd.format(**locals())
  #print cmd
  try:
    stdout_res, stderr_res = "",""
    stdout_res, stderr_res = run_job(cmd,
                                     job_name = "cufflinks",
                                     job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                     job_other_options    = "-w n -S /bin/bash -l h_rt=04:00:00 -w n -l mem=4G -l tmpfs=60G -pe smp 8 -wd /home/sejjctj/Scratch -j yes ",
                                     #job_environment      = { 'BASH_ENV' : '/home/sejjctj/.bashrc' } ,
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
#         stringtie generates gtf files, abundance files and ctab files from sorted hisat/star output
#_______________________________________________________________________________________________________


@active_if(stringtie_check)
@transform([hisat2,star],formatter("(?P<basedir>[/.].+)/(?P<sample>[a-zA-Z0-9_\-\.]+)/(?P<replicate>replicate_[0-9])/(?P<bam_dir>bam)/(?P<bam_file>[a-zA-Z0-9_\-\.]+$)"),
                                  "{basedir[0]}/{sample[0]}/{replicate[0]}/stringtie/{sample[0]}.gtf",
                                  "{basedir[0]}/{sample[0]}/{replicate[0]}/stringtie/{sample[0]}.gene_abund.tab",
                                  "{basedir[0]}/{sample[0]}/{replicate[0]}/qc",
                                   gtf,logger,logger_mutex)
def stringtie(input_file, output_file,abundance_file,qc_path,gtf,logger, logger_mutex):
  bam=os.path.basename(input_file)
  cmd = ( "source ~/.bashrc \n"
          "cd $TMPDIR \n"
          "mkdir reference \n"
          "cp {input_file} . \n"
          "cp {gtf} ./reference/gencode.gtf \n"
          "stringtie -p 8 -G ./reference/gencode.gtf -A {abundance_file} -o {output_file} -B -e -v {bam} \\\n"
          "2>{qc_path}/stringtie.log \n" )
  cmd = cmd.format(**locals())
  #print cmd
  try:
    stdout_res, stderr_res = "",""
    stdout_res, stderr_res = run_job(cmd,
                                     job_name = "stringtie",
                                     job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                     job_other_options    = "-w n -S /bin/bash -l h_rt=04:00:00 -w n -l mem=4G -l tmpfs=60G -pe smp 8 -wd /home/sejjctj/Scratch -j yes ",
                                     #job_environment      = { 'BASH_ENV' : '/home/sejjctj/.bashrc' } ,
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
    logger.debug("stringtie worked")


#_______________________________________________________________________________________________________
# 
#              run QoRTs
#_______________________________________________________________________________________________________
@active_if(star_check or hisat_check)
@collate([hisat2,star],formatter("(?P<basedir>[/.].+)/(?P<sample>[a-zA-Z0-9_\-\.]+)/(?P<replicate>replicate_[0-9])/(?P<bam_dir>bam)/(?P<bam_file>[a-zA-Z0-9_\-\.]+$)"),
                            "{basedir[0]}/{sample[0]}/{replicate[0]}/qc/QC.summary.txt", 
                            "{basedir[0]}/{sample[0]}/{replicate[0]}/qc/qorts.log",
                             gtf,
                             logger, logger_mutex )
def qorts(input_file, output_file, log_file, gtf, logger, logger_mutex):
    bam=os.path.basename(input_file[0])
    cmd = (" source ~/.bashrc \n"
           " cd $TMPDIR; mkdir tmp \n"
           " cp {input_file[0]} ./ \n"
           " samtools sort -n -m 12G -T prefix -O bam {bam} > namesort.bam \n"
           " java -Xmx48G -Djava.io.tmpdir=./tmp \\\n"
           " -jar ~/applications/QoRTs/QoRTs.jar QC \\\n" 
           " --nameSorted \\\n"
           " --minMAPQ 60 \\\n"
           " --maxReadLength 100 \\\n"
           " namesort.bam \\\n"
           " {gtf} \\\n"
           " {output_file} \\\n"
           " 2>{log_file} " )
    cmd = cmd.format(**locals())
    #print cmd
    try:
      stdout_res, stderr_res = "",""
      stdout_res, stderr_res = run_job(cmd,
                                     job_name = "qorts",
                                     job_script_directory = "/home/sejjctj/Scratch/test_dir",
                                     job_other_options    = "-w n -S /bin/bash  -l h_rt=08:00:00 -w n -l mem=48G -l tmpfs=30G -wd /home/sejjctj/Scratch -j yes ",
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

    cmd = ("source ~/.bashrc \n"
           "cd $TMPDIR \n"
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
                                     job_other_options    = "-w n -S /bin/bash -l h_rt=04:00:00 -l mem=8G -w n -pe smp 4 -l tmpfs=60G -wd /home/sejjctj/Scratch -j yes ",
                                     #job_environment      = { 'BASH_ENV' : '/home/sejjctj/.bashrc' } ,
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
  pipeline_printout_graph ("bulk_rna-seq.jpg", "jpg", [trim_fastq,hisat2,star,kallisto,cufflinks,qorts],
                          no_key_legend=True,
                          ignore_upstream_of_target=True,
                          pipeline_name="bulk RNA-seq",
                          user_colour_scheme = {
                                                "colour_scheme_index" :2,
                                                "Bulk RNA-seq"      :{"fontcolor" : '"#FF3232"' },
                                                "Task to run"       :{"linecolor" : '"#0044A0"' },
                                                "Final target"      :{"fillcolor" : '"#EFA03B"',
                                                                       "fontcolor" : "black",
                                                                       "dashed"    : 0           }
                                               })
  pipeline_printout()


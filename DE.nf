#! /usr/bin/env nextflow

params.help = null

log.info ""
log.info "--------------------------------------------------------------------"
log.info "  DE 1.0 : Pipeline RNAseq for the Differential Expression analysis."
log.info "--------------------------------------------------------------------"

if (params.help) {
    log.info "------------------------------------------------------------------------------------------------------------------------------"
    log.info "  USAGE : nextflow run Lipinski-B/DE-nf --input /data/ --GTF /data/fichier.gtf --FNA /data/fichier.fna --output /output/ "
    log.info "------------------------------------------------------------------------------------------------------------------------------"
    log.info ""
    log.info "nextflow run Lipinski-B/DE-nf [-r vX.X -profile docker/singularity] [OPTIONS]"
    log.info ""
    log.info "Mandatory arguments:"
    log.info ""
    log.info "--input                       FOLDER                      Folder where you can find your data (fasta/fastq files)."
    log.info "--output                      FOLDER                      Folder where you want to find your result."
    log.info "--GTF                         FILE                        Path where you can find the annotation to use."
    log.info ""
    log.info "Optional arguments:"
    log.info "--STAR_Index                  FOLDER                      Folder where you can find the STAR index. If this option is not used, please make sure to provide the --FNA option in addition to the --GTF option to perform the STAR index"
    log.info "--FNA                         FILE                        Path where you can find the FNA file to use for the STAR index."
    log.info "--R                           STRING                      on/off : Chose to use or not the standard R analyses from the pipeline."
    log.info "--metadata                    FILE                        Path where you can find the XLS file to use as metadata for the R analyse. Mandatory is the option --R in on."
    log.info "--thread                      INT                         Number of thread to use."
  
    exit 0
} else {
    /* Software information */
    log.info "help:                               ${params.help}"
}

// -- Path :
params.input = null
params.output = null
params.GTF = null

// -- Option :
params.R = "off"
params.thread = 1
params.STAR_Index = null
params.FNA = params.STAR_Index
params.metadata = null
//params.metadata = "!{baseDir}/data/Metadata.xls"

// -- Pipeline :
process MultiQC{ 
  publishDir params.output+'/QC/', mode: 'copy'
  
  input:
  file data from Channel.fromPath(params.input+'*').collect()

  //output:
  //file "*fastqc.html" into result_QC1
  //file "multiqc*" into result_QC2
  
  shell:
  '''
  #Multi QC analysis
  files=(*)
  #for file in *; do
    #fastqc $file
  #done
  #multiqc .
  '''}



process Mapping{ 
  publishDir params.output+'/mapping/', mode: 'copy'
  cpus params.thread
  
  input:
  file data from Channel.fromPath(params.input+'*').collect()
  file GTF from Channel.fromPath(params.GTF).collect()
  file FNA from Channel.fromPath(params.FNA).collect()

  output:
  file "*Aligned.sortedByCoord.out.bam" into Mapping_bam
  file "other/" into Mapping_Log
  
  shell:
  if(params.STAR_Index==null) {
    '''
    ## -- Index construction -------------------------------------------------------------- ##
    mkdir STARIndex_last/
    STAR --runThreadN !{params.thread} \
      --runMode genomeGenerate --genomeDir STARIndex_last/ --genomeFastaFiles !{FNA} \
      --sjdbGTFfile !{GTF} --sjdbOverhang 149 --genomeSAsparseD 2

    mkdir data/ ; mv *gz data/ ; cd data/ 

    ## -- Mapping analyse ----------------------------------------------------------------- ## --limitBAMsortRAM 127598325157
    #ulimit -v 127598325157
    list=`ls -1 | sed 's/_R.*//' | uniq`
    for file in $list; do
      STAR --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN !{params.thread} \
      --runThreadN !{params.thread} --genomeDir ../STARIndex_last --readFilesCommand gunzip -c \
      --readFilesIn $file'_R1.fastq.gz' $file'_R2.fastq.gz' \
      --outFileNamePrefix $file --outSAMunmapped Within
    done

    mv * ../. ; cd .. ; rm -r data/
    mkdir other ; mv STARIndex_last/ other/ ; mv *Log* other/
    '''
  } else {
    '''
    mkdir data/ ; mv *gz data/ ; cd data/ 

    ## -- Mapping analyse ----------------------------------------------------------------- ## --limitBAMsortRAM 127598325157
    #ulimit -v 17598325157
    list=`ls -1 | sed 's/_R.*//' | uniq`
    for file in $list; do
      STAR --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN !{params.thread} \
      --runThreadN !{params.thread} --genomeDir ../!{FNA} --readFilesCommand gunzip -c \
      --readFilesIn $file'_R1.fastq.gz' $file'_R2.fastq.gz' \
      --outFileNamePrefix $file --outSAMunmapped Within
    done
    
    mv * ../. ; cd .. ; rm -r data/
    mkdir other ; mv *Log* other/
    '''
    }}


process Intersection{ 
  publishDir params.output+'/intersect/', mode: 'copy'
  cpus params.thread

  input:
  file data from Mapping_bam
  file GTF from Channel.fromPath(params.GTF).collect()
  
  output:
  file "*.txt" into Intersect
  
  shell:
  '''
  ## -- Intersection analyse ----------------------------------------------- ##
  for file in *.bam; do
    htseq-count -f bam -n !{params.thread} $file !{GTF} > ${file}_intersect.txt
  done
  '''}

process Merge_result{ 
  publishDir params.output+'/merge/', mode: 'copy'
  
  input:
  file data from Intersect
  
  output:
  file "finale.txt" into Result
  
  shell:
  '''
  ## -- Differancial expression analyse ------------------------------------------------- ##
  files=(*)
  awk '{print $1}' ${files[0]} > AAAA.txt

  for file in *_intersect.txt; do
    awk '{print $2}' $file > ${file}.tmp
    rm $file
    mv ${file}.tmp $file 
    tail -n +2 "$file" > "$file.tmp" && mv "$file.tmp" "$file"
    echo "${file%%.*}" > ${file}.name
    cat ${file}.name ${file} > ${file}.tmp && mv ${file}.tmp ${file}
    rm ${file}.name
  done

  paste -d "\t" * > finale.txt
  rm AAAA.txt
  '''}


if(params.R=="on"){
  process DEA{ 
    publishDir params.output+'/R/', mode: 'copy'
    
    input:
    file data from Result.collect()
    file metadata from Channel.fromPath(params.metadata).collect()
    
    output:
    file "*.pdf" into Result_DE
    
    shell:
    '''
    Rscript !{baseDir}/bin/DE.r finale.txt !{metadata} !{baseDir}/bin/grid_0.7-4.tar.gz !{baseDir}/bin/gridExtra_2.3.tar.gz
    '''
    }}

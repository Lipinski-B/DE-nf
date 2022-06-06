 process Mapping_STAR{ 
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
    '''
    list=`ls *fast[q.a]* -1 | sed 's/_R.*//' | uniq`

    ## -- Index construction -------------------------------------------------------------- ##
    if [ "!{params.index}" == "null" ] ; then
      mkdir STARIndex_last/
      STAR --runThreadN !{params.thread} \
        --runMode genomeGenerate --genomeDir STARIndex_last/ --genomeFastaFiles !{FNA} \
        --sjdbGTFfile !{GTF} --sjdbOverhang 149 --genomeSAsparseD 2
      IDX=STARIndex_last
    else
      IDX=!{FNA}
    fi

    ## -- Mapping analyse ----------------------------------------------------------------- ##
    for file in $list; do
      STAR --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN !{params.thread} \
      --runThreadN !{params.thread} --genomeDir $IDX --readFilesCommand gunzip -c \
      --readFilesIn $file'_R1.fastq.gz' $file'_R2.fastq.gz' \
      --outFileNamePrefix $file --outSAMunmapped Within
    done

    mkdir other

    if [ "!{params.index}" == "null" ] ; then
        mv STARIndex_last/ other/
    fi

    mv *Log* other/

    '''
  }
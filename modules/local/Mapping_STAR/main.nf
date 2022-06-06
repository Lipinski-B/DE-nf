process Mapping_STAR {
    publishDir params.output+'/mapping/', mode: 'copy'
    cpus params.thread
    
    input:
        tuple val(ID), file(R)

    output:
        tuple file "*Aligned.sortedByCoord.out.bam", file "other/", emit: Mapping_STAR
    
    script:
    """
    list=`ls *fast[q.a]* -1 | sed 's/_R.*//' | uniq`

    ## -- Index construction -------------------------------------------------------------- ##
    if [ "${params.index}" == "null" ] ; then
      mkdir STARIndex_last/
      STAR --runThreadN ${params.thread} \
        --runMode genomeGenerate --genomeDir STARIndex_last/ --genomeFastaFiles ${params.FNA} \
        --sjdbGTFfile ${params.GTF} --sjdbOverhang 149 --genomeSAsparseD 2
      IDX=STARIndex_last
    else
      IDX=${params.FNA}
    fi

    ## -- Mapping analyse ----------------------------------------------------------------- ##
    for file in $list; do
      STAR --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN ${params.thread} \
      --runThreadN ${params.thread} --genomeDir $IDX --readFilesCommand gunzip -c \
      --readFilesIn $file'_R1.fastq.gz' $file'_R2.fastq.gz' \
      --outFileNamePrefix $file --outSAMunmapped Within
    done

    mkdir other

    if [ "${params.index}" == "null" ] ; then
        mv STARIndex_last/ other/
    fi

    mv *Log* other/

    """
}
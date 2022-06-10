process Mapping_STAR {
    publishDir params.output+'/mapping/', mode: 'copy'
    cpus params.thread
    
    input:
      tuple val(ID), file(R)

    output:
      tuple file("*Aligned.sortedByCoord.out.bam"), emit: result_STAR
      tuple file("other/"), emit: other_STAR
    
    script:
    """
    mkdir STARIndex_last/
    ## -- Index construction -------------------------------------------------------------- ##
    if [ ${params.index} == null ] ; then
      mkdir STARIndex_last/
      STAR --runThreadN ${params.thread} \
        --runMode genomeGenerate --genomeDir STARIndex_last/ --genomeFastaFiles ${params.FNA} \
        --sjdbGTFfile ${params.GTF} --sjdbOverhang 149 --genomeSAsparseD 2
      IDX=STARIndex_last
    else
      IDX=${params.FNA}
    fi

    ## -- Mapping analyse ----------------------------------------------------------------- ##
    STAR --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN ${params.thread} \
      --runThreadN ${params.thread} --genomeDir "\$IDX" --readFilesCommand gunzip -c \
      --readFilesIn ${R} --outFileNamePrefix ${ID} --outSAMunmapped Within
    
    mkdir other

    if [ "${params.index}" == "null" ] ; then
        mv STARIndex_last/ other/
    fi

    mv *Log* other/

    """
}
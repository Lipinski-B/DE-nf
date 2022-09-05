process Mapping_BWA{ 
    publishDir params.output+'/mapping/', mode: 'copy'
    cpus params.thread
    
    input:
    tuple val(ID), file(R)
    //file data from Channel.fromPath(params.input+'*').collect()
    //file data from Channel.fromPath(params.index+'/*').collect()

    output:
    tuple file("*.sorted.bam*"), emit: result_BAW
    tuple file("other/"), emit: other_BAW
    //file "*.sorted.bam*" into Mapping_bam
    //file "other/" into Mapping_Log
    
    shell:
    '''
    list=`ls *fast[q.a]* -1 | sed 's/_R.*//' | uniq`

    ## -- Index construction -------------------------------------------------------------- ##
    if [ "!{params.index}" == "null" ] ; then
      bwa index !{params.FNA}
      IDX=!{params.FNA}
    else
      IDX="!{params.index}/*.fa"
    fi

    mkdir other/
    
    ## -- Mapping analyse ----------------------------------------------------------------- ##
    for file in $list; do
      bwa mem -o $file.sam -t !{params.thread} $IDX ${R}
      samtools view -@ !{params.thread} -b -O BAM -o $file.bam $file.sam
      samtools sort -@ !{params.thread} $file.bam -o $file.sorted.bam --reference $IDX 
      samtools index -@ !{params.thread} -b $file.sorted.bam
      mv $file.sam $file.bam other/
    done

    
    if [ "!{params.index}" == "null" ] ; then
      mkdir other/BWAIndex_last
      mv !{params.FNA}.* other/BWAIndex_last/
      cp !{params.FNA} other/BWAIndex_last/
    fi

    '''
  }

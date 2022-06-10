process Intersection { 
  publishDir params.output+'/intersect/', mode: 'copy'
  cpus params.thread

  input:
    file bam
    //file GTF from Channel.fromPath(params.GTF).collect()
  
  output:
    tuple file('*.txt'), emit: Intersect
  
  script:
  """
  ## -- Intersection analyse ----------------------------------------------- ##
  for file in *.bam; do
    htseq-count -f bam -n ${params.thread} \$file ${params.GTF} > \${file}_intersect.txt
  done
  """
}
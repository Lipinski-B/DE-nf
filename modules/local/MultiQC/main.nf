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
  #files=(*)
  #for file in *; do
    #fastqc $file
  #done
  #multiqc .
  '''
}

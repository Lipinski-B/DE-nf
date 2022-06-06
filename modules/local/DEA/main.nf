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
}
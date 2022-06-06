#! /usr/bin/env nextflow


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DE pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/Lipinski-B/DE-nf
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2
def nextflowMessage() {
    log.info "N E X T F L O W  ~  DSL 2  ~  version ${workflow.nextflow.version} ${workflow.nextflow.build}"
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETERS VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/




// -- RaiseError :
//if(!file(params.input).exists())      error "ERROR: --input must refer to an existing directory"
//if(!file(params.output).exists())     error "ERROR: --output must refer to an existing directory"




/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOG INFO
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
log.info "--------------------------------------------------------------------------------------"
log.info ""
log.info "          DE 1.0 : Pipeline RNAseq for the Differential Expression analysis."
log.info ""
log.info ""


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
    log.info "--index                       FOLDER                      Folder where you can find the STAR/BWA index. If this option is not used, please make sure to provide the --FNA option in addition to the --GTF option to perform the index"
    log.info "--mapper                      STRING                      [STAR/BWA] Choose the mapper to use between STAR and BWA MEME (Default : BWA MEM)"
    log.info "--FNA                         FILE                        Path where you can find the FNA file to use for the STAR index."
    log.info "--R                           STRING                      [on/off] : Chose to use or not the standard R analyses from the pipeline."
    log.info "--metadata                    FILE                        Path where you can find the XLS file to use as metadata for the R analyse. Mandatory is the option --R in on."
    log.info "--thread                      INT                         Number of thread to use."
  
    exit 0
} else {
    log.info "-------------------------------- Nextflow parameters ---------------------------------"
    log.info ""
    log.info "Project              : $workflow.projectDir"
    log.info "Git repository       : $workflow.repository"
    log.info "Release [Commit ID]  : $workflow.revision [$workflow.commitId]"
    log.info "User Name            : $workflow.userName"
    log.info "Run Name             : $workflow.runName"
    log.info "Resume               : $workflow.resume"
    log.info "Script Name          : $workflow.scriptName"
    log.info "Script File          : $workflow.scriptFile"
    log.info "Home Directory       : $workflow.homeDir"
    log.info "Work Directory       : $workflow.workDir"
    log.info "Launch Directory     : $workflow.launchDir"
    log.info "Command line         : $workflow.commandLine"
    log.info "Config Files         : $workflow.configFiles"
    log.info "Config Profile       : $workflow.profile"
    log.info "Container Engine     : $workflow.containerEngine"
    log.info "Container            : $workflow.container"
    log.info "Session ID           : $workflow.sessionId"
    log.info "Script ID            : $workflow.scriptId"
    log.info ""
    log.info "-------------------------------- Workflow parameters ---------------------------------"
    log.info ""
    log.info "date                 : ${params.date}"
    log.info "input                : ${params.input}"
    log.info "output               : ${params.output}"
    log.info "n                    : ${params.n}"
    log.info ""
    log.info "--------------------------------------------------------------------------------------"
    log.info ""
}





/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// -- Pipeline :

if(params.mapper=="STAR"){
 
} else {
  }

if(params.R=="on"){

}



// -- Modules :
include { MultiQC }       from "${baseDir}/modules/local/MultiQC/main"
include { Mapping_STAR }  from "${baseDir}/modules/local/Mapping_STAR/main"
include { Mapping_BWA }   from "${baseDir}/modules/local/Mapping_BWA/main"
include { Intersection }  from "${baseDir}/modules/local/Intersection/main"
include { Merge_result }  from "${baseDir}/modules/local/Merge_result/main"
include { DEA }           from "${baseDir}/modules/local/DEA/main"

// -- Pipeline :
workflow DE_nf {
    take: fastq
    main:
        //MultiQC           (fastq)
        Mapping_STAR      (fastq)
        //Mapping_BWA       (fastq)
        Intersection      (Mapping_STAR.out.Mapping_STAR.collect())
        Merge_result      (Intersection.out.Intersection.collect())
        //DEA               (Merge_result.out.Merge_result)
}



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow {
    DE_nf( channel.fromFilePairs("${params.input}*_R{1,2}.fastq.gz") )
}




/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow.onComplete {
	this.nextflowMessage()
	log.info "Completed at: " + workflow.complete
	log.info "Duration    : " + workflow.duration
	log.info "Success     : " + workflow.success
	log.info "Exit status : " + workflow.exitStatus
	log.info "Error report: " + (workflow.errorReport ?: '-')}

workflow.onError {
	this.nextflowMessage()
	log.info "Workflow execution stopped with the following message:"
	log.info "  " + workflow.errorMessage}
process seurat_doublet {

    cpus { 26 }
    errorStrategy 'retry'
    maxRetries 6
    memory { 340.GB }
    
    label 'seurat'

    //publishDir 'mut_compiled', mode: 'copy', overwrite: true, pattern: '*snp_summarise.txt.gz'

    input: 
    tuple val(sample), file("filtered_seurat.RData"), file("${sample}_initial_QC_plots")
    
    output:
    file("doublet_seurat.RData")

    script:
    """
    
    #!/bin/bash
    Rscript ${projectDir}/Rscripts/seurat/doublet_finder.R \
	filtered_seurat.RData \
	. \
	${task.cpus} \
	${params.cellcycle_markers} \
	TRUE \
	${projectDir}

 	
    mv outdata/doublet_seurat.RData .
    
    """

}



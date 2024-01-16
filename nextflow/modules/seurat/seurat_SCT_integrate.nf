process seurat_SCT_integrate {

    cpus { 26 }
    errorStrategy 'retry'
    maxRetries 6
    memory { 340.GB }
    
    label 'seurat'

    publishDir 'seurat_objects', mode: 'copy', overwrite: true, pattern: 'integrated_seurat.RData'

    input: 
    file("doublet_seurat.RData")

    output:
    file("integrated_seurat.RData")

    script:
    """
    #!/bin/bash
    echo ${task.memory}
    Rscript ${projectDir}/Rscripts/seurat/SCT_integrate.R \
	.

 	
    mv outdata/integrated_seurat.RData .

    """
}



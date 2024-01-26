process R_var_filtering {

    label 'tidyverse'

    cpus { 4 * task.attempt } 
    errorStrategy 'retry'
    maxRetries 6
    memory { 512.GB * task.attempt }


    publishDir 'filtered_var', mode: 'copy', overwrite: true, pattern: '*snp_summarise.txt.gz'

    input: 
    tuple val(species), val(sample), val(contig), file("${species}_${sample}_${contig}_fin_snp_read.txt.gz")

    output:
    tuple val(sample), file("${sample}_${contig}_snp_summarise.txt.gz")
    script:
    """
    #!/bin/bash

    Rscript ${projectDir}/Rscripts/var_filtering.R \
        ${species}_${sample}_${contig}_fin_snp_read.txt.gz \
	$sample \
	$contig

    """
}

process ws_var_call_HC {

    label 'gatk'
    errorStrategy 'retry'

    cpus { 16 * task.attempt }
    errorStrategy 'retry'
    memory { 128.GB * task.attempt }

    input:
    tuple val(species), val(sample), file("${sample}_dup_NCR.bam"), val(ref)

    output:
    tuple val(species), val(sample), file("${sample}_dup_NCR.bam"), file("${sample}.vcf.gz"), val(ref)

    script:
    """
    #!/bin/bash
    gatk --java-options "-Xmx4g" HaplotypeCaller -R ${params.fasta_dir}/${ref}.fna -I ${sample}_dup_NCR.bam -O ${sample}.vcf.gz -ERC NONE
    """
}


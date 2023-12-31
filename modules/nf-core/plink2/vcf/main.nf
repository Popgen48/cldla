process PLINK2_VCF {
    tag "$meta"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a2.3--h712d239_1' :
        'biocontainers/plink2:2.00a2.3--h712d239_1' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.bed")    , emit: bed
    tuple val(meta), path("*.bim")    , emit: bim
    tuple val(meta), path("*.fam"), emit: fam
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    def mem_mb = task.memory.toMega()
    def max_chrom = params.num_autosomes

    """
    plink2 \\
        --threads $task.cpus \\
        --memory $mem_mb \\
        --chr-set ${max_chrom} \\
        --double-id \\
        $args \\
        --vcf $vcf \\
        --make-bed \\
        --out ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//' )
    END_VERSIONS
    """
}

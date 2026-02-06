#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

// ================================================================
// PARAMETERS (User Could Change on the basis of their requirments)
// ================================================================

params.genome = "hg38" // name
params.genome_fasta = "genome/hg38.fa"        //Refrence genome
params.outdir = "results"               //Output folder
params.email = null                     // Email for completion notification

// Genome sizes for normalization
params.genome_sizes = [
    'hg38': 2700000000,
    'hg19': 2700000000,
    'mm10': 1870000000,
    'mm9': 1870000000
]

// Tool paths (set these to your installations)
params.spp_script = "$HOME/tools/phantompeakqualtools/run_spp.R"

// Processing options
params.skip_trimming = false
params.skip_fastqc = false
params.skip_multiqc = false

// Peak calling options
params.macs_gsize = params.genome_sizes[params.genome] ?: 2.7e9
params.macs_qvalue = 0.05
params.broad_peaks = false                 // Use for histone marks

// ============================================
// VALIDATE PARAMETERS
// ============================================

if (!params.input) {
    exit 1, "Please provide a samplesheet with --input"
}

log.info """
========================================================================================
    ChIP-seq Analysis Pipeline v1.0
========================================================================================
    Input samplesheet : ${params.input}
    Genome           : ${params.genome}
    Output directory : ${params.outdir}
    Genome size      : ${params.macs_gsize}
    Broad peaks      : ${params.broad_peaks}
========================================================================================
"""

// ================================================================
//PROCESS FASTQC - Check Qulaity
// ================================================================

process FASTQC {
    tag "$sample_id"
    publishDir "${params.outdir}/01_fastqc", mode: 'copy'

    when:
    !params.skip_fastqc

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*.{html,zip}", emit: reports

    script:
    """
    fastqc -q  -t 2 ${reads}
    """
}

// ============================================
// PROCESS TRIM_GALORE - Adapter Trimming
// ============================================

process TRIM_GALORE {
    tag "$sample_id"
    publishDir "${params.outdir}/02_trimmed", mode: 'copy',
        saveAs: { filename -> filename.endsWith('.txt') ? filename : null }
    
    when:
    !params.skip_trimming
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("*_trimmed.fq.gz"), emit: reads
    path "*trimming_report.txt", emit: reports
    
    script:
    def single_end = reads instanceof Path
    if (single_end) {
        """
        trim_galore \
            --gzip \
            --cores 2 \
            ${reads}
        
        mv ${reads.simpleName}_trimmed.fq.gz ${sample_id}_trimmed.fq.gz
        """
    } else {
        """
        trim_galore \
            --paired \
            --gzip \
            --cores 2 \
            ${reads[0]} ${reads[1]}
        
        mv ${reads[0].simpleName}_val_1.fq.gz ${sample_id}_R1_trimmed.fq.gz
        mv ${reads[1].simpleName}_val_2.fq.gz ${sample_id}_R2_trimmed.fq.gz
        """
    }
}

// ================================================================
//PROCESS BWA - Align to genome
// ================================================================

process BWA_INDEX {
    tag "$fasta"
    
    input:
    path fasta
    
    output:
    path "${fasta}*", emit: index
    
    script:
    """
    bwa index ${fasta}
    samtools faidx ${fasta}
    """
}

process BWA_ALIGN {
    tag "$sample_id"
    publishDir "${params.outdir}/03_alignment", mode: 'copy',
        saveAs: { filename -> filename.endsWith('.txt') ? filename : null }
    
    input:
    tuple val(sample_id), path(reads)
    path genome_index
    
    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), emit: bam
    path "${sample_id}.flagstat.txt", emit: flagstat
    
    script:
    def single_end = reads instanceof Path
    def read_input = single_end ? reads : "${reads[0]} ${reads[1]}"
    def genome_name = genome_index[0].toString().replaceAll(/\.(amb|ann|bwt|pac|sa|fai)$/, '')
    """
    bwa mem \
        -t ${task.cpus} \
        -M \
        ${genome_name} \
        ${read_input} \
        | samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam -
    
    samtools index ${sample_id}.sorted.bam
    samtools flagstat ${sample_id}.sorted.bam > ${sample_id}.flagstat.txt
    """
}

// ================================================================
// PROCESS Remove Duplicates (Picard only)
// Primary goal is to identify and remove artificial DNA copies created 
// during the laboratory PCR amplification process
// ================================================================

process MARK_DUPLICATES {
    tag "$sample_id"
    publishDir "${params.outdir}/04_dedup", mode: 'copy',
        saveAs: { filename -> filename.endsWith('.txt') ? filename : null }
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    
    output:
    tuple val(sample_id), path("${sample_id}.dedup.bam"), emit: bam
    path "${sample_id}.dup_metrics.txt", emit: metrics
    
    script:
    """
    picard MarkDuplicates \
        INPUT=${bam} \
        OUTPUT=${sample_id}.dedup.bam \
        METRICS_FILE=${sample_id}.dup_metrics.txt \
        REMOVE_DUPLICATES=true \
        ASSUME_SORTED=true \
        VALIDATION_STRINGENCY=LENIENT
    """
}

// ================================================================
// PROCESS Index Deduplicated BAM (Samtools only)
// Separate process to use different container with samtools
// ================================================================

process INDEX_DEDUP_BAM {
    tag "$sample_id"
    publishDir "${params.outdir}/04_dedup", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam)
    
    output:
    tuple val(sample_id), path(bam), path("${bam}.bai"), emit: bam
    path "${sample_id}.dedup.flagstat.txt", emit: flagstat
    
    script:
    """
    samtools index ${bam}
    samtools flagstat ${bam} > ${sample_id}.dedup.flagstat.txt
    """
}

// ================================================================
// PROCESS PhantomPeakQualTools - Quality
// ================================================================

process PHANTOMPEAKQUALTOOLS {
    tag "$sample_id"
    publishDir "${params.outdir}/05_qc", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    
    output:
    tuple val(sample_id), path("${sample_id}.spp.out"), emit: metrics
    path "${sample_id}.spp.pdf", emit: plot
    
    script:
    """
    Rscript ${params.spp_script} \
        -c=${bam} \
        -savp=${sample_id}.spp.pdf \
        -out=${sample_id}.spp.out \
        -rf
    """
}

// ============================================
// PROCESS Extract Fragment Length
// ============================================

process EXTRACT_FRAGMENT_LENGTH {
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(spp_out)
    
    output:
    tuple val(sample_id), env(FRAG_LEN), emit: frag_len
    
    shell:
    '''
    FRAG_LEN=$(awk '{print $3}' !{spp_out} | cut -d',' -f1)
    echo "Fragment length for !{sample_id}: $FRAG_LEN"
    '''
}

// ================================================================
// PROCESS Create bigWig
// ================================================================

process CREATE_BIGWIG {
    tag "$sample_id"
    publishDir "${params.outdir}/bigwig", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai), val(frag_len)
    
    output:
    path "${sample_id}.bw"
    
    script:
    """
    bamCoverage \
        -b ${bam} \
        -o ${sample_id}.bw \
        --binSize 10 \
        --normalizeUsing RPGC \
        --effectiveGenomeSize ${params.macs_gsize} \
        --extendReads ${frag_len} \
        --numberOfProcessors ${task.cpus}
    """
}

// ================================================================
// PROCESS MACS2 - Call peaks
// ================================================================

process MACS2_CALLPEAK {
    tag "$chip_id"
    publishDir "${params.outdir}/07_peaks", mode: 'copy'
    
    input:
    tuple val(chip_id), path(chip_bam), path(chip_bai), val(chip_frag), 
          val(input_id), path(input_bam), path(input_bai)
    
    output:
    tuple val(chip_id), path("${chip_id}_peaks.narrowPeak"), emit: peaks
    path "${chip_id}_peaks.xls", emit: xls
    path "${chip_id}_summits.bed", emit: summits
    path "${chip_id}_model.r", optional: true, emit: model
    
    script:
    def broad = params.broad_peaks ? "--broad --broad-cutoff 0.1" : ""
    """
    macs2 callpeak \
        -t ${chip_bam} \
        -c ${input_bam} \
        -f BAM \
        -g ${params.macs_gsize} \
        -n ${chip_id} \
        -q ${params.macs_qvalue} \
        --extsize ${chip_frag} \
        --keep-dup all \
        ${broad}
    """
}

// ============================================
// PROCESS Annotate Peaks with HOMER
// ============================================

process ANNOTATE_PEAKS {
    tag "$sample_id"
    publishDir "${params.outdir}/08_annotation", mode: 'copy'
    
    input:
    tuple val(sample_id), path(peaks)
    
    output:
    path "${sample_id}.annotated.txt"
    
    script:
    """
    annotatePeaks.pl \
        ${peaks} \
        ${params.genome} \
        > ${sample_id}.annotated.txt
    """
}

// ============================================
// PROCESS MultiQC - Aggregate Reports
// ============================================

process MULTIQC {
    publishDir "${params.outdir}/09_multiqc", mode: 'copy'
    
    when:
    !params.skip_multiqc
    
    input:
    path '*'
    
    output:
    path "multiqc_report.html"
    path "multiqc_data"
    
    script:
    """
    multiqc . \
        --title "ChIP-seq Analysis Report" \
        --filename multiqc_report.html
    """
}

// ============================================
// WORKFLOW - Connect All Processes
// ============================================

workflow {
    
    // Read sample sheet
    Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            def chip_id = row.chip_sample
            def chip_r1 = file(row.chip_fastq_r1, checkIfExists: true)
            def chip_r2 = row.chip_fastq_r2 ? file(row.chip_fastq_r2, checkIfExists: true) : null
            def input_id = row.input_sample
            def input_r1 = file(row.input_fastq_r1, checkIfExists: true)
            def input_r2 = row.input_fastq_r2 ? file(row.input_fastq_r2, checkIfExists: true) : null
            
            return [chip_id, chip_r1, chip_r2, input_id, input_r1, input_r2]
        }
        .set { samples_ch }
        
    // Separate ChIP and Input samples
    samples_ch
        .flatMap { chip_id, chip_r1, chip_r2, input_id, input_r1, input_r2 ->
            def chip_reads = chip_r2 ? [chip_r1, chip_r2] : chip_r1
            def input_reads = input_r2 ? [input_r1, input_r2] : input_r1
            [
                [chip_id, chip_reads, 'chip'],
                [input_id, input_reads, 'input']
            ]
        }
        .unique()
        .map { sample_id, reads, type -> [sample_id, reads] }
        .set { all_samples_ch }
    
    // Get genome
    genome_ch = Channel.fromPath(params.genome_fasta ?: 
        "$HOME/genomes/${params.genome}/${params.genome}.fa")
    
    // Step 1: FastQC
    FASTQC(all_samples_ch)
    
    // Step 2: Trim adapters
    if (!params.skip_trimming) {
        TRIM_GALORE(all_samples_ch)
        reads_for_alignment = TRIM_GALORE.out.reads
    } else {
        reads_for_alignment = all_samples_ch
    }
    
    // Index genome and align
    BWA_INDEX(genome_ch)
    BWA_ALIGN(reads_for_alignment, BWA_INDEX.out.index.collect())
    
    // Step 4: Remove duplicates (now two separate processes)
    MARK_DUPLICATES(BWA_ALIGN.out.bam)
    INDEX_DEDUP_BAM(MARK_DUPLICATES.out.bam)
    
    // Quality metrics (using indexed BAMs)
    PHANTOMPEAKQUALTOOLS(INDEX_DEDUP_BAM.out.bam)
    EXTRACT_FRAGMENT_LENGTH(PHANTOMPEAKQUALTOOLS.out.metrics)
    
    // Create bigWig files
    INDEX_DEDUP_BAM.out.bam
        .join(EXTRACT_FRAGMENT_LENGTH.out.frag_len)
        .set { bam_with_fraglen }
    
    CREATE_BIGWIG(bam_with_fraglen)
    
    // Peak calling - pair ChIP with Input
    samples_ch
        .map { chip_id, chip_r1, chip_r2, input_id, input_r1, input_r2 ->
            [chip_id, input_id]
        }
        .combine(INDEX_DEDUP_BAM.out.bam, by: 0)
        .map { chip_id, input_id, chip_bam, chip_bai ->
            [input_id, chip_id, chip_bam, chip_bai]
        }
        .combine(INDEX_DEDUP_BAM.out.bam, by: 0)
        .map { input_id, chip_id, chip_bam, chip_bai, input_bam, input_bai ->
            [chip_id, input_id, chip_bam, chip_bai, input_bam, input_bai]
        }
        .combine(EXTRACT_FRAGMENT_LENGTH.out.frag_len, by: 0)
        .map { chip_id, input_id, chip_bam, chip_bai, input_bam, input_bai, frag_len ->
            [chip_id, chip_bam, chip_bai, frag_len, input_id, input_bam, input_bai]
        }
        .set { peak_calling_input }
    
    MACS2_CALLPEAK(peak_calling_input)
    
    // Annotate peaks
    ANNOTATE_PEAKS(MACS2_CALLPEAK.out.peaks)
    
    // Collect all QC files
    qc_files = Channel.empty()
        .mix(
            FASTQC.out.reports.collect().ifEmpty([]),
            TRIM_GALORE.out.reports.collect().ifEmpty([]),
            BWA_ALIGN.out.flagstat.collect().ifEmpty([]),
            MARK_DUPLICATES.out.metrics.collect().ifEmpty([]),
            INDEX_DEDUP_BAM.out.flagstat.collect().ifEmpty([]),
            PHANTOMPEAKQUALTOOLS.out.metrics.collect().ifEmpty([])
        )
    
    // MultiQC report
    MULTIQC(qc_files.collect())
}

// ============================================
// COMPLETION NOTIFICATION
// ============================================

workflow.onComplete {
    log.info """
    ========================================================================================
    Pipeline completed at: $workflow.complete
    Execution status: ${ workflow.success ? 'SUCCESS' : 'FAILED' }
    Duration: $workflow.duration
    Results directory: ${params.outdir}
    ========================================================================================
    """
}
// Creates a parameters file and a summary file to 
// be added to later
process Setup {
    input:
        // The name of the reference supplied (for use
        // in the parameters file)
        val refName
        // A value denoting whether PCR deduplication was
        // enabled.
        val pcrDedupe
        // A value containing which alignment
        // mode was used during the analysis.
        val alignmentMode
        // The minimum read length allowed post trimmming (for
        // use in the parameters file)
        val minLen
        // The minimum coverage depth used for masking
        // the consensus genome.
        val minCov
        // The minimum base call quality for a site to be considered
        // in both variant calling and depth masking
        val minBQ
        // The minimum mapping quality for a site to be considered in
        // both variant calling and depth masking
        val minMapQ
        // The name of the host reference file or index
        // if supplied. If not supplied, the value will 
        // be NONE.
        val hostRefName
        // The header to write to the summary file.
        val summaryHeader
        // The output directory to be used.
        val outDir
        
    output:
        // The parameters file created.
        file "analysis-parameters.txt"
        // The blank summary file to be added to.
        file "stats-summary.csv"

    publishDir "${outDir}", mode: 'copy'

    script:
    /*
    Creates a parameters file (in case the user needs to look back at how they ran the analysis)
    as well as a blank summary file, which will later be populated with data from each
    sample.

    The parameters file contains:
        1. The minimum read length allowed after trimmming.
        2. Whether PCR deduplication was enabled.
        3. The name of host reference genome (if supplied)
        4. The name of the reference genome supplied
        5. The alignment mode used for reference-based assembly
        6. The minimum coverage threshold used for masking
        7. The minimum base call quality used for variant calling and masking
        8. The minimum mapping quality used for variant calling and masking

    The summary file will always contain:
        1. The sample
        2. Raw Reads
        3. Reads after trimming
        4. Contigs Generated
        5. Scaffolds Generated
        6. Contigs aligned to reference
        7. Variants applied to the reference
        8. Reads aligned to the corrected reference
        9. Average Read Depth
        10. SNPs passing filtering
        11. Indels passing filtering
        12. sites masked
        13. coverage
    As well, depending on the user's options, the summary file may contain:
        - Reads after PCR Deduplication
        - Reads after host removal
    */
    """
    #!/bin/bash

    touch analysis-parameters.txt

    echo "Minimum Read Length Post-Trimming : ${minLen} bp" >> analysis-parameters.txt
    echo "PCR Deduplication : ${pcrDedupe}" >> analysis-parameters.txt
    echo "Host Removal : ${hostRefName}" >> analysis-parameters.txt
    echo "Reference Supplied : ${refName}" >> analysis-parameters.txt
    echo "Alignment Mode: ${alignmentMode}" >> analysis-parameters.txt
    echo "Minimum Coverage Allowed : ${minCov}" >> analysis-parameters.txt
    echo "Minimum Base Call Quality Allowed: ${minBQ}" >> analysis-parameters.txt
    echo "Minimum Read Mapping Quality Allowed: ${minMapQ}" >> analysis-parameters.txt

    touch stats-summary.csv

    echo "${summaryHeader}" > stats-summary.csv
    """
}

// Builds a bowtie2 index for a provided reference file
process Index_Host_Reference {
    input:
        // Tuple contains the reference name and reference file.
        tuple val(refName), file(ref)
        // The output directory
        val outDir
        // The number of threads provided
        val threads

    output:
        // Tuple containing the reference index directory and the index file name.
        tuple file("${refName}-idx/"), val(refName)

    publishDir "${outDir}", mode: 'copy'

    /*
    Creates an index Directory and then runs bowtie2-build to create the index
    */
    script:
    """
    #!/bin/bash
    
    mkdir ${refName}-idx/

    bowtie2-build --threads ${threads} ${ref} ${refName}-idx/${refName}
    """
}

// Creates a fastqc report for a set of paired-end reads provided.
process QC_Report {
    input:
        // Tuple contains the file basename as well as the paired-end read files
        tuple val(base), file(F1), file(F2)
        // The output directory
        val outDir
        // The name of the directory to place the fastqc output into (allows
        // for this command to be used multiple times and to separate the output
        // i.e. pre-processed-reads vs trimmed-reads.)
        val dirName

    output:
        // The directory cotaining the fastqc files produced.
        file("${base}")
    
    publishDir "${outDir}/${dirName}/", mode: 'copy'

    // Creates a Directory and then runs fastqc on a set of reads.
    script:
    """
    #!/bin/bash

    mkdir ${base}

    fastqc ${F1} ${F2} -o ${base}
    """
}

// Performs quality and adapter trimming on a set of paired-end reads.
process Trimming {
    input: 
        // Tuple cotains the file basename as well as the paired-end read files.
        tuple val(base), file(R1), file(R2)
        // The output directory
        val outDir
        // The adapter file in fasta format.
        file adapters
        //The minimum read lenths to be allowed post trimming
        val minLen

    output:
        // Tuple containing the file basename, the trimmed forward/reverse reads, and 
        tuple val(base), file("${base}_1.trimmed.fastq"), file("${base}_2.trimmed.fastq")
        // Tuple containing the unpaired read files.
        tuple file("${base}_1.unpaired.fastq.gz"), file("${base}_2.unpaired.fastq.gz") 
        // The summary string containing the sample and raw/trimmed read counts.
        env summary

    publishDir "${outDir}/${base}-Intermediate-Files/Processed-Reads", mode: 'copy'

    script:
    /*
    The number of raw reads in the forward and reverse files are grabbed
    and used to calculate the total raw reads.

    Uses trimmomatic to trim the read files.
    Trimmomatic performs trimming steps in the order provided in the 
    command line. Our steps:
    1. ILLUMINACLIP: Removes illumina adapters
    2. LEADING: Begins at the start of the read and trims bases with quality less than
        5 until it hits a base above that threshold.
    3. TRAILING: Begins at the end of the read and trims bases with quality less than 5
        until it hits a base above that threshold.
    4. SLIDINGWINDOW: Begins at the 5' end and scans in 4 base windows, trimming when it hits
        an average quality less than 20.
    5. MINLEN: Removes reads less than 75 bp long.

    Prinseq (used in the following process) cannot take gzipped reads as input. Thus,
    only the unpaired reads are gzipped to save space. 

    Finally, the forward and reverse reads post trimming are grabbed and used
    to calculate the total trimmed reads.

    The sample, raw reads, and trimmed reads are added to the summary string.
    */
    """
    #!/bin/bash

    raw_reads_1=\$((\$(gunzip -c ${R1} | wc -l)/4))
    raw_reads_2=\$((\$(gunzip -c ${R1} | wc -l)/4))

    total_raw=\$((\$raw_reads_1 + \$raw_reads_2))

    trimmomatic PE ${R1} ${R2} ${base}_1.trimmed.fastq ${base}_1.unpaired.fastq \
    ${base}_2.trimmed.fastq ${base}_2.unpaired.fastq ILLUMINACLIP:${adapters}:2:30:10:1:true \
    LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20 MINLEN:${minLen}

    gzip ${base}_1.unpaired.fastq ${base}_2.unpaired.fastq
    
    trimmed_reads_1=\$((\$(cat ${base}_1.trimmed.fastq | wc -l)/4))
    trimmed_reads_2=\$((\$(cat ${base}_2.trimmed.fastq | wc -l)/4))

    total_trimmed=\$((\$trimmed_reads_1 + \$trimmed_reads_2))

    summary="${base},\$total_raw,\$total_trimmed"
    """
        
}

// Removes PCR Duplicates from a set of reads using prinseq.
process Remove_PCR_Duplicates {
    input:
        // Tuple contains the file basename and the two read files.
        tuple val(base), file(R1), file(R2)
        // The output directory
        val outDir
        // The existing summary file
        val existingSummary
    output:
        // Tuple containing the file basename and the deduped, paired-end reads.
        tuple val(base), file("${base}_deduped_1.fastq.gz"), file("${base}_deduped_2.fastq.gz")
        // The log file produced by prinseq for troubleshooting.
        file "${base}-prinseq-log.txt"
        // The summary string containing the number of reads after deduplication
        env summary

    publishDir "${outDir}/${base}-Intermediate-Files/Processed-Reads", mode: 'copy'

    /*
    Calls prinseq on the paired-end reads. The --derep command
    tells prinseq to deduplicate the reads. Prinseq has 5 
    derep modes:

    "1 (exact duplicate), 2 (5' duplicate), 3 (3' duplicate), 4
    (reverse complement exact duplicate), 5 (reverse complement
    5'/3' duplicate))"

    This pipelines only removes exact duplicates. (--derep 1)

    The fastq files produced are gzipped to save space.

    The forward and reverse reads post deduplication are grabbed
    and summed to calculate the total deduplicated reads.

    This value is added to the summary string.
    */
    script:
    """
    #!/bin/bash

    prinseq-lite.pl -fastq ${R1} -fastq2 ${R2} --out_good ${base}_deduped --derep 1 -log ${base}-prinseq-log.txt

    gzip ${R1} ${R2} ${base}_deduped_1.fastq ${base}_deduped_2.fastq

    deduped_reads_1=\$((\$(gunzip -c ${base}_deduped_1.fastq | wc -l)/4))
    deduped_reads_2=\$((\$(gunzip -c ${base}_deduped_2.fastq | wc -l)/4))

    total_deduped=\$((\$deduped_reads_1 + \$deduped_reads_2))

    summary="${existingSummary},\$total_deduped"
    """
}

// Removes host reads by aligning to a host genome using bowtie2.
process Host_Read_Removal {
    input:
        // Tuple contains the file basename and paired-end reads
        tuple val(base), file(R1), file(R2)
        // The output directory
        val outDir
        // Tuple contains the bt2 index directory and basename of the index files.
        tuple file(refDir), val(refName)\
        
        val alignmentMode
        // The number of threads provided.
        val threads
        // The existing summary string
        val existingSummary
    output:
        // Tuple contains the file basename and the paired-end read files with host reads removed.
        tuple val(base), file("${base}_host_removed_1.fq.gz"), file("${base}_host_removed_2.fq.gz")
        // A directory containing the alignment to the host file.
        file "host-reads/${base}-host.sam"
        // The summary string containing the number of reads after host removal
        env summary
    
    publishDir "${outDir}/${base}-Intermediate-Files/Processed-Reads", mode: 'copy'

    script:
    /*
    Aligns the reads to a host reference genome using bowtie2. Local alignment is used (--local)
    to ensure the reads are aligne without gaps. The unaligned reads are sent back ot into a
    fastq files (--un-conc).

    The unaligned read files are then renamed to give them the typical paired end
    read name scheme.

    Additionally, the reads are gzipped to preserve space.

    Finally, the number of forward and reverse reads are grabbed and used to calculate
    the total number of reads post host removal. This value is added to the summary string.
    */
    """
    #!/bin/bash
    mkdir host-reads
    bowtie2 --threads ${threads} -x ${refDir}/${refName} -1 ${R1} -2 ${R2} ${alignmentMode} -S host-reads/${base}-host.sam --un-conc ${base}_host_removed

    mv ${base}_host_removed.1 ${base}_host_removed_1.fq
    mv ${base}_host_removed.2 ${base}_host_removed_2.fq

    gzip ${base}_host_removed_1.fq ${base}_host_removed_2.fq

    nonHost_reads_1=\$((\$(gunzip -c ${base}_host_removed_1.fastq | wc -l)/4))
    nonHost_reads_2=\$((\$(gunzip -c ${base}_host_removed_2.fastq | wc -l)/4))

    total_nonHost=\$((\$nonHost_reads_1 + \$nonHost_reads_2))

    summary="${existingSummary},\$total_nonHost"
    """
}

// Uses spades to produce a de novo assembly.
process Spades_Assembly {
    input:
        // Tuple contains the file basename and paired-end reads.
        tuple val(base), file(R1), file(R2)
        // The output directory
        val outDir
        // The number of threads provided.
        val threads
        // The existing summary string
        val existingSummary
    output:
        // Tuple contains the basename of the sample and the assembled contigs 
        // produced by spades.
        tuple val(base), file(R1), file(R2), file("${base}-contigs.fasta")
        // The scaffolds produced by spades
        file "${base}-scaffolds.fasta"
        // The summary string containing the number of contigs and scaffolds
        env summary
    
    publishDir "${outDir}", mode: 'copy', pattern: "${base}-contigs.fasta"
    publishDir "${outDir}/${base}-Intermediate-Files", mode: 'copy', pattern: "${base}-scaffolds.fasta"

    script:
    /*
    Runs spades using the provided paired-end reads.

    The contigs and scaffolds are renamed and moved, and
    the number of contigs/scaffolds are recorded.
    
    The number of contigs and scaffolds are added to the summary string.
    */
    """
    #!/bin/bash

    spades.py --threads ${threads} -1 ${R1} -2 ${R2} -o ${base}-Assembly

    mv ${base}-Assembly/contigs.fasta ./${base}-contigs.fasta

    num_contigs=\$(grep ">" ${base}-contigs.fasta | wc -l)

    mv ${base}-Assembly/scaffolds.fasta ./${base}-scaffolds.fasta

    num_scaffolds=\$(grep ">" ${base}-scaffolds.fasta | wc -l)

    summary="${existingSummary},\$num_contigs, \$num_scaffolds"
    """
}

// Generates an alignment of the asembly contigs to a reference genome
// using minimap.
process Contig_Alignment {
    input:
        // Tuple contains the sample basename
        // and the assembly fasta to align
        tuple val(base), file(R1), file(R2), file(assembly)
        // The output directory
        val outDir
        // the reference fasta file to be aligned to.
        tuple val(refName), file(ref)
        // The existing summary string
        val existingSummary

    output:
        // Tuple contains the file basename and the alignment bam file.
        tuple val(base), file(R1), file(R2), file("${base}-contig-align.bam")
        // The summary string containing the number of mapped contigs
        env summary

    publishDir "${outDir}/${base}-Intermediate-Files", mode: 'copy', pattern: "${base}-contig-align.bam"
    
    script:
    /*
    Uses minimap2 to align the contigs to the reference fasta. An option to ignore secondary
    alignments is provided to minimap2.

    Then uses samtools to convert the alignment sam into bam format
    (samtools view). Additionally, any non-primary alignments (i.e. supplementary alignments)
    are filetered out durign this step (using the -F 2048 option).
    The output then sorted and stored in a bam file (samtools sort).
    */
    """
    #!/bin/bash

    minimap2 --secondary=no -a ${ref} ${assembly} > align.sam

    samtools view -F 2048 -b align.sam | samtools sort > ${base}-contig-align.bam

    mapped_contigs=\$(samtools view -F 0x04 -c ${base}-contig-align.bam)

    summary="${existingSummary},\$mapped_contigs"
    """
}

process Correct_Reference_with_Contigs {
    input: 
        tuple val(base), file(R1), file(R2), file(bam)
        // Tuple contains the reference name and reference fasta file
        tuple val(refName), file(ref)
        // The output directory
        val outDir
        // The existing summary string
        val existingSummary

    output:
        // Tuple contains the sample base name, raw reads, and the corrected reference genome.
        tuple val(base), file(R1), file(R2), file("${base}-corrected-ref.fasta")
        // The VCF File containing filtered variants from the contigs that were applied
        // to the reference.
        file "${base}-contig-align-filtered.vcf.gz"
        // The summary string containing the number of variants applied to the 
        // reference
        env summary

    publishDir "${outDir}/${base}-Intermediate-Files", mode: "copy", pattern: "${base}-contig-align-filtered.vcf.gz"
    publishDir "${outDir}", mode: "copy", pattern: "${base}-corrected-ref.fasta"

    script:
    """
    #!/bin/bash

    samtools index ${bam}

    freebayes -f ${ref} -C 1 -m 0 ${bam} > ${base}-contig-align.vcf
    
    python3 ${baseDir}/scripts/remove_genotype_vcf.py -i ${base}-contig-align.vcf -o ${base}-contig-align-no-genotype.vcf

    python3 ${baseDir}/scripts/fix_multi_allelic.py -i ${base}-contig-align-no-genotype.vcf -o ${base}-contig-align-biallelic.vcf

    bcftools view -i "((INFO/AO /  INFO/DP) > (INFO/RO / INFO/DP)) & (REF != 'N')" ${base}-contig-align-biallelic.vcf > ${base}-contig-align-filtered.vcf

    correction_variants=\$(grep -v "^#" ${base}-contig-align-filtered.vcf | wc -l)

    bgzip ${base}-contig-align-filtered.vcf

    tabix ${base}-contig-align-filtered.vcf.gz

    bcftools consensus -f ${ref} ${base}-contig-align-filtered.vcf.gz > ${base}-corrected-ref.fasta

    summary="${existingSummary},\$correction_variants"
    """
}

process Align_Reads_to_Corrected_Reference {
    input:
        tuple val(base), file(R1), file(R2), file(correctedRef)

        val alignmentMode

        val outDir

        val existingSummary
    
    output:

        tuple val(base), file("${base}-corrected-ref-align.bam"), file(correctedRef)

        env summary

    publishDir "${outDir}/${base}-Intermediate-Files", mode: "copy", pattern: "${base}-corrected-ref-align.bam"

    script: 
    """
    #!/bin/bash

    bowtie2-build ${correctedRef} ${base}-corrected-idx

    bowtie2 -x ${base}-corrected-idx -1 ${R1} -2 ${R2} ${alignmentMode} -S ${base}-align.sam

    samtools view -b ${base}-align.sam | samtools sort > ${base}-corrected-ref-align.bam

    mapped_reads=\$(samtools view -F 0x04 -c ${base}-corrected-ref-align.bam)

    avg_read_depth=\$(samtools depth -a -J -q 0 -Q 0 ${base}-corrected-ref-align.bam | awk -F'\t' 'BEGIN{totalCov=0} {totalCov +=\$3} END{print totalCov/NR}')

    summary="${existingSummary},\$mapped_reads,\$avg_read_depth"
    """

}

// Performs variant calling and filtering.
process Call_Variants {
    input:
        // Tuple contains the file basename and alignment bam file
        tuple val(base), file(bam), file(correctedRef)
        // The script base directory name (to call python scripts)
        val baseDir
        // The output directory name
        val outDir
        // The minimum coverage cutoff
        val minCov
        // The minimum base call quality for a site to be
        // used in variant calling
        val minBQ
        // The minimum mapping quality for a site to be used in variant
        // calling
        val minMapQ
        // The existing summary string.
        val existingSummary
    output:
        // Tuple contains the file basename, alignment bamfile, filtered snp vcf, and filtered indel vcf
        tuple val(base), file(bam), file(correctedRef), file("${base}-snps-filtered.vcf"), file("${base}-indels-filtered.vcf")
        // Tuple contains the unfiltered vcf and multiallelic filtered vcf files.
        tuple file("${base}.vcf"), file("${base}-biallelic.vcf")
        // The summary string with the number of snps and indels added.
        env summary

    publishDir "${outDir}/${base}-Intermediate-Files", mode: 'copy', pattern: "*.vcf"

    script:
    /*
    Beings by indexing the bam file using samtools.

    Next, variants are called using freebayes. Freebayes can produce multiallelic
    variants (variants were multiple alternative alleles are provided). Thus, a 
    custom python script is used to select only the most prevalent allele at these 
    multiallelic positions. The INFO fields are corrected to account for this
    filtering step. (NOTE: this script also removes genotype information as
    this data is not used for downstreadm processing).

    VCFtools is used to split indels and snps apart, and each set of variatns are
    compressed with bgzip and tab indexed. (Required for filtering)

    BCFtools is then used to remove variants at position with less than the
    minimum coverage and that are not majority. The alternative and reference
    allele frequencies are calculated by dividing the occurrances (AO = alternative occurrance
    and RO = reference occurance) by the depth.

    The number of snps and indels is determined by grabbing all lines that do not begin with an '#'
    character (only variant entries do not begin with # in a VCF file). These values are
    added to the summary string.
    */
    """
    #!/bin/bash

    samtools index ${bam}
    freebayes -q ${minBQ} -m ${minMapQ} -f ${correctedRef} ${bam} > ${base}.vcf

    python3 ${baseDir}/scripts/remove_genotype_vcf.py -i ${base}.vcf -o ${base}-no-genotype.vcf

    python3 ${baseDir}/scripts/fix_multi_allelic.py -i ${base}-no-genotype.vcf -o ${base}-biallelic.vcf

    vcftools --keep-only-indels --vcf ${base}-biallelic.vcf --recode --recode-INFO-all --stdout > ${base}-indels.vcf
    bgzip ${base}-indels.vcf
    tabix ${base}-indels.vcf.gz

    vcftools --remove-indels --vcf ${base}-biallelic.vcf --recode --recode-INFO-all --stdout > ${base}-snps.vcf
    bgzip ${base}-snps.vcf
    tabix ${base}-snps.vcf.gz

    bcftools view -i "(INFO/DP >= ${minCov}) && ((INFO/AO / INFO/DP) > (INFO/RO / INFO/DP))" ${base}-indels.vcf.gz > ${base}-indels-filtered.vcf
    num_indels=\$(grep -v "^#" ${base}-indels-filtered.vcf | wc -l)

    bcftools view -i "(INFO/DP >= ${minCov}) && ((INFO/AO / INFO/DP) > (INFO/RO / INFO/DP))" ${base}-snps.vcf.gz > ${base}-snps-filtered.vcf
    num_snps=\$(grep -v "^#" ${base}-snps-filtered.vcf | wc -l)

    summary="${existingSummary},\$num_snps,\$num_indels"
    """
}

// Compiles a consensus by masking the reference and applying the variants.
process Generate_Consensus {
    input:
        // Tuple contains the file basename, the alignment bam, the snp vcf file
        // and the indel vcf file
        tuple val(base), file(bam), file(correctedRef), file(snps), file(indels)
        // The name of the base directory
        val baseDir
        // The name of the base directory
        val outDir
        // The minimum coverage threshold
        val minCov
        // The minimum base call quality for a site to be considered in
        // depth masking
        val minBQ
        // The minimum mapping quality for a site to be used in variant
        // calling
        val minMapQ
        // The existing summary string.
        val existingSummary
    output:
        // The consensus file created.
        file "${base}-consensus.fasta"
        // A file containing the sites that were masked in bed format
        file "${base}-mask-sites.bed"
        // The summary string with the number of masked positions and coverage
        // added.
        env summary

    publishDir "${outDir}", mode: 'copy', pattern: "${base}-consensus.fasta"
    publishDir "${outDir}/${base}-Intermediate-Files", mode: "copy", pattern: "${base}-mask-sites.bed"

    script:
    /*
    The script first computes sites to mask by identifying sites that have less than
    the minimum coverage provided. However, there are formatting issues with the pileup
    format that make this difficult. Samtools mpileup's output has the depth in 0-based
    format, which makes it impossible to distinguish between a site with 0 and 1 coverage.

    Thus, the pipeline instead makes use of bedtools subtract. It first creates a pileup for only sites with
    coverage, uses an in-house script to filter sites with less than the minimum coverage
    and converts these into a bed file.

    Next, the pipeline creates a pileup containing every site, and converts this into a bed file using
    the in-house script. 

    Finally, the sites we want to keep (those above the minimum coverage threshold) are substracted
    from the bed file with every site, to give us the low-coverage sites.

    Additionally, there is an interesting case when the sites that fall within a deletion are marked as masked.
    Because masking is applied first, this will cause an error when applying the variants (as the deletion site will contain 
    an N character and will not match the VCF reference). Thus, the pipeline uses bcftools to create a bed file for all indel sites,
    and then subtracts these from the low coverage sites. Now, this results in the sites to mask.

    The bedtools maskfasta command is then used to mask the reference at these positions.

    Then the variants are applied to the mask fasta. The reason this is done after masking is 
    because the pileup (and therefore masking) positions do not account for indels, which would
    shift the genomic coordinates (we would end up masking things we did not want to).

    The fasta is then wrapped using bioawk to make the sequence one line.

    Finally, the coverage is calculated by grabbing the sequence length using bioawk,
    and subtracting the masked positions divided by the length from 1 using awk.

    The number of masked sites and coverage are then added to the summary string. 
    */
    """
    #!/bin/bash

    samtools mpileup --no-BAQ -d 100000 -x -A -q ${minMapQ} -Q ${minBQ} -f ${correctedRef} ${bam} > ${base}.pileup
    python3 ${baseDir}/scripts/pileup_to_bed.py -i ${base}.pileup -o passed-sites.bed --minCov ${minCov}

    bioawk -c fastx '{print \$name"\t0\t"length(\$seq)}' ${correctedRef} > all-sites.bed

    if [[ -s all-sites.bed ]]; then

        bedtools subtract -a all-sites.bed -b passed-sites.bed > ${base}-low-cov-sites.bed

        bcftools query -f'%CHROM\t%POS0\t%END\n' ${indels} > indel-sites.bed

        bedtools subtract -a ${base}-low-cov-sites.bed -b indel-sites.bed > ${base}-mask-sites.bed
    else
        bioawk -c fastx '{print \$name"\t0\t"length(\$seq)}' ${correctedRef} > ${base}-mask-sites.bed
    fi

    num_mask=\$(bioawk -c bed 'BEGIN{SITES=0} {SITES+=\$end-\$start } END{print SITES}' ${base}-mask-sites.bed)

    bedtools maskfasta -fi ${correctedRef} -bed ${base}-mask-sites.bed -fo masked.fasta

    bgzip ${snps}
    tabix ${snps}.gz

    bgzip ${indels}
    tabix ${indels}.gz

    bcftools consensus -f masked.fasta ${snps}.gz > with-snps.fasta

    bcftools consensus -f with-snps.fasta ${indels}.gz > with-indels-snps.fasta
    
    bioawk -c fastx '{ gsub(/\\n/,"",seq); print ">${base}-"\$name; print \$seq }' with-indels-snps.fasta > ${base}-consensus.fasta

    seq_len=\$(bioawk -c fastx '{ print length(\$seq) }' < ${base}-consensus.fasta)

    coverage=\$(python3 ${baseDir}/scripts/calculate_genome_coverage.py -i ${base}-consensus.fasta)

    summary="${existingSummary},\$num_mask,\$coverage"
    """
}

// Writes a line to the summary file for the sample.
process Write_Summary {
    input:
        // The summary string containing statistics collected as
        // the pipeline ran.
        val summary
        // The output directory.
        val outDir

    script:
    /*
    The summary statistics are written to the summary file.
    */
    """
    #!/bin/bash

    echo "${summary}" >> ${outDir}/stats-summary.csv
    """  

}
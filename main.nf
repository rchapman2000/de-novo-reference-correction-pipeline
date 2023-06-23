#!/usr/bin/env/ nextflow

nextflow.enable.dsl=2

def helpMessage() {
    log.info"""
De Novo Reference Correction Assembly Pipeline:
    
Takes an input of paired-end fastq files from illumina
sequencers, processes the reads, removes host contamination,
and performs a de novo assembly.

USAGE: nextflow run main.nf [options] --input INPUT_DIR --output OUTPUT_DIR [--host_fasta HOST_FASTA | --host_bt2_index INDEX_DIR]

OPTIONS:

--input INPUT_DIR - [Required] A directory containing paired-end fastq files

--output OUTPUT_DIR - [Required] A directory to place output files (If not existing, pipeline will create)

--reference REFERENCE_FASTA - [Required] A FASTA file containing a reference genome to use for assembly.

OPTIONAL:
    
    --dedupe - If supplied, the pipeline will use prinseq to remove PCR duplicates from the reads [Default = OFF]

    --localAlignment - Supplying this option enables bowtie2's local alignment mode (will be used for host removal is that option is supplied as well) [Default= end-to-end mode]

    --host_fasta HOST_FASTA - A host fasta file that the reads will be aligned to to remove host contamination
            
    --host_bt2_index INDEX_DIRECTORY - To save time, an existing bowtie2 index can be supplied. Must be in its own directory

    --minCov INT - The minimum coverage below which a position will be masked [Default = 20]

    --minBQ INT - The minimum base call quality for a site to be considered in variant calling and depth-masking [Default = 10]

    --minMapQ INT - The minimum mapping quality for a site to be considered in variant calling and depth-masking [Default = 0]

    --minLen INT - the minimum length of a read to keep post trimming [Default = 75bp]

    --threads INT - the number of threads that can be use to run pipeline tools in parallel
    """
}

// Dynamically creates the header for the summary file based on
// the options supplied by the user.
def createSummaryHeader(dedupe, hostRef, hostIdx) {
    
    // Statistics regarding the sample, number of raw reads
    // and number of trimmed reads will always be collected.
    FinalHeader = "Sample,Raw Reads,Trimmed Reads,"

    // If the user supplied the --dedupe option, this 
    // process will occur next and a statistic of the number
    // of reads post-deduplication will be collected.
    if (dedupe != false) {
        FinalHeader = FinalHeader + "Deduped Reads,"
    }

    // If the user supplied an option for host removal (either a
    // FASTA file or bowtie2 index), then this will occur next, and the number
    // of reads post host removal will be collected.
    if (hostRef != false || hostIdx != false) {
        Finalheader = FinalHeader + "Non-Host Reads"
    }

    // Finally, the following statistics will be collected no matter what options the
    // user supplies.
    FinalHeader = FinalHeader + "Contigs,Scaffolds,Contigs Aligned to Reference,Variants Applid to Reference,Reads Mapped to Corrected Reference,Average Read Depth,SNPs,Indels,Masked Sites,Coverage"

    // Return the assembled header.
    return FinalHeader
}

// If the help parameter is supplied, link display the help message
// and quit the pipeline
params.help = false
if (params.help) {
    helpMessage()
    exit(0)
}

// Defines input parameters. Setting to false by default
// allows us to check that these have been set by the user.
params.input = false
params.reference = false
params.output = false
params.dedupe = false
params.localAlignment = false
params.minCov = 20
params.minLen = 75
params.minBQ = 10
params.minMapQ = 0
params.host_bt2_index = false
params.host_reference = false
params.threads = 1

// Grab the Illumina sequencing adapters from a
// fasta file included with the pipeline.
adapters = file("${baseDir}/data/adapters.fa")

// Import workflow modules from the modules.nf file.
include { Setup } from "./modules.nf"
include { QC_Report } from "./modules.nf"
include { Trimming } from "./modules.nf"
include { QC_Report as QC_Report_Trimmed } from "./modules.nf"
include { Remove_PCR_Duplicates } from "./modules.nf"
include { QC_Report as QC_Report_Deduped } from "./modules.nf"
include { Host_Read_Removal } from "./modules.nf"
include { QC_Report as QC_Report_Host_Removed } from "./modules.nf"
include { Spades_Assembly } from "./modules.nf"
include { Contig_Alignment } from "./modules.nf"
include { Correct_Reference_with_Contigs } from "./modules.nf"
include { Align_Reads_to_Corrected_Reference } from "./modules.nf"
include { Call_Variants } from "./modules.nf"
include { Generate_Consensus } from "./modules.nf"
include { Write_Summary } from "./modules.nf"



// Checks the input parameter
if (params.input == false) {
    // If the parameter is not set, notify the user and exit.
    println "ERROR: No input directory provided. Pipeline requires an input directory."
    exit(1)
}
else if (!(file(params.input).isDirectory())) {
    // If the input directory is not set, notify the user and exit.
    println "ERROR: ${params.input} is not an existing directory."
    exit(1)
}

// Create a channel for hte input files.
inputFiles_ch = Channel
    // Pull from pairs of files (illumina fastq files denoted by having R1 or R2 in
    // the file name).
    .fromFilePairs("${params.input}*_R{1,2}*.fastq.gz")
    // The .fromFilePairs() function spits out a list where the first 
    // item is the base file name, and the second is a list of the files.
    // This command creates a tuple with the base file name and two files.
    .map { it -> [it[0], it[1][0], it[1][1]]}

// Checks the output parameter.
outDir = ''
if (params.output == false) {
    // If the parameter is not set, notify the user and exit.
    println "ERROR: No output directory provided. Pipeline requires an output directory."
    exit(1)
}
else {
    // If the parameter is set, convert the value provided to a file type
    // to get the absolute path, and then convert back to a string to be
    // used in the pipeline.
    outDir = file(params.output).toString()
    println "Files will be written to ${outDir}"
}


// Checks the reference parameter. For this, we cannot use an
// input channel like was used for the input files. Using an input channel
// will cause Nextflow to only iterate once as the reference 
// channel would only only have 1 file in it. Thus, we manually parse
// the reference file into a tuple.
refData = ''
refName = ''
if (params.reference == false) {
    // If the parameter is not set, notify the user and exit.
    println "ERROR: no reference file proivded. Pipeline requires a reference file."
    exit(1)
}
else if (!(file(params.reference).exists())) {
    // If the reference file provided does not exist, notify the user and exit.
    println "ERROR: ${params.reference} does not exist."
    exit(1)
}
else {
    // Process the reference file to be supplied to the index step.
    
    // Parse the file provided into a file object.
    ref = file(params.reference)

    // Grab the basename of the file.
    refName = ref.getBaseName()

    // Place the file basename and file object into
    // a tuple.
    refData = tuple(refName, ref)
}

// Handles the --dedupe parameter. The dedupeParamValue
// is passed to the analysis summary file to denote 
// whether or not the option was supplied.
dedupeParamValue = "DISABLED"
if (params.dedupe != false) {
    dedupeParamValue = "ENABLED"
}

// Parses the --localAlignment Parameter
// By default, the pipeline will use end-to-end
// alignment mode from bowtie2. Thus, we will
// set these values by default. One value stores the 
// actual parameter to be supplied to bowtie2, and the other
// stores the value to be written in the summary file.
alignmentParam = "--end-to-end"
alignmentModeSummary = "End-to-End Alignment Mode"
if (params.localAlignment != false) {
    // If the user supplied the --localAlignment option,
    // the parameter and summary value will be upated
    // to enable and reflect this option.
    alignmentParam = "--local"
    alignmentModeSummary = "Local Alignment Mode"
}

// Parses the host options (--host_reference and --host_bt2_index).
// Again, we cannot use a channel like we did for the input files
// because then Nextfow would only run other modules once.
// Thus, we need to manually create a tuple of input data to pass to the indexing
// step and the alignment step.
hostRefData = ''
hostRefIdxData = ''
hostRefName = 'NONE'
if ((params.host_reference != false) && (params.host_bt2_index != false)) {
    // If both options are supplied, notify the user and exit.
    println "ERROR: you have specified both a host fasta file and bowtie2 index. Please only supply one."
    exit(1)
}
else {
    // If the user supplied the --host_reference option
    if (params.host_reference != false) {
        if (!(file(params.host_reference).exists())) {
            // If the file supplied does not exist, notify the user and exit.
            println "ERROR: ${params.host_reference} does not exist."
            exit(1)
        }
        else {
            // Parse the file into a file object
            hostRef = file(params.host_reference)
            // Use the getBaseName() function to 
            // get the reference name. This will be
            // used to name the bowtie2 index.
            hostRefName = hostRef.getBaseName()
            // Place these both into a tuple.
            hostRefData = tuple(hostRefName, hostRef)
        }
    }
    // If the user supplied the --host_bt2_index
    else if (params.host_bt2_index != false) {
        if (!(file(params.host_bt2_index).exists())) {
            // If the index provided does not exist, notify the user and exit.
            println "Error: ${params.host_bt2_index} does not exist."
            exit(1)
        }
        else {
            // Parse the directory into a file object
            hostRefDir = file(params.host_bt2_index)
            println hostRefDir
            // Grab a list of file objects from the directory
            // ending in .bt2
            indexFiles = file("${hostRefDir}/*.bt2")
            if (indexFiles.size() == 0) {
                // If there are no file in the directory ending in bt2, notify the user and exit
                println "Index Directory provided (${params.host_bt2_index}) does not contain any bt2 files"
                exit(1)
            }
            else {
                // Use the getSimpleName() function to grab the base name
                // of the index files (getSimpleName() removes anything following
                // the last . in a file name.)
                hostRefName = indexFiles[0].getSimpleName()
                println hostRefName
                // Place the index dir and name into a tuple.
                hostRefIdxData = tuple(hostRefDir, hostRefName)
            }
        }
    }
}

// Creates the summary header based on the options provided.
summaryHeader = createSummaryHeader(params.dedupe, params.host_reference, params.host_bt2_index)

workflow {

    Setup( refName, dedupeParamValue, alignmentModeSummary, params.minLen, params.minCov, params.minBQ, params.minMapQ, hostRefName, summaryHeader, outDir )

    if (params.host_reference != false) {
        Index_Host_Reference( hostRefData, outDir, params.threads )
    }

    QC_Report( inputFiles_ch, outDir, "FASTQC-Pre-Processing")

    // Perform adapter and quality trimming with trimmomatic.
    Trimming( inputFiles_ch, outDir, adapters, params.minLen )

    // Use FASTQC to perform a QC check on the trimmed reads
    QC_Report_Trimmed( Trimming.out[0], outDir, "FASTQC-Trimmed")

    if (params.dedupe != false) {
        // Perform PCR Duplicate removal using prinseq.
        Remove_PCR_Duplicates( Trimming.out[0], outDir, Trimming.out[2] )

        QC_Report_Deduped( Remove_PCR_Duplicates.out[0], outDir, "FASTQC-Deduped")

        if (params.host_reference != false) {
            Host_Read_Removal( Remove_PCR_Duplicates.out[0], outDir, Index_Host_Reference.out, alignmentParam, params.threads, Remove_PCR_Duplicates.out[2] )

            QC_Report_Host_Removed( Host_Read_Removal.out[0], outDir, "FASTQ-Host-Removed")

            Spades_Assembly( Host_Read_Removal.out[0], outDir, params.threads, Host_Read_Removal.out[2] )
        }
        else if (params.host_bt2_index != false) {

            Host_Read_Removal( Remove_PCR_Duplicates.out[0], outDir, hostRefIdxData, alignmentParam, params.threads, Remove_PCR_Duplicates.out[2] )

            QC_Report_Host_Removed( Host_Read_Removal.out[0], outDir, "FASTQ-Host-Removed")

            Spades_Assembly( Host_Read_Removal.out[0], outDir, params.threads, Host_Read_Removal.out[2] )
        }
        else {
            Spades_Assembly( Remove_PCR_Duplicates.out[0], outDir, params.threads, Remove_PCR_Duplicates.out[2] )
        }
    }
    else {
        if (params.host_reference != false) {
            Host_Read_Removal( Trimming.out[0], outDir, Index_Host_Reference.out, alignmentParam, params.threads, Trimming.out[2] )

            QC_Report_Host_Removed( Host_Read_Removal.out[0], outDir, "FASTQ-Host-Removed")

            Spades_Assembly( Host_Read_Removal.out[0], outDir, params.threads, Host_Read_Removal.out[2] )
        }
        else if (params.host_bt2_index != false) {

            Host_Read_Removal( Trimming.out[0], outDir, hostRefIdxData, alignmentParam, params.threads, Trimming.out[2] )

            QC_Report_Host_Removed( Host_Read_Removal.out[0], outDir, "FASTQ-Host-Removed")

            Spades_Assembly( Host_Read_Removal.out[0], outDir, params.threads, Host_Read_Removal.out[2] )
        }
        else {
            Spades_Assembly( Trimming.out[0], outDir, params.threads, Trimming.out[2] )
        }
    }

    Contig_Alignment( Spades_Assembly.out[0], outDir, refData, Spades_Assembly.out[2] )

    Correct_Reference_with_Contigs( Contig_Alignment.out[0], refData, outDir, Contig_Alignment.out[1] )

    Align_Reads_to_Corrected_Reference( Correct_Reference_with_Contigs.out[0], alignmentParam, outDir, Correct_Reference_with_Contigs.out[2] )

    Call_Variants( Align_Reads_to_Corrected_Reference.out[0], baseDir, outDir, params.minCov, params.minBQ, params.minMapQ, Align_Reads_to_Corrected_Reference.out[1] )

    Generate_Consensus( Call_Variants.out[0], baseDir, outDir, params.minCov, params.minBQ, params.minMapQ, Call_Variants.out[2] )

    Write_Summary( Generate_Consensus.out[2], outDir )

}
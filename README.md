# Nextflow De Novo Reference Correction Pipeline
This pipeline performs a combination between a de novo and reference-based assembly from Illumina Data. The pipeline works by performing the following steps:

1. Adapter and Quality Trimming
2. (Optional) Removing PCR Duplicates
3. (Optional) Host Read Removal
4. De Novo Assembly using SPAdes
5. Align the contigs produced to a Reference Genome
6. Correct the Reference using the contigs
7. Align processed reads to the Corrected Reference Genome
8. Variant Calling
9. Consensus Generation

This method is useful because it is more sensitive to insertions and deletions than a simple reference based assembly, where these mutations will often result in alignment errors. The de novo assembly will capture these complex mutations and they will be applied to the reference genome. That way alignment of reads to these regions will be much smoother.

## Technical Considerations

Some aspects of this analysis are important to consider and understand when interpreting the output.

### Variant Calling
This pipeline calls majority variants. Thus, the variant allele must be applied in the final consensus so long as it is more prevalent than the reference allele.

As well, the pipeline uses an in-house script to filter multiallelic sites from the variant call files. A multiallelic site is one where multiple alternative alleles are present. See the "CCA,GCA,ACA,ACG,ATA" in the following example.
```
K03455.1        2911    .       ATG     CCA,GCA,ACA,ACG,ATA     10482.7 .       AB=0,0,0,0,0;ABP=0,0,0,0,0;AC=0,0,2,0,0;AF=0,0,1,0,0;AN=2;AO=1,2,310,1,2;CIGAR=3X,3X,1M2X,1M1X1M,2M1X;DP=316;DPB=316;DPRA=0,0,0,0,0;EPP=5.18177,3.0103,5.8122,5.18177,3.0103;EPPR=0;GTI=0;LEN=3,3,2,1,1;MEANALT=5,5,5,5,5;MQM=60,60,59.9032,60,60;MQMR=0;NS=1;NUMALT=5;ODDS=413.347;PAIRED=1,1,1,1,1;PAIREDR=0;PAO=0,0,0,0,0;PQA=0,0,0,0,0;PQR=0;PRO=0;QA=27,76,11724,39,72;QR=0;RO=0;RPL=1,2,165,1,2;RPP=5.18177,7.35324,5.8122,5.18177,7.35324;RPPR=0;RPR=0,0,145,0,0;RUN=1,1,1,1,1;SAF=0,1,178,1,1;SAP=5.18177,3.0103,17.8324,5.18177,3.0103;SAR=1,1,132,0,1;SRF=0;SRP=0;SRR=0;TYPE=complex,complex,mnp,snp,snp      GT:DP:AD:RO:QR:AO:QA:GL 3/3:316:0,1,2,310,1,2:0:0:1,2,310,1,2:27,76,11724,39,72:-1054.34,-1052.22,-1051.91,-1048.11,-1045.98,-1047.51,-93.3193,-91.2078,-87.0746,0,-1051.14,-1049.01,-1044.9,-90.1068,-1050.84,-1048.47,-1046.34,-1042.24,-87.4429,-1045.26,-1047.87
```

## Installation

To install this pipeline, enter the following commands:

1. Clone the repository
```
git clone https://github.com/rchapman2000/de-novo-reference-correction-pipeline
```

2. Create a conda environment using the provided environment.yml file
```
conda env create -f environment.yml
```

3. Activate the conda environment
```
conda activate DeNovoRefCorrect
```
### Updating the Pipeline
If you already have the pipeline installed, you can update it using the following commands:

1. Navigate to your installation directory
```
# Navigate to your installation directory
cd ref-based-assembly-pipeline
```

2. Use ```git pull``` to get the latest update
```
git pull
```
3. Activate the conda environment and use the environment.yml file to download updates
```
conda activate DeNovoRefCorrect

conda env update --file environment.yml --prune
```

## Usage
To run the pipeline, use the following command:
```
# You must either be in the same directory as the main.nf file or reference the file location.
nextflow run main.nf [options] --input INPUT_DIR --output OUTPUT_DIR --reference REFERENCE_FASTA
```

### Options
The pipeline also supports the following optional arguments:

| Option | Type | Description |
|---|---|---|
| --dedupe | *None* | If supplied, the pipeline will use prinseq to remove PCR duplicates from the reads [Default = OFF] |
| --localAlignment | *None* | Supplying this option switches Bowtie2's alignment mode to local alignment (by default, it uses end-to-end mode). For more on this see the bowtie2 manual. [Default = OFF (End-to-end mode)] |
| --host_reference | *FASTA File* | If provided a host fasta file, the pipeline will perform host read removal via alignment (Supplying this option will require a bowtie2 index to be built). |
| --host_bt2_index | *Directory* | Alternative to providing a fasta file, a pre-built bowtie2 index can be provided for host removal. The index must be in its own directory. |
| --minCov | *int* | The minimum depth of coverage, below which a a position will be masked. [Default = 20] |
| --minLen | *int* | The minimum length of a read (in base pairs) to keep post trimming. [Default = 75] |
| --minBQ | *int* |  The minimum base call quality for a site to be considered in variant calling and depth-masking [Default = 10] |
| --minMapQ | *int* | The minimum mapping quality for a site to be considered in variant calling and depth-masking [Default = 0] |
| --threads | *int* | The number of threads that can be use to run pipeline tools in parallel. [Default = 1] |

To view the list of options from the command line, use the following command:
```
nextflow run main.nf --help
```
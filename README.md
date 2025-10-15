[![Static Badge](https://img.shields.io/badge/DOI-10.1261%2Frna.061945.117-blue?style=for-the-badge)](https://rnajournal.cshlp.org/content/early/2017/11/07/rna.061945.117)

ShapeMapper automates the calculation of RNA chemical probing reactivities from mutational profiling (MaP) experiments, in which chemical adducts on RNA are detected as internal mutations in cDNA through reverse transcription and read out by massively parallel sequencing. While originally built for analyzing SHAPE structure probing data, ShapeMapper is broadly useful for other types probing experiments, and includes a DMS mode for analyzing DMS experiments. ShapeMapper performs

* Reference sequence correction
* Read basecall quality trimming
* Paired read merging (using BBmerge)
* Alignment to reference sequences (using bowtie2 or STAR)
* Enforcement of read location requirements and primer trimming
* Multinucleotide and ambiguously aligned mutation handling
* Post-alignment basecall quality filtering
* Chemical adduct location inference from detected mutations
* Mutation rates calculation from mutation counts and effective read depths
* Reactivity profile calculation and normalization
* Heuristic quality control checks



### Inputs
**Required**
target
* **Description:** _FASTA file or list of files (.fa or .fasta) containing one or more target DNA sequences ('T' not 'U'). Lowercase positions will be excluded from reactivity profile, and should be used to indicate primer-binding sites if using directed primers. If multiple primer pairs were used, provide the primer sequences in a separate file with '--primers' ._
* **Format:** _FASTA_
* **Default:** -

sample:
modified_fastq 
* **Description:** _SHAPE Modified RNA. Please provide a folder directory for the modified fastq files._
* **Format:** _FASTQ_
* **Default:** TRUE

untreated_fastq
* **Description:** _Untreated/Unmodified RNA. Please provide a folder directory for the unmodified fastq files._
* **Format:** _FASTQ DIRECTORY_
* **Default:** TRUE

**Optional**
denatured_fastq
* **Description:** _Denatured RNA. Please provide a folder directory for the denatured fastq files._
* **Format:** _FASTQ_
* **Default:** TRUE

dms
* **Description:** _Run using DMS mode. Data will be normalized on a per-nucleotide basis. Note that this option is optimized for DMS-MaP data collected using Marathon and/or TGIRT enzmyes._
* **Format:** _CHECKBOX_
* **Default:** TRUE

aligner
* **Description:** _Use Bowtie2 or STAR for sequence alignment. STAR is recommended for sequences longer than several thousand nucleotides. STAR may also be slightly less sensitive than Bowtie2 (fewer aligned reads)._
* **Format:** _DROPDOWN_
* **Options:** _"BowTie2", "STAR"_
* **Default:** "BowTie2"

amplicon
* **Description:** _Require reads to align near expected primer pair locations, and intelligently trim primer sites. If a single pair of primers on the ends of the RNA sequence is used, simply set primer sequences to lowercase in the '--target' fasta file. If multiple pairs or internal locations are needed, specify primers with a '--primers' file._
* **Format:** _CHECKBOX_
* **Default:** FALSE

min_depth
* **Description:** _Minimum effective sequencing depth for including data (threshold must be met for all provided samples). Default=5000_
* **Format:** _INPUT_
* **Default:** 5000

Other Parameters
* **Description:** _Please see the other main parameters that can be entered here from ShapeMapper2 documentation._
* **Format:** _INPUT_
* **Default:** ""
* **Details:**
```
--log:  Location of output log file. Default="<name>_shapemapper_log.txt"
--verbose: Display full commands for each executed process, and display more process output messages in the event of an error. Default=False
--random-primer-len <n>: Length of random primers used (if any). Mutations within (length+1) of the 3-prime end of a read will be excluded, as will read depths over this region. Unused if '--amplicon' and/or '--primers' are provided. Default=0
--preserve-order: Preserve the order of input reads through all analysis stages. May slow down execution, but can be useful for debugging. Default=False
--max-paired-fragment-length <n>: Maximum distance between aligned ends of non-overlapping mate pairs to be merged into a single read (analogous to bowtie2 '--maxins'). Default=800
--max-search-depth <n>: Set bowtie2 '-D' parameter. If negative, shapemapper calls bowtie2 with a default -D 15. Unused with --star-aligner. Default=-1
--max-reseed <n>: Set bowtie2 '-R' parameter. If negative, shapemapper calls bowtie2 with a default -R 2. Unused with --star-aligner. Default=-1
--min-mapq <n>: Minimum aligner-reported mapping quality for included reads. Default=10 Note: When using Bowtie2, mutations contribute to lower mapping quality. Therefore, raising this threshold will have the side effect of excluding highly mutated reads.Note: This option does not apply to sequence correction, which uses a threshold of 10 regardless of this option
--min-qual-to-trim <n>: Minimum phred score in initial basecall quality trimming. Default=20
--window-to-trim <n>: Window size in initial basecall quality trimming. Default=5
--min-length-to-trim: Minimum trimmed read length in initial basecall quality trimming. Default=25
--min-qual-to-count <n>: Only count mutations with all basecall quality scores meeting this minimum score (including the immediate upstream and downstream  basecalls). This threshold is also used when computing the effective read depth. Default=30
--min-seq-depth <n>: Minimum sequencing depth for making a sequence correction (with '--correct-seq'). Default=50
--min-freq <n>: Minimum mutation frequency for making a sequence correction (with '--correct-seq'). Default=0.6
--disable-soft-clipping: Disable soft-clipping (i.e. perform end-to-end rather than local alignment). Default=False Note: this does not apply to sequence correction, which uses soft-clipping regardless.
--right-align-ambig: Realign ambiguous deletions/insertions to their rightmost valid position instead of leftmost. Not recommended, since left-realignment produces empirically better reactivity profiles than right-realignment. Default=False
--min-mutation-separation <n>: For two mutations to be treated as distinct, they must be separated by at least this many unchanged reference sequence nucleotides. Otherwise, they will be merged and treated as a single mutation. Does not apply to sequence correction. Default=6
--max-pages <n>: Maximum pages to render for '--render-mutations'. Default=100
--render-must-span <n>-<n>: Only render reads that cover a given nucleotide range. Disabled by default
--per-read-histograms: Output read length and per-read mutation frequency histogram tables in log file.
--serial: Run pipeline components one at a time and write all intermediate files to disk. Useful for debugging, but not generally recommended, as this will use large amounts of disk space. Default=False
--N7: Add N7 information to data visualization. Adds a graph of mutation rates and reactivities specific to N7 data in profiles.pdf. Prior to usage, ensure proper protocol was followed to generate valid N7 data.
--output-temp: Preserves temp files. Default=False.
--pernt-norm-factor-threshold: Set the number of NTs needed for effective per-nt normalization factor calculation. May need to change in the case of short RNAs. Default=20
--ignore_low_N7: Bypass N7 quality control filters.
--bypass_filters: Bypass N7 quality control filters and set threshold for NTs needed for effective per-nt normalization factor calculation to 1. (Equivalent to "--ignore_low_N7 --theshold 1")
```
### Usage examples
(Note: commandline argument examples only; will not produce output. 
 For a runnable example, execute `run_example.sh`)

Three-sample experiment, input FASTQ files:

  shapemapper --name example --target TPP.fa --out TPP_shapemap --amplicon --modified --R1 TPPplus_R1.fastq.gz --R2 TPPplus_R2.fastq.gz --untreated --R1 TPPminus_R1.fastq.gz --R2 TPPminus_R2.fastq.gz --denatured --R1 TPPdenat_R1.fastq.gz --R2 TPPdenat_R2.fastq.gz

Two-sample experiment, input from folders:

  shapemapper --name example2 --target TPP.fa --out TPP_shapemap --amplicon --modified --folder TPPplus --untreated --folder TPPminus

Only generate corrected sequence:

  shapemapper --name example3 --target TPP.fa --out TPP_mutant --amplicon --correct-seq --folder sequence_variant/A100

Generate corrected sequence using untreated sample,
then perform SHAPE-MaP analysis:

  shapemapper --name example4 --target TPP.fa --out --amplicon TPP_mutant --correct-seq --folder TPPminus --modified --folder TPPplus --untreated --folder TPPminus --denatured --folder TPPdenat

Multiple RNAs, randomly-primed experiment, STAR aligner:

  shapemapper --name example5 --target 16S.fa 23S.fa --out ribosome --random-primer-len 9 --star-aligner --modified --folder ribosome_plus --untreated --folder ribosome_minus --denatured --folder ribosome_denat

Process single DMS modified sample using DMS mode:

shapemapper --name example6 --target add.fa --out add_dms --dms --amplicon --modified --folder ribosome_plus

### Outputs
**Reported**
pdf
* **Description:** _PDF reports of drawn plots_
* **Format:** _PDF_
* **Location:** "shapemapper_plots/shapemapper_out/"

modified parsed mut
* **Description:** _modified fastq output data_
* **Format:** _MUT_
* **Location:** "modified_mut/shapemapper_out/"

untreated parsed mut
* **Description:** _untreated fastq output data_
* **Format:** _MUT_
* **Location:** "untreated_mut/shapemapper_out/"
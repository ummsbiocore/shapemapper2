$HOSTNAME = ""
params.outdir = 'results'  


if (!params.untreated_fastq){params.untreated_fastq = ""} 
if (!params.denatured_fastq){params.denatured_fastq = ""} 
if (!params.target_fasta){params.target_fasta = ""} 
if (!params.modified_fastq){params.modified_fastq = ""} 
if (!params.named){params.named = ""} 
if (!params.mate){params.mate = ""} 
// Stage empty file to be used as an optional input where required
ch_empty_file_1 = file("$baseDir/.emptyfiles/NO_FILE_1", hidden:true)
ch_empty_file_2 = file("$baseDir/.emptyfiles/NO_FILE_2", hidden:true)

g_15_1_g_2 = file(params.untreated_fastq, type: 'any')
g_16_2_g_2 = params.denatured_fastq && file(params.denatured_fastq, type: 'any').exists() ? file(params.denatured_fastq, type: 'any') : ch_empty_file_1
g_17_3_g_2 = file(params.target_fasta, type: 'any')
g_18_0_g_2 = file(params.modified_fastq, type: 'any')
Channel.value(params.named).set{g_19_4_g_2}
Channel.value(params.mate).set{g_20_5_g_2}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 4
    $MEMORY = 40
}
//* platform
//* platform
//* autofill

process shapemapper {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*$/) "shapemapper_out/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_Modified_.*.mut$/) "modified_mut/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_Untreated_.*.mut$/) "untreated_mut/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_profile.txt$/) "profile/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_flowchart.svg$/) "flowchart/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.pdf$/) "ShapeMapper_Plots/$filename"}
input:
 path modified_fastq
 path untreated_fastq
 path denatured_fastq
 path target_fasta
 val named
 val mate

output:
 path "*"  ,emit:g_2_directory00 
 path "*_Modified_*.mut" ,optional:true  ,emit:g_2_mutation11 
 path "*_Untreated_*.mut" ,optional:true  ,emit:g_2_mutation22 
 path "*_profile.txt"  ,emit:g_2_outputFileTxt33 
 path "*_flowchart.svg" ,optional:true  ,emit:g_2_outputFile44 
 path "*.pdf"  ,emit:g_2_outputFilePdf55 

container 'quay.io/ummsbiocore/shapemapper2:1.1.2'

when:
!params.run_shapedance || (params.run_shapdance == "yes")

script:
threads = task.cpus

mate_folder = (mate == "pair") ? "--folder" : "--unpaired-folder"

denatured_input = denatured_fastq.toString().startsWith("NO_FILE") ? "" : "--denatured ${mate_folder} ${denatured_fastq}"
name_input = (named == "") ? "" : "--name ${named}"

dms = params.shapemapper.dms
dms_option = (dms == "true") ? "--dms" : ""

indiv_norm = params.shapemapper.indiv_norm
indiv_norm_option = (indiv_norm == "true") ? "--indiv-norm" : ""

max_bg = params.shapemapper.max_bg

aligner = params.shapemapper.aligner
genomeSAindexNbase = params.shapemapper.genomeSAindexNbase
rerun_on_star_segfault = params.shapemapper.rerun_on_star_segfault
rerun_genomeSAindexNbase = params.shapemapper.rerun_genomeSAindexNbase
star_aligner = (aligner == "STAR") ? "--star-aligner --genomeSAindexNbase ${genomeSAindexNbase} --rerun-on-star-segfault ${rerun_on_star_segfault} --rerun-genomeSAindexNbase ${rerun_genomeSAindexNbase}" : ""

sh_amplicon = params.shapemapper.sh_amplicon
amp_primers = params.shapemapper.amp_primers
max_primer_offset = params.shapemapper.max_primer_offset
amplicon_option = (sh_amplicon == "true") ? "--amplicon --primers ${amp_primers} --max-primer-offset ${max_primer_offset}" : ""
amplicon_option = (sh_amplicon == "true") ? (amp_primers.getName().startsWith('NO_FILE')) ? "" : "--amplicon --primers ${amp_primers} --max-primer-offset ${max_primer_offset}" : ""
min_depth = params.shapemapper.min_depth

render_fc = params.shapemapper.render_fc
rfc = (render_fc == "true") ? "--render-flowchart" : ""
processed_reads = params.shapemapper.processed_reads
prs = (processed_reads == "true") ? "--output-processed-reads" : ""
aligned_reads = params.shapemapper.aligned_reads
ars = (aligned_reads == "true") ? "--output-aligned-reads" : ""
parsed_muts = params.shapemapper.parsed_muts
pmut = (parsed_muts == "true") ? "--output-parsed-mutations" : ""
counted_muts = params.shapemapper.counted_muts
cmut = (counted_muts == "true") ? "--output-counted-mutations" : ""

other_parameters = params.shapemapper.other_parameters
//* @style @condition:{aligner="STAR", genomeSAindexNbase, rerun_on_star_segfault, rerun_genomeSAindexNbase},{aligner="Bowtie2"},{sh_amplicon="true", amp_primers, max_primer_offset},{sh_amplicon="false"} @multicolumn:{dms, max_bg},{aligner, genomeSAindexNbase, rerun_on_star_segfault, rerun_genomeSAindexNbase},{sh_amplicon, amp_primers, max_primer_offset},{processed_reads, aligned_reads, parsed_muts, counted_muts}

"""
mkdir -p shapemapper_out/

shapemapper ${name_input} --target ${target_fasta} \
	${amplicon_option} \
	--modified ${mate_folder} ${modified_fastq} \
	--untreated ${mate_folder} ${untreated_fastq} \
	${denatured_input} \
	--out ./ \
	--min-depth ${min_depth} \
	--max-bg ${max_bg} \
	--nproc ${threads} \
	${star_aligner} ${dms_option} ${indiv_norm_option} \
	${rfc} ${prs} ${ars} ${pmut} ${cmut} \
	${other_parameters}
"""

}


workflow {

g_19_4_g_2= g_19_4_g_2.ifEmpty("") 


if (!(!params.run_shapedance || (params.run_shapdance == "yes"))){
g_18_0_g_2.set{g_2_directory00}
g_15_1_g_2.set{g_2_directory00}
g_16_2_g_2.set{g_2_directory00}
g_17_3_g_2.set{g_2_directory00}
g_2_mutation11 = Channel.empty()
g_2_mutation22 = Channel.empty()
g_2_outputFileTxt33 = Channel.empty()
g_2_outputFile44 = Channel.empty()
g_2_outputFilePdf55 = Channel.empty()
} else {

shapemapper(g_18_0_g_2,g_15_1_g_2,g_16_2_g_2,g_17_3_g_2,g_19_4_g_2,g_20_5_g_2)
g_2_directory00 = shapemapper.out.g_2_directory00
g_2_mutation11 = shapemapper.out.g_2_mutation11
g_2_mutation22 = shapemapper.out.g_2_mutation22
g_2_outputFileTxt33 = shapemapper.out.g_2_outputFileTxt33
g_2_outputFile44 = shapemapper.out.g_2_outputFile44
g_2_outputFilePdf55 = shapemapper.out.g_2_outputFilePdf55
}


}

workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}

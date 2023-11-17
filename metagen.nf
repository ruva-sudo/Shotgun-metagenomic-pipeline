nextflow.enable.dsl = 2

params.reads = "$projectDir/metagen/*_R{1,2}.fq.gz"
params.outdir = "results"
params.cpus = 2
params.mem = 80
params.fastas = "contigs.fasta"
params.blast_db = "$projectDir/16S_ribosomal_RNA/*"
params.threads = 24
params.quastref = "$projectDir/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/97/silva_132_97_16S.fna"

reads_pair_ch = Channel.fromFilePairs(params.reads).dump()
blastribosomal_ch = Channel.fromPath(params.blast_db)
reference_ch = Channel.fromPath(params.quastref)

process ASSEMBLY {

	publishDir "${params.outdir}/assembly", mode: 'copy'

	input:
	tuple val(pair_id), path(reads)

	output:
	path pair_id,										emit: full_folder
	tuple val(pair_id), path("${pair_id}.contigs.fasta"),					emit: contig_file
	tuple val(pair_id), path("${pair_id}/corrected/*_R1.fq.00.0_0.cor.fastq.gz"),	emit: first_corrected_reads
	tuple val(pair_id), path("${pair_id}/corrected/*_R2.fq.00.0_0.cor.fastq.gz"),	emit: second_corrected_reads

	script:
	"""
	mkdir -pv assembly

	metaspades.py --threads $task.cpus -m ${task.memory.toGiga()} -1 ${reads[0]} -2 ${reads[1]} -o $pair_id

	cp $pair_id/contigs.fasta ${pair_id}.contigs.fasta
	"""
}

process BINNING_MAXBIN {

	publishDir "${params.outdir}/maxbin2", mode: 'copy'

	input:
	tuple val(pair_id), path(contigs)
	tuple val(pair_id), path(first_corr_reads)
	tuple val(pair_id), path(second_corr_reads)

	output:
	tuple val(pair_id), path("${pair_id}.*.fasta"), 		emit: binned_fastas		  
	tuple val(pair_id), path("*.summary"), 				emit: binning_summary	  
	tuple val(pair_id), path("*.log"), 				emit: bin_log
	tuple val(pair_id), path("*.noclass"), 				emit: unbinned_fastas	  
	tuple val(pair_id), path("*.tooshort"),				emit: short_contigs
	tuple val(pair_id), path("*.marker_of_each_bin.tar.gz"),	emit: gene_file

	script:
	"""
	run_MaxBin.pl -contig ${contigs} \
		      -reads ${first_corr_reads}\
		      -reads2 ${second_corr_reads} \
		      -thread $task.cpus \
		      -out ${pair_id}
	"""
}

process FILTER_BINS {

	publishDir "${params.outdir}/filter_bins", mode: 'copy'

	input:
	tuple val(pair_id), path(summary)
	tuple val(pair_id), path(bins)

	output:
	tuple val(pair_id), path("passed_bins/${pair_id}.*.fasta"),	emit: passed_fastas
	tuple val(pair_id), path("failed_bins/"),			emit: failed_fastas

	script:
	"""
	bins_filter.py 
	"""	
}

process PROKKA {

	publishDir "${params.outdir}/prokka", mode: 'copy'

	input:
	tuple val(pair_id), path(binned_contigs)
	
	output:
	tuple val(pair_id), path("${prefix}/*.gff"), emit: gff
	tuple val(pair_id), path("${prefix}/*.gbk"), emit: gbk
	tuple val(pair_id), path("${prefix}/*.fna"), emit: fna
	tuple val(pair_id), path("${prefix}/*.faa"), emit: faa
	tuple val(pair_id), path("${prefix}/*.ffn"), emit: ffn
	tuple val(pair_id), path("${prefix}/*.sqn"), emit: sqn
	tuple val(pair_id), path("${prefix}/*.fsa"), emit: fsa
	tuple val(pair_id), path("${prefix}/*.tbl"), emit: tbl
	tuple val(pair_id), path("${prefix}/*.err"), emit: err
	tuple val(pair_id), path("${prefix}/*.log"), emit: log
	tuple val(pair_id), path("${prefix}/*.txt"), emit: txt
	tuple val(pair_id), path("${prefix}/*.tsv"), emit: tsv
	
	script:

	prefix=binned_contigs.getBaseName()
	"""
	prokka ${binned_contigs} \
		--prefix ${prefix} \
		--kingdom Bacteria \
		--cpus $task.cpus \
		--mincontiglen 200 \
		--centre XXX
	"""
	
}
 
/*
process MAKEBLAST_DB {

	input:
	path fasta

	output:
	path "blast_db", emit: database

	script:
	"""
	mkdir blast_db
	makeblastdb -dtype nucl \
		    -in ${fasta}
	
	mv ${fasta}* blast_db/
	"""
}
*/

process BLASTN {

	publishDir "${params.outdir}/blast_outdir", mode: 'copy'

	input:
	tuple val(pair_id), path(cds)
	path db

	output:
	tuple val(pair_id), path("${prefix}.blastn.txt"), emit: blast_out

	script:

	prefix=cds.getBaseName()
	"""
	blastn -db 16S_ribosomal_RNA \
		-query ${cds} \
		-num_threads $task.cpus \
		-task blastn \
		-dust no \
		-outfmt "7 delim=, qacc sacc evalue bitscore qcovus pident" \
		-max_target_seqs 1 \
		-out ${prefix}.blastn.txt
	"""
}

process DEEPARG_DOWNLOAD {

	output:
	path "db/", emit: deep_db

	script:
	"""
	deeparg \
		download_data \
		-o db/
	"""
}

process DEEPARG_PREDICT {

	publishDir "${params.outdir}/deeparg", mode: 'copy'

	input:
	tuple val(pair_id), path(orfs)
	path db

	output:
	tuple val(pair_id), path("*.align.daa"),	emit: daa
	tuple val(pair_id), path("*.align.daa.tsv"),	emit: daa_tsv
	tuple val(pair_id), path("*.ARG"),		emit: mapping
	tuple val(pair_id), path("*potential.ARG"),	emit: probability_arg

	script:
	prefix=orfs.getBaseName()
	"""
	deeparg predict \
		--model LS \
		-i ${orfs} \
		-o ${prefix} \
		-d ${db}
	"""
}

process QUAST {

	publishDir "${params.outdir}/quast", mode: 'copy'

	input:
	tuple val(pair_id), path(contigs)
	tuple val(pair_id), path(first_trimmed_reads)
	tuple val(pair_id), path(second_trimmed_reads)
	path reference 

	output:
	tuple val(pair_id), path("${pair_id}/combined_reference/"),	emit: combined_references 
	tuple val(pair_id), path("${pair_id}/icarus_viewers/"),		emit: icarus
	tuple val(pair_id), path("${pair_id}/not_aligned/"),		emit: not_aligned
	tuple val(pair_id), path("${pair_id}/icarus.html"),		emit: icarus_html
	tuple val(pair_id), path("${pair_id}/runs_per_reference/"),	emit: runs_per_reference
	tuple val(pair_id), path("${pair_id}/report.html"),		emit: report_hmtl
	tuple val(pair_id), path("${pair_id}/metaquast.log"),		emit: log_file
	tuple val(pair_id), path("${pair_id}/summary/"),		emit: summary

	script:
	"""
	metaquast.py ${contigs} \
			-r ${reference} \
			-1 ${first_trimmed_reads} \
			-2 ${second_trimmed_reads} \
			-t $task.cpus \
			-l SPAdes \
			-o ${pair_id}
	"""
}	
		

workflow {

	ASSEMBLY(reads_pair_ch)

	ASSEMBLY.out.contig_file.collect().dump(tag: "ASSEMBLY")
	
	BINNING_MAXBIN(ASSEMBLY.out.contig_file, ASSEMBLY.out.first_corrected_reads, ASSEMBLY.out.second_corrected_reads)

	FILTER_BINS(BINNING_MAXBIN.out.binning_summary, BINNING_MAXBIN.out.binned_fastas)

	bins_ch = FILTER_BINS.out.passed_fastas.transpose()

	PROKKA(bins_ch)

	//MAKEBLAST_DB(db_ch)

	blastdb_ch = blastribosomal_ch.collect()

	BLASTN(PROKKA.out.ffn, blastdb_ch)

	DEEPARG_DOWNLOAD()

	DEEPARG_PREDICT(PROKKA.out.ffn, DEEPARG_DOWNLOAD.out.deep_db)

	QUAST(ASSEMBLY.out.contig_file, ASSEMBLY.out.first_corrected_reads, ASSEMBLY.out.second_corrected_reads, reference_ch)
}
 

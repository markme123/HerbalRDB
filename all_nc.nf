
//--queueq
//--genome_ddir

Channel
    .fromFilePairs( "${params.genome_ddir}/*.{fa,fasta}", size: 1 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.genome_ddir}."}
    .set {pass_genome2}

Channel
    .fromFilePairs( "${params.genome_ddir}/*.{fa,fasta}", size: 1 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.genome_ddir}."}
    .set {pass_genome3}

process run_tRNA{
	publishDir "Result/01.nc/${gg}", mode: 'link'
	tag "${gg}"
	executor "slurm"
	errorStrategy "ignore"
	queue "${params.queueq}"
	cpus 25
	input:
		set val(gg),path(genome) from pass_genome2
	output:
		path "${gg}.tRNA.gff"
	script:
		"""
		tRNAscan-SE -E -o ${genome}.tRNA -f ${genome}.structure --thread 16 ${genome}
		perl /work/home/shuziqiang/Pipe/02.AreYouOK/01.nc/bin/tRNAscan/tRNAscan_to_gff3.pl ${genome}.tRNA ${genome}.structure > ${gg}.tRNA.gff
		"""
}

process run_cmscan{
	publishDir "Result/01.nc/${gg}", mode: 'link'
	tag "${gg}"
	executor "slurm"
	errorStrategy "ignore"
	queue "${params.queueq}"
	cpus 25
	input:
		set val(gg),path(genome) from pass_genome3
	output:
		path "*.rRNA_cmscan.gff"
        path "*.snRNA.gff"
        path "*.miRNA.gff"
        path "*.tblout"
	script:
		"""
		python3 ${projectDir}/bin/cmscan/Run_cmscan.py -r ${genome} -o ./ -c cmscan -Rfam_cm /work/home/shuziqiang/Pipe/02.AreYouOK/database/rfam/14.10/Rfam_part.cm -Rfam_clanin /work/home/shuziqiang/Pipe/02.AreYouOK/database/rfam/14.10/Rfam_part.clanin -s ${genome} -cm_tpye cut_ga -cpu 16
		python3 ${projectDir}/bin/cmscan/cmscan_ncRNA2gff.py -i ${genome}.tblout -t /work/home/shuziqiang/Pipe/02.AreYouOK/database/rfam/14.10/rfam_anno.xls -s ${gg}
		"""
}

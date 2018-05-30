import itertools
from itertools import chain

# build metadata
# run once outside of snakemake 
# bash ~/git/NGS_db/build_metadata.sh Huf > ~/git/ipsc_rpe_atac/metadata.csv

# read metadata into dictionary for snakemake to link samples to lane bams 
# one sample can have 1 to n lane bams
SAMPLE_PATH = dict()
metadata = open(config['metadata_file'])
for line in metadata:
	path = line.split(',')[3][:-1]
	sample = line.split(',')[1]
	# skip header
	if sample == 'Sample':
		continue
	if sample not in SAMPLE_PATH:
		SAMPLE_PATH[sample] = [path]
	else:
		old_path = SAMPLE_PATH[sample]
		old_path.append(path)
		SAMPLE_PATH[sample] = old_path

localrules: pull_lane_bams_from_nisc, retrieve_and_process_black_list, black_list, \
	ucsc_view, total_reads, union_peaks, merge_peaks, bootstrap_peaks, \
	peak_fasta_bootstrap, peak_fasta, remove_tss_promoters, build_tss_regions, \
	build_cisbp_master_file
	#find_closest_TSS, find_closest_TSS_bootstrap

wildcard_constraints:
	sample='|'.join(list(SAMPLE_PATH.keys())),
	lane_bam='|'.join([x.split('/')[-1].split('.bam')[0] for x in list(itertools.chain.from_iterable(SAMPLE_PATH.values()))])

rule all:
	input:
		expand('macs_peak/{sample}_peaks.blackListed.hg19.narrowPeak', sample = list(SAMPLE_PATH.keys())),
		'fastqc/multiqc/multiqc_report.html',
		'deeptools/multiBamSummary.tsv',
		'metrics/reads_by_sample.txt',
		expand('/data/mcgaugheyd/datashare/hufnagel/hg19/{sample}.bw', sample = list(SAMPLE_PATH.keys())),
		expand('/data/mcgaugheyd/datashare/hufnagel/hg19/{sample}_peaks.blackListed.hg19.narrowPeak.bb', sample = list(SAMPLE_PATH.keys())),
		expand('closest_TSS_motifs/{sample}.fimo.closestTSS.dat.gz', sample = list(SAMPLE_PATH.keys())),
		expand('closest_TSS_motifs/bootstrapping/{sample}.bootstrap_{bootstrap_num}.fimo.closestTSS.dat.gz', \
			sample = list(SAMPLE_PATH.keys()), \
			bootstrap_num = [str(x).zfill(3) for x in list(range(1,config['motif_bootstrap_num']))]),
		expand('fimo_motifs/bootstrapping_stats/{sample}.bootstrap_{bootstrap_num}.fimo_counts.dat.gz', \
			sample = list(SAMPLE_PATH.keys()), \
			bootstrap_num = [str(x).zfill(3) for x in list(range(1,config['motif_bootstrap_num']))]),
		expand('closest_TSS_motifs/processed/{sample}.fimo.closestTSS.processed.dat.gz',
			sample = list(SAMPLE_PATH.keys()))

rule pull_lane_bams_from_nisc:
	output:
		'original/{lane_file}.bam'
	run:
		lane_files = [x for x in list(itertools.chain.from_iterable(SAMPLE_PATH.values()))]
		for bam in lane_files:
			command = 'rsync -av trek.nhgri.nih.gov:' + bam + ' original/'
			echo_command = 'echo ' + command
			shell('echo ' + str(output))
			shell(echo_command)
			shell(command)
		

# bam to fastq to bwa to make bam, again
rule realign:
	input:
		'original/{lane_file}.bam'
	output:
		temp('realigned/{lane_file}.bam')
	threads: 8
	shell:
		"""
		module load {config[samtools_version]}
		module load {config[bwa_version]}
		tmp_name=$( echo {input} | sed 's/\///g' )
		samtools collate -uOn 128 {input} /tmp/TMP_$tmp_name | \
			samtools fastq - | \
			bwa mem -t {threads} -B 4 -O 6 -E 1 -M -p {config[bwa_genome]} - | \
			samtools view -1 - > \
			{output}
		"""

# if sample has more than one lane bam, then merge into one bam
rule merge_bam:
	input:
		lambda wildcards: expand('realigned/{lane_file}.bam', lane_file = [x.split('/')[-1][:-4] for x in SAMPLE_PATH[wildcards.sample]])
	output:
		bam = 'merged_bam/{sample}.bam',
		bai = 'merged_bam/{sample}.bam.bai'
	shell:
		"""
		module load {config[samtools_version]}
		module load {config[picard_version]}
		picard_i=""
		for bam in {input}; do
			picard_i+=" I=$bam"
		done
		java -Xmx8g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
			MergeSamFiles \
			TMP_DIR=/scratch/$SLURM_JOB_ID \
			$picard_i \
			O={output.bam}
		samtools index {output.bam}
		"""

# remove dups, only keep reads with over q 5
rule filter_bam:
	input:
		'merged_bam/{sample}.bam'
	output:
		metrics = 'metrics/{sample}.picard.metrics',
		bam = 'merged_bam_HQ/{sample}.q5.rmdup.bam',
		bai = 'merged_bam_HQ/{sample}.q5.rmdup.bam.bai'
	threads: 5
	shell:
		"""
		module load {config[samtools_version]}
		module load {config[picard_version]}
		samtools view -O bam -q 5 -b {input} | \
			samtools sort - -O bam -o {output.bam}TEMP
		samtools index {output.bam}TEMP
		java -Xmx4g -XX:ParallelGCThreads={threads} -jar $PICARDJARPATH/picard.jar \
			MarkDuplicates \
			INPUT={output.bam}TEMP \
			OUTPUT={output.bam} \
			REMOVE_DUPLICATES=true \
			METRICS_FILE={output.metrics}
		samtools index {output.bam}
		rm {output.bam}TEMP*
		"""

# downsample to n reads
rule downsample:
	input:
		bam = 'merged_bam_HQ/{sample}.q5.rmdup.bam',
		bai = 'merged_bam_HQ/{sample}.q5.rmdup.bam.bai'
	output:
		bam = 'downsample_bam/{sample}.q5.rmdup.ds.bam',
		bai = 'downsample_bam/{sample}.q5.rmdup.ds.bam.bai'
	shell:
		"""
		module load {config[samtools_version]}
		frac=$( samtools idxstats {input.bam} | cut -f3 | awk 'BEGIN {{total=0}} {{total += $1}} END {{frac=15000000/total;if (frac > 1) {{print 1}} else {{print frac}}}}' )
		samtools view -s $frac -b {input.bam} > {output.bam}
		samtools index {output.bam}
		"""

# bam to bigwig
# for UCSC visualization
# normalize by CPM
rule bam_to_bigWig:
	input:
		'downsample_bam/{sample}.q5.rmdup.ds.bam'
	output:
		bedgraph = temp('bigWig/{sample}.bG'),
		bw = 'bigWig/{sample}.bw'
	threads: 
		16
	shell:
		"""
		module load {config[deeptools_version]}
		module load ucsc
		bamCoverage --bam {input} -o {output.bedgraph} \
			--numberOfProcessors {threads} \
			--binSize 10 \
			--normalizeUsing RPGC \
			--effectiveGenomeSize 2864785220 \
			--ignoreForNormalization MT \
			--outFileFormat bedgraph

		 /home/mcgaugheyd/git/ChromosomeMappings/convert_notation.py \
			-c /home/mcgaugheyd/git/ChromosomeMappings/GRCh37_ensembl2UCSC.txt \
			-f <( grep ^[0-9XY] {output.bedgraph} ) | sort -k1,1 -k2,2n > {output.bedgraph}TEMP
		bedGraphToBigWig {output.bedgraph}TEMP /data/mcgaugheyd/genomes/hg19/hg19.chrom.sizes {output.bw}
		rm {output.bedgraph}TEMP
		"""

# basic stats on the bam file
# run on the bam file before q5, rmdup, and downsample happen
rule fastqc:
	input:
		'merged_bam/{sample}.bam'
	output:
		'fastqc/{sample}'
	threads: 8
	shell:
		"""
		module load fastqc
		mkdir -p {output}
		fastqc -t {threads} -o {output} {input}
		"""

# multiqc
rule multiqc_fastqc:
	input:
		expand('fastqc/{sample}', sample = list(SAMPLE_PATH.keys()))
	output:
		'fastqc/multiqc/multiqc_report.html'
	shell:
		"""
		module load multiqc
		multiqc fastqc/ -o fastqc/multiqc
		"""

# post filtering metrics 
rule total_reads:
	input:
		expand('merged_bam_HQ/{sample}.q5.rmdup.bam', sample = list(SAMPLE_PATH.keys()))
	output:
		'metrics/reads_by_sample.txt'
	shell:
		"""
		for i in {input}; do samtools idxstats $i | cut -f3 | awk 'BEGIN {{total=0}} {{total += $1}} END {{print total}}'; done > metrics/TEMP2
		for i in {input}; do echo $i; done > metrics/TEMP1
		paste metrics/TEMP1 metrics/TEMP2 > {output}
		rm metrics/TEMP1; rm metrics/TEMP2
		"""

# macs2
# ATAC-seq doesn't have control (no input sample). So run in single mode
# https://groups.google.com/forum/#!topic/macs-announcement/4OCE59gkpKY
rule peak_calling:
	input:
		'downsample_bam/{sample}.q5.rmdup.ds.bam'
	output:
		peaks = 'macs_peak/{sample}_peaks.xls',
		narrow_peaks = 'macs_peak/{sample}_peaks.narrowPeak',
		summits = 'macs_peak/{sample}_summits.bed'
	shell:
		"""
		module load {config[macs2_version]}
		macs2 callpeak -f BAM -g "hs" -t {input} -q 0.01 \
			--keep-dup all --nomodel --extsize 200 --shift -100 \
			-n {wildcards.sample} \
			--outdir macs_peak
		"""

# blacklist
# remove peaks called in ENCODE black list regions
# https://sites.google.com/site/anshulkundaje/projects/blacklists
rule retrieve_and_process_black_list:
	output:
		'/data/mcgaugheyd/genomes/GRCh37/ENCFF001TDO.bed.gz'
	shell:
		"""
		wget https://www.encodeproject.org/files/ENCFF001TDO/@@download/ENCFF001TDO.bed.gz
		mv ENCFF001TDO.bed.gz /data/mcgaugheyd/genomes/hg19/
		/home/mcgaugheyd/git/ChromosomeMappings/convert_notation.py \
			-c /home/mcgaugheyd/git/ChromosomeMappings/GRCh37_UCSC2ensembl.txt \
			-f /data/mcgaugheyd/genomes/hg19/ENCFF001TDO.bed.gz | gzip > \
			{output}
		"""

rule black_list:
	input:
		peaks = 'macs_peak/{sample}_peaks.narrowPeak',
		blacklist = '/data/mcgaugheyd/genomes/GRCh37/ENCFF001TDO.bed.gz'
	output:
		'macs_peak/{sample}_peaks.blackListed.narrowPeak'
	shell:
		"""
		module load {config[bedtools_version]}
		intersectBed -a {input.peaks} -b {input.blacklist} -v > {output}
		"""

# create 1000bp exclusion regions at the TSS
# these are promoters (which are also open)
# we want enhancers
# only use appris prinicpal transcripts
# map chr notation to num notation
# bed output
rule build_tss_regions:
	input:
		config['gtf_file']	
	output:
		gene = temp('annotation/gene.bed'),
		tss = temp('annotation/tss_hg19.bed'),
		tssT2 = temp('annotation/tss_hg19fix2.bed'),
		tssT3 = temp('annotation/tss_hg19fix3.bed'),
		tssG = 'annotation/tss_GRCh37.bed'
	shell:
		"""
		module load ucsc
		zcat {input} | \
			grep appris_principal | \
			awk '{{if ($3 != "gene") print $0;}}' | \
			grep -v "^#" | \
			gtfToGenePred /dev/stdin /dev/stdout | \
			genePredToBed stdin {output.gene}
		awk '{{if($6 == "-") {{$2 = $3 - 1}} else {{$3 = $2 + 1}} print}}' {output.gene} | \
			awk -v OFS='\t' '{{if($6 == "-") {{$3 = $2 + 1000}} else {{$2 = $3 - 1000}} print}}' > \
			{output.tss}
		awk '{{if($2 < 1) {{$2 = 1}} print}}' {output.tss} > {output.tssT2}
		awk '{{if($3 < 1) {{$3 = 1}} print}}' {output.tssT2} > {output.tssT3}
		/home/mcgaugheyd/git/ChromosomeMappings/convert_notation.py -f {output.tssT3} -c /home/mcgaugheyd/git/ChromosomeMappings/GRCh37_gencode2ensembl.txt > \
			{output.tssG}
		"""

# subtract the peaks in the tss regions
rule remove_tss_promoters:
	input:
		peaks = 'macs_peak/{sample}_peaks.blackListed.narrowPeak',
		tss = 'annotation/tss_GRCh37.bed'
	output:
		'macs_peak/{sample}_peaks.blackListed.tss_subtract.narrowPeak'
	shell:
		"""
		#module load {config[bedtools_version]}
		bedtools subtract -a {input.peaks} -b {input.tss} > {output} 
		"""

# create n bootstraps for each peak file set
# as background distribution for the motifs found
# with fimo
rule bootstrap_peaks:
	input:
		'macs_peak/{sample}_peaks.blackListed.tss_subtract.narrowPeak'
	output:
		'macs_peak/bootstrapping/{sample}_peaks.bootstrap_{bootstrap_num}.blackListed.tss_subtract.narrowPeak'
	shell:
		"""
		#module load {config[bedtools_version]}
		bedtools shuffle -excl annotation/tss_GRCh37.bed -g {config[bwa_genome_sizes]} -i {input} > {output}
		"""

# turn the peak data (minus tss) into fasta to ID motifs with fimo
rule peak_fasta:
	input:
		'macs_peak/{sample}_peaks.blackListed.tss_subtract.narrowPeak'
	output:
		'macs_peak/fasta/{sample}_peaks.blackListed.tss_subtract.narrowPeak.fasta'
	shell:
		"""
		#module load {config[bedtools_version]}
		bedtools getfasta -fi {config[bwa_genome]} -bed {input} -fo {output}
		"""

# identify TF which have differential expression
# with RNA-seq data
# will only use TF motifs which have a abs(log2FC) > 1 between iPSC and RPE (GFP / RFP)
# diff expression data derived from 01_simple_analysis.Rmd from ipsc_ rpe_RNA-seq repo
rule build_cisbp_master_file:
	input:
		config['diff_expression_csv']
	output:
		motifs = temp('tf.motifs'),
		diff_gene = temp('diff.gene'),
		join = temp('result.join'),
		meme = 'master_motifs.meme'
	params:
		diff_exp = config['diff_expression_value']
	run:
		shell("awk '{{print $7, $4}}' {config[cisbp_metadata_file]} /data/mcgaugheyd/motifs/TF_Information.txt | sort | uniq | grep -v '\.$'  > {output.motifs}")
		shell("awk -F, '{{if($3 > {params.diff_exp}  || $3 < -{params.diff_exp}) {{print $1}}}}'  {input}  > diff.gene")
		shell("join <( sort -k1,1 diff.gene ) <( sort -k1,1 tf.motifs) > {output.join}")
		join = open(output.join)
		meme_base_path = config['meme_cisbp_path'] #'/data/mcgaugheyd/motifs/meme_motifs'
		header = 'MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies (from uniform background):\nA 0.25000 C 0.25000 G 0.25000 T 0.25000\n\n'
		output = open(output.meme, 'w')
		output.write(header)
		for line in join:
			line = line.split()
			gene = line[0]
			motif = line[1]
			try:
				meme_file = open(meme_base_path + '/' + motif + '.meme')
				contents = meme_file.readlines()
				index = [ind for ind, line in enumerate(contents) if 'letter-prob' in line][0]
				if index != '':
					output.write('MOTIF ' + gene + ' ' + motif + '\n\n')
					output.write(''.join(contents[index:]))
			except:
				print(motif + ' not present')
		output.close()
			
		
		

# fimo calls motifs over a ~1e-5 threshold (based on background A-G-C-T rates)
rule call_motifs:
	input:
		fasta = 'macs_peak/fasta/{sample}_peaks.blackListed.tss_subtract.narrowPeak.fasta',
		meme_file = 'master_motifs.meme'
	output:
		'fimo_motifs/{sample}.fimo.dat.gz'
	shell:
		"""
		module load {config[meme_version]}
		fimo --no-qvalue --max-stored-scores 1000000 --text --parse-genomic-coord {input.meme_file} {input.fasta} | gzip -f > {output}
		"""

# turn the peak data (minus tss) into fasta to ID motifs with fimo
# runs on the bootstraps
rule peak_fasta_bootstrap:
	input:
		'macs_peak/bootstrapping/{sample}_peaks.bootstrap_{bootstrap_num}.blackListed.tss_subtract.narrowPeak'
	output:
		'macs_peak/bootstrapping/{sample}_peaks.bootstrap_{bootstrap_num}.blackListed.tss_subtract.narrowPeak.fasta'
	shell:
		"""
		#module load {config[bedtools_version]}
		bedtools getfasta -fi {config[bwa_genome]} -bed {input} -fo {output}
		"""

# fimo calls motifs over a ~1e-5 threshold (based on background A-G-C-T rates)
# runs on the bootstraps
rule call_motifs_bootstrap:
	input:
		fasta = 'macs_peak/bootstrapping/{sample}_peaks.bootstrap_{bootstrap_num}.blackListed.tss_subtract.narrowPeak.fasta',
		meme_file = 'master_motifs.meme'
	output:
		'fimo_motifs/bootstrapping/{sample}.bootstrap_{bootstrap_num}.fimo.dat.gz'
	shell:
		"""
		module load {config[meme_version]}
		fimo --no-qvalue --max-stored-scores 1000000 --text --parse-genomic-coord {input.meme_file} {input.fasta} | gzip -f > {output}
		"""

# get stats per bootstrap file
rule bootstrap_stats:
	input:
		'fimo_motifs/bootstrapping/{sample}.bootstrap_{bootstrap_num}.fimo.dat.gz'
	output:
		'fimo_motifs/bootstrapping_stats/{sample}.bootstrap_{bootstrap_num}.fimo_counts.dat.gz'
	shell:
		"""
		module load {config[R_version]}
		Rscript ~/git/ipsc_rpe_atac/src/fimo_counter.R {input} {output}
		"""


# find closest TSS/transcript to each motif
# finds 10 closest (-k 10)
# why so high? because I'm thinking that a bunch of similar transcripts 
# will artificially inflate the hits
# I want the two closest genes
# I'll collapse the tx into genes later in R
rule find_closest_TSS:
	input:
		'fimo_motifs/{sample}.fimo.dat.gz'
	output:
		'closest_TSS_motifs/{sample}.fimo.closestTSS.dat.gz'
	shell:
		"""
		module load {config[bedtools_version]}
		zcat {input} | \
			grep -v '^#' | \
			awk -v OFS='\t' '{{print $3, $4, $5, $2, $8, $6}}' | \
			sort -k1,1 -k2,2n | \
					closestBed -g <( sort -k1,1 -k2,2n {config[bwa_genome_sizes]} ) \
						-k 10 -d -a - -b <( sort -k1,1 -k2,2n annotation/tss_GRCh37.bed ) | \
			gzip -f > {output}
		"""

# bootstrap version
rule find_closest_TSS_bootstrap:
	input:
		'fimo_motifs/bootstrapping/{sample}.bootstrap_{bootstrap_num}.fimo.dat.gz'
	output:
		'closest_TSS_motifs/bootstrapping/{sample}.bootstrap_{bootstrap_num}.fimo.closestTSS.dat.gz'
	shell:
		"""
		module load {config[bedtools_version]}
		zcat {input} | \
			grep -v '^#' | \
			awk -v OFS='\t' '{{print $3, $4, $5, $2, $8, $6}}' | \
			sort -k1,1 -k2,2n | \
			closestBed -g <( sort -k1,1 -k2,2n {config[bwa_genome_sizes]} ) \
				-k 10 -d -a - -b <( sort -k1,1 -k2,2n annotation/tss_GRCh37.bed ) | \
			gzip -f > {output}
		"""

# parses find_closet_TSS with R
# pretty slow to do on local comp
# so let's make the cluster do it
rule process_closest_TSS_data:
	input:
		'closest_TSS_motifs/{sample}.fimo.closestTSS.dat.gz'
	output:
		'closest_TSS_motifs/processed/{sample}.fimo.closestTSS.processed.dat.gz'
	shell:
		"""
		module load {config[R_version]}
		Rscript ~/git/ipsc_rpe_atac/src/tss_motif_parser.R {input} {config[gtf_metadata]} {output} 
		"""

rule convert_peaks_to_hg19:
	input:
		'macs_peak/{sample}_peaks.blackListed.narrowPeak'
	output:
		'macs_peak/{sample}_peaks.blackListed.hg19.narrowPeak'
	shell:
		"""
		/home/mcgaugheyd/git/ChromosomeMappings/convert_notation.py \
			-c /home/mcgaugheyd/git/ChromosomeMappings/GRCh37_ensembl2UCSC.txt \
			-f <( grep -v hs37d5 {input} )  > \
			{output}
		"""

# create superset of all peak file
rule union_peaks:
	input:
		expand('macs_peak/{sample}_peaks.blackListed.narrowPeak', sample = list(SAMPLE_PATH.keys()))
	output:
		'macs_peak/union_peaks.blackListed.narrowPeak'
	shell:
		"""
		cat {input} | sort -k1,1 -k2,2n > {output}
		"""

# merge overlapping peaks together
rule merge_peaks:
	input:
		'macs_peak/union_peaks.blackListed.narrowPeak'
	output:
		'macs_peak/master_peaks.blackListed.narrowPeak'
	shell:
		"""
		module load {config[bedtools_version]}
		bedtools merge -i {input} > {output}
		"""

# compute read coverage across master peaks 
rule multiBamSummary:
	input:
		peaks = 'macs_peak/master_peaks.blackListed.narrowPeak',
		bam = expand('merged_bam/{sample}.bam', sample = list(SAMPLE_PATH.keys())),
		bai = expand('merged_bam/{sample}.bam.bai', sample = list(SAMPLE_PATH.keys()))
	output:
		default = 'deeptools/multiBamSummary.npz',
		tsv = 'deeptools/multiBamSummary.tsv'
	threads:
		16
	shell:
		"""
		module load {config[deeptools_version]}
		multiBamSummary BED-file \
			--bamfiles {input.bam} \
			--BED {input.peaks} \
			-o {output.default} \
			-p {threads} \
			--outRawCounts {output.tsv}
		"""

# set up data for ucsc viewing
# soft link bigWig to datashare
# create bigBed and soft link to datashare
rule ucsc_view:
	input:
		bigWig = 'bigWig/{sample}.bw', 
		bed = 'macs_peak/{sample}_peaks.blackListed.hg19.narrowPeak'
	params:
		base_path = '/data/mcgaugheyd/projects/nei/hufnagel/iPSC_RPE_ATAC_Seq/'
	output:
		bigWig = '/data/mcgaugheyd/datashare/hufnagel/hg19/{sample}.bw',
		bigBed = '/data/mcgaugheyd/datashare/hufnagel/hg19/{sample}_peaks.blackListed.hg19.narrowPeak.bb'
	shell:
		"""
		module load ucsc
		cut -f1,2,3,4 {input.bed} | sort -k1,1 -k2,2n > {input.bed}TEMP 
		bedToBigBed {input.bed}TEMP /data/mcgaugheyd/genomes/hg19/hg19.chrom.sizes {input.bed}.bb
		ln -s  {params.base_path}{input.bed}.bb {output.bigBed}
		ln -s {params.base_path}{input.bigWig} {output.bigWig}
		rm {input.bed}TEMP
		"""
	
	

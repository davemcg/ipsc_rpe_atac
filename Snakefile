import itertools
from itertools import chain

# build metadata
# run once outside of snakemake 
# bash ~/git/NGS_db/build_metadata.sh Huf > ~/git/ipsc_rpe_atac/metadata.csv

# read metadata into dictionary for snakemake to link samples to lane fastqs 
# one sample can have 1 to n lane fastq
# also making a TYPE <-> SAMPLE dict (where type is RFP / GFP / iPSC)
SAMPLE_PATH = dict() 
TYPE_SAMPLE = dict()
metadata = open(config['metadata_file'])
for line in metadata:
	path = line.split(',')[4][:-1]
	sample = line.split(',')[1]
	cell_type = sample.split('_')[0].upper() + '_' + line.split(',')[2]
	# skip header
	if sample == 'Sample':
		continue
	# skip SRA samples, do internal samples
#	if path != 'SRA':
	if sample not in SAMPLE_PATH:
		SAMPLE_PATH[sample] = [path]
	else:
		old_path = SAMPLE_PATH[sample]
		old_path.append(path)
		SAMPLE_PATH[sample] = old_path
	# SRA samples, build dict
#	elif path == 'SRA':
#		if sample not in SAMPLE_RUN:
#			SAMPLE_RUN[sample] = [path] 
#		else:
#			old_run = SAMPLE_PATH[sample]
#			old_run.append(path)
#			SAMPLE_PATH[sample] = old_run
#	else:
#		print(line)
#		print(path)
#		break

	if cell_type not in TYPE_SAMPLE:
		TYPE_SAMPLE[cell_type] = [sample] 
	else:
		old_sample = TYPE_SAMPLE[cell_type]
		old_sample = list(set(old_sample))
		old_sample.append(sample)
		TYPE_SAMPLE[cell_type] = old_sample
#
def fastq_pair_maker(wildcards):
	if 'SRR' in wildcards.lane_sample:
		out = 'fastq/' + wildcards.lane_sample
	else:
		out = ['fastq/' + wildcards.lane_sample + '_R1_001.fastq.gz', \
				'fastq/' + wildcards.lane_sample + '_R2_001.fastq.gz']
	return(out)

# given homer folder, go into homerResults and knownResults subfolders and extract all motif names
# and write file names for homer annotation for each motif
def yank_all_motifs(wildcards):
	import glob
	path = 'homer_unique_peaks_' + wildcards.comparison + '/' + wildcards.peak_type + '/'
	known_motifs = glob.glob(path + 'knownResults/*motif')
	novel_motifs = glob.glob(path + 'homerResults/*motif')
	all_motifs = known_motifs + novel_motifs
	motif_names = [i.split('/')[-1] for i in all_motifs]
	return [path + 'peak_by_motif/' + i + '.bed' for i in motif_names]
	#return [path + 'peak_by_motif/' + i + '.bed' for i in motif_names][0:5] # use for config_pretty.yaml for making dag

localrules: pull_lane_fastq_from_nisc_or_sra, retrieve_and_process_black_list, black_list, \
	total_reads, build_tss_regions, \
	download_HOCOMOCO_meme, reformat_motifs, \
	union_TFBS, union_TFBS_pretty_ucsc, prettify_union_TFBS, \
	common_peaks_across_all, find_closest_TSS_against_unique_peaks, \
	common_peaks_by_type, find_closest_TSS_to_all_common, \
	label_with_TF_RNA_seq_expression, \
	find_closest_TSS, union_summits_by_type, \
	summits_overlapping_common_peaks, \
	TF_to_target,  highlighted_homer_motif, cat_homer_annotate_peak, \
	intersect_homer_motifs_with_all_common_peaks, clean_all_common_peaks, \
	ucsc_view_homer_motifs, ucsc_view_bigWig, ucsc_view_common_peaks, \
	closest_TSS_each_cell_type, unique_peaks, union_HINT_by_type, \
	filter_down_all_common_peaks, intersect_homer_motifs_with_comparison_peaks

wildcard_constraints:
	sample = '|'.join(list(SAMPLE_PATH.keys())),
	cell_type = '|'.join(list(TYPE_SAMPLE.keys())),
	lane_fastq = '|'.join([x.split('/')[-1].split('.bam')[0] for x in list(itertools.chain.from_iterable(SAMPLE_PATH.values()))]),
	comparison = '|'.join(config['peak_comparison_pair'])

rule all:
	input:
		'network_reports/analysis.html',
		expand('punctate/punctate_analysis_{num}.txt', num = 1 ),
		expand('HINT_{comparison}/stats.txt', comparison = config['peak_comparison_pair']),
		'/data/mcgaugheyd/datashare/hufnagel/hg19/all_common_HINT_footprints.bb',
	#	expand('merged_bam_HQ/{cell_type}.q5.rmdup.bam', cell_type = ['GFP_ATAC-Seq', 'RFP_ATAC-Seq', 'IPSC_ATAC-Seq']),
	##	expand('macs_peak/{comparison}_common_peaks.blackListed.narrowPeak', comparison = config['peak_comparison_pair']),
	##	expand('homer/{cell_type}_{num}_homer_TFBS_run/knownResults_{num}.html', cell_type = ['GFP_ATAC-Seq'], num = [1,2,3,4,5,6,7,8]),
	#	expand('macs_peak/{cell_type}_common_peaks.blackListed.closestTSS.narrowPeak', cell_type = ['GFP_ATAC-Seq']),
	#	expand('computeMatrix/{cell_type}.matrix.gz', cell_type = ['IPSC_ChIP_h3k27ac','GFP_ATAC-Seq', 'RFP_ATAC-Seq', 'IPSC_ATAC-Seq']),
	#	expand('CGM/{cell_type}_common_peaks.blackListed.narrowPeak.CGM_score.tsv', cell_type = ['GFP_ATAC-Seq', 'RFP_ATAC-Seq', 'IPSC_ATAC-Seq']),
	#	expand('macs_peak/{cell_type}_common_peaks.h3k27ac_intersect.blackListed.narrowPeak', cell_type = ['GFP_ATAC-Seq', 'RFP_ATAC-Seq', 'IPSC_ATAC-Seq'])
	#	expand('fastq/{SRA_runs}_pass.fastq.gz', SRA_runs = [x for x in list(itertools.chain.from_iterable(SAMPLE_RUN.values()))]),
		#expand('/data/mcgaugheyd/datashare/hufnagel/hg19/{sample}.bw', sample = list(SAMPLE_PATH.keys())),
	##	expand('/data/mcgaugheyd/datashare/hufnagel/hg19/{sample}_peaks.blackListed.hg19.narrowPeak.bb', sample = list(SAMPLE_PATH.keys())),
		#'deeptools/multiBamSummary.npz',
		#'deeptools/multiBamSummary.tsv',
		#'metrics/reads_by_sample.txt',
	#	'fastqc/multiqc/multiqc_report.html',
	#	expand('downsample_bam/{sample}.q5.rmdup.ds.bam', sample = list(SAMPLE_PATH.keys())),
		#expand('msCentipede/closest_TSS/{cell_type}.{motif}.closestTSS.dat.gz', cell_type = list(TYPE_SAMPLE.keys()), motif = config['motif_IDs']),
		#expand('/data/mcgaugheyd/datashare/hufnagel/hg19/{motif}.union.HQ.pretty.bb', motif = config['motif_IDs']),
		#'/data/mcgaugheyd/datashare/hufnagel/hg19/all_common_peaks.blackListed.narrowPeak.bb',
		#'/data/mcgaugheyd/datashare/hufnagel/hg19/interesting_homer_motif.bb',
	##	expand('network_reports/{comparison}_{peak_type}/{comparison}_{peak_type}_networkAnalysis.html', comparison = config['peak_comparison_pair'], peak_type = ['all','enhancers','promoters']),
	#	expand('CGM/{cell_type}_common_peaks.blackListed.narrowPeak.CGM.tsv', cell_type = list(TYPE_SAMPLE.keys()))

rule pull_lane_fastq_from_nisc_or_sra:
	output:
		'fastq/{lane_files}'
	run:
		lane_files_full = [x for x in list(itertools.chain.from_iterable(SAMPLE_PATH.values()))]
		lane_files = [x.split('/')[-1] for x in list(itertools.chain.from_iterable(SAMPLE_PATH.values()))]
		for fastq in lane_files_full:
			# SRA files have "pass" in their fastq name
			if 'pass' not in fastq: 
				command = 'rsync -av trek.nhgri.nih.gov:' + fastq + ' fastq/'
				echo_command = 'echo ' + command
				shell('echo ' + str(output))
				shell(echo_command)
				shell(command)
			# SRA files
			else:
				srr = fastq.split('_')[0]
				command = 'module load sratoolkit; fastq-dump ' + srr + ' -O fastq \
							--gzip --skip-technical --readids --read-filter pass \
							--dumpbase --split-3 --clip'
				echo_command = 'echo ' + command
				shell('echo ' + str(output))
				shell('echo ' + command)
				shell(echo_command)
				shell(command)

localrules: build_hg19	
# define build custom hg19
# the TYR enhancer is a MOUSE enhancer
rule build_hg19:
	input:
		fa = config['bwa_genome']
	output:
		'/data/mcgaugheyd/genomes/hg19/hg19_bharti_TYR_enhancer.fa'
	shell:
		"""
		module load {config[samtools_version]}
		cat {input.fa} ~/git/ipsc_rpe_atac/data/mouse_tyr_enhancer.fa > {output}
		samtools faidx {output}
		"""

# build the custom bwa index
rule build_bwa_index:
	input:
		'/data/mcgaugheyd/genomes/hg19/hg19_bharti_TYR_enhancer.fa'
	output:
		'/data/mcgaugheyd/genomes/hg19/hg19_bharti_TYR_enhancer.fa.bwt'
	shell:
		"""
		module load {config[bwa_version]}
		bwa index {input}
		"""

# fastq to bwa for alignment
rule align:
	input:
		fa = fastq_pair_maker,
		bwa_genome = '/data/mcgaugheyd/genomes/hg19/hg19_bharti_TYR_enhancer.fa',
		bwa_index = '/data/mcgaugheyd/genomes/hg19/hg19_bharti_TYR_enhancer.fa.bwt'
	output:
		temp('realigned/{lane_sample}.bam')
	threads: 12 
	run:
		if 'SRR' not in input[0]:
			shell('module load {config[bwa_version]}; \
						bwa mem -t {threads} -B 4 -O 6 -E 1 -M {input.bwa_genome} {input.fa} | \
						samtools view -1 - >  {output}')
		else:
			shell('module load {config[bwa_version]}; \
						bwa mem -t {threads} -B 4 -O 6 -E 1 -M {input.bwa_genome} fastq/{wildcards.lane_sample} | \
						samtools view -1 - >  {output}')
		

# if sample has more than one lane bam, then merge into one bam
rule merge_bam:
	input:
		lambda wildcards: expand('realigned/{lane_sample}.bam', lane_sample = list(set([re.split(r'_R[1|2]_',x.split('/')[-1])[0] for x in SAMPLE_PATH[wildcards.sample]])))
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

# remove dups, only keep mapped (-F 4) reads with over q 5
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
		samtools view -O bam -F 4 -q 5 -b {input} | \
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

rule cell_type_bam:
	input:
		lambda wildcards: expand('merged_bam_HQ/{sample}.q5.rmdup.bam', sample = TYPE_SAMPLE[wildcards.cell_type])
	output:
		'merged_bam_HQ/{cell_type}.q5.rmdup.bam'
	shell:
		"""
		module load {config[samtools_version]}
		samtools merge -O bam {output} {input}
		samtools index {output}
		"""	

# downsample to 50e6 reads
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
		frac=$( samtools idxstats {input.bam} | cut -f3 | awk 'BEGIN {{total=0}} {{total += $1}} END {{frac=25000000/total;if (frac > 1) {{print 1}} else {{print frac}}}}' )
		samtools view -s $frac -b {input.bam} > {output.bam}
		samtools index {output.bam}
		"""

# bam to bigwig
# for UCSC visualization
rule bam_to_bigWig:
	input:
		'downsample_bam/{sample}.q5.rmdup.ds.bam'
	output:
		no_TYR = 'downsample_bam/{sample}.q5.rmdup.ds.noTYRchr.bam',
		bedgraph = temp('bigWig/{sample}.bG'),
		bw = 'bigWig/{sample}.bw'
	threads: 
		16
	shell:
		"""
		module load {config[samtools_version]}
		module load {config[deeptools_version]}
		module load ucsc
		# remove the mm10_TYR_enhancer_chr7_87538324_87543016
		samtools view -h {input} | grep -v "mm10_TYR_enhancer_chr7_87538324_87543016" | samtools view -b > {output.no_TYR}
		samtools index {output.no_TYR}
		bamCoverage --bam {output.no_TYR} -o {output.bedgraph} \
			--numberOfProcessors {threads} \
			--binSize 10 \
			--normalizeUsing RPGC \
			--effectiveGenomeSize 2864785220 \
			--ignoreForNormalization MT mm10_TYR_enhancer_chr7_87538324_87543016 \
			--outFileFormat bedgraph
		
		sort -k1,1 -k2,2n {output.bedgraph} > {output.bedgraph}TEMP 
		bedGraphToBigWig {output.bedgraph}TEMP /data/mcgaugheyd/genomes/hg19/hg19.chrom.sizes {output.bw}
		rm {output.bedgraph}TEMP
		"""

# basic stats on the bam file
# run on the bam file before q5, rmdup, and downsample happen
rule fastqc:
	input:
		'merged_bam/{sample}.bam'
	output:
		directory('fastqc/{sample}')
	threads: 8
	shell:
		"""
		module load fastqc
		mkdir -p {output}
		fastqc -t {threads} -o {output} {input}
		"""

# multiqc:
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

# Peak calling
# macs2
# ATAC-seq doesn't have control (no input sample). So run in single mode
# https://groups.google.com/forum/#!topic/macs-announcement/4OCE59gkpKY
rule peak_calling:
	input:
		'merged_bam_HQ/{sample}.q5.rmdup.bam'
	output:
		peaks = 'macs_peak/{sample}_peaks.xls',
		narrow_peaks = 'macs_peak/{sample}_peaks.narrowPeak',
		summits = 'macs_peak/{sample}_summits.bed'
	run:
		if 'SRS' in wildcards.sample:
			type = 'BAM'
		else:
			type = 'BAMPE'
		shell("module load {config[macs2_version]}; \
				macs2 callpeak -f " + type + " -g \"hs\" -t {input} -q 0.01 \
				--keep-dup all \
				-n {wildcards.sample} \
				--outdir macs_peak")

# Peak calling for GFP vs RFP and RFP vs iPSC
# macs2
rule peak_calling_comparator:
	input:
		t = lambda wildcards: expand('merged_bam_HQ/{sample}.q5.rmdup.bam', sample = TYPE_SAMPLE[wildcards.comparison.split('__not__')[0].upper()]),
		c = lambda wildcards: expand('merged_bam_HQ/{sample}.q5.rmdup.bam', sample = TYPE_SAMPLE[wildcards.comparison.split('__not__')[1].upper()])
	output:
		peaks = 'macs_peak/{comparison}_direct_peaks.xls',
		narrow_peaks = 'macs_peak/{comparison}_direct_peaks.narrowPeak',
		summits = 'macs_peak/{comparison}_direct_summits.bed'
	shell:
		"""
		module load {config[macs2_version]}
		macs2 callpeak -f BAMPE -g "hs" -t {input.t} -c {input.c} -q 0.01 \
			--keep-dup all \
			-n {wildcards.comparison}_direct \
			--outdir macs_peak
		"""

# blacklist
# remove peaks called in ENCODE black list regions
# https://sites.google.com/site/anshulkundaje/projects/blacklists
rule retrieve_and_process_black_list:
	output:
		'/data/mcgaugheyd/genomes/hg19/ENCFF001TDO.bed.gz'
	shell:
		"""
		wget https://www.encodeproject.org/files/ENCFF001TDO/@@download/ENCFF001TDO.bed.gz
		mv ENCFF001TDO.bed.gz /data/mcgaugheyd/genomes/hg19/
		"""

rule black_list:
	input:
		peaks = 'macs_peak/{sample}_peaks.narrowPeak',
		blacklist = '/data/mcgaugheyd/genomes/hg19/ENCFF001TDO.bed.gz'
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
		tss = temp('annotation/tss_hg19fix1.bed'),
		tssT2 = temp('annotation/tss_hg19fix2.bed'),
		tssG = 'annotation/tss_hg19.bed'
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
		awk '{{if($3 < 1) {{$3 = 1}} print}}' {output.tssT2} > {output.tssG}
		"""

		
rule union_peaks_by_type:
	input:
		lambda wildcards: expand('macs_peak/{sample}_peaks.blackListed.narrowPeak', sample = TYPE_SAMPLE[wildcards.cell_type])
	output:
		'macs_peak/{cell_type}_peaks.blackListed.narrowPeak'
	shell:
		"""
		cat {input} | sort -k1,1 -k2,2n > {output}
		"""

rule union_summits_by_type:
	input:
		lambda wildcards: expand('macs_peak/{sample}_summits.bed', sample = TYPE_SAMPLE[wildcards.cell_type])
	output:
		'macs_peak/{cell_type}_summits.bed'
	shell:
		"""
		cat {input} | sort -k1,1 -k2,2n > {output}
		"""

# merge overlapping peaks together
rule merge_peaks_by_type:
	input:
		'macs_peak/{cell_type}_peaks.blackListed.narrowPeak'
	output:
		'macs_peak/{cell_type}_peaks.blackListed.narrowPeak.merged'
	shell:
		"""
		module load {config[bedtools_version]}
		bedtools merge -d 150 -i {input} -c 5,7,10 -o mean > {output}
		"""

# only keep peaks which appear in two or individual peak calls and are 40% or more overlapping
# also add bed 9+ rgb coloring by cell type
rule common_peaks_by_type:
	input:
		union = 'macs_peak/{cell_type}_peaks.blackListed.narrowPeak',
		merged = 'macs_peak/{cell_type}_peaks.blackListed.narrowPeak.merged'
	output:
		'macs_peak/{cell_type}_common_peaks.blackListed.narrowPeak'
	run:
		shell("module load {config[bedtools_version]}; \
			bedtools intersect -a {input.merged} -b {input.union} -c -f 0.4 | awk '$7>1 {{print $0}}' > {output}T")
		tsv = open(output[0] + 'T')
		out = open(output[0], 'w')
		if wildcards.cell_type == 'RFP_ATAC-Seq':
			color = '255,0,0'
		elif wildcards.cell_type == 'GFP_ATAC-Seq':
			color = '0,255,0'
		elif wildcards.cell_type == 'IPSC_ATAC-Seq':
			color = '0,0,255'
		else:
			color = '255,255,255'
		for line in tsv:
			line = line.split()
			line[3] = str(min(999, round(float(line[3])))) # round for bigBed, no more than 999
			new_line = '\t'.join(line) + '\t1\t.\t' + line[1] + '\t' + line[2] + '\t' + color + '\n'	
			out.write(new_line)
		tsv.close()
		out.close()

# run HINT to find footprints within peaks
# run for each individual sample
# remove if peak in blacklist
rule HINT:
	input:
		blacklist = '/data/mcgaugheyd/genomes/hg19/ENCFF001TDO.bed.gz',
		peaks = 'macs_peak/{sample}_peaks.narrowPeak',
		bam = 'merged_bam_HQ/{sample}.q5.rmdup.bam'
	output:
		'HINT/{sample}.bed'
	shell:
		"""
		module unload python
		module load {config[rgt_version]}
		module load {config[samtools_version]}
		module load {config[bedtools_version]}
		grep -v "gl0\|hap\|TYR" {input.peaks} | intersectBed -a - -b {input.blacklist} -v > {input.peaks}GREP
		samtools view -bL chrs.bed {input.bam} > {input.bam}CHRS
		samtools index {input.bam}CHRS
		rgt-hint footprinting \
			--atac-seq \
			--paired-end \
			--organism=hg19 \
			--output-prefix {wildcards.sample} \
			--output-location=HINT \
			{input.bam}CHRS \
			{input.peaks}GREP
		sed -i 's/[[:blank:]]*$//' {output}
		rm {input.bam}CHRS
		rm {input.peaks}GREP
		"""

# keep footprints which are seen twice or more within a cell type
# also overlap peaks seen in two or more cell types
rule union_HINT_by_type:
	input:
		hint = lambda wildcards: expand('HINT/{sample}.bed', sample = TYPE_SAMPLE[wildcards.cell_type]),
		common_peaks = 'macs_peak/{cell_type}_common_peaks.blackListed.narrowPeak'
	output:
		union = 'HINT/{cell_type}.bed',
		merge = 'HINT/{cell_type}.merge.bed',
		intersect = 'HINT/{cell_type}.intersect.bed',
		processed = 'HINT/{cell_type}.intersect.colors.bed'
	run:
		shell("cat {input.hint} | sort -k1,1 -k2,2n > {output.union}")
		shell("module load {config[bedtools_version]}; \
			bedtools merge -i {output.union} -c 5 -o mean > {output.merge}")
		shell("module load {config[bedtools_version]}; \
			bedtools intersect -a {output.merge} -b {output.union} -c -f 0.4 | \
			awk '$5>1 {{print $0}}' | \
			bedtools intersect -a - -b {input.common_peaks} > {output.intersect}")
		tsv = open(output[2])
		out = open(output[3], 'w')
		if wildcards.cell_type == 'RFP_ATAC-Seq':
			color = '255,0,0'
		elif wildcards.cell_type == 'GFP_ATAC-Seq':
			color = '0,255,0'
		elif wildcards.cell_type == 'IPSC_ATAC-Seq':
			color = '0,0,255'
		else:
			color = '255,255,255'
		for line in tsv:
			line = line.split()
			line[3] = str(min(999, round(float(line[3])))) # round for bigBed, no more than 999
			new_line = '\t'.join(line) + '\t1\t.\t' + line[1] + '\t' + line[2] + '\t' + color + '\n'	
			out.write(new_line)
		tsv.close()
		out.close()

# scan for HOCOMOCO motifs in HINT footprint
rule HINT_motif_scan:
	input:
		'HINT/{cell_type}.intersect.colors.bed'
	output:
		'HINT/{cell_type}.intersect.colors_mpbs.bed'
	run:
		shell("mkdir -p HINT_motif_scan")
		shell("module load {config[rgt_version]}; \
			rgt-motifanalysis matching --motif-dbs $RGTDATA/motifs/hocomoco \
				--organism=hg19 \
				--input-files {input} \
				--output-location HINT")
		shell("sort -k1,1 -k2,2n {output} > {output}TEMP")
		if wildcards.cell_type == 'RFP_ATAC-Seq':
			color = '255,0,0'
		elif wildcards.cell_type == 'GFP_ATAC-Seq':
			color = '0,255,0'
		elif wildcards.cell_type == 'IPSC_ATAC-Seq':
			color = '0,0,255'
		else:
			color = '255,255,255'
		awk_call = "awk -v OFS='\t' '$5<0 {{$5=0}} {{print $1,$2,$3,$4,int($5),$6,$2,$3,\"" + color + "\"}}' {output}TEMP > {output}"
		shell(awk_call)
		shell("rm {output}TEMP")

localrules: make_chunks
# split into 13 lines each creates 115 files
rule make_chunks:
	input:
		'HINT/{cell_type}.intersect.colors_mpbs.bed'
	output:
		'HINT/chunks/{cell_type}_TF'
	shell:
		"""
		mkdir -p HINT/chunks
		cut -f4 {input} | sort | uniq  > {output} 
		split -l 13 --additional-suffix {wildcards.cell_type} {output}
		"""

localrules: HINT_motif_scan_chunker
# break HINT scanned motif bed into
# 100 chunks by motif
rule HINT_motif_scan_chunker:
	input:
		footprint = 'HINT/{cell_type}.intersect.colors_mpbs.bed',
		chunks = 'HINT/chunks/{cell_type}_TF'
	output:
		'HINT/chunks/{cell_type}.{chunk}.intersect.colors_mpbs.bed'
	shell:
		"""
		LC_ALL=C
		grep -f {wildcards.chunk}{cell_type} {input.footprint} > {output}
		"""
	
# differential footprint testing
rule HINT_differential:
	input:
		motif_A = lambda wildcards: ancient('HINT/chunks/' + wildcards.comparison.split('__not__')[0] + '.{chunk}.intersect.colors_mpbs.bed'),
		motif_B = lambda wildcards: ancient('HINT/chunks/' + wildcards.comparison.split('__not__')[1] + '.{chunk}.intersect.colors_mpbs.bed'),
		bam_A = lambda wildcards: 'merged_bam_HQ/' + wildcards.comparison.split('__not__')[0] + '.q5.rmdup.bam',
		bam_B = lambda wildcards: 'merged_bam_HQ/' + wildcards.comparison.split('__not__')[1]+ '.q5.rmdup.bam'
	output:
		directory = directory('HINT_{comparison}/{chunk}/'),
#		stats = 'HINT_{comparison}/{chunk}/' + wildcards.comparison.split('__not__')[0].split('_')[0] + '_' + \
#						wildcards.comparison.split('__not__')[1].split('_')[0] + '_factor.txt'
	threads: 1 
	run:
		shell("mkdir -p {output.directory}")
		condition1 = wildcards.comparison.split('__not__')[0].split('_')[0]
		condition2 = wildcards.comparison.split('__not__')[1].split('_')[0]
		shell("module load {config[rgt_version]}; module load tex; \
		rgt-hint differential --organism=hg19 --bc --nc {threads} \
			--mpbs-file1={input.motif_A} \
			--mpbs-file2={input.motif_B} \
			--reads-file1={input.bam_A} \
			--reads-file2={input.bam_B} \
			--output-location {output.directory} \
			--condition1=" + condition1 + " --condition2=" + condition2)
		shell("mkdir -p HINT_{wildcards.comparison}/Lineplots")
		shell("cp {output.directory}/Lineplots/* HINT_{wildcards.comparison}/Lineplots/")

localrules: merge_HINT_diff_results
rule merge_HINT_diff_results:
	input:
		expand('HINT_{{comparison}}/{chunk}/',
			chunk = ['xaa','xab','xac','xad','xae','xaf','xag','xah','xai','xaj','xak','xal','xam','xan','xao','xap','xaq','xar','xas','xat','xau','xav','xaw','xax','xay','xaz','xba','xbb','xbc','xbd','xbe','xbf','xbg','xbh','xbi','xbj','xbk','xbl','xbm','xbn','xbo','xbp','xbq','xbr','xbs','xbt','xbu','xbv','xbw','xbx','xby','xbz','xca','xcb','xcc','xcd','xce','xcf','xcg','xch','xci','xcj','xck','xcl','xcm','xcn','xco','xcp','xcq','xcr','xcs','xct','xcu','xcv','xcw','xcx','xcy','xcz','xda','xdb','xdc','xdd','xde','xdf','xdg','xdh','xdi','xdj','xdk','xdl','xdm','xdn','xdo','xdp','xdq','xdr','xds','xdt','xdu','xdv','xdw','xdx','xdy','xdz','xea','xeb','xec','xed','xee'])
	output:
		all = 'HINT_{comparison}/stats.txt',
		sig = 'HINT_{comparison}/sig_TF.txt'
	shell:
		"""
		cat HINT_{wildcards.comparison}/*/*_statistics.txt | grep "^Motif" | uniq > {output.all}
		grep -hv "^Motif" HINT_{wildcards.comparison}/*/*_statistics.txt >> {output.all}
		awk '$NF < 0.1 {{print $0}}' {output.all} | grep -v "^Motif" | cut -f1 > {output.sig} 
		"""

# create tss peaks
rule peaks_in_tss:
	input:
		peaks = 'macs_peak/{cell_type}_common_peaks.blackListed.narrowPeak',
		tss = 'annotation/tss_hg19.bed'
	output:
		'macs_peak/{cell_type}_common_peaks.blackListed.narrowPeak.tss_overlap.bed'
	shell:
		"""
		module load {config[bedtools_version]}
		bedtools intersect -a {input.peaks} -b {input.tss} -f 0.2 -wa -wb > {output}
		"""


# merge peaks
rule common_peaks_across_all:
	input:
		hint = expand('HINT/{cell_type}.intersect.colors_mpbs.bed', cell_type  = ['GFP_ATAC-Seq', 'RFP_ATAC-Seq','IPSC_ATAC-Seq']),
		macs = expand('macs_peak/{cell_type}_common_peaks.blackListed.narrowPeak', cell_type  = ['GFP_ATAC-Seq', 'RFP_ATAC-Seq','IPSC_ATAC-Seq'])
	output:
		hint = 'HINT/all_common_footprints.intersect.colors_mpbs.bed',
		macs = 'macs_peak/all_common_peaks.blackListed.narrowPeak.bed'
	shell:
		"""
		cat {input.macs} | sort -k1,1 -k2,2n -T /scratch/ | awk -v OFS='\t' '{{print $1,$2,$3,".",$4,".",$10,$11,$12}}' > {output.macs}
		cat {input.hint} | sort -k1,1 -k2,2n -T /scratch/ > {output.hint}
		"""

# filter down HINT to diff sig TFBS footprints
localrules: HINT_interesting_footprints
rule HINT_interesting_footprints:
	input:
		footprints = 'HINT/all_common_footprints.intersect.colors_mpbs.bed',
		sig = expand('HINT_{comparison}/sig_TF.txt', comparison = config['peak_comparison_pair']) 
	output:
		sig = 'HINT/interesting_tf.txt', 
		bed = 'HINT/all_common_footprints.intersect.colors_mpbs.sig_footprints.bed'
	shell:
		"""
		cat {input.sig} | sort | uniq  > {output.sig}
		LC_ALL=C
		grep -hFf {output.sig} {input} | grep "^chr"  > {output.bed}
		"""
		
	
# finds 10 closest (-k 10)
# also reports distance
rule find_closest_TSS_to_all_common:
	input:
		footprint = 'HINT/all_common_footprints.intersect.colors_mpbs.bed',
		bed = 'macs_peak/all_common_peaks.blackListed.narrowPeak.bed',
		tss = 'annotation/tss_hg19.bed'
	output:
		hint = 'HINT/all_common_footprints.intersect.colors_mpbs.closestTSS.bed.gz',
		macs = 'macs_peak/all_common_peaks.blackListed.narrowPeak.closestTSS.bed.gz'
	shell:
		"""
		module load samtools
		module load {config[bedtools_version]}
		cat {input.bed} | grep -v mm10 | \
			sort -k1,1 -k2,2n | \
			closestBed -g <( sort -k1,1 -k2,2n {config[bwa_genome_sizes]} ) \
				-k 4 -d -a - -b <( sort -k1,1 -k2,2n {input.tss} ) | \
			bgzip -f > {output.macs}

		cat {input.footprint} | \
			closestBed -g <( sort -k1,1 -k2,2n {config[bwa_genome_sizes]} ) \
				-k 4 -d -a - -b <( sort -k1,1 -k2,2n {input.tss} ) | \
			bgzip -f > {output.hint}
		"""

# filter out unique peaks
# with user given pair
# the first in the pair given in config['peak_comparison_pair'] are the peaks kept against the second of the pair
rule unique_peaks:
	input:
		one = lambda wildcards: 'macs_peak/' + wildcards.comparison.split('__not__')[0] + '_common_peaks.h3k27ac_intersect.blackListed.narrowPeak',
		two = lambda wildcards: 'macs_peak/' + wildcards.comparison.split('__not__')[1] + '_common_peaks.h3k27ac_intersect.blackListed.narrowPeak'
	output:
		'macs_peak/{comparison}_common_peaks.blackListed.narrowPeak'
	shell:
		"""
		echo wildcards.comparison
		module load {config[bedtools_version]}
		bedtools intersect -v -f 0.1 -a {input.one} -b {input.two} > {output}
		"""

# process 'macs_peak/all_common_peaks.blackListed.narrowPeak.closestTSS.bed'
# only keep two closest genes under 500,00bp away
# also clean up the column fields a bit
rule clean_all_common_peaks:
	input:
		hint_TF = 'HINT/interesting_tf.txt',
		hint = 'HINT/all_common_footprints.intersect.colors_mpbs.closestTSS.bed.gz',
		macs = 'macs_peak/all_common_peaks.blackListed.narrowPeak.closestTSS.bed.gz'
	output:
		hint = 'HINT/all_common_footprints.intersect.colors_mpbs.closestTSS.cleaned.bed.gz',
		macs = 'macs_peak/all_common_peaks.blackListed.narrowPeak.closestTSS.cleaned.bed.gz'
	shell:
		"""
		module load samtools
		module load {config[R_version]}
		zgrep -hFf {input.hint_TF} {input.hint} | bgzip  > {input.hint}GREP
		Rscript /home/mcgaugheyd/git/ipsc_rpe_atac/src/merge_peaks_homer_motifs.R {input.hint}GREP 2 {output.hint}
		rm {input.hint}GREP
		Rscript /home/mcgaugheyd/git/ipsc_rpe_atac/src/merge_peaks_homer_motifs.R {input.macs} 2 {output.macs}
		"""

# create network analysis R report
rule TF_gene_network_R:
	input:
		full = 'HINT/all_common_footprints.intersect.colors_mpbs.closestTSS.cleaned.bed.gz'
	output:
		file = 'network_reports/analysis.html'
	run:
		if wildcards.comparison == 'GFP_ATAC-Seq__not__IPSC_ATAC-Seq':
			shell("module load {config[R_version]}; \
				cp ~/git/ipsc_rpe_atac/src/network_analysis.Rmd {params.folder}/network_analysis.Rmd; \
				cp {params.gfp_vs_ipsc} {params.folder}/; \
				Rscript -e \"rmarkdown::render('{params.folder}/network_analysis.Rmd', output_file = '{params.file}', output_dir = '{params.folder}', params = list(comparison_1 = 'GFP', comparison_2 = 'iPSC', expression = '{params.gfp_vs_ipsc_file}', datafile = '../../{input}'))\" ")
		elif wildcards.comparison == 'GFP_ATAC-Seq__not__RFP_ATAC-Seq':
			shell("module load {config[R_version]}; \
				cp ~/git/ipsc_rpe_atac/src/network_analysis.Rmd {params.folder}/network_analysis.Rmd; \
				cp {params.gfp_vs_rfp} {params.folder}/; \
				Rscript -e \"rmarkdown::render('{params.folder}/network_analysis.Rmd', output_file = '{params.file}', output_dir = '{params.folder}', params = list(comparison_1 = 'GFP', comparison_2 = 'RFP', expression = '{params.gfp_vs_rfp_file}', datafile = '../../{input}'))\" ")
		elif wildcards.comparison == 'RFP_ATAC-Seq__not__GFP_ATAC-Seq':
			shell("module load {config[R_version]}; \
				cp ~/git/ipsc_rpe_atac/src/network_analysis.Rmd {params.folder}/network_analysis.Rmd; \
				cp {params.gfp_vs_rfp} {params.folder}/; \
				Rscript -e \"rmarkdown::render('{params.folder}/network_analysis.Rmd', output_file = '{params.file}', output_dir = '{params.folder}', params = list(comparison_1 = 'RFP', comparison_2 = 'GFP', expression = '{params.gfp_vs_rfp_file}', datafile = '../../{input}'))\" ")


# create CGM for GFP, RFP, iPSC
rule make_CGM:
	input:
		atac = 'macs_peak/{cell_type}_common_peaks.blackListed.narrowPeak',
		atac_h3k27ac = 'macs_peak/{cell_type}_common_peaks.h3k27ac_intersect.blackListed.narrowPeak',
	output:
		atac_score = 'CGM/{cell_type}_common_peaks.blackListed.narrowPeak.CGM_score.tsv',
		atac_h3k27ac_score = 'CGM/{cell_type}_common_peaks.h3k27ac_intersect.blackListed.narrowPeak.CGM_score.tsv',
		atac_count = 'CGM/{cell_type}_common_peaks.blackListed.narrowPeak.CGM_count.tsv',
		atac_h3k27ac_count = 'CGM/{cell_type}_common_peaks.h3k27ac_intersect.blackListed.narrowPeak.CGM_count.tsv'
	params:
		hint_value_column = 5,
		atac_value_column = 4
	shell:
		"""
		module load R
		module load bedtools
		bash ~/git/CGM/./matrix_maker.sh {config[gtf_file]} {input.atac} {output.atac_score} {params.atac_value_column}
		bash ~/git/CGM/./matrix_maker.sh {config[gtf_file]} {input.atac} {output.atac_count}
		bash ~/git/CGM/./matrix_maker.sh {config[gtf_file]} {input.atac_h3k27ac} {output.atac_h3k27ac_score} {params.atac_value_column}
		bash ~/git/CGM/./matrix_maker.sh {config[gtf_file]} {input.atac_h3k27ac} {output.atac_h3k27ac_count}
		"""

# create punctate chromatin differences report
# TF_gene_network.R does a broad (1mb) search for diff open regions around a gene
localrules: TF_gene_punctate_R
rule TF_gene_punctate_R:
	input:
		expand('CGM/{cell_type}.intersect.colors_mpbs.CGM_score.tsv', cell_type = ['GFP_ATAC-Seq', 'RFP_ATAC-Seq', 'IPSC_ATAC-Seq']),
		'/home/mcgaugheyd/git/ipsc_rpe_RNA-seq/data/lsTPM_by_Line.tsv', # expression
		expand('CGM/{cell_type}_common_peaks.blackListed.narrowPeak.CGM_score.tsv', cell_type = ['GFP_ATAC-Seq', 'RFP_ATAC-Seq', 'IPSC_ATAC-Seq'])# CGM
	output:
		'punctate/punctate_analysis_{num}.txt'
	shell:
		"""
		module load {config[pandoc_version]}
		module load {config[R_version]}
		Rscript -e "rmarkdown::render('/home/mcgaugheyd/git/ipsc_rpe_atac/src/punctate_CGM_analysis.Rmd')"
		"""

# compute read coverage across reference point
rule computeMatrix:
	input:
		bigWig = lambda wildcards: expand('bigWig/{sample}.bw', sample = TYPE_SAMPLE[wildcards.cell_type]), 
		region = config['gtf_file']
	output:
		'computeMatrix/{cell_type}.matrix.gz'
	threads: 16
	shell:
		"""
		module load {config[deeptools_version]}
		computeMatrix reference-point -p {threads} -S {input.bigWig} -R {input.region} -b 10000 -a 10000 -o {output}
		"""

# compute read coverage across common  peaks 
rule multiBamSummary:
	input:
		peaks = 'macs_peak/all_common_peaks.blackListed.narrowPeak',
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

# make bigBed for UCSC genome browse
localrules: HINT_bigBed
rule HINT_bigBed:
	input:
		all = 'HINT/all_common_footprints.intersect.colors_mpbs.bed',
		sig = 'HINT/all_common_footprints.intersect.colors_mpbs.sig_footprints.bed'
	output:
		all = '/data/mcgaugheyd/datashare/hufnagel/hg19/all_common_HINT_footprints.bb',
		sig = '/data/mcgaugheyd/datashare/hufnagel/hg19/sig_HINT_footprints.bb'
	params:
		base_path = '/data/mcgaugheyd/projects/nei/hufnagel/iPSC_RPE_ATAC_Seq/'
	shell:
		"""
		module load ucsc
		bedToBigBed {input.all} /data/mcgaugheyd/genomes/hg19/hg19.chrom.sizes {input.all}.bb
		ln -s {params.base_path}{input.all}.bb {output.all}
		
		bedToBigBed {input.sig} /data/mcgaugheyd/genomes/hg19/hg19.chrom.sizes {input.sig}.bb
		ln -s {params.base_path}{input.sig}.bb {output.sig}
		"""

# set up data for ucsc viewing
# soft link bigWig to datashare
# create bigBed and soft link to datashare
rule ucsc_view_bigWig:
	input:
		bigWig = 'bigWig/{sample}.bw', 
		bed = 'macs_peak/{sample}_peaks.blackListed.narrowPeak',
	params:
		base_path = '/data/mcgaugheyd/projects/nei/hufnagel/iPSC_RPE_ATAC_Seq/'
	output:
		bigWig = '/data/mcgaugheyd/datashare/hufnagel/hg19/{sample}.bw',
		bigBed = '/data/mcgaugheyd/datashare/hufnagel/hg19/{sample}_peaks.blackListed.hg19.narrowPeak.bb',
	shell:
		"""
		module load ucsc
		cut -f1,2,3,4 {input.bed} | grep -v mm10_TYR_enhancer_chr7_87538324_87543016 | sort -k1,1 -k2,2n > {input.bed}TEMP 
		bedToBigBed {input.bed}TEMP /data/mcgaugheyd/genomes/hg19/hg19.chrom.sizes {input.bed}.bb
		ln -s  {params.base_path}{input.bed}.bb {output.bigBed}
		ln -s {params.base_path}{input.bigWig} {output.bigWig}
		rm {input.bed}TEMP
		"""

rule ucsc_view_common_peaks:
	input:
		master = 'macs_peak/all_common_peaks.blackListed.narrowPeak.bed'
	params:
		 base_path = '/data/mcgaugheyd/projects/nei/hufnagel/iPSC_RPE_ATAC_Seq/'
	output:
		bigBedMasterPeaks = '/data/mcgaugheyd/datashare/hufnagel/hg19/all_common_peaks.blackListed.narrowPeak.bb'
	shell:
		"""
		module load ucsc
		bedToBigBed {input.master} /data/mcgaugheyd/genomes/hg19/hg19.chrom.sizes {input.master}.bb
		ln -s {params.base_path}{input.master}.bb {output.bigBedMasterPeaks}
		"""	

# get homer motif cat bed file into bigBed for online viewing
rule ucsc_view_homer_motifs:
	input:
		bed = 'peak_full/homer/interesting_homer_motif.bed.gz'
	output:
		'/data/mcgaugheyd/datashare/hufnagel/hg19/interesting_homer_motif.bb'
	params:
		base_path = '/data/mcgaugheyd/projects/nei/hufnagel/iPSC_RPE_ATAC_Seq/'
	shell:
		"""
		module load ucsc
		zcat {input} | cut -f1,2,3,4 > {input}TEMP
		bedToBigBed {input}TEMP /data/mcgaugheyd/genomes/hg19/hg19.chrom.sizes {input}.bb
		cp -f {params.base_path}{input}.bb {output}
		rm {input}TEMP
		"""


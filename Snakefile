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
		expand('merged_bam_HQ/{cell_type}.q5.rmdup.bam', cell_type = ['GFP_ATAC-Seq', 'RFP_ATAC-Seq', 'IPSC_ATAC-Seq']),
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
		expand('homer_unique_peaks_{comparison}/{peak_type}/homerResults.html', comparison = config['peak_comparison_pair'], peak_type = ['all','enhancers','promoters']), 
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

# further filter down ATAC common_peaks_by_type by overlapping with h3k27Ac ChIP-Seq
rule overlap_with_h3k27ac:
	input:
		atac = 'macs_peak/{cell_type}_common_peaks.blackListed.narrowPeak',
		ipsc_h3k27ac = 'macs_peak/IPSC_ChIP_h3k27ac_common_peaks.blackListed.narrowPeak',
		rpe1 = '/home/mcgaugheyd/git/ipsc_rpe_atac/data/ipsc_rpe_h3k27ac_frazer_1.bed.gz',
		rpe2 = '/home/mcgaugheyd/git/ipsc_rpe_atac/data/ipsc_rpe_h3k27ac_frazer_2.bed.gz',
		rpe3 = '/home/mcgaugheyd/git/ipsc_rpe_atac/data/ipsc_rpe_h3k27ac_frazer_3.bed.gz',
		rpe4 = '/home/mcgaugheyd/git/ipsc_rpe_atac/data/ipsc_rpe_h3k27ac_frazer_4.bed.gz',
		rpe5 = '/home/mcgaugheyd/git/ipsc_rpe_atac/data/ipsc_rpe_h3k27ac_frazer_5.bed.gz',
		rpe6 = '/home/mcgaugheyd/git/ipsc_rpe_atac/data/ipsc_rpe_h3k27ac_frazer_6.bed.gz',
		rpe7 = '/home/mcgaugheyd/git/ipsc_rpe_atac/data/ipsc_rpe_h3k27ac_frazer_7.bed.gz',
		rpe8 = '/home/mcgaugheyd/git/ipsc_rpe_atac/data/ipsc_rpe_h3k27ac_frazer_8.bed.gz',
		rpe9 = '/home/mcgaugheyd/git/ipsc_rpe_atac/data/ipsc_rpe_h3k27ac_frazer_9.bed.gz',
		ips1 = '/home/mcgaugheyd/git/ipsc_rpe_atac/data/ENCFF111WIP.bed.gz',
		ips2 = '/home/mcgaugheyd/git/ipsc_rpe_atac/data/ENCFF369AMU.bed.gz',
		ips3 = '/home/mcgaugheyd/git/ipsc_rpe_atac/data/ENCFF680AQK.bed.gz',
	output:
		'macs_peak/{cell_type}_common_peaks.h3k27ac_intersect.blackListed.narrowPeak'
	run:
		if wildcards.cell_type == 'IPSC_ATAC-Seq':
			#shell("module load {config[bedtools_version]}; \
			#		bedtools intersect -a {input.atac} -b <(zcat {input.ips1} {input.ips2} {input.ips3} | sort -k1,1 -k2,2n)  \
			#		> {output}")
			shell("module load {config[bedtools_version]}; \
					bedtools intersect -a {input.atac} -b {input.ipsc_h3k27ac} -c | \
					awk '$13 > 2 {{print $0}}' > {output}")
		elif wildcards.cell_type == 'GFP_ATAC-Seq' or wildcards.cell_type == 'RFP_ATAC-Seq':
			shell("module load {config[bedtools_version]}; \
					bedtools intersect -a {input.atac} -b <(zcat {input.rpe1} {input.rpe2} {input.rpe3} \
					{input.rpe4} {input.rpe5} {input.rpe6} {input.rpe7} {input.rpe8} {input.rpe9} | sort -k1,1 -k2,2n) -c | \
					awk '$13 > 2 {{print $0}}' > {output}")

# run HINT to find footprints within peaks
# run for each individual sample
rule HINT:
	input:
		peaks = 'macs_peak/{sample}_peaks.narrowPeak',
		bam = 'merged_bam_HQ/{sample}.q5.rmdup.bam'
	output:
		'HINT/{sample}.bed'
	shell:
		"""
		module unload python
		module load {config[rgt_version]}
		module load {config[samtools_version]}
		grep -v "gl0\|hap\|TYR" {input.peaks} > {input.peaks}GREP
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
rule union_HINT_by_type:
	input:
		lambda wildcards: expand('HINT/{sample}.bed', sample = TYPE_SAMPLE[wildcards.cell_type])
	output:
		union = 'HINT/{cell_type}.bed',
		merge = 'HINT/{cell_type}.merge.bed',
		intersect = 'HINT/{cell_type}.intersect.bed',
		processed = 'HINT/{cell_type}.intersect.colors.bed'
	run:
		shell("cat {input} | sort -k1,1 -k2,2n > {output.union}")
		shell("module load {config[bedtools_version]}; \
			bedtools merge -i {output.union} -c 5 -o mean > {output.merge}")
		shell("module load {config[bedtools_version]}; \
			bedtools intersect -a {output.merge} -b {output.union} -c -f 0.4 | awk '$5>1 {{print $0}}' > {output.intersect}")
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


# reformat union_TFBS to turn into proper bed and prep for UCSC viewing
rule prettify_peaks:
	input:
		'macs_peak/{cell_type}_common_peaks.blackListed.narrowPeak'
	output:
		'macs_peak/{cell_type}_common_peaks.blackListed.narrowPeak.bed'
	run:
		tsv = open(input[0])
		out = open(output[0], 'w')
		out.write("track name=\"" + wildcards.motif + "\" description=\"" + wildcards.motif + " msCentipede TFBS\" visibility=2 itemRgb=\"On\"\n")
		for line in tsv:
			line = line.split()
			new_line = line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[5] + '\t' + line[4] + '\t' + line[3] + '\t' + line[1] + '\t' + line[2]  
			if 'IPSC' in line[5]:
				out.write(new_line + '\t0,0,255\n')
			elif 'GFP' in line[5]:
				out.write(new_line + '\t0,255,0\n')
			else:
				out.write(new_line + '\t255,0,0\n')
		tsv.close()
		out.close()

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
		expand('macs_peak/{cell_type}_common_peaks.blackListed.narrowPeak', cell_type  = list(TYPE_SAMPLE.keys()))
	output:
		'macs_peak/all_common_peaks.blackListed.narrowPeak.bed'
	shell:
		"""
		cat {input} | sort -k1,1 -k2,2n | awk -v OFS='\t' '{{print $1,$2,$3,".",$4,".",$10,$11,$12}}' > {output}
		"""

# get summits in the 'macs_peak/{cell_type}_common_peaks.blackListed.narrowPeak
# for use in homer motif finding
rule summits_overlapping_common_peaks:
	input:
		peak_bed = 'macs_peak/{cell_type}_common_peaks.blackListed.narrowPeak',
		summit_bed = 'macs_peak/{cell_type}_summits.bed'
	output:
		'macs_peak/{cell_type}_common_peaks.blackListed.summits_intersect.bed'
	shell:
		"""
		module load {config[bedtools_version]}
		bedtools intersect -a {input.summit_bed} -b {input.peak_bed} > {output}
		"""
 
# find closest TSS/transcript to each peak
# finds 10 closest (-k 10)
# also reports distance
rule find_closest_TSS_to_all_common:
	input:
		bed = 'macs_peak/all_common_peaks.blackListed.narrowPeak.bed',
		tss = 'annotation/tss_hg19.bed'
	output:
		'macs_peak/all_common_peaks.blackListed.narrowPeak.closestTSS.bed'
	shell:
		"""
		module load {config[bedtools_version]}
		cat {input.bed} | grep -v mm10 | \
			sort -k1,1 -k2,2n | \
					closestBed -g <( sort -k1,1 -k2,2n {config[bwa_genome_sizes]} ) \
						-k 10 -d -a - -b <( sort -k1,1 -k2,2n {input.tss} ) | \
			gzip -f > {output}
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

# overlap summits in the unique peaks for homer
rule summits_in_unique_peaks:
	input:
		unique_peaks = 'macs_peak/{comparison}_common_peaks.blackListed.narrowPeak',
		summit_bed = lambda wildcards: 'macs_peak/' + wildcards.comparison.split('__not__')[0] + '_summits.bed'
	output:
		'macs_peak/{comparison}_common_peaks.blackListed.summits.bed'
	shell:
		"""
		module load {config[bedtools_version]}
		bedtools intersect -a {input.summit_bed} -b {input.unique_peaks} > {output}
		"""
	
# find closest TSS/transcript to each unique peak, as above for TFBS
# finds 10 closest (-k 10)
rule find_closest_TSS_against_unique_peaks:
	input:
		peak = 'macs_peak/{comparison}_common_peaks.blackListed.narrowPeak',
		tss = 'annotation/tss_hg19.bed'
	output:
		 'macs_peak/{comparison}_common_peaks.closestTSS.blackListed.narrowPeak'
	shell:
		"""
		module load {config[bedtools_version]}
		cat {input.peak} | \ ]
			sort -k1,1 -k2,2n | \
					closestBed -g <( sort -k1,1 -k2,2n {config[bwa_genome_sizes]} ) \
						-k 10 -d -a - -b <( sort -k1,1 -k2,2n {input.tss} ) | \
			gzip -f > {output}
		"""

localrules: prep_direct_summits_for_homer
# prep summits for homer
# remove in ENCODE black list regions
# h3k27ac overlap with smith frazer data
# split on tss / peak type
rule prep_direct_summits_for_homer:
	input:
		summits = 'macs_peak/{comparison}_direct_summits.bed',
		peaks = 'macs_peak/{comparison}_direct_peaks.narrowPeak',
		blackList = '/data/mcgaugheyd/genomes/hg19/ENCFF001TDO.bed.gz', 
		rpe1 = '/home/mcgaugheyd/git/ipsc_rpe_atac/data/ipsc_rpe_h3k27ac_frazer_1.bed.gz',
		rpe2 = '/home/mcgaugheyd/git/ipsc_rpe_atac/data/ipsc_rpe_h3k27ac_frazer_2.bed.gz',
		rpe3 = '/home/mcgaugheyd/git/ipsc_rpe_atac/data/ipsc_rpe_h3k27ac_frazer_3.bed.gz',
		rpe4 = '/home/mcgaugheyd/git/ipsc_rpe_atac/data/ipsc_rpe_h3k27ac_frazer_4.bed.gz',
		rpe5 = '/home/mcgaugheyd/git/ipsc_rpe_atac/data/ipsc_rpe_h3k27ac_frazer_5.bed.gz',
		rpe6 = '/home/mcgaugheyd/git/ipsc_rpe_atac/data/ipsc_rpe_h3k27ac_frazer_6.bed.gz',
		rpe7 = '/home/mcgaugheyd/git/ipsc_rpe_atac/data/ipsc_rpe_h3k27ac_frazer_7.bed.gz',
		rpe8 = '/home/mcgaugheyd/git/ipsc_rpe_atac/data/ipsc_rpe_h3k27ac_frazer_8.bed.gz',
		rpe9 = '/home/mcgaugheyd/git/ipsc_rpe_atac/data/ipsc_rpe_h3k27ac_frazer_9.bed.gz',
		ips1 = '/home/mcgaugheyd/git/ipsc_rpe_atac/data/ENCFF111WIP.bed.gz',
		ips2 = '/home/mcgaugheyd/git/ipsc_rpe_atac/data/ENCFF369AMU.bed.gz',
		ips3 = '/home/mcgaugheyd/git/ipsc_rpe_atac/data/ENCFF680AQK.bed.gz',
	output:
		peaks = 'macs_peak/{comparison}_direct.h3k27ac_intersect.blackListed.peaks.bed',
		summits = 'macs_peak/{comparison}_direct.h3k27ac_intersect.blackListed.summits.bed'
	shell:
		"""
		module load {config[bedtools_version]}
		# awk removes peaks with less than 1.8 fold enrichment over background
		# first intersect removes in blacklist
		# second intersect (big one) only keeps overlaps with 3 or more smith frazer h3k27ac 
		# third intersect keeps summits that overlap peaks that survive the first two intersects above
		awk '$7>1.8 {{print $0}}' {input.peaks} | \
			bedtools intersect -a - -b {input.blackList} -v | \
			bedtools intersect -a - -b <(zcat {input.rpe1} {input.rpe2} {input.rpe3} \
				{input.rpe4} {input.rpe5} {input.rpe6} {input.rpe7} {input.rpe8} {input.rpe9} | sort -k1,1 -k2,2n) -c | \
				 awk '$11>2 {{print $0}}' - > {output.peaks} 
		bedtools intersect -a {input.summits} -b {output.peaks} > {output.summits}
		"""	

# run homer on unique_peaks
rule homer_find_motifs_unique_peaks:
	input:
		hint = expand('HINT/{cell_type}.intersect.colors.bed', cell_type = ['GFP_ATAC-Seq', 'RFP_ATAC-Seq', 'IPSC_ATAC-Seq']),
		peak = 'macs_peak/{comparison}_common_peaks.blackListed.narrowPeak',
		tss = 'annotation/tss_hg19.bed'
	output:
		known = 'homer_unique_peaks_{comparison}/{peak_type}/knownResults.html',
		novel = 'homer_unique_peaks_{comparison}/{peak_type}/homerResults.html'
	params:
		out_dir = 'homer_unique_peaks_{comparison}/{peak_type}/'
	threads:
		8
	run:
		shell('mkdir -p {params.out_dir}')
		hint_peak = 'HINT/' + wildcards.comparison.split('__')[0] + '.intersect.colors.bed'
		if wildcards.peak_type == 'all':
			shell("module load {config[homer_version]}; module load {config[bedtools_version]}; \
					findMotifsGenome.pl <(intersectBed -a " + hint_peak + " -b {input.peak} -wa | awk '$4>50 {{print $0}}') \
						hg19 {params.out_dir} -p {threads} -size given -preparsedDir homer_preparsed/")
		elif wildcards.peak_type == 'enhancers':
			shell("module load {config[homer_version]};  module load {config[bedtools_version]}; \
					findMotifsGenome.pl <(intersectBed -a " + hint_peak + " -b {input.peak} | \
							awk '$4>50 {{print $0}}' | intersectBed -a - -b {input.tss} -v ) \
						hg19 {params.out_dir} -p {threads} -size given -preparsedDir homer_preparsed/")
		else: # promoter
			shell("module load {config[homer_version]};  module load {config[bedtools_version]}; \
					findMotifsGenome.pl <(intersectBed -a " + hint_peak + " -b {input.peak} | \
							awk '$4>50 {{print $0}}' | intersectBed -a - -b {input.tss} ) \
						hg19 {params.out_dir} -p {threads} -size given -preparsedDir homer_preparsed/")

localrules: homer_annotate_peaks
# label unique_peaks with motifs that homer finds
rule homer_annotate_peaks:
	input:
		homer_known = 'homer_unique_peaks_{comparison}/{peak_type}/knownResults.html',
		homer_novel = 'homer_unique_peaks_{comparison}/{peak_type}/homerResults.html',
		peak_file = 'macs_peak/all_common_peaks.blackListed.narrowPeak.bed'
	output:
		motif_bed = 'homer_unique_peaks_{comparison}/{peak_type}/peak_by_motif/{homer_motif}.bed',
		peak_bed = 'homer_unique_peaks_{comparison}/{peak_type}/peak_info/{homer_motif}.DELETE_LATER.xls'
	run:
		import glob
		motif_path = glob.glob(str(input.homer_known).split('/')[0] + '/' + wildcards.peak_type + '/*/' + wildcards.homer_motif)[0]
		shell_call = 'module load {config[homer_version]}; \
						annotatePeaks.pl {input.peak_file} hg19 -mbed {output.motif_bed} -m ' + \
						motif_path + \
						' > {output.peak_bed}'
		print(shell_call)
		shell(shell_call)

# label homer ID'ed TF with differential RNA-seq expression
rule label_with_TF_RNA_seq_expression:
	input:
		homer_known = 'homer_unique_peaks_{comparison}/{peak_type}/knownResults.html',
		homer_novel = 'homer_unique_peaks_{comparison}/{peak_type}/homerResults.html'
	output:
		known = 'homer_unique_peaks_{comparison}/{peak_type}/knownResults_matched_with_RNAseq.tsv',
		novel = 'homer_unique_peaks_{comparison}/{peak_type}/homerResults_matched_with_RNAseq.tsv'
	params:
		gfp_vs_ipsc = '/home/mcgaugheyd/git/ipsc_rpe_RNA-seq/data/GFP_vs_iPSC.results.csv',
		gfp_vs_rfp = '/home/mcgaugheyd/git/ipsc_rpe_RNA-seq/data/GFP_vs_RFP.results.csv',
		rfp_vs_ipsc = '/home/mcgaugheyd/git/ipsc_rpe_RNA-seq/data/RFP_vs_iPSC.results.csv'
	run:
		if wildcards.comparison == 'GFP_ATAC-Seq__not__IPSC_ATAC-Seq':
			shell("module load {config[R_version]}; \
			Rscript ~/git/ipsc_rpe_atac/src/match_TF_to_expression.R {input.homer_known} 'Name' {output.known} {params.gfp_vs_ipsc}; \
			Rscript ~/git/ipsc_rpe_atac/src/match_TF_to_expression.R {input.homer_novel} 'Best Match/Details' {output.novel} {params.gfp_vs_ipsc}")
		elif wildcards.comparison == 'RFP_ATAC-Seq__not__IPSC_ATAC-Seq':
			shell("module load {config[R_version]}; \
			Rscript ~/git/ipsc_rpe_atac/src/match_TF_to_expression.R {input.homer_known} 'Name' {output.known} {params.gfp_vs_rfp}; \
			Rscript ~/git/ipsc_rpe_atac/src/match_TF_to_expression.R {input.homer_novel} 'Best Match/Details' {output.novel} {params.gfp_vs_rfp}")
		else:	
			shell("module load {config[R_version]}; \
			Rscript ~/git/ipsc_rpe_atac/src/match_TF_to_expression.R {input.homer_known} 'Name' {output.known} {params.gfp_vs_rfp}; \
			Rscript ~/git/ipsc_rpe_atac/src/match_TF_to_expression.R {input.homer_novel} 'Best Match/Details' {output.novel} {params.gfp_vs_rfp}")

# filter label_with_TF_RNA_seq_expression to create
# list of targets
rule TF_to_target:
	input:
		known = 'homer_unique_peaks_{comparison}/{peak_type}/knownResults_matched_with_RNAseq.tsv',
		novel = 'homer_unique_peaks_{comparison}/{peak_type}/homerResults_matched_with_RNAseq.tsv'
	output:
		targets = 'homer_unique_peaks_{comparison}/{peak_type}/targets.txt'
	shell:
		"""
		module load {config[R_version]}
		Rscript ~/git/ipsc_rpe_atac/src/filter_match_results.R {input.known} {input.novel} {output}
		"""

# concatenate homer_annotate_peaks
# homer output is a bed {output.motif_bed} with the coordinates of the TFBS
rule cat_homer_annotate_peaks:
	input:
		motifs = yank_all_motifs,
		homer = 'homer_unique_peaks_{comparison}/{peak_type}/homerResults.html'
	output:
		'peak_full/homer_{comparison}/{peak_type}/all_homer_motif.bed.gz'
	shell:
		"""
		module load {config[samtools_version]}
		cat {input.motifs} | grep -v "track name" | uniq | sort -k1,1 -k2,2n | uniq | bgzip -cf > {output}
		"""

# process 'macs_peak/all_common_peaks.blackListed.narrowPeak.closestTSS.bed'
# only keep two closest genes under 500,00bp away
# also clean up the column fields a bit
rule clean_all_common_peaks:
	input:
		'macs_peak/all_common_peaks.blackListed.narrowPeak.closestTSS.bed'
	output:
		'macs_peak/all_common_peaks.blackListed.narrowPeak.closestTSS.cleaned.bed'
	shell:
		"""
		module load {config[R_version]}
		Rscript /home/mcgaugheyd/git/ipsc_rpe_atac/src/merge_peaks_homer_motifs.R {input} 2 {output}
		"""

# intersect homer motif bed 'peak_full/homer/interesting_homer_motif.bed.gz' 
# with all 'macs_peak/all_common_peaks.blackListed.narrowPeak.closestTSS.bed'
rule intersect_homer_motifs_with_comparison_peaks:
	input:
		peak_gene = 'macs_peak/all_common_peaks.blackListed.narrowPeak.closestTSS.cleaned.bed',
		a = 'macs_peak/{comparison}_common_peaks.blackListed.narrowPeak',
		b = 'peak_full/homer_{comparison}/{peak_type}/all_homer_motif.bed.gz'
	output:
		'peak_full/homer_{comparison}/{peak_type}/common_peaks.blackListed.narrowPeak__all_homer_motif.bed.gz'
	shell:
		"""
		module load {config[bedtools_version]}
		module load {config[samtools_version]}
		LC_ALL=C
		intersectBed -a {input.a} -b {input.peak_gene} -wb | sort -k1,1 -k2,2n | \
			intersectBed -a stdin -b {input.b} -wa -wb -sorted | bgzip -cf > {output}
		"""

# create network analysis R report
rule TF_gene_network_R:
	input:
		full = 'peak_full/homer_{comparison}/{peak_type}/common_peaks.blackListed.narrowPeak__all_homer_motif.bed.gz'
	output:
		file = 'network_reports/{comparison}_{peak_type}/{comparison}_{peak_type}_networkAnalysis.html'
	params:
		file = '{comparison}_{peak_type}_networkAnalysis.html',
		folder = 'network_reports/{comparison}_{peak_type}',
		gfp_vs_ipsc = '/home/mcgaugheyd/git/ipsc_rpe_RNA-seq/data/GFP_vs_iPSC.results.csv',
		gfp_vs_rfp = '/home/mcgaugheyd/git/ipsc_rpe_RNA-seq/data/GFP_vs_RFP.results.csv',
		gfp_vs_ipsc_file = 'GFP_vs_iPSC.results.csv',
		gfp_vs_rfp_file = 'GFP_vs_RFP.results.csv'
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


# closest 10 TSS (within 1e6) for each peak
rule closest_TSS_each_cell_type:
	input:
		bed = 'macs_peak/{cell_type}_common_peaks.blackListed.narrowPeak',
		tss = 'annotation/tss_hg19.bed'
	output:
		'macs_peak/{cell_type}_common_peaks.blackListed.closestTSS.narrowPeak'
	shell:
		"""
		module load {config[bedtools_version]}
		cat {input.bed} | grep -v mm10 | \
			sort -k1,1 -k2,2n | \
					closestBed -g <( sort -k1,1 -k2,2n {config[bwa_genome_sizes]} ) \
						-k 10 -d -a - -b <( sort -k1,1 -k2,2n {input.tss} ) | \
			gzip -f > {output}
		"""

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
		value_column = 4
	shell:
		"""
		module load R
		module load bedtools
		bash ~/git/CGM/./matrix_maker.sh {config[gtf_file]} {input.atac} {output.atac_score} {params.value_column}
		bash ~/git/CGM/./matrix_maker.sh {config[gtf_file]} {input.atac} {output.atac_count}
		bash ~/git/CGM/./matrix_maker.sh {config[gtf_file]} {input.atac_h3k27ac} {output.atac_h3k27ac_score} {params.value_column}
		bash ~/git/CGM/./matrix_maker.sh {config[gtf_file]} {input.atac_h3k27ac} {output.atac_h3k27ac_count}
		"""

# create punctate chromatin differences report
# TF_gene_network.R does a broad (1mb) search for diff open regions around a gene
localrules: TF_gene_punctate_R
rule TF_gene_punctate_R:
	input:
		#'peak_full/homer_{comparison}/{peak_type}/all_common_peaks.blackListed.narrowPeak.closestTSS__interesting_homer_motif.bed.gz', 
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

# take gene lists from TF_gene_punctate_R
# and do HOMER TFBS enrichment tests
#localrules: homer_TFBS
rule homer_TFBS:
	input:
		gene_list = 'punctate/punctate_analysis_{num}.txt',
		bed = 'macs_peak/{cell_type}_common_peaks.blackListed.closestTSS.narrowPeak',
		summit = 'macs_peak/{cell_type}_summits.bed'
	output:
		homer_in = '{cell_type}_{num}_common_peaks.blackListed.closestTSS.narrowPeak.summits',
		#homer_out = directory('homer/{cell_type}_{num}_homer_TFBS_run'),
		known = 'homer/{cell_type}_{num}_homer_TFBS_run/knownResults_{num}.html'
	params:
		out_dir = 'homer/{cell_type}_{num}_homer_TFBS_run',
		known = 'homer/{cell_type}_{num}_homer_TFBS_run/knownResults.html'
	threads:
		8
	shell:
		"""
		module load {config[bedtools_version]}
		module load {config[homer_version]}
		
		# peaks within 1000bp of target genes
		zcat {input.bed} | awk '$NF < 1000 {{print $0}}' | grep -f {input.gene_list} - | \
			cut -f1,2,3 | sort -k1,1 -k2,2n | uniq | \
			bedtools intersect -a - -b {input.summit} -wb > {output.homer_in}
		findMotifsGenome.pl {output.homer_in} hg19 {params.out_dir} -p {threads} -size 100 -preparsedDir homer_preparsed/
		cp {params.known} {output.known}
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


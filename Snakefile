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
	cell_type = sample.split('_')[0].upper()
	# skip header
	if sample == 'Sample':
		continue
	if sample not in SAMPLE_PATH:
		SAMPLE_PATH[sample] = [path]
	else:
		old_path = SAMPLE_PATH[sample]
		old_path.append(path)
		SAMPLE_PATH[sample] = old_path
	if cell_type not in TYPE_SAMPLE:
		TYPE_SAMPLE[cell_type] = [sample] 
	else:
		old_sample = TYPE_SAMPLE[cell_type]
		old_sample = list(set(old_sample))
		old_sample.append(sample)
		TYPE_SAMPLE[cell_type] = old_sample


localrules: pull_lane_fastq_from_nisc, retrieve_and_process_black_list, black_list, \
	ucsc_view, total_reads, union_peaks, merge_peaks, bootstrap_peaks, \
	peak_fasta, remove_tss_promoters, build_tss_regions, \
	build_cisbp_master_file, download_HOCOMOCO_meme, common_peaks, reformat_motifs, \
	merge_HOCOMOCO_cisbp, union_TFBS_pretty_ucsc, prettify_union_TFBS, \
	ucsc_view_master, common_peaks_across_all, find_closest_TSS_against_unique_peaks, \
	common_peaks_by_type
	#find_closest_TSS, find_closest_TSS_bootstrap

wildcard_constraints:
	sample = '|'.join(list(SAMPLE_PATH.keys())),
	motif = '|'.join(list(config['motif_IDs'])),
	cell_type = '|'.join(list(TYPE_SAMPLE.keys())),
	lane_fastq = '|'.join([x.split('/')[-1].split('.bam')[0] for x in list(itertools.chain.from_iterable(SAMPLE_PATH.values()))])

rule all:
	input:
		expand('/data/mcgaugheyd/datashare/hufnagel/hg19/{sample}.bw', sample = list(SAMPLE_PATH.keys())),
		expand('/data/mcgaugheyd/datashare/hufnagel/hg19/{sample}_peaks.blackListed.hg19.narrowPeak.bb', sample = list(SAMPLE_PATH.keys())),
		'deeptools/multiBamSummary.npz',
		'deeptools/multiBamSummary.tsv',
		'metrics/reads_by_sample.txt',
		'fastqc/multiqc/multiqc_report.html',
		#'macs_peak/all_common_peaks.blackListed.narrowPeak',
		expand('msCentipede/closest_TSS/{cell_type}.{motif}.closestTSS.dat.gz', cell_type = list(TYPE_SAMPLE.keys()), motif = config['motif_IDs']),
		expand('/data/mcgaugheyd/datashare/hufnagel/hg19/{motif}.union.HQ.pretty.bb', motif = config['motif_IDs']),
		'/data/mcgaugheyd/datashare/hufnagel/hg19/all_common_peaks.blackListed.narrowPeak.bb',
		 'macs_peak/' + config['peak_comparison_pair'][0] + '_vs_' + \
			config['peak_comparison_pair'][1] + '_common_peaks.closestTSS.blackListed.narrowPeak',
		'homer/'

rule pull_lane_fastq_from_nisc:
	output:
		'fastq/{lane_files}'
	run:
		lane_files = [x for x in list(itertools.chain.from_iterable(SAMPLE_PATH.values()))]
		for fastq in lane_files:
			command = 'rsync -av trek.nhgri.nih.gov:' + fastq + ' fastq/'
			echo_command = 'echo ' + command
			shell('echo ' + str(output))
			shell(echo_command)
			shell(command)
		

# fastq to bwa for alignment
rule align:
	input:
		expand('fastq/{{lane_sample}}{pair}.fastq.gz', pair = ['_R1_001', '_R2_001']) 
	output:
		temp('realigned/{lane_sample}.bam')
	threads: 12 
	shell:
		"""
		module load {config[bwa_version]}
		bwa mem -t {threads} -B 4 -O 6 -E 1 -M {config[bwa_genome]} {input} | \
			samtools view -1 - > \
			{output}
		"""

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
	shell:
		"""
		module load {config[macs2_version]}
		macs2 callpeak -f BAMPE -g "hs" -t {input} -q 0.01 \
			--keep-dup all \
			-n {wildcards.sample} \
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


rule download_HOCOMOCO_meme:
	output:
		'meme_pwm/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme'
	shell:
		"""
		wget http://hocomoco11.autosome.ru/final_bundle/hocomoco11/core/HUMAN/mono/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme -P meme_pwm
		"""

rule merge_HOCOMOCO_cisbp_jolma:
	input:
		'meme_pwm/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme'
	output:
		'meme_pwm/HOCOMOCOv11_core_HUMAN_mono_meme_format__cisBP__jolma.meme'
	run:
		shell("tail -n +10 /fdb/meme/motif_databases/EUKARYOTE/jolma2013.meme > jolma")
		shell("tail -n +10 /fdb/meme/motif_databases/CIS-BP/Homo_sapiens.meme > cisbp")
		file = open('cisbp')
		out_cis = open('out_cis', 'w')
		for line in file:
			if 'MOTIF' in line[0:10]:
				line = line.split()
				out_cis.write(line[0] + '\t' + line[1] + '__' + line[2] + '\n')
			else:
				out_cis.write(line)
		file.close()
		out_cis.close()
		shell("cat {input} out_cis jolma > {output}")
		shell("rm cisbp; rm out_cis; rm jolma")

# fimo calls motifs over a ~1e-5 threshold (based on background A-G-C-T rates)	
rule call_motifs:
	input:
		fasta = config['genome'],
		meme_file = 'meme_pwm/HOCOMOCOv11_core_HUMAN_mono_meme_format__cisBP__jolma.meme'
	output:
		'fimo_motifs/{motif}/meme.fimo.dat.gz'
	shell:
		"""
		module load {config[meme_version]}
		fimo --thresh 0.001 --no-qvalue  --text --parse-genomic-coord --motif {wildcards.motif} {input.meme_file} {input.fasta} | gzip -f > {output}
		"""	

rule reformat_motifs:
	input:
		'fimo_motifs/{motif}/meme.fimo.dat.gz'
	output:
		header = ('fimo_motifs/{motif}/meme.fimo.reformatted_for_msCentipede.HEADER'),
		body = ('fimo_motifs/{motif}/meme.fimo.reformatted_for_msCentipede.BODY'),
		intermediate = ('fimo_motifs/{motif}/meme.fimo.reformatted_for_msCentipede.INTERMEDIATE'),
		ready = 'fimo_motifs/{motif}/meme.fimo.reformatted_for_msCentipede.dat.gz'
	shell:
		"""
		module load samtools; # || true
		# only keep top 10,000 scoring motifs for msCentipede
		printf "Chr\tStart\tStop\tStrand\tPwmScore\n" > {output.header}
		zcat {input} | tail -n +2 | sort -k6,6nr | cut -f3,4,5,6,7 | head -n 10005 | sort -k1,1 -k2,2n > {output.body} || true
		cat {output.header} {output.body} > {output.intermediate} || true
		bgzip -fc {output.intermediate} > {output.ready} || true
		"""

# select parameters to optimize msCentipede ID of TFBS
# use replicates for each cell type
rule msCentipede_learn:
	input:
		motifs = 'fimo_motifs/{motif}/meme.fimo.reformatted_for_msCentipede.dat.gz',
		atac_bams = lambda wildcards: expand('merged_bam_HQ/{sample}.q5.rmdup.bam', sample = TYPE_SAMPLE[wildcards.cell_type]),
	output:
		model = 'msCentipede/{cell_type}/{motif}.model'
	shell:
		"""
		module load {config[msCentipede_version]}
		call_binding.py --task learn {input.motifs} {input.atac_bams} \
			--protocol ATAC_seq \
			--model_file {output.model} \
			--mintol 1e-2 \
			--log_file msCentipede/{wildcards.cell_type}/{wildcards.motif}.learn.logfile 
		"""

# ID TFBS
rule msCentipede_infer:
	input:
		motifs = 'fimo_motifs/{motif}/meme.fimo.reformatted_for_msCentipede.dat.gz',
		atac_bams = lambda wildcards: expand('merged_bam_HQ/{sample}.q5.rmdup.bam', sample = TYPE_SAMPLE[wildcards.cell_type]),
		model = 'msCentipede/{cell_type}/{motif}.model'
	output:
		posterior = 'msCentipede/{cell_type}/{motif}.posterior'
	shell:
		"""
		module load {config[msCentipede_version]}
		call_binding.py --task infer {input.motifs} {input.atac_bams} \
			--protocol ATAC_seq \
			--model_file {input.model} \
			--posterior_file {output.posterior} \
			--log_file msCentipede/{wildcards.cell_type}/{wildcards.motif}.infer.logfile 
		"""

# retain highest quality hits, and return as bed
# also only keep one that overlap macs peaks
# https://github.com/rajanil/msCentipede/issues/11
# One suggested set of conditions to select the most likely bound sites is
# LogPosOdds > 2 AND MultLikeRatio > 1 AND NegBinLikeRatio > 1 .
# LogPosOdds is the log of the posterior odds that the TF is bound. (LogPosOdds > 2 is equivalent to the binding posterior > 0.99)
# MultLikeRatio is the log likelihood ratio for the cleavage profile model.
# NegBinLikeRatio is the log likelihood ratio for the total chromatin accessibility model.
# (all logs are base 10)
rule msCentipede_motif_HQ:
	input:
		peaks = 'macs_peak/{cell_type}_common_peaks.blackListed.narrowPeak',
		tfbs = 'msCentipede/{cell_type}/{motif}.posterior'
	output:
		'msCentipede/{cell_type}/{motif}.posterior.HQ.tsv'
	shell:
		"""
		module load bedtools
		zcat {input.tfbs} | \
			awk -v s="{wildcards.cell_type}" -v m="{wildcards.motif}" '$5>2 && $7>1 && $8>1 {{print $0, "\t"s"_"m}}'| \
			cut -f1,2,3,4,5,9 | \
			tail -n +2 | \
			bedtools intersect -a - -b {input.peaks} > {output} 
		"""

# merge TFBS found by motifs
rule msCentipede_cat_by_motif:
	input:
		expand('msCentipede/{cell_type}/{{motif}}.posterior.HQ.tsv', cell_type = list(TYPE_SAMPLE.keys())) 
	output:
		'msCentipede/motif_cat/{motif}.posterior.HQ.tsv'
	shell:
		"""
		cat {input} | sort -k1,1 -k2,2n > {output}
		"""

# find closest TSS/transcript to each TFBS
# finds 10 closest (-k 10)
# why so high? because I'm thinking that a bunch of similar transcripts 
# will artificially inflate the hits
# I want the two closest genes
# I'll collapse the tx into genes later in R
rule find_closest_TSS:
	input:
		'msCentipede/{cell_type}/{motif}.posterior.HQ.tsv'	
	output:
		'msCentipede/closest_TSS/{cell_type}.{motif}.closestTSS.dat.gz'
	shell:
		"""
		module load {config[bedtools_version]}
		cat {input} | \
			sort -k1,1 -k2,2n | \
					closestBed -g <( sort -k1,1 -k2,2n {config[bwa_genome_sizes]} ) \
						-k 10 -d -a - -b <( sort -k1,1 -k2,2n annotation/tss_hg19.bed ) | \
			gzip -f > {output}
		"""

# create superset of all TFBS grouped by TFBS
rule union_TFBS:
	input:
		expand('msCentipede/{cell_type}/{{motif}}.posterior.HQ.tsv', cell_type  = list(TYPE_SAMPLE.keys()))
	output:
		'msCentipede_TFBS/{motif}.union.HQ.tsv'
	shell:
		"""
		cat {input} | sort -k1,1 -k2,2n > {output}
		"""

# reformat union_TFBS to turn into proper bed and prep for UCSC viewing
rule prettify_union_TFBS:
	input:
		'msCentipede_TFBS/{motif}.union.HQ.tsv'
	output:
		'msCentipede_TFBS/{motif}.union.HQ.pretty.bed'
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

# make bigBed for ucsc
# turn score into int and cap at 1000
rule union_TFBS_pretty_ucsc:
	input:
		'msCentipede_TFBS/{motif}.union.HQ.pretty.bed'
	output:
		nohead = temp('msCentipede_TFBS/{motif}.union.HQ.pretty.bedNOHEAD'),
		scorefix = temp('msCentipede_TFBS/{motif}.union.HQ.pretty.bedSCOREFIX'),
		bb = '/data/mcgaugheyd/datashare/hufnagel/hg19/{motif}.union.HQ.pretty.bb'
	params:
		base_path = '/data/mcgaugheyd/projects/nei/hufnagel/iPSC_RPE_ATAC_Seq/'
	run:
		shell("tail -n +2 {input} > {output.nohead}")
		data = open(output.nohead)
		out = open(output.scorefix, 'w')
		for line in data:
			line = line.split()
			score = int(float(line[4]))
			if score >= 1000:
				score = 999
			line[4] = str(score)
			out.write('\t'.join(line) + '\n')
		data.close()
		out.close()
		shell("module load ucsc; \
			bedToBigBed {output.scorefix} /data/mcgaugheyd/genomes/hg19/hg19.chrom.sizes {output.bb}")
		
rule union_peaks_by_type:
	input:
		lambda wildcards: expand('macs_peak/{sample}_peaks.blackListed.narrowPeak', sample = TYPE_SAMPLE[wildcards.cell_type])
	output:
		'macs_peak/{cell_type}_peaks.blackListed.narrowPeak'
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
		bedtools merge -i {input} -c 5,7,10 -o mean > {output}
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
			bedtools intersect -a {input.merged} -b {input.union} -c -f 0.4 | awk '$4>1 {{print $0}}' > {output}T")
		tsv = open(output[0] + 'T')
		out = open(output[0], 'w')
		if wildcards.cell_type == 'RFP':
			color = '255,0,0'
		elif wildcards.cell_type == 'GFP':
			color = '0,255,0'
		else:
			color = '0,0,255'
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

# filter out unique peaks
# with user given pair
# the first in the pair given in config['peak_comparison_pair'] are the peaks kept against the second of the pair
rule unique_peaks:
	input:
		one = 'macs_peak/' + config['peak_comparison_pair'][0] + '_common_peaks.blackListed.narrowPeak',
		two = 'macs_peak/' + config['peak_comparison_pair'][1] + '_common_peaks.blackListed.narrowPeak'
	output:
		'macs_peak/' + config['peak_comparison_pair'][0] + '_vs_' + \
			config['peak_comparison_pair'][1] + '_common_peaks.blackListed.narrowPeak'
	shell:
		"""
		module load {config[bedtools_version]}
		bedtools intersect -v -f 0.1 -a {input.one} -b {input.two} > {output}
		"""

# find closest TSS/transcript to each unique peak, as above for TFBS
# finds 10 closest (-k 10)
rule find_closest_TSS_against_unique_peaks:
	input:
		'macs_peak/' + config['peak_comparison_pair'][0] + '_vs_' + \
			config['peak_comparison_pair'][1] + '_common_peaks.blackListed.narrowPeak'
	output:
		 'macs_peak/' + config['peak_comparison_pair'][0] + '_vs_' + \
			config['peak_comparison_pair'][1] + '_common_peaks.closestTSS.blackListed.narrowPeak'
	shell:
		"""
		module load {config[bedtools_version]}
		cat {input} | \
			sort -k1,1 -k2,2n | \
					closestBed -g <( sort -k1,1 -k2,2n {config[bwa_genome_sizes]} ) \
						-k 10 -d -a - -b <( sort -k1,1 -k2,2n annotation/tss_hg19.bed ) | \
			gzip -f > {output}
		"""

# extract fasta
rule get_fasta_from_unique_peaks:
	input:
		'macs_peak/' + config['peak_comparison_pair'][0] + '_vs_' + \
			config['peak_comparison_pair'][1] + '_common_peaks.blackListed.narrowPeak'
	output:
		'macs_peak/' + config['peak_comparison_pair'][0] + '_vs_' + \
			config['peak_comparison_pair'][1] + '_common_peaks.blackListed.fasta'
	shell:
		"""
		module load {config[bedtools_version]}
		bedtools getfasta -fi {config[genome]} -bed {input} > {output}
		"""

# scramble fasta
rule homer_scramble_fasta:
	input:
		'macs_peak/' + config['peak_comparison_pair'][0] + '_vs_' + \
			config['peak_comparison_pair'][1] + '_common_peaks.blackListed.fasta'
	output:
		'macs_peak/' + config['peak_comparison_pair'][0] + '_vs_' + \
			config['peak_comparison_pair'][1] + '_common_peaks.blackListed.scramble.fasta'
	shell:
		"""
		module load {config[homer_version]}
		scrambleFasta.pl {input} > {output}
		"""

# run homer on unique_peaks
rule homer_find_motifs:
	input:
		exp = 'macs_peak/' + config['peak_comparison_pair'][0] + '_vs_' + \
			config['peak_comparison_pair'][1] + '_common_peaks.blackListed.fasta',
		control = 'macs_peak/' + config['peak_comparison_pair'][0] + '_vs_' + \
			config['peak_comparison_pair'][1] + '_common_peaks.blackListed.scramble.fasta'
	output:
		'homer/'
	shell:
		"""
		mkdir -p {output}
		module load {config[homer_version]}
		findMotifs.pl {input.exp} human {output} -fasta {input.control} -minlp -3 -humanGO
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
rule ucsc_view:
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
		cut -f1,2,3,4 {input.bed} | sort -k1,1 -k2,2n > {input.bed}TEMP 
		bedToBigBed {input.bed}TEMP /data/mcgaugheyd/genomes/hg19/hg19.chrom.sizes {input.bed}.bb
		ln -s  {params.base_path}{input.bed}.bb {output.bigBed}
		ln -s {params.base_path}{input.bigWig} {output.bigWig}
		rm {input.bed}TEMP
		"""

rule ucsc_view_master:
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
	

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

localrules: pull_lane_bams_from_nisc, retrieve_and_process_black_list, black_list, ucsc_view

wildcard_constraints:
	sample='|'.join(list(SAMPLE_PATH.keys())),
	lane_bam='|'.join([x.split('/')[-1].split('.bam')[0] for x in list(itertools.chain.from_iterable(SAMPLE_PATH.values()))])

rule all:
	input:
		expand('macs_peak/{sample}_peaks.blackListed.hg19.narrowPeak', sample = list(SAMPLE_PATH.keys())),
		'fastqc/multiqc/multiqc_report.html',
		'deeptools/multiBamSummary.tsv',
		expand('/data/mcgaugheyd/datashare/hufnagel/{sample}.bw', sample = list(SAMPLE_PATH.keys())),
		expand('/data/mcgaugheyd/datashare/hufnagel/{sample}_peaks.blackListed.hg19.narrowPeak.bb', sample = list(SAMPLE_PATH.keys()))


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
		java -Xmx20g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
			MergeSamFiles \
			TMP_DIR=/lscratch/$SLURM_JOB_ID \
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
		'merged_bam_HQ/{sample}.q5.rmdup.ds.bam'
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
			-f <( grep ^[0-9] {output.bedgraph} ) | sort -k1,1 -k2,2n > {output.bedgraph}TEMP
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
rule multiqc:
	input:
		expand('fastqc/{sample}', sample = list(SAMPLE_PATH.keys()))
	output:
		'fastqc/multiqc/multiqc_report.html'
	shell:
		"""
		module load multiqc
		multiqc fastqc/ -o fastqc/multiqc
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
		bigWig = '/data/mcgaugheyd/datashare/hufnagel/{sample}.bw',
		bigBed = '/data/mcgaugheyd/datashare/hufnagel/{sample}_peaks.blackListed.hg19.narrowPeak.bb'
	shell:
		"""
		module load ucsc
		cut -f1,2,3,4 {input.bed} | sort -k1,1 -k2,2n > {input.bed}TEMP 
		bedToBigBed {input.bed}TEMP /data/mcgaugheyd/genomes/hg19/hg19.chrom.sizes {input.bed}.bb
		ln -s  {params.base_path}{input.bed}.bb {output.bigBed}
		ln -s {params.base_path}{input.bigWig} {output.bigWig}
		rm {input.bed}TEMP
		"""
	
	

### snakemake --local-cores 1 --jobs 50 --use-envmodules --latency-wait 30 --restart-times 10 --cluster "sbatch -c {threads} --time={resources.maxtime} --mem={resources.mem} -e ./logdump/slurm-%j.err -o ./logdump/slurm-%j.out" 

configfile: "snake_config.yaml"
localrules: all, directories, update_alignments_table, extract_relevant_alignment_info, stage_in, multiqc_pre, quantify_genes_all_noexon_pairedend, fix_names_and_stage_out_featureCounts, setup_DE, stage_out_DE #, kallisto_rehead, kallisto_stage_out

import os, sys, re
import pandas as pd
import subprocess as sp
from itertools import compress


# System modules
sysenv22 = "build-env/f2022"
py = "python/3.9.6-gcccore-11.2.0"
scipy = "scipy-bundle/2021.10-foss-2021b"
rbase = "r/4.1.2-foss-2021b"
rlibs = "r-bundle-bioconductor/3.14-foss-2021b-r-4.1.2"
pdoc = "pandoc/2.18"
samt = "samtools/1.15-gcc-11.2.0"
kal = "kallisto/0.48.0-foss-2021b"

sysenv21 = "build-env/f2021"
featc = "subread/2.0.2-gcc-10.2.0"
fqc="fastqc/0.11.9-java-11"
mqc= "multiqc/1.11-foss-2020b-python-3.8.6"

# I'll have to type these many times, so simplify them.
# Snakemake recommends not using absolute paths, for portability, but it is causing problems for cluster jobs.
HOME = config["projectHome"]
SCRATCH = config["projectScratch"]
pairedend = config["pairedend"]
sys.path.append(config["codeDir"])

THREADS=8


## Parse sample IDs and names from description file
samples_raw = None
dates_raw = None
samples_temp = {}
dates_temp = {}
with open(config["sampleinfo"], 'r') as fin:
	next(fin) # skip heqader line
	# Get sample identifiers as tuples.
	samples_raw = [ (x[config["idcols"]["id"]].rstrip("\n"), \
	                 x[config["idcols"]["name"]].rstrip("\n")) \
	               for x in [line.split("\t") for line in fin] ]
with open(config["sampleinfo"], 'r') as fin:
	next(fin) # skip heqader line
	dates_raw = [ (x[config["idcols"]["id"]].rstrip("\n"), \
	                 x[config["idcols"]["aligndate"]].rstrip("\n")) \
	               for x in [line.split("\t") for line in fin] ]
# Convert tuples to dictionary.
def Convert(tup, di):
    di = dict(tup)
    return di
samples = Convert(samples_raw, samples_temp)
dates = Convert(dates_raw, dates_temp)


# Create file names based on assgned name and original sample ID.
uscore = config["idjoin"]
COMBONAME = [uscore.join([samples[x], x]) for x in samples.keys()]

FLAVOURS = ["exon_genecounts", "intron_genecounts", "spliced_genecounts","all_genecounts","all_noexon_genecounts"]


## Parse contrasts and covariates, to work out what output file names to expect.
contexts = re.sub(' |\t', '', config['DE']['contexts']).split(',')
context_names = [x.split('-')[0] for x in contexts]
context_levs = [int(x.split('-')[1]) for x in contexts]
contexts = [x for x in zip(context_names, context_levs)]
comparisons = re.sub(' |\t', '', config['DE']['comparisons']).split(',')
comparison_names = [x.split('-')[0] for x in comparisons]
comparison_vs = [x.split('-')[1] for x in comparisons]
comparison_treat = [int(x.split('v')[0]) for x in comparison_vs]
comparison_ref = [int(x.split('v')[1]) for x in comparison_vs]
comparisons = [x for x in zip(comparison_names, comparison_treat, comparison_ref)]

covars = pd.read_table(config['DE']['covars'], header=0, index_col=None, dtype=str)

def unique(x):     # R style set of unique values, preserving order of occurence.
	return(list(dict.fromkeys(x)))
contexts = [x + '_' + unique(covars.loc[:, x])[y - 1] for x,y in contexts]
comparisons = [x + '_' + unique(covars.loc[:, x])[y - 1] + '_vs_' + unique(covars.loc[:, x])[z - 1] for x,y,z in comparisons]

EXPECTED_DE = [x + '.' + y for x,y in zip(contexts, comparisons)]



## One chunk to rule them all
if config["do_DE"]:

	rule all:
		input:
			HOME + "/process/multiqc_pre/multiqc_report.html",
			HOME + "/process/duprates/duprates.txt",
			expand(HOME + "/process/featureCounts/{countsfile}.txt", countsfile=FLAVOURS),
			expand(HOME + "/results/PCA/{countsfile}.default.pca.html", countsfile=FLAVOURS),
			expand(HOME + "/results/DE/{countsfile}.all.custom.deseq2.spotfire.1.txt.gz", countsfile=FLAVOURS),
			expand(HOME + "/results/DE/{countsfile}.all.custom.deseq2.spotfire.2.txt.gz", countsfile=FLAVOURS),
			expand(HOME + "/results/DE/{countsfile}.all.custom.deseq2.spotfire.3.txt.gz", countsfile=FLAVOURS),
			expand(HOME + "/results/DE/{countsfile}.{dename}.deseq2.tsv", countsfile=FLAVOURS, dename=EXPECTED_DE)
		# params:
		# 	scratch=SCRATCH
		# shell:
		# 	"rm -r {params.scratch}"   # Clean up temporary/intermediate files

else:

	rule all:
		input:
			HOME + "/process/multiqc_pre/multiqc_report.html",
			HOME + "/process/duprates/duprates.txt",
			expand(HOME + "/process/featureCounts/{countsfile}.txt", countsfile=FLAVOURS),
			expand(HOME + "/results/PCA/{countsfile}.default.pca.html", countsfile=FLAVOURS)
		params:
			scratch=SCRATCH
		# shell:
		# 	"rm -r {params.scratch}"   # Clean up temporary/intermediate files
##




# Set up project home.
rule directories:
	output:
		"dirsdone"
	shell:
		"""
		mkdir -p aux data logdump process results 
		touch dirsdone
		"""
##

# Obtain data from VBCF

rule update_alignments_table:
	input:
		"dirsdone"
	output:
		HOME + "/aux/all_alignments_table.txt"
	params:
		url=config["vbcf_alignments_table_url"],
		un=config["un"],
		pw=config["pw"]
	shell:
		"""
		/usr/bin/wget -c --no-check-certificate --auth-no-challenge --user $(cat {params.un}) --password $(cat {params.pw}) -O {output} {params.url}
		rm dirsdone
		"""


rule extract_relevant_alignment_info:
	input:
		HOME + "/aux/all_alignments_table.txt"
	output:
		HOME + "/aux/sample_vbcf_metadata.txt"
	run:
		r = re.compile('|'.join(samples.keys()))
		q = re.compile('|'.join(dates.values()))
		with open(str(input), 'r') as fin:
			with open(str(output), 'w') as fout:
				for line in fin:
					l = line.split("\t")
					matched = r.match(l[11])  # 11 sample id, 47 date, 50 paired, 52 bam path, 54 bam url. 
					special_matched = r.match(l[0])  # Lazy data management means sample IDs end up in column 0 instead, with a lot of metadata left undefined.
					dated = q.match(l[47])
					if ((matched or special_matched) and dated) or l[0] == "flowcellId":
						fout.write(line)


if config["do_UMI"]:
    
	rule fetch_umi_alignments:
		input:
			HOME + "/" + config["sampleinfo"]
		output:
			HOME + "/data/{comboname}.aligned.bam"
		params:
			sep=uscore
		threads: 1
		resources:
			maxtime=10,
			mem='1G'
		shell:
			""" 
   			ncol=$(head -n 1 {input}  | awk '{{ print NF }}') &&
      		cp $(grep $(perl -e '$a="{output}";$a=~s/.*?{params.sep}(\d+)\.aligned\.bam/$1/;print $a') {input} | cut -f $ncol) {output}
			"""

else:
	rule fetch_alignments:
		input:
			HOME + "/aux/sample_vbcf_metadata.txt"
		output:
			HOME + "/data/{comboname}.aligned.bam"
		params:
			sep=uscore,
			un=config["un"],
			pw=config["pw"]
		threads: 1
		resources:
			maxtime=10,
			mem='1G'
		shell:
			"""
   			# added ".sorted.bam" to distinguish from umi_dedup.sorted.bam
			/usr/bin/wget -c --no-check-certificate --auth-no-challenge --user $(cat {params.un}) --password $(cat {params.pw}) -O {output} $(grep $(perl -e '$a="{output}";$a=~s/.*?{params.sep}(\d+)\.aligned\.bam/$1/;print $a').sorted.bam {input} | cut -f 55)
			"""
##




## Prepare

rule stage_in:
	input:
		HOME + "/data/{comboname}.aligned.bam"
	output:
		SCRATCH + "/data/{comboname}.aligned.bam"
	params:
		project=SCRATCH
	shell:
		"cp {input} {output}"


# Chr names, paired-end duplicated
# Duplicates already marked by facility
if pairedend:

	rule fix_bam:      
		input:
			SCRATCH + "/data/{comboname}.aligned.bam"
		output:
			SCRATCH + "/data_fixed/{comboname}.aligned.bam"
		envmodules:
			sysenv22,
			samt
		threads: 1
		resources:
			maxtime=30,
			mem='10G'
		shell:
			"samtools view -h -G 1024 {input} | sed 's/chrMT/chrM/g' | samtools view -Sb - > {output}"


	if config.get("kallisto_index", None) is not None:   # Providing a genome index or not controls whether to also run Kallisto or not.
		rule split_pair_bam:     # for Kallisto
			input:
				SCRATCH + "/data/{comboname}.aligned.bam"
			output:
				r1=SCRATCH + "/data_fastq/{comboname}.1.fastq",
				r2=SCRATCH + "/data_fastq/{comboname}.2.fastq"
			envmodules:
				sysenv22,
				samt
			threads: THREADS - 1
			resources:
				maxtime=5,
				mem='4G'
			shell:
				"samtools collate -@ $(({threads} - 1)) -O -u {input} | samtools fastq -s /dev/null -0 /dev/null -1 {output.r1} -2 {output.r2} -n"

else:

	rule fix_bam:      
		input:
			SCRATCH + "/data/{comboname}.aligned.bam"
		output:
			SCRATCH + "/data_fixed/{comboname}.aligned.bam"
		envmodules:
			sysenv22,
			samt
		threads: 1
		resources:
			maxtime=30,
			mem='1G'
		shell:
			"samtools view -h {input} | sed 's/chrMT/chrM/g' | samtools view -Sb - > {output}"


	if config.get("kallisto_index", None) is not None:   # Providing a genome index or not controls whether to also run Kallisto or not.
		rule split_pair_bam:     # for Kallisto
			input:
				SCRATCH + "/data/{comboname}.aligned.bam"
			output:
				SCRATCH + "/data_fastq/{comboname}.fastq"
			envmodules:
				sysenv22,
				samt
			threads: THREADS - 1
			resources:
				maxtime=5,
				mem='4G'
			shell:
				"samtools collate -@ $(({threads} - 1)) -O -u {input} | samtools fastq -s /dev/null -0 /dev/null -o {output} -n"
##




## Quality control


rule flagstats:
	input:
		SCRATCH + "/data_fixed/{comboname}.aligned.bam"
	output:
		SCRATCH + "/qc_pre/{comboname}.flagstats.txt"
	envmodules:
		sysenv22,
		samt
	threads: 1
	resources:
		maxtime=5,
		mem='1G'
	shell:
		"samtools flagstat -@ $(({threads} - 1)) {input} > {output}"


rule bamstats:
	input:
		SCRATCH + "/data_fixed/{comboname}.aligned.bam"
	output:
		SCRATCH + "/qc_pre/{comboname}.bamstats.txt"
	params:
		ref=config["reference"]
	envmodules:
		sysenv22,
		samt
	threads: 1
	resources:
		maxtime=5,
		mem='1G'
	shell:     # The facility's reference has contig names that are not in our reference fasta and it breaks the -r option.
		#"samtools stats -@ $(({threads} - 1)) -p -r {params.ref} {input} > {output}"
		"samtools stats -@ $(({threads} - 1)) -p {input} > {output}"


rule fastqc_mapped:
	input:
		SCRATCH + "/data_fixed/{comboname}.aligned.bam"
	output:
		SCRATCH + "/qc_pre/{comboname}.aligned_fastqc.html",
		SCRATCH + "/qc_pre/{comboname}.aligned_fastqc.zip"
	params:
		project=SCRATCH
	threads: 1 # multiple threads not applicable when each input has its own instance
	resources:
		maxtime=20,
		mem='1G'
	envmodules:
		sysenv21,
		fqc
	shell:
		"fastqc -f bam_mapped -t {threads} -o {params.project}/qc_pre {input}"


if pairedend:

	rule dupRadar:
		input:
			SCRATCH + "/data/{comboname}.aligned.bam"
		output:
			table=SCRATCH + "/duprate/{comboname}.duprate.txt",
			figure=SCRATCH + "/duprate/{comboname}.duprate.pdf"
		params:
			code=config["codeDir"],
			ann=config["annotation"]["exons"],
			stranding=config["stranding"],
			home=HOME
		threads: 4
		resources:
			maxtime=30,
			mem='4G'
		envmodules:
			sysenv22,
			rbase,
			rlibs
		shell:
			"""
			Rscript {params.code}/dupRadar_run.R -b {input} -g {params.ann} -o {output.table} -p -s {params.stranding} -t {threads}
			"""
			
else:	

	rule dupRadar:
		input:
			SCRATCH + "/data/{comboname}.aligned.bam"
		output:
			table=SCRATCH + "/duprate/{comboname}.duprate.txt",
			figure=SCRATCH + "/duprate/{comboname}.duprate.pdf"
		params:
			code=config["codeDir"],
			ann=config["annotation"]["exons"],
			stranding=config["stranding"],
			home=HOME
		threads: 4
		resources:
			maxtime=30,
			mem='4G'
		envmodules:
			sysenv22,
			rbase,
			rlibs
		shell:
			"""
			Rscript {params.code}/dupRadar_run.R -b {input} -g {params.ann} -o {output.table} -s {params.stranding} -t {threads}
			"""
#

rule multiqc_pre:
	input:						# don't os.path.dirname here, because it won't check for folder contents
		flagstats=expand(SCRATCH + "/qc_pre/{comboname}.flagstats.txt", comboname=COMBONAME),
		bamstats=expand(SCRATCH + "/qc_pre/{comboname}.bamstats.txt", comboname=COMBONAME),
		fastqc=expand(SCRATCH + "/qc_pre/{comboname}.aligned_fastqc.zip", comboname=COMBONAME)
	output:
		HOME + "/process/multiqc_pre/multiqc_report.html"
	params:
		project=HOME
	envmodules:
		sysenv21,
		mqc
	shell:
		"multiqc -f -o {params.project}/process/multiqc_pre $(dirname {input.fastqc[0]}) $(dirname {input.flagstats[0]}) $(dirname {input.bamstats[0]})"


rule merge_dup:
	input:
		pdf=expand(SCRATCH + "/duprate/{comboname}.duprate.pdf", comboname=COMBONAME),
		table=expand(SCRATCH + "/duprate/{comboname}.duprate.txt", comboname=COMBONAME)
	output:
		pdf=HOME + "/process/duprates/duprates.pdf",
		table=HOME + "/process/duprates/duprates.txt"
	params:
		code=config["codeDir"],
		scratch=SCRATCH,
		home=HOME
	envmodules:
		sysenv22,
		py,
		scipy,
		rbase,
		rlibs
	resources:
		maxtime=5,
		mem='1G'
	shell:
		"""
		mkdir -p {params.home}/process/duprates &&

		pdfmerge() {{ gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/default -dNOPAUSE -dQUIET -dBATCH -dDetectDuplicateImages -dCompressFonts=true -r150 -sOutputFile=$@ ; }}
		
		pdfmerge {output.pdf} {input.pdf}

		{params.code}/fileutilities.py T {params.scratch}/duprate/*duprate.txt -i -l --cols 10 12 | perl -pe '~s/_\|10//g;~s/\.duprate_\|12/.RPKb/g' > {output.table}
		"""
##




## Quantify


if pairedend:
	# Omitting -P -d -D because the VBCF reported fragment size range seems unreliable, nothing gets counted if I use it.
	# Overlapping features from transcript isoforms make it necessary to use -O --fraction . Summarised at the gene level, the fractional assignments of such a read should add up to 1, each read counted once per gene regardless of isoform ambiguity.

	rule quantify_genes_exons_pairedend:
		input:
			expand(SCRATCH + "/data_fixed/{comboname}.aligned.bam", comboname=COMBONAME)
		output:
			SCRATCH + "/featureCounts/exon_genecounts.txt"
		params:
			stranding=config["stranding"],
			ann=config["annotation"]["exons"]
		envmodules:
			sysenv21,
			featc
		group: "featureCounts"
		threads: THREADS
		resources:
			maxtime=30,
			mem='4G'
		shell:
			"""
			mkdir -p $(dirname {output})
			featureCounts -T {threads} -t exon -O --fraction -p --countReadPairs -B -C -a {params.ann} -o {output} -s {params.stranding} {input}
			"""


	rule quantify_genes_introns_pairedend:
		input:
			expand(SCRATCH + "/data_fixed/{comboname}.aligned.bam", comboname=COMBONAME)
		output:
			SCRATCH + "/featureCounts/intron_genecounts.txt"
		params:
			stranding=config["stranding"],
			ann=config["annotation"]["introns"]
		envmodules:
			sysenv21,
			featc
		group: "featureCounts"
		threads: THREADS
		resources:
			maxtime=30,
			mem='4G'
		shell:
			"""
			mkdir -p $(dirname {output})
			featureCounts -T {threads} -t intron -O --fraction -p --countReadPairs -B -C -a {params.ann} -o {output} -s {params.stranding} {input}
			"""


	rule quantify_genes_all_pairedend:
		input:
			expand(SCRATCH + "/data_fixed/{comboname}.aligned.bam", comboname=COMBONAME)
		output:
			SCRATCH + "/featureCounts/all_genecounts.txt"
		params:
			stranding=config["stranding"],
			ann=config["annotation"]["genes"],
			code=config["codeDir"]
		envmodules:
			sysenv21,
			featc
		group: "featureCounts"
		threads: THREADS
		resources:
			maxtime=30,
			mem='4G'
		shell:
			"""
			mkdir -p $(dirname {output})
			featureCounts -T {threads} -t gene -O --fraction -p --countReadPairs -B -C -a {params.ann} -o {output} -s {params.stranding} {input}
			"""

	rule quantify_genes_exon_split_pairedend:
			input:
				expand(SCRATCH + "/data_fixed/{comboname}.aligned.bam", comboname=COMBONAME)
			output:
				SCRATCH + "/featureCounts/spliced_genecounts.txt"
			params:
				stranding=config["stranding"],
				ann=config["annotation"]["exons"]
			envmodules:
				sysenv21,
				featc
			group: "featureCounts"
			threads: THREADS
			resources:
				maxtime=30,
				mem='4G'
			shell:
				"featureCounts -T {threads} -t exon --splitOnly -O --fraction -p --countReadPairs -B -C -a {params.ann} -o {output} -s {params.stranding} {input}"

	if config.get("kallisto_index", None) is not None:   # Providing a genome index or not controls whether to also run Kallisto or not.
		rule kallisto_pairedend:
			input:
				r1=SCRATCH + "/data_fastq/{comboname}.1.fastq",
				r2=SCRATCH + "/data_fastq/{comboname}.2.fastq"
			output:
				txt=SCRATCH + "/kallisto/{comboname}/abundance.tsv",
				boot=SCRATCH + "/kallisto/{comboname}/abundance.h5"
			params:
				idx=config["kallisto_index"]
			envmodules:
				sysenv22,
				kal
			group: "kallisto"
			threads: THREADS
			resources:
				maxtime=120,
				mem="48G"
			shell:
				"kallisto quant -t {threads} -i {params.idx} -o $(dirname {output.txt}) --bias -b 100 --rf-stranded {input.r1} {input.r2}"

else:

	rule quantify_genes_exons_singleend:
		input:
			expand(SCRATCH + "/data_fixed/{comboname}.aligned.bam", comboname=COMBONAME)
		output:
			SCRATCH + "/featureCounts/exon_genecounts.txt"
		params:
			stranding=config["stranding"],
			ann=config["annotation"]["exons"]
		envmodules:
			sysenv21,
			featc
		group: "featureCounts"
		threads: THREADS
		resources:
			maxtime=30,
			mem='4G'
		shell:
			"""
			mkdir -p $(dirname {output})
			featureCounts -T {threads} -t exon -O --fraction -a {params.ann} -o {output} -s {params.stranding} {input}
			"""


	rule quantify_genes_introns_singleend:
		input:
			expand(SCRATCH + "/data_fixed/{comboname}.aligned.bam", comboname=COMBONAME)
		output:
			SCRATCH + "/featureCounts/intron_genecounts.txt"
		params:
			stranding=config["stranding"],
			ann=config["annotation"]["introns"]
		envmodules:
			sysenv21,
			featc
		group: "featureCounts"
		threads: THREADS
		resources:
			maxtime=30,
			mem='4G'
		shell:
			"""
			mkdir -p $(dirname {output})
			featureCounts -T {threads} -t intron -O --fraction -a {params.ann} -o {output} -s {params.stranding} {input}
			"""


	rule quantify_genes_exons_split_singleend:
			input:
				expand(SCRATCH + "/data_fixed/{comboname}.aligned.bam", comboname=COMBONAME)
			output:
				SCRATCH + "/featureCounts/spliced_genecounts.txt"
			params:
				stranding=config["stranding"],
				ann=config["annotation"]["exons"]
			envmodules:
				sysenv21,
				featc
			group: "featureCounts"
			threads: THREADS
			resources:
				maxtime=30,
				mem='4G'
			shell:
				"featureCounts -T {threads} -t exon --splitOnly -O --fraction -a {params.ann} -o {output} -s {params.stranding} {input}"


	if config.get("kallisto_index", None) is not None:   # Providing a genome index or not controls whether to also run Kallisto or not.
		rule kallisto_singleend:
				input:
					SCRATCH + "/data_fastq/{comboname}.fastq"
				output:
					txt=SCRATCH + "/kallisto/{comboname}/abundance.tsv",
					boot=SCRATCH + "/kallisto/{comboname}/abundance.h5"
				params:
					idx=config["kallisto_index"],
					fragLen=config["fragLen"],
					fragLenSD=config["fragLenSD"]
				envmodules:
					sysenv21,
					kal
				group: "kallisto"
				threads: THREADS
				resources:
					maxtime=120,
					mem="48G"
				shell:
					"kallisto quant -t {threads} -i {params.idx} -o $(dirname {output.txt}) --bias -b 100 --rf-stranded --single -l {params.fragLen} -s {fragLenSD} {input}"
##


rule quantify_genes_all_noexon_pairedend: # run after featureCount
	input:
		exon_counts=SCRATCH + "/featureCounts/exon_genecounts.txt",
		all_counts=SCRATCH + "/featureCounts/all_genecounts.txt"
	output:
		all_noexon_counts=SCRATCH + "/featureCounts/all_noexon_genecounts.txt"
	params:
		ann=config["annotation"]["exons"],
		code=config["codeDir"]
	envmodules:
		sysenv22,
		rbase,
		rlibs
	shell:
		"""
		Rscript {params.code}/subtract_featureCount.v2.R -a {input.exon_counts} -b {input.all_counts} -c {params.ann} -m subtract -o {output.all_noexon_counts}
		"""


rule fix_names_and_stage_out_featureCounts: # featureCounts uses full input paths as sample names, which would suck for everything downstream.
	input:
		SCRATCH + "/featureCounts/{countsfile}.txt"
	output:
		tmp=SCRATCH + "/process/featureCounts_fixed/{countsfile}.txt",
		staged=HOME + "/process/featureCounts/{countsfile}.txt"
	params:
		prefix=(SCRATCH + "/data_fixed/").replace('/','\/')
	shell:  # Skip the top line comment and strip the paths and file extensions from the sample names.
		"""
		mkdir -p $(dirname {output.tmp}) &&
		head -n2 {input} | tail -n1 | perl -ne 's/.aligned.bam//g;s/{params.prefix}//g;print' > {output.tmp} &&
		tail -n +3 {input} >> {output.tmp} &&

		cp {output.tmp} {output.staged}  &&
		cp {input}.summary {output.staged}.summary 
		"""
##


if config.get("kallisto_index", None) is not None:   # Providing a genome index or not controls whether to also run Kallisto or not.

	rule kallisto_rehead:
		input:
			SCRATCH + "/kallisto/{comboname}/abundance.tsv"
		output:
			SCRATCH + "/kallisto/{comboname}/reheaded.tsv"
		shell:
			"""
				head -n 1 {input} | perl -e 'foreach $a (<STDIN>){{$a=~s/tpm/$ARGV[0]/;$a=~s/est_counts/$ARGV[0]/;print $a}}' $(basename $(dirname {input})) > {output} &&
				tail -n +2 {input} >> {output}
			"""


	rule kallisto_collate_separate:
		input:
			expand(SCRATCH + "/kallisto/{comboname}/reheaded.tsv", comboname=COMBONAME)
		output:
			counts=SCRATCH + "/process/kallisto/est_counts.txt",
			tpms=SCRATCH + "/process/kallisto/tpm.txt",
			counts_intron=SCRATCH + "/process/kallisto/est_counts.only_premRNA.txt",
			tpms_intron=SCRATCH + "/process/kallisto/tpm.only_premRNA.txt",
			counts_exon=SCRATCH + "/process/kallisto/est_counts.only_tx.txt",
			tpms_exon=SCRATCH + "/process/kallisto/tpm.only_tx.txt"
		envmodules:
			sysenv22,
			py,
			scipy
		resources:
			maxtime=5,
			mem="4G"
		shell:
			"""
			code/fileutilities.py T {input} -r -i --cols 3 > {output.counts} &&
			code/fileutilities.py T {input} -r -i --cols 4 > {output.tpms} &&

			head -n 1 {output.counts} > {output.counts_intron} &&
			grep _pre {output.counts} >> {output.counts_intron} &&
			head -n 1 {output.tpms} > {output.tpms_intron} &&
			grep _pre {output.tpms} >> {output.tpms_intron} &&
			
			head -n 1 {output.counts} > {output.counts_exon} &&
			grep _no_pre {output.counts} | perl -ne 's/_no_pre[+-]//;print' >> {output.counts_exon} &&
			grep -E -v _pre {output.counts} >> {output.counts_exon} &&
			head -n 1 {output.tpms} > {output.tpms_exon} &&
			grep _no_pre {output.tpms} | perl -ne 's/_no_pre[+-]//;print' >> {output.tpms_exon} &&
			grep -E -v _pre {output.tpms} >> {output.tpms_exon}
			"""


	rule kallisto_stage_out:
		input:
			SCRATCH + "/process/kallisto/est_counts.txt",
			SCRATCH + "/process/kallisto/tpm.txt",
			SCRATCH + "/process/kallisto/est_counts.only_premRNA.txt",
			SCRATCH + "/process/kallisto/tpm.only_premRNA.txt",
			SCRATCH + "/process/kallisto/est_counts.only_tx.txt",
			SCRATCH + "/process/kallisto/tpm.only_tx.txt"
		output:
			HOME + "/process/kallisto/est_counts.txt",
			HOME + "/process/kallisto/tpm.txt",
			HOME + "/process/kallisto/est_counts.only_premRNA.txt",
			HOME + "/process/kallisto/tpm.only_premRNA.txt",
			HOME + "/process/kallisto/est_counts.only_tx.txt",
			HOME + "/process/kallisto/tpm.only_tx.txt"
		params:
			destdir=HOME + "/process/kallisto/"
		shell:
			"cp {input} {params.destdir}"
##




## PCA
rule PCA:
	input:
		HOME + "/process/featureCounts/{countsfile}.txt"
	output:
		HOME + "/results/PCA/{countsfile}.default.pca.html",
		HOME + "/process/featureCounts/{countsfile}.scaled.txt"
	envmodules:
		sysenv22,
		rbase,
		rlibs,
		pdoc
	params:
		home=HOME,
		minMean=config["PCA"]["minMean"],
		minSingle=config["PCA"]["minSingle"],
		ntop=config["PCA"]["ntop"],
		covars=config["PCA"]["covars"],
		code=config["codeDir"],
		pattern=config["PCA"]["specialnorm"],
		forvar=config["PCA"]["loopgroup"],
		plotPDF=config["PCA"]["plotPDF"]
	group: "pca"
	threads: 1
	resources:
		maxtime=10,
		mem='4G'
	shell:
		"""
		Rscript {params.code}/pca_run.R -F {params.forvar} -b {params.home} -f process/featureCounts/$(basename {input}) -G {params.pattern} -I 6 -i 1 -M {params.minMean} -m {params.minSingle} -n {params.ntop} -o results/PCA -r process/PCA -s {params.covars} -T {params.code}/pca_report_template.Rmd -W 6 -P {params.plotPDF}
		"""
##




## DE

if config["do_DE"]:

	rule setup_DE:
		input:
			covars=os.path.join(HOME, config["DE"]["covars"])
		output:
			covars=os.path.join(SCRATCH, config["DE"]["covars"])
		params:
			home=HOME,
			scratch=SCRATCH
		shell:
			"""
			mkdir -p $(dirname {output.covars}) &&
			cp {input.covars} {output.covars} &&

			mkdir -p {params.home}/results/DE
			"""


	rule DE:
		input:
			counts=SCRATCH + "/process/featureCounts_fixed/{countsfile}.txt",
			samples=os.path.join(SCRATCH, config["DE"]["covars"])       # listed only to trigger the setup rule
		output:
			vanilla=expand(SCRATCH + "/results/DE/{countsfile}/{countsfile}.{dename}.deseq2.tsv", dename=EXPECTED_DE, allow_missing=True),  # allow_missing: {dename} takes all the specified values in each run of this rule, while the rule is ran for every value of {countsfile} specified later.
			custom=expand(SCRATCH + "/results/DE/{countsfile}/{countsfile}.{dename}.custom.deseq2.tsv", dename=EXPECTED_DE, allow_missing=True)
			# One table per contrast, but contrasts may get grouped into fewer HTML reports
		wildcard_constraints:
			countsfile="[^.]+",
			dename=".+_vs_.+"
		envmodules:
			sysenv22,
			rbase,
			rlibs,
			pdoc
		params:
			samples=config["DE"]["covars"],     # not taking it from {input} because I don't want the full path value, nor the basename only
			scratch=SCRATCH,
			form=config["DE"]["formula"],
			reform=config['DE']['reduced'],
			minCnt=config["DE"]["minCnt"],
			minTPM=config["DE"]["minTPM"],
			ntop=config["DE"]["ntop"],
			comps=re.sub(' |\t', '', config['DE']["comparisons"]),
			bgr=re.sub(' |\t', '', config['DE']['contexts']),
			lfc=config["DE"]["lfc_thresh"],
			pcut=config["DE"]["p_cut"],
			code=config["codeDir"],
			pattern=config["DE"]["specialnorm"]
		group: "deseq"
		threads: 1
		resources:
			maxtime=60,
			mem='5G'
		shell:
			"""
			Rscript {params.code}/deseq2_run_LR.R -b {params.scratch} -c {params.lfc} -f process/featureCounts_fixed/$(basename {input.counts}) -F {params.form} -G {params.pattern} -i 1 -I 6 -l -m {params.minCnt} -M {params.minTPM} -n {params.ntop} -o results/DE -p {params.pcut} -P $(basename {input.counts} .txt) -r process/DE -R {params.reform} -s {params.samples} -T {params.code}/deseq2_report_template_LR.Rmd -W 6 -x {params.comps} -X {params.bgr}
			"""


	rule merge_DE_table:
			input:
				expand(SCRATCH + "/results/DE/{countsfile}/{countsfile}.{dename}.custom.deseq2.tsv", dename=EXPECTED_DE, allow_missing=True)   # allow_missing: {dename} takes all the specified values in each run of this rule, while the rule is ran for every value of {countsfile} specified later.
			output:
				SCRATCH + "/results/DE/{countsfile}/{countsfile}.all.custom.deseq2.tsv"
			wildcard_constraints:
				countsfile="[^.]+"
			envmodules:
				sysenv22,
				py,
				scipy
			params:
				code=config["codeDir"],
				scratch=SCRATCH,
				home=HOME
			group: "stage_out"
			threads: 1
			resources:
				maxtime=10,
				mem='4G'
			shell:
				"""
				{params.code}/fileutilities.py T {input} -i -r --appnd > {output}
				"""


	rule prepare_for_Spotfire:
		input:
			de=SCRATCH + "/results/DE/{countsfile}/{countsfile}.all.custom.deseq2.tsv",
			cnt=SCRATCH + "/featureCounts/{countsfile}.txt"
		output:
			SCRATCH + "/results/DE/{countsfile}/{countsfile}.all.custom.deseq2.spotfire.1.txt.gz",
			SCRATCH + "/results/DE/{countsfile}/{countsfile}.all.custom.deseq2.spotfire.2.txt.gz",
			SCRATCH + "/results/DE/{countsfile}/{countsfile}.all.custom.deseq2.spotfire.3.txt.gz"
			# the third one is a dummy and need to be ignored in further processing
		wildcard_constraints:
				countsfile="[^.]+",
				denamenl="all\.deseq2|.+\.deseq2\.nolab"
		envmodules:
			sysenv22,
			rbase,
			rlibs
		params:
			code=config["codeDir"],
			scratch=SCRATCH,
			home=HOME,
			countCap=config["Spotfire"]["countCap"],
			tpmCap=config["Spotfire"]["tpmCap"],
			fc=config["Spotfire"]["fc_thresh"],
			pcut=config["Spotfire"]["p_cut"],
			urlbase=config["ucsc_session"]
		group: "spotfire"
		threads: 1
		resources:
			maxtime=20,
			mem='16G'
		shell:
			"""
			Rscript {params.code}/deseq2_prepare_for_spotfire.R -C {input.cnt} -c {params.countCap} -D {input.de} -f {params.fc} -p {params.pcut} -s -t {params.tpmCap} -u {params.urlbase}
			"""


	rule stage_out_DE:
		input:
			spot1=expand(SCRATCH + "/results/DE/{countsfile}/{countsfile}.all.custom.deseq2.spotfire.1.txt.gz", countsfile=FLAVOURS),
			spot2=expand(SCRATCH + "/results/DE/{countsfile}/{countsfile}.all.custom.deseq2.spotfire.2.txt.gz", countsfile=FLAVOURS),
			spot3=expand(SCRATCH + "/results/DE/{countsfile}/{countsfile}.all.custom.deseq2.spotfire.3.txt.gz", countsfile=FLAVOURS),
			vanilla=expand(SCRATCH + "/results/DE/{countsfile}/{countsfile}.{dename}.deseq2.tsv", countsfile=FLAVOURS, dename=EXPECTED_DE),
		output:
			spot1=expand(HOME + "/results/DE/{countsfile}.all.custom.deseq2.spotfire.1.txt.gz", countsfile=FLAVOURS),
			spot2=expand(HOME + "/results/DE/{countsfile}.all.custom.deseq2.spotfire.2.txt.gz", countsfile=FLAVOURS),
			spot3=expand(HOME + "/results/DE/{countsfile}.all.custom.deseq2.spotfire.3.txt.gz", countsfile=FLAVOURS),
			vanilla=expand(HOME + "/results/DE/{countsfile}.{dename}.deseq2.tsv", countsfile=FLAVOURS, dename=EXPECTED_DE)
		params:
			outdir=HOME + "/results/DE/",
			indir=SCRATCH + "/results/DE/"
		shell:
			"""
			cp {input.spot1} {params.outdir}
			cp {input.spot2} {params.outdir}
			cp {input.spot3} {params.outdir}
			cp {input.vanilla} {params.outdir}
			cp {params.indir}/*/*.html {params.outdir} 
			cp {params.indir}/*/*.tpms.tsv {params.outdir} 
			cp {params.indir}/*/*.vsts.tsv {params.outdir} 
   			cp {params.indir}/*/*.norm.tsv {params.outdir} 
			"""
##


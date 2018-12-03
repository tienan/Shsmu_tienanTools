#==============================================================================#
#  RNA-seq Protocol: Raw Fastq Data in Multiple Datasets to Expression Matrix  #
#    Using: HISAT + StringTie + Ballgown                                       #
#                                                    Designed by: Wen-Xing Li  #
#                                                    Last Update: 2016-12-27   #
#==============================================================================#

# Require software and data:
# 1. Human Genome Reference: hg37_index, Homo_sapiens.GRCh37.75.gtf
# 2. FASTX-Toolkit: fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
# 3. SRAtoolkit: sratoolkit.current-centos_linux64.tar.gz
# 4. HISAT2: hisat2-2.0.4-Linux_x86_64.zip
# 5. StringTie: stringtie-1.2.3.Linux_x86_64.tar.gz
# 6. Ballgown: 	ballgown_2.6.0.tar.gz

# Step 1: 
# Copy and unzip human genome reference data and software 2-5 into Linux
# Install "ballgown" in R:
# source("https://bioconductor.org/biocLite.R");
# biocLite("ballgown");

# Step 2: Add software path to environment variable "PATH"
# PATH=$PATH:/home/liwenxing/Software/fastx_toolkit_0.0.13
# PATH=$PATH:/home/liwenxing/Software/sratoolkit.2.6.2-centos_linux64/bin
# PATH=$PATH:/home/liwenxing/Software/hisat2-2.0.4
# PATH=$PATH:/home/liwenxing/Software/stringtie-1.2.3.Linux_x86_64

# Step 3: Data preparation
# In NABI-SRA: .sra to .fastq
# Single
# fastq-dump --fasta /home/liwenxing/AD_RNAseq/RawData/SRR1931812.sra  -O /home/liwenxing/AD_RNAseq/Fastq
# Paired
# fastq-dump --split-3  ../SRR002018.sra -O out_path ==> SRR002018_1.fastq, SRR002018_2.fastq
# For .tar, .rar, .gz, .zip, ... files, please decompress to your data directory.
# Directory format (example):
# ../Alzheimer
#    |-GSE53697
#    | |-SRR2422918.fq
#    | |-SRR2422919.fq
#    | |-...
#    |-GSE67333
#      |-SRR1931812.fq
#      |-SRR1931813.fq
#      |-...

# Step 4: Set parameters (in === parameters === )
# dataPath: your RNA-seq data path, the data should be .fastq or .fq format
# outPath: your output results path
# alignment_genome_ref: ../hg37_index/hg37
# assembly_genome_ref: ../Homo_sapiens.GRCh37.75.gtf
# seq_type: 1: Single-read, 2: Paired-end

# Step 5: Run and wait
# In R environment:
# source("../RNAseq_Protocol.R");

library(ballgown);

# === parameters ===
dataPath <- "/home/liwenxing/SLE/RawData";
outPath <- "/home/liwenxing/SLE";
alignment_genome_ref <- "/home/liwenxing/genome_ref/hg37_index/hg37";
assembly_genome_ref <- "/home/liwenxing/genome_ref/Homo_sapiens.GRCh37.75.gtf";
seq_type <- 2;
# === parameters ===

resultsPath <- paste(outPath, "Results", sep="/");
if (!file.exists(resultsPath)) {
	dir.create(resultsPath);
}

qualityFilterPath <- paste(resultsPath, "1_Fastq_Quality_Filter", sep="/");
if (!file.exists(qualityFilterPath)) {
	dir.create(qualityFilterPath);
}
alignmentPath <- paste(resultsPath, "2_HISAT_Alignment", sep="/");
if (!file.exists(alignmentPath)) {
	dir.create(alignmentPath);
}
convertPath <- paste(resultsPath, "3_SAMtools_Convert", sep="/");
if (!file.exists(convertPath)) {
	dir.create(convertPath);
}
sortPath <- paste(resultsPath, "4_SAMtools_Sort", sep="/");
if (!file.exists(sortPath)) {
	dir.create(sortPath);
}
correctPath <- paste(resultsPath, "5_SAMtools_Correct", sep="/");
if (!file.exists(correctPath)) {
	dir.create(correctPath);
}
assemblyPath <- paste(resultsPath, "6_StringTie_Assembly", sep="/");
if (!file.exists(assemblyPath)) {
	dir.create(assemblyPath);
}
expressionPath <- paste(resultsPath, "7_Ballgown_Expression", sep="/");
if (!file.exists(expressionPath)) {
	dir.create(expressionPath);
}

perlCode <- "perl -ne 'if(/HI:i:(\\d+)/) { $m=$m1-1; $_ =~ s/HI:i:(\\d+)/HI:i:$m/} print $_;'";

matchFilesExtension <- function(filelist, extension) {
	for (i in c(1:length(filelist))) {
		isMatched <- FALSE;
		temp <- unlist(strsplit(filelist[i], ".", fixed=TRUE));
		if (extension == tolower(temp[length(temp)])) {
			isMatched <- TRUE;
		}
	}
	return(isMatched);
}

studyIDs <- dir(dataPath);
studyPaths <- paste(dataPath, studyIDs, sep="/");
for (i in c(1:length(studyIDs))) {
	tempFileNames <- dir(studyPaths[i]);
	filenum <- length(tempFileNames);
	isFastq <- FALSE;
	if (matchFilesExtension(tempFileNames, "fq") | matchFilesExtension(tempFileNames, "fastq")) {
		isFastq <- TRUE;
	}
	fastqFileNames <- tempFileNames;
	if (isFastq == FALSE) {
		for (j in c(1:filenum)) {
			fastqFileNames[j] <- paste(unlist(strsplit(tempFileNames[j], ".", fixed=TRUE))[1], ".fq", sep="")
		}
	}
	tempFilePaths <- paste(studyPaths[i], tempFileNames, sep="/");
	fastqFilePaths <- paste(studyPaths[i], fastqFileNames, sep="/");
	qualityFilterStudyPath <- paste(qualityFilterPath, studyIDs[i], sep="/");
	if (!file.exists(qualityFilterStudyPath)) {
		dir.create(qualityFilterStudyPath);
	}
	alignmentStudyPath <- paste(alignmentPath, studyIDs[i], sep="/");
	if (!file.exists(alignmentStudyPath)) {
		dir.create(alignmentStudyPath);
	}
	convertStudyPath <- paste(convertPath, studyIDs[i], sep="/");
	if (!file.exists(convertStudyPath)) {
		dir.create(convertStudyPath);
	}
	sortStudyPath <- paste(sortPath, studyIDs[i], sep="/");
	if (!file.exists(sortStudyPath)) {
		dir.create(sortStudyPath);
	}
	correctStudyPath <- paste(correctPath, studyIDs[i], sep="/");
	if (!file.exists(correctStudyPath)) {
		dir.create(correctStudyPath);
	}
	assemblyStudyPath <- paste(assemblyPath, studyIDs[i], sep="/");
	if (!file.exists(assemblyStudyPath)) {
		dir.create(assemblyStudyPath);
	}
	j <- 1;
	while (j <= filenum) {
		sampleName <- unlist(strsplit(fastqFileNames[j], ".", fixed=TRUE))[1];
		if (seq_type == 1) {
			if (isFastq == FALSE) {
				systemCode <- paste("gunzip -c \"", tempFilePaths[j], "\" > \"", fastqFilePaths[j], "\"", sep="");
				system(systemCode);
			}
			qualityFilterFilePath <- paste(qualityFilterStudyPath, fastqFileNames[j], sep="/");
			systemCode <- paste("fastq_quality_filter -Q 33 -q 20 -p 100 -i \"", fastqFilePaths[j], "\" -o \"", qualityFilterFilePath, "\"", sep="");
			system(systemCode);
			samFileName <- paste(sampleName, ".sam", sep="");
			alignmentFilePath <- paste(alignmentStudyPath, samFileName, sep="/");
			systemCode <- paste("hisat2 -x \"", alignment_genome_ref, "\" -U \"", qualityFilterFilePath, "\" -S \"", alignmentFilePath, "\"", sep="");
			system(systemCode);
		}
		if (seq_type == 2) {
			if (isFastq == FALSE) {
				systemCode <- paste("gunzip -c \"", tempFilePaths[j], "\" > \"", fastqFilePaths[j], "\"", sep="");
				system(systemCode);
				systemCode <- paste("gunzip -c \"", tempFilePaths[j+1], "\" > \"", fastqFilePaths[j+1], "\"", sep="");
				system(systemCode);
			}
			sampleName <- gsub("_[12]", "", sampleName);
			qualityFilterFilePath_1 <- paste(qualityFilterStudyPath, fastqFileNames[j], sep="/");
			systemCode <- paste("fastq_quality_filter -Q 33 -q 20 -p 100 -i \"", fastqFilePaths[j], "\" -o \"", qualityFilterFilePath_1, "\"", sep="");
			system(systemCode);
			qualityFilterFilePath_2 <- paste(qualityFilterStudyPath, fastqFileNames[j+1], sep="/");
			systemCode <- paste("fastq_quality_filter -Q 33 -q 20 -p 100 -i \"", fastqFilePaths[j+1], "\" -o \"", qualityFilterFilePath_2, "\"", sep="");
			system(systemCode);
			samFileName <- paste(sampleName, ".sam", sep="");
			alignmentFilePath <- paste(alignmentStudyPath, samFileName, sep="/");
			systemCode <- paste("hisat2 -x \"", alignment_genome_ref, "\" -1 \"", qualityFilterFilePath_1, "\" -2 \"", qualityFilterFilePath_2, "\" -S \"", alignmentFilePath, "\"", sep="");
			system(systemCode);
			j <- j + 1;
		}
		bamFileName <- paste(sampleName, ".bam", sep="");
		convertFilePath <- paste(convertStudyPath, bamFileName, sep="/");
		systemCode <- paste("samtools view -bS \"", alignmentFilePath, "\" > \"", convertFilePath, "\"", sep="");
		system(systemCode);
		sortFilePath <- paste(sortStudyPath, sampleName, sep="/");
		systemCode <- paste("samtools sort \"", convertFilePath, "\" \"", sortFilePath, "\"", sep="");
		system(systemCode);
		sortFilePath <- paste(sortStudyPath, bamFileName, sep="/");
		correctFilePath <- paste(correctStudyPath, bamFileName, sep="/");
		systemCode <- paste("samtools view -h \"", sortFilePath, "\" | ", perlCode, " | samtools view -bS - > \"", correctFilePath, "\"", sep="");
		system(systemCode);
		txtFileName <- paste(sampleName, ".txt", sep="");
		assemblySamplePath <- paste(assemblyStudyPath, sampleName, sep="/");
		if (!file.exists(assemblySamplePath)) {
			dir.create(assemblySamplePath);
		}
		assemblyFilePath <- paste(assemblySamplePath, txtFileName, sep="/");
		systemCode <- paste("stringtie \"", correctFilePath, "\" -p 25 -B -G \"", assembly_genome_ref, "\" -o \"", assemblyFilePath, "\"", sep="");
		system(systemCode);
		if (seq_type == 1) {
			file.remove(qualityFilterFilePath);
			if (isFastq == FALSE) {
				file.remove(fastqFilePaths[j]);
			}
		}
		if (seq_type == 2) {
			file.remove(qualityFilterFilePath_1);
			file.remove(qualityFilterFilePath_2);
			if (isFastq == FALSE) {
				file.remove(fastqFilePaths[j-1]);
				file.remove(fastqFilePaths[j]);
			}
		}
		file.remove(alignmentFilePath)
		file.remove(convertFilePath)
		file.remove(sortFilePath)
		file.remove(correctFilePath)
		j <- j + 1;
	}
	result <- ballgown(samples=dir(assemblyStudyPath, full.name=TRUE), meas="all");
	transcript_fpkm <- texpr(result, "FPKM");
	transcript_cov <- texpr(result, "cov");
	whole_tx_table <- texpr(result, "all");
	exon_mcov <- eexpr(result, "mcov");
	junction_rcount <- iexpr(result);
	whole_intron_table <- iexpr(result, "all");
	gene_expression <- gexpr(result);
	resultStudyPath <- paste(expressionPath, studyIDs[i], sep="/");
	if (!file.exists(resultStudyPath)) {
		dir.create(resultStudyPath);
	}
	write.table(transcript_fpkm, file=paste(resultStudyPath, "transcript_fpkm.txt", sep="/"), quote=FALSE, sep="\t");
	write.table(transcript_cov, file=paste(resultStudyPath, "transcript_cov.txt", sep="/"), quote=FALSE, sep="\t");
	write.table(whole_tx_table, file=paste(resultStudyPath, "whole_tx_table.txt", sep="/"), quote=FALSE, sep="\t");
	write.table(exon_mcov, file=paste(resultStudyPath, "exon_mcov.txt", sep="/"), quote=FALSE, sep="\t");
	write.table(junction_rcount, file=paste(resultStudyPath, "junction_rcount.txt", sep="/"), quote=FALSE, sep="\t");
	write.table(whole_intron_table, file=paste(resultStudyPath, "whole_intron_table.txt", sep="/"), quote=FALSE, sep="\t");
	write.table(gene_expression, file=paste(resultStudyPath, "gene_expression.txt", sep="/"), quote=FALSE, sep="\t");
}

# ÍøÖ·
# templist <- paste(resultsPath, "templist.txt", sep="/");
# system(paste("ls -l \"", fastqPaths[i], "\" > \"", templist, "\"", sep=""));
# system(paste("sort -k9 \"", templist, "\" -o \"", templist, "\"", sep=""));
# http://blog.sciencenet.cn/blog-759995-990471.html
# http://blog.biochen.com/archives/337
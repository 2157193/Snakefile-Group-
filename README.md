# Snakefile-Group-
Snakefile Group  2023
# Snakefile

# Define input samples
samples = ["SRR17237659.fastq", "SRR17237663.fastq"]

# Define rule for quality control with FastQC for raw fastq files
rule raw_fastqc:
    input:
        fastq1 = "SRR17237663.fastq",
        fastq2 = "SRR17237659.fastq"
    output:
        html1 = "results/fastqc/SRR17237663_fastqc.html",
        zip1 = "results/fastqc/SRR17237663_fastqc.zip",
        html2 = "results/fastqc/SRR17237659_fastqc.html",
        zip2 = "results/fastqc/SRR17237659_fastqc.zip"
    shell:
        """
        fastqc {input.fastq1} -o results/fastqc
        fastqc {input.fastq2} -o results/fastqc
        """

# Define rule for trimming the input fastq files
rule trim_fastq:
    input:
        fastq1 = "SRR17237663.fastq",
        fastq2 = "SRR17237659.fastq"
    output:
        trimmed1 = "data/fastq/{sample}_1.trimmed.fastq",
        trimmed2 = "data/fastq/{sample}_2.trimmed.fastq"
    shell:
        """
        # Trimming commands here
        trim_fastq_script {input.fastq1} {input.fastq2} {output.trimmed1} {output.trimmed2}
        """

# Define rule for quality control with FastQC for trimmed fastq files
rule trimmed_fastqc:
    input:
        fastq1 = "data/fastq/{sample}_1.trimmed.fastq",
        fastq2 = "data/fastq/{sample}_2.trimmed.fastq"
    output:
        html1 = "results/fastqc/trimmed_SRR17237663_fastqc.html",
        zip1 = "results/fastqc/trimmed_SRR17237663_fastqc.zip",
        html2 = "results/fastqc/trimmed_SRR17237659_fastqc.html",
        zip2 = "results/fastqc/trimmed_SRR17237659_fastqc.zip"
    shell:
        """
        fastqc {input.fastq1} -o results/fastqc
        fastqc {input.fastq2} -o results/fastqc
        """

# Define rule for generating a MultiQC report for raw fastq files
rule raw_multiqc:
    input:
        fastqc_reports = expand("results/fastqc/raw_{sample}_fastqc.html", sample=samples)
    output:
        html = "results/multiqc/raw_multiqc_report.html"
    shell:
        """
        multiqc results/fastqc -o results/multiqc
        """

# Define rule for generating a MultiQC report for trimmed fastq files
rule trimmed_multiqc:
    input:
        fastqc_reports = expand("results/fastqc/trimmed_{sample}_fastqc.html", sample=samples)
    output:
        html = "results/multiqc/trimmed_multiqc_report.html"
    shell:
        """
        multiqc results/fastqc -o results/multiqc
        """

# Define rule for aligning reads with HISAT2
rule hisat2:
    input:
        fastq1 = "data/fastq/{sample}_1.trimmed.fastq",
        fastq2 = "data/fastq/{sample}_2.trimmed.fastq"
    output:
        sam = "results/hisat2/{sample}.sam",
        summary = "results/hisat2/{sample}_summary.txt"
    shell:
        """
        hisat2 -x /path/to/genome_index -1 {input.fastq1} -2 {input.fastq2} -S {output.sam} 2> {output.summary}
        """

# Define rule for sorting and converting SAM to BAM
rule sam_to_bam:
    input:
        sam = "results/hisat2/{sample}.sam"
    output:
        bam = "results/hisat2/{sample}.bam",
        sorted_bam = "results/hisat2/{sample}_sorted.bam"
    shell:
        """
        samtools view -bS {input.sam} > {output.bam}
        samtools sort -o {output.sorted_bam} {input.bam}
        """

# Define rule for indexing BAM files
rule index_bam:
    input:
        sorted_bam = "results/hisat2/{sample}_sorted.bam"
    output:
        bai = "results/hisat2/{sample}_sorted.bam.bai"
    shell:
        """
        samtools index {input.sorted_bam}
        """

# Define rule for counting reads with featureCounts
rule featurecounts:
    input:
        sorted_bam = "results/hisat2/{sample}_sorted.bam"
    output:
        counts = "results/counts/{sample}.counts"
    shell:
        """
        featureCounts -a /path/to/annotation.gtf -o {output.counts} {input.sorted_bam}
        """

# Define rule for pathway enrichment analysis
rule pathway_enrichment:
    input:
        differential_expression="GSE190972_Normalised_smallrna_counts_matrix_PLAC.csv"
    output:
        enrichment_results="enrichment_results_file"
    shell:
        """
        # Pathway enrichment analysis commands here
        Rscript pathway_enrichment.R {input.differential_expression} {output.enrichment_results}
        """

# Define the final rule for all output files
rule all:
    input:
        expand("results/fastqc/raw_{sample}_fastqc.html", sample=samples),
        expand("results/fastqc/trimmed_{sample}_fastqc.html", sample=samples),
        "results/multiqc/raw_multiqc_report.html",
        "results/multiqc/trimmed_multiqc_report.html",
        expand("results/hisat2/{sample}.sam", sample=samples),
        expand("results/hisat2/{sample}_summary.txt", sample=samples),
        expand("results/hisat2/{sample}_sorted.bam", sample=samples),
        expand("results/hisat2/{sample}_sorted.bam.bai", sample=samples),
        expand("results/counts/{sample}.counts", sample=samples),
        "Result.csv",
        "enrichment_results_file"


       VISUALIZATION
       Volcano plot
# Create the volcano plot
ggplot(result_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(size = 1, color = ifelse(result_df$padj < 0.01, "red", "black")) +
  xlim(c(-3, 3)) + ylim(c(0, 10)) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "red") +
  xlab("log2 Fold Change") + ylab("-log10 p-value") +
  ggtitle("Volcano Plot") +
  theme_bw()
## Warning: Removed 19400 rows containing missing values (`geom_point()`).
 
Heatmap
de_genes <- rownames(resSig)
counts_subset <- countsN[de_genes, ]
scaled_counts <- scale(counts_subset)
pheatmap(scaled_counts, 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         show_rownames = TRUE, show_colnames = TRUE)
 
Hierarchical Clustering
dist_mat <- dist(t(countsN))
hclust_res <- hclust(dist_mat)
clusters <- cutree(hclust_res, k = 3)
plot(hclust_res, hang = -1, main = "Hierarchical Clustering")
rect.hclust(hclust_res, k = 3, border = "red")
 
Correlation Heatmap
# Calculate the correlation matrix
correlation_matrix <- cor(countsN)

# Create the correlation heatmap
heatmap(correlation_matrix,
        col = colorRampPalette(c("#4575b4", "#ffffbf", "#d73027"))(100),
        main = "Correlation Heatmap")
 
Downregulated Genes
resSigrep <- read.csv("Result.csv")
resSigrep <- resSigrep[resSigrep$log2FoldChange < 0, ]
resSigrep

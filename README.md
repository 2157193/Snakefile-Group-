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
##                  X    baseMean log2FoldChange     lfcSE      stat       pvalue
## 17  hsa-miR-210-3p   584.26381      -1.509530 0.3583774 -4.212125 2.529795e-05
## 21    hsa-miR-4454   407.30221      -1.551082 0.4329283 -3.582770 3.399703e-04
## 25 hsa-miR-200a-3p    30.03784      -2.506626 0.7141131 -3.510126 4.478952e-04
## 38     hsa-miR-375    90.14094      -2.721609 0.7571055 -3.594756 3.246962e-04
## 42  hsa_piR_016856    14.25866      -3.113135 0.8354959 -3.726093 1.944707e-04
## 45  hsa_piR_012129    54.09109      -2.193192 0.6189310 -3.543517 3.948286e-04
## 48  hsa_piR_021242    29.59546      -4.848674 0.8709721 -5.566968 2.592104e-08
## 49  hsa_piR_002468   113.88339      -1.744133 0.4910129 -3.552112 3.821515e-04
## 53  hsa_piR_000295    16.71216      -6.538386 1.3623765 -4.799251 1.592604e-06
## 54  hsa_piR_008488 14724.77345      -1.406834 0.3909533 -3.598471 3.200930e-04
## 56  hsa_piR_011324  2144.28509      -1.332360 0.3794271 -3.511505 4.455775e-04
##            padj
## 17 1.514900e-03
## 21 8.038710e-03
## 25 8.472828e-03
## 38 8.038710e-03
## 42 5.999127e-03
## 45 8.038710e-03
## 48 1.319381e-05
## 49 8.038710e-03
## 53 2.026589e-04
## 54 8.038710e-03
## 56 8.472828e-03
Up regulated Genes
resSigind <- read.csv("Result.csv")
resSigind <- resSigind[resSigind$log2FoldChange > 0, ]
resSigind
##                  X    baseMean log2FoldChange     lfcSE     stat       pvalue
## 1    hsa-miR-16-5p   818.19728       1.721343 0.4436858 3.879645 1.046091e-04
## 2     hsa-miR-1-3p   222.65495       1.698357 0.4457030 3.810512 1.386790e-04
## 3  hsa-miR-519b-3p   249.39676       1.486884 0.4108911 3.618682 2.961071e-04
## 4   hsa-miR-369-3p   142.26985       1.558947 0.4492458 3.470144 5.201800e-04
## 5   hsa-miR-501-5p    26.85552       2.323915 0.6202109 3.746975 1.789798e-04
## 6   hsa-miR-652-3p   484.76456       2.192305 0.4522098 4.847981 1.247243e-06
## 7   hsa-miR-214-3p   836.15615       1.746879 0.4299645 4.062844 4.847836e-05
## 8  hsa-miR-125a-3p    94.25059       1.603238 0.4356421 3.680173 2.330757e-04
## 9   hsa-miR-504-5p   580.16480       1.912891 0.3963389 4.826404 1.390205e-06
## 10  hsa-miR-377-5p    68.86600       1.751699 0.4907622 3.569345 3.578752e-04
## 11 hsa-miR-6501-5p    58.85009       1.862377 0.5222685 3.565937 3.625581e-04
## 12  hsa-miR-431-3p   377.89726       1.680164 0.3681010 4.564410 5.009003e-06
## 13 hsa-miR-3614-5p    98.20993       2.253915 0.5449438 4.136050 3.533355e-05
## 14 hsa-miR-3120-5p    23.82424       2.108697 0.6016557 3.504823 4.569107e-04
## 15  hsa-miR-324-5p   135.50728       2.657915 0.5678689 4.680508 2.861647e-06
## 16  hsa-miR-337-3p    64.16120       1.603852 0.4570198 3.509370 4.491689e-04
## 18  hsa-miR-373-3p   237.04618       1.364138 0.3812908 3.577686 3.466501e-04
## 19  hsa-miR-382-5p   254.20358       1.456353 0.4088771 3.561835 3.682717e-04
## 20    hsa-miR-3182  6545.90766       2.455828 0.6281070 3.909889 9.233863e-05
## 22  hsa-miR-485-3p   191.96939       1.924108 0.3858395 4.986809 6.138466e-07
## 23  hsa-miR-145-5p  5260.74172       2.018997 0.4451156 4.535894 5.735992e-06
## 24 hsa-miR-487a-5p    28.23904       2.062909 0.5796233 3.559051 3.721967e-04
## 26 hsa-miR-500a-5p    22.38850       2.437643 0.6956098 3.504325 4.577658e-04
## 27  hsa-miR-15b-5p   105.94131       1.716911 0.4557444 3.767268 1.650440e-04
## 28  hsa-miR-139-5p   586.82433       2.106377 0.4417906 4.767818 1.862321e-06
## 29  hsa-miR-432-5p   401.04872       2.253810 0.4003662 5.629372 1.808674e-08
## 30 hsa-miR-146a-5p   686.56110       1.856862 0.4699721 3.951005 7.782378e-05
## 31   hsa-let-7e-5p   814.10701       1.791356 0.3696820 4.845667 1.261869e-06
## 32  hsa-miR-139-3p    26.81505       1.959365 0.5638201 3.475160 5.105479e-04
## 33  hsa-miR-654-5p   106.01432       1.778917 0.4387122 4.054860 5.016432e-05
## 34  hsa-miR-299-5p    55.47315       1.817474 0.5127656 3.544453 3.934280e-04
## 35  hsa-miR-584-5p  1147.20260       1.863352 0.4972320 3.747450 1.786417e-04
## 36  hsa-miR-197-3p   819.75722       1.544270 0.4180087 3.694348 2.204517e-04
## 37 hsa-miR-1301-3p   171.07555       1.870091 0.4001209 4.673816 2.956549e-06
## 39  hsa-miR-335-3p   411.00583       1.729935 0.4717064 3.667398 2.450313e-04
## 40  hsa_piR_019752   347.22620       1.761797 0.4875568 3.613522 3.020661e-04
## 41  hsa_piR_020499 12145.61922       1.591482 0.4157974 3.827541 1.294296e-04
## 43  hsa_piR_006950    12.72298       3.316538 0.8749502 3.790545 1.503169e-04
## 44  hsa_piR_016866    82.12237       2.123981 0.5018900 4.231964 2.316590e-05
## 46  hsa_piR_009295  1942.41180       2.405469 0.5617274 4.282271 1.849957e-05
## 47  hsa_piR_000823  2272.44709       2.865879 0.5943535 4.821842 1.422385e-06
## 50  hsa_piR_019914   199.42598       2.676568 0.6396664 4.184318 2.860234e-05
## 51  hsa_piR_019825   635.37666       3.062360 0.6690491 4.577183 4.712795e-06
## 52  hsa_piR_004153   136.15021       1.514243 0.4183697 3.619390 2.952978e-04
## 55  hsa_piR_019822    52.48511       3.429339 0.8241060 4.161284 3.164626e-05
## 57         SNORD48  1266.86872       1.503999 0.3785848 3.972689 7.106592e-05
##            padj
## 1  4.095850e-03
## 2  5.041974e-03
## 3  7.884701e-03
## 4  9.290232e-03
## 5  5.693796e-03
## 6  2.026589e-04
## 7  2.321240e-03
## 8  6.779174e-03
## 9  2.026589e-04
## 10 8.038710e-03
## 11 8.038710e-03
## 12 3.922435e-04
## 13 1.798478e-03
## 14 8.472828e-03
## 15 2.736152e-04
## 16 8.472828e-03
## 18 8.038710e-03
## 19 8.038710e-03
## 20 3.760029e-03
## 22 2.026589e-04
## 23 4.170886e-04
## 24 8.038710e-03
## 26 8.472828e-03
## 27 5.600495e-03
## 28 2.106492e-04
## 29 1.319381e-05
## 30 3.301025e-03
## 31 2.026589e-04
## 32 9.281032e-03
## 33 2.321240e-03
## 34 8.038710e-03
## 35 5.693796e-03
## 36 6.600582e-03
## 37 2.736152e-04
## 39 6.928941e-03
## 40 7.884701e-03
## 41 4.879977e-03
## 43 5.276641e-03
## 44 1.473931e-03
## 46 1.255504e-03
## 47 2.026589e-04
## 50 1.617621e-03
## 51 3.922435e-04
## 52 7.884701e-03
## 55 1.695573e-03
## 57 3.145440e-03
![image](https://github.com/2157193/Snakefile-Group-/assets/117467941/fbb289b1-08ac-4ea7-ba0d-699d720f4d3a)


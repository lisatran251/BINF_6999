# Evaluating Primer Efficiency and Antimicrobial Resistance Gene Distribution Across Publicly Available Metagenomic Datasets

## Introduction
Antimicrobial resistance (AMR) has emerged as a critical global public health crisis, with projections indicating up to 10 million deaths annually by 2050 (Block M & Blanchard DL, 2023; Christaki et al., 2020). This alarming trend is largely due to the indiscriminate use of antibiotics in human health, veterinary medicine, and agriculture, leading to an increase in resistant bacterial strains. Recognizing the magnitude and potential impacts of this issue is essential for developing effective countermeasures (Harbarth et al., 2015).

In response to this challenge, the Diversity of Antibiotic Resistance genes and Transfer Elements-Quantitative Monitoring (DARTE-QM) system has been established. Utilizing the power of high-throughput sequencing, DARTE-QM quantifies antibiotic resistance genes (ARGs) in environmental samples. Surpassing traditional detection methods, DARTE-QM offers a combination of specificity, sensitivity, and cost-effectiveness, representing a significant enhancement in diagnostic capabilities (Smith et al., 2022).

The initial deployment of DARTE-QM targeted a spectrum of 662 ARGs, examining a variety of environmental samples such as manure, soil, and livestock feces. These preliminary investigations have laid the groundwork for subsequent evaluations, aimed at assessing the performance of the initial primer set used in DARTE-QM. The primary objective is to evaluate the effectiveness of each primer pair and to estimate the prevalence of targeted genes in diverse metagenomic sequences available in public databases.

The research objectives of this study are threefold. First, to develop an in-silico PCR analysis tool using existing primers for predicting PCR product abundance. Second, to map the distribution of AMR genes within publicly available metagenomic datasets and identify key AMR targets for further investigation. Third, to compare the outcomes of the in-silico PCR analysis with results from the ABRicate program, which screens contigs for antimicrobial resistance and virulence genes using the Resfinder and CARD databases.

This paper aims to highlight the severity of AMR and outline the methodological approach of DARTE-QM in addressing this critical issue. Our study contributes to the scientific community’s efforts in monitoring, understanding, and combating the spread of antibiotic resistance. The One Health approach, recognizing the interconnectedness of human, animal, and environmental health, is particularly relevant. By adopting this holistic framework, our research not only targets specific ARGs but also considers the broader ecological and societal contexts in which these genes circulate. This comprehensive perspective is essential for devising more effective and sustainable strategies to tackle AMR (Velazquez-Meza et al., 2022). Finally, ABRicate will be utilized for the mass screening of contigs for antimicrobial resistance, employing the Resfinder and CARD databases as tools to cross-validate program's final results (Seemann T, 2023).

## Methods
### Bioinformatics Toolkit and Environment 
The analysis was carried out in Python (Python Software Foundation, Python Language Reference, version 3.7.7) using various libraries including pandas (for data manipulation), BioPython (for handling biological computation), argparse (for command-line scripting), os (for interacting with the operating system), re (for regular expression operations), and time (for tracking execution time).

### Program Workflow
### 'run_program.sh': 
The present method is designed to enable efficient and distributed analysis of DNA sequences through the Slurm Workload Manager. We chose the Bash Shell script for this task due to its ability to handle extensive workloads and to address the challenges posed by processing large datasets, such as those containing millions of sequences. The Python script used in previous experiments required approximately 45 minutes to process 10,000 sequences, making the direct analysis of vast datasets impractical. The proposed solution utilizes Shell scripting and the Slurm Workload Manager to partition contigs into smaller chunks, allowing parallel processing and reducing overall execution time. This method enhances the reproducibility and consistency of the analysis by encapsulating the entire workflow in a single user command, streamlining the analysis, and facilitating error recovery.

The method involves a series of scripts that are automatically submitted as separate jobs to the Slurm Workload Manager. The primary Bash Shell script, "run_program.sh," serves as the main orchestrator for the analysis. The user initiates the entire workflow with a single command:
```sbatch run_program.sh contig_fasta_file primer_file chunk_size email_address```

### Input Parameters:
•	contig_fasta_file: Input FASTA file containing DNA sequences to be analyzed (contigs).

•	primer_file: CSV file containing primer sequences. 

•	chunk_size: Number of sequences to be processed per chunk for parallelization.

•	email_address: Email address of the user to fetch sequences from NCBI.

## Workflow Execution:
1.	Chunking of Contigs: The input FASTA file is split into smaller chunks based on the specified chunk_size. Each chunk is saved in the "chunks" directory.
2.	Parallel Processing: A job is submitted for each chunk to process the DNA sequences in parallel. The "process1.py" script is used to extract products between primer pairs within each chunk. Jobs are executed on the Slurm cluster nodes with specific resource requirements.
3.	Results Aggregation: After the successful completion of all "process1.py" jobs, a second script, "process2.py," is executed. This script processes the results generated by the individual chunks and consolidates them into a single "raw_results.csv" file.
4.	Final Analysis: Upon the successful completion of "process2.py," a third script, "process3.py," is executed. This script performs any final analysis steps required and generates the "final_result.csv" file.
   
### 'process1.py': 
This script extracts product sequences between primers from FASTA-formatted DNA sequences. It cleans and preprocesses the data, ensuring data integrity. The script utilizes BioPython to import primer information, including reverse complements, from a CSV file. It then identifies primer locations in the DNA sequences using regular expressions and extracts product sequences. The results are stored in a CSV file, providing valuable information for further analyses.
1. Data Cleaning and Preprocessing: In this step, any missing or inconsistent entries in the dataset were removed to ensure the integrity of the data. The script was executed separately using the command:
```python3 process1.py fasta_file primer_file```
If there were any inconsistencies in the naming of the primer file, the script automatically handled them by removing, replacing, or adjusting the entries to the appropriate format. This process ensured the smooth execution of subsequent analyses.
2. Primer Information Extraction: The product extraction process commenced with the import of primer sequences from a CSV file. The primer information was organized in a tabular format using pandas, a popular data manipulation library in Python. The unique forward and reverse primer sequences were then isolated from the data frame. Additionally, the reverse complements of each primer were generated using the 'Seq' module from BioPython, a comprehensive library for bioinformatics tasks.
3. FASTA File Parsing: To process the FASTA file containing the target DNA sequences, the 'SeqIO' module from BioPython was employed. This module facilitated the efficient parsing of the sequences from the FASTA file, making them available for subsequent analysis.
4. Primer Sequence Search and Product Extraction: The script combined all the forward and reverse primers from both files, along with their corresponding reverse complements, into a single dictionary. Utilizing regular expressions, the script then searched for the location of each primer within the DNA sequences from the FASTA file. If a sequence contained two or more primers, the script identified and sorted the positions where the primers were found. The first identified primer location served as the start position, while the subsequent ones were considered end positions. Based on this information, the product sequence between the start and end primers was extracted.
5. Origin Identification of Primers: The origin of each primer, whether it belonged to the forward primer group ('F_primers'), the reverse primer group ('R_primers'), or one of their respective reverse complements, was traced back and recorded in the 'Combination' column. This step was performed after finding all the primers, as it was computationally less expensive and more efficient.
6. Results Storage and Reporting: All the relevant information associated with the product, including the 'Record ID', 'Start Primer', 'End Primer', 'Start Position', 'End Position', 'Length', 'Combination', and 'Product', was recorded. The script organized this information and returned it as a CSV file named 'raw_results.csv'. This file contained the detailed results of the product extraction process, making it easily accessible for further analysis and interpretation.
7. Execution Time Measurement: Due to the potentially time-consuming nature of the product extraction process, the script also measured the total execution time. The execution time was printed out, providing valuable information and references for assessing the performance of the script during each run.

### 'process2.py': 
The purpose of this script is to analyze primer-targeted DNA sequences, identify qualified products and execute the pipeline and produce final results in CSV format. The pipeline involves data loading, data preprocessing, primer-target matching, product filtering, and results storage. The script was executed separately using the command:
```python3 process2.py primer_file raw_results.csv```
1. Primer-Target Matching: The script iterates through the 'results' DataFrame, analyzing the 'Combination' column to determine the primer orientation in each product. Based on the combination, the script converts the product sequence, forward primer, and reverse primer to their original direction using the 'complement_dict'. The script then filters the results DataFrame to include only valid combinations: 'F_primers-R_primers', 'R_primers-F_primers', 'R_reverse_complement-F_reverse_complement', and 'F_reverse_complement-R_reverse_complement'.
2. Product Filtering: To meet the requirements for qualified products, the script filters the results based on product length. Only products with lengths between 150 and 400 nucleotides are considered for further analysis.
3. Primer-Target Association: The script searches for matching forward and reverse primers in the 'primers' DataFrame for each qualified product. If a matching pair of primers is found, the corresponding target gene is assigned to the product in the 'Target' column.
4. Qualification Determination: The script evaluates whether a product qualifies as a target gene match based on its length and primer association. If a qualified target gene is assigned to a product, the 'Qualified?' column is marked as 'Yes'; otherwise, it is marked as 'No'.
5. Results Storage: The final results, including product information, target gene associations, and qualification status, are saved to a new CSV file named 'final_result.csv'.
6. Additional Output: For better visualization and further analysis, the pipeline generates two additional CSV files: 'match.csv' containing the qualified products and 'nonmatch.csv' containing the non-qualified products.

### 'process3.py': 
This script automates the process of gene identification and sequence analysis including ABRicate integration, sequence fetching from NCBI, gene name processing, scenario-based analysis and BLAST analysis. The script was executed separately using the command:
```python3 process3.py email_address original_fasta_file primer_file result_file```
1. Run ABRicate at 10% Coverage: Processes the contig FASTA file using ABRicate and prepares results for further analysis.
2. Data Preparation and Cleaning: Loads primer and result data, and cleans it for processing. This includes mapping original contig sequences and filtering results.
3. Gene Name Conversion: Applies functions to standardize gene names across different datasets for accurate comparison.
4. Data Merging: Merges program results with ABRicate findings, applying logic to handle various scenarios and identify genes found by both or exclusively by one method.
5. BLAST Analysis: Runs BLAST for contigs found exclusively by ABRicate to check for potential mismatches, and analyzes the BLAST report to update mismatch information.
6. Final Data Compilation: Combines and finalizes data from different sources into comprehensive CSV reports, including reasons for mismatches and suggestions for potential fixes.
7. Output Files: Generates multiple output files, such as final_report.csv and abricate_only.fasta, for downstream analysis.

### 'process4.py': 
This script is designed for analyzing and visualizing genomic data, specifically focusing on the results of antimicrobial resistance gene detection. This script processes data from various sources, performs statistical analysis, and generates insightful visualizations. The script was executed separately using the command:
```python3 process4.py final_report_file primer_file```
1. Data Loading and Preprocessing: 
    - Reads in genomic result data and primer information from provided files.
    - Preprocesses data by removing unnecessary columns, handling individual assemblies, and setting up for further analysis.
2. Statistical Analysis: 
    - Calculates total products found, along with categorization based on the detection method.
    - Analyzes primer usage and frequency, identifying primers that were not found in the analysis.
3. Data Visualization: 
    - Creates bar plots to compare the number of products found in each category.
    - Generates a heatmap for primer combinations leading to correct target genes.
    - Plots target gene frequencies and identifies the top 30 most found target genes.
4. Output File Generation:
    - Saves all analysis results, including figures and data summaries, to an output text file (analysis_result.txt).
    - Exports various plots as image files for visual insights.
5. Output Files
product_count.png: Bar chart showing the comparison of product findings in each category.
target_frq.png: Bar plot of target gene counts.
target_genes.png: Scatter plot showing the top 30 target genes found.
heatmap.png: Heatmap of primer combinations leading to correct target genes.
analysis_result.txt: A comprehensive text file including all analysis summaries and references to the figures.

## Citation
Achard, A., Villers, C., Pichereau, V., & Leclercq, R. (2005). New lnu (C) Gene Conferring Resistance to Lincomycin by Nucleotidylation in Streptococcus agalactiae UCN36. Antimicrobial Agents and Chemotherapy, 49(7), 2716–2719. https://doi.org/10.1128/AAC.49.7.2716-2719.2005

Amini, F., Krimpour, H. A., Ghaderi, M., Vaziri, S., Ferdowsi, S., Azizi, M., & Amini, S. (2018). Prevalence of Aminoglycoside Resistance Genes in Enterococcus Strains in Kermanshah, Iran. Iranian Journal of Medical Sciences, 43(5), 487–493.

Block M, & Blanchard DL. (2023). StatPearls. StatPearls.

Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 34(17), i884–i890. https://doi.org/10.1093/bioinformatics/bty560

Christaki, E., Marcou, M., & Tofarides, A. (2020). Antimicrobial Resistance in Bacteria: Mechanisms, Evolution, and Persistence. Journal of Molecular Evolution, 88(1), 26–40. https://doi.org/10.1007/s00239-019-09914-3

Flasche, S., & Atkins, K. E. (2018). Balancing Benefits and Risks of Antibiotic Use. The Journal of Infectious Diseases, 218(9), 1351–1353. https://doi.org/10.1093/infdis/jiy344

Green, K. D., Punetha, A., Hou, C., Garneau-Tsodikova, S., & Tsodikov, O. V. (2019). Probing the Robustness of Inhibitors of Tuberculosis Aminoglycoside Resistance Enzyme Eis by Mutagenesis. ACS Infectious Diseases, 5(10), 1772–1778. https://doi.org/10.1021/acsinfecdis.9b00228

Harbarth, S., Balkhy, H. H., Goossens, H., Jarlier, V., Kluytmans, J., Laxminarayan, R., Saam, M., Van Belkum, A., & Pittet, D. (2015). Antimicrobial resistance: one world, one fight! Antimicrobial Resistance and Infection Control, 4(1), 49. https://doi.org/10.1186/s13756-015-0091-2

Holman, D. B., Gzyl, K. E., & Kommadath, A. (2023). The gut microbiome and resistome of conventionally vs. pasture-raised pigs. Microbial Genomics, 9(7). https://doi.org/10.1099/mgen.0.001061

Li, D., Liu, C.-M., Luo, R., Sadakane, K., & Lam, T.-W. (2015). MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. Bioinformatics, 31(10), 1674–1676. https://doi.org/10.1093/bioinformatics/btv033

Seemann T. (2023). Abricate. Github. https://github.com/tseemann/abricate

Smith, S. D., Choi, J., Ricker, N., Yang, F., Hinsa-Leasure, S., Soupir, M. L., Allen, H. K., & Howe, A. (2022). Diversity of Antibiotic Resistance genes and Transfer Elements-Quantitative Monitoring (DARTE-QM): a method for detection of antimicrobial resistance in environmental samples. Communications Biology, 5(1), 216. https://doi.org/10.1038/s42003-022-03155-9

Sydenham, T. V., Overballe-Petersen, S., Hasman, H., Wexler, H., Kemp, M., & Justesen, U. S. (2019). Complete hybrid genome assembly of clinical multidrug-resistant Bacteroides fragilis isolates enables comprehensive identification of antimicrobial-resistance genes and plasmids. Microbial Genomics, 5(11). https://doi.org/10.1099/mgen.0.000312

Tang, K. W. K., Millar, B. C., & Moore, J. E. (2023). Antimicrobial Resistance (AMR). British Journal of Biomedical Science, 80. https://doi.org/10.3389/bjbs.2023.11387

UZUN, B., GÜNGÖR, S., PEKTAŞ, B., AKSOY GÖKMEN, A., YULA, E., KOÇAL, F., & KAYA, S. (2014). Macrolide-Lincosamide-Streptogramin B (MLSB) Resistance Phenotypes in Clinical Staphylococcus Isolates and Investigation of Telithromycin Activity. Mikrobiyoloji Bulteni, 48(3), 469–476. https://doi.org/10.5578/mb.7748

Velazquez-Meza, M. E., Galarde-López, M., Carrillo-Quiróz, B., & Alpuche-Aranda, C. M. (2022). Antimicrobial resistance: One Health approach. Veterinary World, 743–749. https://doi.org/10.14202/vetworld.2022.743-749

Zankari, E., Hasman, H., Cosentino, S., Vestergaard, M., Rasmussen, S., Lund, O., Aarestrup, F. M., & Larsen, M. V. (2012). Identification of acquired antimicrobial resistance genes. Journal of Antimicrobial Chemotherapy, 67(11), 2640–2644. https://doi.org/10.1093/jac/dks261

Zhan, Y., Liu, Y., Lin, J., Fu, X., Zhuang, C., Liu, L., Xu, W., Li, J., Chen, M., Zhao, G., Huang, W., & Cai, Z. (2015). Synthetic Tet-inducible artificial microRNAs targeting β-catenin or HIF-1α inhibit malignant phenotypes of bladder cancer cells T24 and 5637. Scientific Reports, 5(1), 16177. https://doi.org/10.1038/srep16177

## Author
Lisa Thuy Duyen Tran | thuyduye@uoguelph.ca

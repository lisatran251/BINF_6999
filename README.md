# Evaluating Primer Efficiency and Antimicrobial Resistance Gene Distribution Across Publicly Available Metagenomic Datasets

## Introduction
In recent years, the rise of AMRs on a global scale has been alarming, spreading between countries at an unprecedented rate. Superbugs and multidrug-resistant bacteria are now prevalent in many regions globally (Christaki et al., 2020). This resistance surge can be largely attributed to the overuse and misuse of antibiotics across various sectors, including human medicine, agriculture, animal farming, and industry (Harbarth et al., 2015). As a result, addressing this issue necessitates a One-Health perspective, recognizing the interconnectedness of human health with animal health and environmental sustainability. Monitoring ARGs in environmental ecosystems is challenging due to the high costs and inefficiencies of current methods. 

In 2019, a new method called DARTE-QM for detecting ARGs in environmental samples was introduced. DARTE-QM employs TruSeq high-throughput sequencing, enabling the simultaneous sequencing of thousands of antibiotic resistance gene targets that span a wide array of resistance classes prevalent in environmental systems (Smith et al., 2022). The methodology involves synthesizing primers to target specific ARGs, proceeding with barcode-multiplexing, and preparing an amplicon library for sequencing. One of DARTE-QM's primary objectives is to overcome the financial constraints inherent in metagenomic approaches to ARG monitoring in environmental samples, which typically require substantial sequencing depth and coverage for accurate identification.

The study showed that DARTE-QM was successful in detecting 662 of diverse ARGs across different types of samples, including soil, manure, water, and mock-community samples. Given that DARTE-QM can enrich and detect ARGs present in low abundance—a common limitation in shotgun metagenomics—this study suggests that DARTE-QM holds substantial promise as a tool for investigating antibiotic resistance in environmental microbiomes. The second phase of DARTE-QM platform development prioritizes primer targets based on risk characterization and prevalence. This process involves evaluating the initial primer set to gauge the effectiveness of each primer pair used in DARTE-QM's inaugural run, as well as the anticipated abundance of the targeted genes across a range of metagenomes available in public databases. Thus, the primary goal of this project is to create a program capable of performing in silico PCR using the ~1500 primers already developed. 

This program will aim to determine the theoretical abundance of the expected PCR products and identify any potential artifacts that may be generated. Additionally, the program's results will be benchmarked against the existing ABRicate program, which conducts comprehensive screenings of contigs for antimicrobial resistance or virulence genes using the Resfinder database (Seemann T, 2023; Zankari et al., 2012). Once developed, the program will analyze the theoretical distribution of AMR genes across various publicly available metagenomic datasets, aiming to identify priority AMR targets within the One Health framework. For our project, we have selected Python and Shell scripting, both known for their adaptability and efficiency. To ensure optimal performance and efficiency, we've chosen the Graham HPC system as our primary platform for execution. This strategic combination ensures that our research is both accurate and efficient, allowing us to delve deep into our datasets with precision.

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
This script was developed to analyze AMGs in contig FASTA files using both ABRicate and a custom Python script. The pipeline involved several steps, including running ABRicate, fetching sequences from NCBI, matching genes between ABRicate and the program, exploring different scenarios, and conducting sequence comparisons. The script was executed separately using the command:
```python3 process3.py email_address primer_file result_file original_fasta_file```
1. ABRicate and Unix Environment: ABRicate can only be executed in a Unix environment. To overcome the need for switching between Python and Unix, a function was defined in the Python script to run Unix commands. Separate Bash Shell scripts were created for each command required to be run.
2. ABRicate Execution and Gene Fetching: ABRicate was run to identify AMGs in the contig FASTA file using the 'resfinder' database. For each gene found by ABRicate, the original fasta sequence was fetched from NCBI using the provided accession number. The sequence was then added as a new column in the ABRicate result DataFrame, named 'originalSEQUENCE', using the fetch_sequence function.
3. Gene Name Conversion: To match the gene names between the original primer file and ABRicate results, a function was defined to shorten the target gene names. This was necessary as the naming systems between the two datasets were not similar.
4. Gene and Variant Matching: Matching between ABRicate and the program results was performed at two levels: gene level and variant level. Genes with similar names and variants found in both ABRicate and the program were combined into a new CSV file called 'both_program_ab_result.csv'.
5. Scenario 1 - Target Genes with Different Variants: For target genes found in both ABRicate and the program but with different variants, the original contigs containing the program-found genes were extracted from the original FASTA file. Sequence similarity between these program-found genes and the corresponding ABRicate-found genes was explored using CD-HIT. However, due to significant differences in sequence lengths, an alternative solution was considered to compare the sequences effectively.
6. Scenario 2 - Genes Found Only by the Program: Genes found only by the program were extracted from the original contigs into a new FASTA file. ABRicate was re-run on this file using a lower coverage percentage (10%) to capture all possible target genes missing from the first run. The results were compared with the program-only results from the first run, and matching genes were saved to 'merged_program_only_ab10_match.csv', while non-matching genes were saved to 'merged_program_only_ab10_nonmatch.csv'.
7. Scenario 3 - Genes Found Only by ABRicate: To identify the reasons why ABRicate found genes that the program did not, a 'Reason' column was created in the ABRicate-only DataFrame. If the gene name was available in the original primer file, the 'Reason' column was filled with 'Mismatches'; otherwise, it was filled with 'Primers not available'. Contigs associated with genes found only by ABRicate were extracted into a new FASTA file named 'abricate_only.fasta'. Primers were then run against these contigs using the BLAST command in Unix to detect potential mismatches and alignment with the extracted contigs.

## Citation
Christaki, E., Marcou, M., & Tofarides, A. (2020). Antimicrobial Resistance in Bacteria: Mechanisms, Evolution, and Persistence. Journal of Molecular Evolution, 88(1), 26–40. https://doi.org/10.1007/s00239-019-09914-3
Harbarth, S., Balkhy, H. H., Goossens, H., Jarlier, V., Kluytmans, J., Laxminarayan, R., Saam, M., Van Belkum, A., & Pittet, D. (2015). Antimicrobial resistance: one world, one fight! Antimicrobial Resistance and Infection Control, 4(1), 49. https://doi.org/10.1186/s13756-015-0091-2
Seemann T. (2023). Abricate. Github. https://github.com/tseemann/abricate
Smith, S. D., Choi, J., Ricker, N., Yang, F., Hinsa-Leasure, S., Soupir, M. L., Allen, H. K., & Howe, A. (2022). Diversity of Antibiotic Resistance genes and Transfer Elements-Quantitative Monitoring (DARTE-QM): a method for detection of antimicrobial resistance in environmental samples. Communications Biology, 5(1), 216. https://doi.org/10.1038/s42003-022-03155-9
Zankari, E., Hasman, H., Cosentino, S., Vestergaard, M., Rasmussen, S., Lund, O., Aarestrup, F. M., & Larsen, M. V. (2012). Identification of acquired antimicrobial resistance genes. Journal of Antimicrobial Chemotherapy, 67(11), 2640–2644. https://doi.org/10.1093/jac/dks261

## Author
Lisa Thuy Duyen Tran

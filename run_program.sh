#!/bin/bash
#SBATCH --account=def-nricker
#SBATCH --time=0-02:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=1G

# How to run: sbatch run_program.sh input_file primer_file chunk_size email address
# Example: sbatch run_program.sh contigs_ex.fasta DARTE-QM_primer_design.csv 10000 lisaa.tran2501@gmail.com

# Variables
input_file=$1   # Input fasta file, contigs 
primer_file=$2  # Primer file
chunk_size=$3   # Number of sequences per chunk
email_address=$4 # Email address of person who run this script


# Fixed variables
chunk_prefix="chunk"
output_dir="chunks"  # Directory to store chunks
result_dir="results"  # Directory to store results

# Create output and result directories if they don't exist
mkdir -p $output_dir
mkdir -p $result_dir

# Split the input fasta file into chunks
awk -v size="$chunk_size" -v prefix="$chunk_prefix" -v dir="$output_dir" 'BEGIN {n_seq=0;} /^>/ {if(n_seq%size==0){file=sprintf("%s/%s_%d.fasta", dir, prefix, n_seq/size);} print >> file; n_seq++; next;} { print >> file; }' < "$input_file"

# Determine number of chunks
num_chunks=$(ls $output_dir/$chunk_prefix*.fasta | wc -l)

# Save job ids to a list for the --dependency option
job_ids=""

# Submit a job for each chunk
for i in $(seq 0 $((num_chunks-1)))
do
    # Variables
    chunk_file="$output_dir/${chunk_prefix}_${i}.fasta"   # Chunk fasta file

    # Run the job and save the job id
    job_id=$(sbatch --parsable --account=def-nricker \
    --time=0-08:00 \
    --nodes=1 \
    --ntasks-per-node=1 \
    --mem=1G \
    --job-name=extract_product_${i} \
    --output=$result_dir/extract_product_${i}_%j.out \
    --wrap="source ~/my_venv/bin/activate && python3 process1.py \"$chunk_file\" \"$primer_file\"")

    # Append this job_id to job_ids with ':' as separator
    job_ids+="${job_id}:"
done

# Remove trailing ':'
job_ids=${job_ids%:}

# Run second script after all chunks have been processed
sbatch --account=def-nricker \
--dependency=afterok:$job_ids \
--time=0-02:00 \
--nodes=1 \
--ntasks-per-node=1 \
--mem=1G \
--job-name=process2 \
--output=$result_dir/process_final_results_%j.out \
--wrap="source ~/my_venv/bin/activate && python3 process2.py \"$primer_file\" raw_results.csv"

# Run third script after the second one is completed 
sbatch --account=def-nricker \
--dependency=afterok:$job_id_process2 \
--time=0-01:00 \ 
--nodes=1 \
--ntasks-per-node=1 \
--mem=1G \
--job-name=process3 \
--output=. 
--wrap="source ~/my_venv/bin/activate && python3 process3.py \"$email_address\" \"$primer_file\" final_result.csv \"$input_file\"

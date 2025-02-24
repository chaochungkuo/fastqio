from fastqio import FASTQReader

# Initialize reader (automatically handles .gz files)
fastq = FASTQReader("/Users/ckuo/Downloads/Pool1_cell_surface_protein_Merged_R2_001.fastq.gz", thread=4,
                    chunk_size=1000000)

# Count reads
# total_reads = fastq.count_reads()
# print(total_reads)

htos = fastq.extract(start=0, end=16, save_parquet=False, parquet_prefix="hto")
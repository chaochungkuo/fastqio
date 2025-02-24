# fastqio

**fastqio** is a FASTQ parser for Python designed for high performance and low memory usage. It loads FASTQ files (plain or gzipped) in chunks using multi-threaded I/O and uses Cython to accelerate core string‐processing functions (such as trimming and quality filtering).

## Features

- **Multi-threaded Chunk Loading:**  
  Uses Python’s `ThreadPoolExecutor` to load file chunks concurrently and overcome I/O limitations.

- **Cython-accelerated String Handling:**  
  Core functions for trimming and quality filtering are implemented in Cython for faster string processing.

- **Lazy Iteration:**  
  Provides a lazy iterator for FASTQ records so that the entire file is not loaded into memory at once.

- **Bulk Operations:**  
  Supports operations like counting reads, trimming sequences, filtering by quality, and sequence extraction (with optional Parquet export via PyArrow).

## Installation

Ensure you have Cython and NumPy installed. Then, install the package from the repository root:

```bash
pip install .
```

This will build the Cython extension and install the Python package.

## Usage Example

```python
from fastqio import FASTQReader

# Initialize reader (automatically handles .gz files)
fastq = FASTQReader("sample.fastq.gz", thread=4, chunk_size=1000000)

# Count reads
total_reads = fastq.count_reads()

# Iterate lazily over records
for read in fastq:
    print(read.info, read.seq, read.quality)

# Trim sequences (remove 5 bases from 5' and 3 bases from 3')
trimmed = fastq.trim(five_prime=5, three_prime=3)

# Filter records by quality threshold (Phred+33)
filtered = fastq.filter_quality(threshold=30)

# Extract a substring from sequences and save as Parquet
fastq.extract(start=0, end=16, save_parquet=True, parquet_prefix="cell_barcode")
```

## Running Tests

Run the tests with pytest from the repository root:

```bash
pytest
```
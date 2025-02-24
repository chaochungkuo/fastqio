import gzip
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
import pyarrow as pa
import pyarrow.parquet as pq

# Import the Cython-accelerated functions.
from .fastq_cython import trim_records_cython, filter_quality_cython

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


class FASTQRecord:
    def __init__(self, info, seq, quality):
        self.info = info
        self.seq = seq
        self.quality = quality

    def __repr__(self):
        return f"FASTQRecord(info={self.info}, seq={self.seq}, quality={self.quality})"


def is_gzipped(file_path):
    return file_path.endswith('.gz')


def open_file(file_path):
    if is_gzipped(file_path):
        return gzip.open(file_path, 'rt')
    else:
        return open(file_path, 'r')


def load_chunk(file_handle, chunk_size):
    """
    Reads chunk_size * 4 lines (each FASTQ record has 4 lines) from file_handle.
    Returns a list of lines.
    """
    lines = []
    try:
        for _ in range(chunk_size * 4):
            lines.append(next(file_handle).rstrip("\n"))
    except StopIteration:
        pass
    return lines


class FASTQReader:
    def __init__(self, file_path, thread=4, chunk_size=1000000):
        """
        Parameters:
          file_path: Path to a FASTQ or FASTQ.gz file.
          thread: Number of threads to use for concurrent chunk loading.
          chunk_size: Number of records to process per chunk.
        """
        self.file_path = file_path
        self.thread = thread
        self.chunk_size = chunk_size
        self._is_gz = is_gzipped(file_path)
        self._file = None

    def _open_file(self):
        self._file = open_file(self.file_path)

    def _reset_file(self):
        if self._file:
            self._file.close()
        self._open_file()

    def __iter__(self):
        self._reset_file()
        return self._record_generator()

    def _record_generator(self):
        """Yield FASTQRecord objects one at a time (lazy iteration)."""
        while True:
            try:
                info = next(self._file).rstrip("\n")
                seq = next(self._file).rstrip("\n")
                next(self._file)  # skip the plus line
                quality = next(self._file).rstrip("\n")
                yield FASTQRecord(info, seq, quality)
            except StopIteration:
                break

    def count_reads(self):
        """Counts the total number of records by processing the file in chunks using multiple threads."""
        self._reset_file()
        count = 0
        executor = ThreadPoolExecutor(max_workers=self.thread)
        futures = []

        while True:
            lines = load_chunk(self._file, self.chunk_size)
            if not lines:
                break
            futures.append(executor.submit(lambda l: len(l) // 4, lines))

        for future in as_completed(futures):
            count += future.result()

        executor.shutdown()
        self._reset_file()
        return count

    def trim(self, five_prime=0, three_prime=0):
        """
        Trims each record's sequence and quality strings by removing bases from
        the 5' and 3' ends using Cython-accelerated code.
        Returns a list of trimmed FASTQRecord objects.
        """
        self._reset_file()
        trimmed_records = []
        executor = ThreadPoolExecutor(max_workers=self.thread)
        futures = []
        while True:
            lines = load_chunk(self._file, self.chunk_size)
            if not lines:
                break
            # Group lines into records (4 lines each)
            records = []
            for i in range(0, len(lines), 4):
                if i + 3 < len(lines):
                    records.append((lines[i], lines[i+1], lines[i+3]))
            futures.append(executor.submit(trim_records_cython, records, five_prime, three_prime))
        for future in as_completed(futures):
            result = future.result()
            for rec in result:
                trimmed_records.append(FASTQRecord(*rec))
        executor.shutdown()
        self._reset_file()
        return trimmed_records

    def filter_quality(self, threshold):
        """
        Filters records based on average quality (Phred+33) using Cython-accelerated code.
        Returns a list of FASTQRecord objects that meet or exceed the threshold.
        """
        self._reset_file()
        filtered_records = []
        executor = ThreadPoolExecutor(max_workers=self.thread)
        futures = []
        while True:
            lines = load_chunk(self._file, self.chunk_size)
            if not lines:
                break
            records = []
            for i in range(0, len(lines), 4):
                if i + 3 < len(lines):
                    records.append((lines[i], lines[i+1], lines[i+3]))
            futures.append(executor.submit(filter_quality_cython, records, threshold))
        for future in as_completed(futures):
            result = future.result()
            for rec in result:
                filtered_records.append(FASTQRecord(*rec))
        executor.shutdown()
        self._reset_file()
        return filtered_records

    def extract(self, start, end, save_parquet=False, parquet_prefix="extracted"):
        """
        Extracts a substring from each record's sequence.
        If save_parquet is True, writes the results to a Parquet file.
        Otherwise, returns a list of extracted substrings.
        """
        self._reset_file()
        extracted = []
        chunk_index = 0
        executor = ThreadPoolExecutor(max_workers=self.thread)
        futures = []
        while True:
            lines = load_chunk(self._file, self.chunk_size)
            if not lines:
                break
            records = []
            for i in range(0, len(lines), 4):
                if i + 1 < len(lines):
                    records.append(lines[i+1])  # sequence only
            # For extraction, simple Python slicing is fast.
            futures.append(executor.submit(lambda recs: [seq[start:end] for seq in recs], records))
        print(f"Processed chunks with {self.chunk_size} sequences each",
              flush=True, end="")
        for future in as_completed(futures):
            extracted.extend(future.result())
            print(".", flush=True, end=".")
        print()
            
        executor.shutdown()
        self._reset_file()
        if save_parquet:
            table = pa.Table.from_pydict({"extracted": extracted})
            filename = f"{parquet_prefix}.parquet"
            pq.write_table(table, filename)
            logger.info(f"Saved parquet file: {filename}")
            return None
        else:
            return extracted
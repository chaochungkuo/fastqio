import os
import tempfile
import pytest
from fastqio import FASTQReader

@pytest.fixture
def sample_fastq(tmp_path):
    content = (
        "@read1\n"
        "ACGTACGTACGT\n"
        "+\n"
        "IIIIIIIIIIII\n"
        "@read2\n"
        "TGCAAGCTTGCA\n"
        "+\n"
        "JJJJJJJJJJJJ\n"
    )
    file_path = tmp_path / "sample.fastq"
    file_path.write_text(content)
    return str(file_path)

def test_count_reads(sample_fastq):
    reader = FASTQReader(sample_fastq, thread=4, chunk_size=1)
    count = reader.count_reads()
    assert count == 2

def test_iteration(sample_fastq):
    reader = FASTQReader(sample_fastq, thread=4, chunk_size=1)
    records = list(reader)
    assert len(records) == 2
    assert records[0].info.startswith("@read1")
    assert records[1].info.startswith("@read2")

def test_trim(sample_fastq):
    reader = FASTQReader(sample_fastq, thread=4, chunk_size=1)
    trimmed = reader.trim(five_prime=2, three_prime=2)
    # For read1: "ACGTACGTACGT" trimmed -> "GTACGTAC"
    # For read2: "TGCAAGCTTGCA" trimmed -> "CAAGCTTG"
    assert trimmed[0].seq == "GTACGTAC"
    assert trimmed[1].seq == "CAAGCTTG"

def test_filter_quality():
    content = (
        "@read1\n"
        "ACGTACGTACGT\n"
        "+\n"
        "!!!!!!!!!!!!\n"  # very low quality (all scores 0)
        "@read2\n"
        "TGCAAGCTTGCA\n"
        "+\n"
        "JJJJJJJJJJJJ\n"  # high quality
    )
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp:
        tmp.write(content)
        tmp_path = tmp.name
    reader = FASTQReader(tmp_path, thread=4, chunk_size=1)
    filtered = reader.filter_quality(30)
    # Only the high-quality record should remain.
    assert len(filtered) == 1
    assert filtered[0].info.startswith("@read2")
    os.unlink(tmp_path)

def test_extract(tmp_path, sample_fastq):
    reader = FASTQReader(sample_fastq, thread=4, chunk_size=1)
    extracted = reader.extract(start=2, end=6, save_parquet=False)
    # For read1: "ACGTACGTACGT" -> substring "GTAC"
    # For read2: "TGCAAGCTTGCA" -> substring "CAAG"
    assert extracted[0] == "GTAC"
    assert extracted[1] == "CAAG"
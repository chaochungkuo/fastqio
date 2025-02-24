# cython: boundscheck=False, wraparound=False, language_level=3

cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)
def trim_records_cython(list records, int five_prime, int three_prime):
    """
    Cython-accelerated function to trim FASTQ records.
    Each record is a tuple: (info, seq, quality)
    Returns a list of tuples with trimmed (info, seq, quality).
    """
    cdef int n = len(records)
    cdef list result = []
    cdef int len_seq, start, end, i
    cdef tuple rec
    cdef str info, seq, quality
    for i in range(n):
        rec = records[i]
        info = rec[0]
        seq = rec[1]
        quality = rec[2]
        len_seq = len(seq)
        start = five_prime if five_prime < len_seq else len_seq
        end = len_seq - three_prime if (len_seq - three_prime) > start else start
        result.append((info, seq[start:end], quality[start:end]))
    return result

@cython.boundscheck(False)
@cython.wraparound(False)
def filter_quality_cython(list records, int threshold):
    """
    Cython-accelerated function to filter FASTQ records based on average quality (Phred+33).
    Each record is a tuple: (info, seq, quality).
    Returns a list of tuples that meet or exceed the quality threshold.
    """
    cdef int n = len(records)
    cdef list result = []
    cdef int i, j, total, q, length
    cdef float avg
    cdef tuple rec
    cdef str info, seq, quality
    for i in range(n):
        rec = records[i]
        info = rec[0]
        seq = rec[1]
        quality = rec[2]
        length = len(quality)
        if length == 0:
            continue
        total = 0
        for j in range(length):
            # Convert character to quality score assuming Phred+33.
            q = ord(quality[j]) - 33
            total += q
        avg = total / length
        if avg >= threshold:
            result.append((info, seq, quality))
    return result
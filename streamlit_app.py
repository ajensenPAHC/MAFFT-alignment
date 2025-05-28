import streamlit as st
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import substitution_matrices
from Bio import Align
from Bio import AlignIO
import tempfile
import requests
import time
from io import StringIO
from matplotlib import cm
import re
import os

st.set_page_config(page_title="Amino Acid Sequence Analyzer", layout="wide")
st.title("🥜 Amino Acid Sequence Analyzer and Classifier")

aligner = Align.PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.mode = 'global'
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5

valid_aa_chars = set("ACDEFGHIKLMNPQRSTVWY-")

def clean_sequence(seq):
    seq = seq.upper()
    return ''.join([aa for aa in seq if aa in valid_aa_chars])

def sanitize_id(base_id, seen):
    safe_id = re.sub(r'[^A-Za-z0-9_]', '_', base_id)[:30]
    orig_id = safe_id
    counter = 1
    while safe_id in seen:
        safe_id = f"{orig_id}_{counter}"
        counter += 1
    seen.add(safe_id)
    return safe_id

def compute_gapped_identity(seq1, seq2):
    pairs = [(a, b) for a, b in zip(seq1, seq2) if a != '-' and b != '-']
    if not pairs:
        return 0.0
    matches = sum(a == b for a, b in pairs)
    return round((matches / len(pairs)) * 100, 2)

def compute_jalview_identity(seq1, seq2):
    matches = 0
    aligned = 0
    for a, b in zip(seq1, seq2):
        if a == '-' and b == '-':
            continue
        aligned += 1
        if a == b:
            matches += 1
    return round((matches / aligned) * 100, 2) if aligned > 0 else 0.0

def clustalo_pairwise_alignment(ref_seq, test_seq):
    fasta_pair = f">ref\n{ref_seq}\n>query\n{test_seq}\n"
    url = "https://www.ebi.ac.uk/Tools/services/rest/clustalo/run"
    payload = {
        'email': 'anonymous@example.com',
        'sequence': fasta_pair,
        'stype': 'protein',
        'outfmt': 'clustal'
    }
    response = requests.post(url, data=payload)
    if not response.ok:
        return None
    job_id = response.text.strip()

    status_url = url.replace("/run", f"/status/{job_id}")
    result_url = url.replace("/run", f"/result/{job_id}/aln-clustal")

    for _ in range(30):
        status = requests.get(status_url).text.strip()
        if status == "FINISHED":
            break
        time.sleep(2)

    alignment_text = requests.get(result_url).text.strip()
    if not alignment_text.startswith("CLUSTAL"):
        return None

    alignments = AlignIO.read(StringIO(alignment_text), "clustal")
    seqs = {rec.id: str(rec.seq) for rec in alignments}
    if 'ref' in seqs and 'query' in seqs:
        return compute_jalview_identity(seqs['ref'], seqs['query'])
    return None

# (the rest of your code remains unchanged except during the FASTA record creation):

        seen_ids = set()
        records = []
        for name, seq in zip(names, sequences):
            safe_id = sanitize_id(name, seen_ids)
            records.append(SeqRecord(Seq(seq), id=safe_id, description=""))

        id_map = {rec.id: name for rec, name in zip(records, names)}

# (rest of code continues using `records` and `id_map` as normal)

import streamlit as st
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import substitution_matrices
from Bio import Align
import tempfile
import requests
import time
from io import StringIO
from matplotlib import cm
import re

st.set_page_config(page_title="Amino Acid Sequence Analyzer", layout="wide")
st.title("🧬 Amino Acid Sequence Analyzer and Classifier")

aligner = Align.PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.mode = 'global'
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5

valid_aa_chars = set("ACDEFGHIKLMNPQRSTVWY-")

def clean_sequence(seq):
    seq = seq.upper()
    return ''.join([aa for aa in seq if aa in valid_aa_chars])

def compute_jalview_identity(seq1, seq2):
    matches = 0
    aligned = 0
    for a, b in zip(seq1, seq2):
        if a == '-' and b == '-':
            continue
        aligned += 1
        if a == b:
            matches += 1
    return round((matches / aligned) * 100, 2) if aligned > 0 else 0.0, matches, aligned

def compute_identity(seq1, seq2):
    matches = sum(a == b for a, b in zip(seq1, seq2))
    aligned_length = len(seq1)
    return round((matches / aligned_length) * 100, 2)

def compute_gapped_identity(seq1, seq2):
    pairs = [(a, b) for a, b in zip(seq1, seq2) if a != '-' and b != '-']
    if not pairs:
        return 0.0
    matches = sum(a == b for a, b in pairs)
    return round((matches / len(pairs)) * 100, 2)

def compute_gap_penalty_identity(seq1, seq2):
    try:
        aligned_score = aligner.align(seq1.replace('-', ''), seq2.replace('-', ''))[0].score
        score_ratio = round(aligned_score / len(seq1.replace('-', '')), 2)
        return score_ratio
    except Exception as e:
        print(f"[compute_gap_penalty_identity] Alignment error: {e}")
        return 0.0

def pairwise_align_and_identity(ref_seq, test_seq):
    try:
        test_record = SeqRecord(Seq(test_seq), id="test", description="")
        ref_record = SeqRecord(Seq(ref_seq), id="ref", description="")

        with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".fasta") as f:
            SeqIO.write([ref_record, test_record], f.name, "fasta")
            fasta_str = open(f.name).read()

        submit_url = "https://www.ebi.ac.uk/Tools/services/rest/clustalo/run"
        response = requests.post(submit_url, data={
            'sequence': fasta_str,
            'stype': 'protein',
            'email': 'anonymous@example.com'
        })

        if response.status_code != 200:
            return 0.0

        job_id = response.text.strip()
        status_url = f"https://www.ebi.ac.uk/Tools/services/rest/clustalo/status/{job_id}"
        result_url = f"https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/{job_id}/aln-clustal_num"

        for _ in range(40):
            status = requests.get(status_url).text.strip()
            if status == "FINISHED":
                break
            time.sleep(3)

        alignment = requests.get(result_url).text
        clustal_io = StringIO(alignment)
        aligned = list(SeqIO.parse(clustal_io, "clustal"))
        ref_aligned = str(aligned[0].seq)
        test_aligned = str(aligned[1].seq)

        identity, _, _ = compute_jalview_identity(ref_aligned, test_aligned)
        return identity

    except Exception as e:
        print(f"[pairwise_align_and_identity] Error: {e}")
        return 0.0

def map_ref_positions(seq):
    mapping = {}
    pos = 0
    for idx, aa in enumerate(seq):
        if aa != '-':
            pos += 1
            mapping[pos] = idx
    return mapping

def color_identity(val):
    try:
        val = float(val)
        norm_val = val / 100
        rgba = cm.Blues(norm_val)
        bg_color = f"rgba({int(255*rgba[0])},{int(255*rgba[1])},{int(255*rgba[2])}, {rgba[3]})"
        text_color = "#FFF" if val >= 99.9 else "#000"
        return f"background-color: {bg_color}; color: {text_color}"
    except Exception as e:
        print(f"[color_identity] Color mapping error: {e}")
        return ""

show_pairwise = st.checkbox("⚠️ Include individual pairwise alignments (slower, more accurate Jalview-style identity)")
st.caption("Note: This will recompute reference vs. each sequence separately, which improves accuracy but takes longer.")

# Existing app logic continues unchanged after alignment dictionary creation
# Add this line in the loop after computing other identity scores:
# row["Pairwise Identity %"] = pairwise_align_and_identity(ref_seq.seq, test_seq) if show_pairwise else "-"

# Add this column to the styled DataFrame visualization subset for highlighting.
# Example: subset=["MSA Identity %", "Gapped Identity %", "Jalview Identity %", "Pairwise Identity %", "Alignment Score / Len"]

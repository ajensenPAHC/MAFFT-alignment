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
st.title("🦜 Amino Acid Sequence Analyzer and Classifier")

aligner = Align.PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.mode = 'global'
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5

valid_aa_chars = set("ACDEFGHIKLMNPQRSTVWY-")

def clean_sequence(seq):
    seq = seq.upper()
    return ''.join([aa for aa in seq if aa in valid_aa_chars])

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
        text_color = "#FFF" if val >= 70 else "#000"
        return f"background-color: {bg_color}; color: {text_color}"
    except:
        return ""

def parse_rows_input(input_str):
    rows = set()
    for part in input_str.split(','):
        part = part.strip()
        if not part:
            continue
        if '-' in part:
            try:
                start, end = map(int, part.split('-'))
                rows.update(range(start, end + 1))
            except:
                continue
        else:
            try:
                rows.add(int(part))
            except:
                continue
    return sorted(rows)

uploaded_file = st.file_uploader("Upload Excel file", type=[".xlsx"])
ref_fasta = st.file_uploader("(Optional) Upload reference sequence (FASTA format)", type=[".fasta"])

if uploaded_file:
    excel = pd.ExcelFile(uploaded_file)
    sheet_name = st.selectbox("Select sheet", excel.sheet_names)
    df = excel.parse(sheet_name)
    df.index += 2

    st.subheader("📋 Excel Preview")
    st.dataframe(df.head(10))

    name_cols = st.multiselect("Select columns to create sequence names", df.columns)
    seq_col = st.selectbox("Select column for amino acid sequence", df.columns)

    selection_type = st.radio("How do you want to select rows?", ["Range", "Specific Rows"])
    if selection_type == "Range":
        row_range = st.slider("Select row range (inclusive)", min_value=2, max_value=int(df.index.max()), value=(2, 5))
        selected_rows = list(range(row_range[0], row_range[1] + 1))
    else:
        selected_rows_input = st.text_input("Enter specific rows or ranges (e.g., 2,4-6,8)", "2,3,4")
        selected_rows = parse_rows_input(selected_rows_input)

    ref_row = None
    use_uploaded_ref = st.checkbox("Use uploaded FASTA as reference instead of selecting from Excel")
    if not use_uploaded_ref:
        ref_row = st.number_input("Enter the row number of the reference sequence (≥2)", min_value=2, step=1)

    aa_pos_input = st.text_input("Enter amino acid positions or ranges (e.g. 5,10-12)", "5,10-12")

    def parse_positions(pos_string):
        pos_set = set()
        for part in pos_string.split(','):
            if '-' in part:
                start, end = part.split('-')
                pos_set.update(range(int(start), int(end)+1))
            else:
                pos_set.add(int(part))
        return sorted(pos_set)

    compute_individual_alignments = st.checkbox("⚠️ Include individual pairwise alignments for each sequence against the reference (slower but more accurate)")

    if st.button("Submit Sequences for Alignment"):
        st.info("Alignment in progress...")

        aa_positions = parse_positions(aa_pos_input)
        names = []
        sequences = []
        invalid_rows = []

        for i in selected_rows:
            raw_value = df.at[i, seq_col] if pd.notna(df.at[i, seq_col]) else ""
            cleaned_seq = clean_sequence(str(raw_value))
            name = "_".join(str(df.at[i, col]) for col in name_cols)
            if cleaned_seq:
                names.append(name[:30].replace(" ", "_"))  # truncate and sanitize ID
                sequences.append(cleaned_seq)
            else:
                invalid_rows.append(i)

        if invalid_rows:
            st.error(f"The following rows have empty or invalid sequences and will be skipped: {invalid_rows}")

        records = [SeqRecord(Seq(seq), id=name, description="") for name, seq in zip(names, sequences)]

        if not records:
            st.error("❌ No valid sequences to process. Please check your input data.")
            st.stop()

        if use_uploaded_ref:
            ref_record = list(SeqIO.parse(ref_fasta, "fasta"))[0]
        else:
            if ref_row not in selected_rows or ref_row in invalid_rows:
                st.error("❌ Reference row is invalid or was skipped. Please choose a valid reference row with a sequence.")
                st.stop()
            ref_idx = selected_rows.index(ref_row)
            ref_record = records[ref_idx]

        # (The rest of the app logic continues unchanged...)

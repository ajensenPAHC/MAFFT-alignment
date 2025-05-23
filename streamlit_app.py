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

def compute_pairwise_identity(ref_seq, test_seq):
    alignment = aligner.align(ref_seq, test_seq)[0]
    aligned_ref = alignment.aligned[0]
    aligned_test = alignment.aligned[1]
    matches = sum(ref_seq[start1:end1] == test_seq[start2:end2] for (start1, end1), (start2, end2) in zip(aligned_ref, aligned_test))
    length = sum(end1 - start1 for (start1, end1) in aligned_ref)
    return round((matches / length) * 100, 2) if length else 0.0

def color_identity(val):
    try:
        val = float(val)
        norm_val = val / 100
        rgba = cm.Blues(norm_val)
        bg_color = f"rgba({int(255*rgba[0])},{int(255*rgba[1])},{int(255*rgba[2])}, {rgba[3]})"
        text_color = "#FFF" if val >= 99.9 else "#000"
        return f"background-color: {bg_color}; color: {text_color}"
    except:
        return ""

uploaded_file = st.file_uploader("Upload Excel file", type=[".xlsx"])
ref_fasta = st.file_uploader("(Optional) Upload reference sequence (FASTA format)", type=[".fasta"])

if uploaded_file:
    excel = pd.ExcelFile(uploaded_file)
    sheet_name = st.selectbox("Select sheet", excel.sheet_names)
    df = excel.parse(sheet_name)
    df = df.astype(str)
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
        selected_rows_input = st.text_input("Enter specific rows (comma-separated)", "2,3,4")
        selected_rows = [int(x.strip()) for x in selected_rows_input.split(",") if x.strip().isdigit()]

    ref_row = None
    use_uploaded_ref = st.checkbox("Use uploaded FASTA as reference instead of selecting from Excel")
    if not use_uploaded_ref:
        ref_row = st.number_input("Enter the row number of the reference sequence (≥2)", min_value=2, step=1)

    aa_pos_input = st.text_input("Enter amino acid positions or ranges (e.g. 5,10-12)", "5,10-12")

    compute_individual_alignments = st.checkbox("⚠️ Include individual pairwise alignments for each sequence against the reference (slower but more accurate)")

    engine = st.selectbox("Select alignment engine", ["Clustal Omega", "MAFFT"])

    if st.button("Submit Sequences for Alignment"):
        st.info("Alignment in progress...")

        names = ["_".join(str(df.loc[i][col]) for col in name_cols) for i in selected_rows]
        sequences = [clean_sequence(str(df.loc[i][seq_col])) for i in selected_rows]
        records = [SeqRecord(Seq(seq), id=name, description="") for name, seq in zip(names, sequences)]

        if use_uploaded_ref:
            ref_record = list(SeqIO.parse(ref_fasta, "fasta"))[0]
        else:
            ref_idx = selected_rows.index(ref_row)
            ref_record = records[ref_idx]

        all_records = [ref_record] + [r for r in records if r.id != ref_record.id]

        with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".fasta") as fasta_file:
            SeqIO.write(all_records, fasta_file, "fasta")
            fasta_path = fasta_file.name

        with open(fasta_path, 'r') as preview:
            st.subheader("🧾 Preview of FASTA File Sent to Alignment Server")
            st.code(preview.read())

        with open(fasta_path, 'r') as f:
            seq_data = f.read()

        if engine == "Clustal Omega":
            url = "https://www.ebi.ac.uk/Tools/services/rest/clustalo/run"
            result_url = "https://www.ebi.ac.uk/Tools/services/rest/clustalo/result"
        else:
            url = "https://www.ebi.ac.uk/Tools/services/rest/mafft/run"
            result_url = "https://www.ebi.ac.uk/Tools/services/rest/mafft/result"

        payload = {'email': 'anonymous@example.com', 'sequence': seq_data}
        response = requests.post(url, data=payload)
        if not response.ok:
            st.error(f"❌ {engine} submission failed.")
            st.stop()
        job_id = response.text.strip()

        status = "RUNNING"
        while status in ["RUNNING", "PENDING"]:
            time.sleep(2)
            status = requests.get(f"{url.replace('/run', '')}/status/{job_id}").text.strip()

        if status != "FINISHED":
            st.error(f"❌ {engine} job failed with status: {status}")
            st.stop()

        result = requests.get(f"{result_url}/{job_id}/aln-clustal")
        aln_text = result.text

        if not aln_text.strip().startswith("CLUSTAL"):
            st.error(f"❌ {engine} result is not a valid Clustal alignment.")
            st.text_area("Raw Output", aln_text)
            st.stop()

        st.code(aln_text, language='clustal')
        alignment = AlignIO.read(StringIO(aln_text), "clustal")

        ref_aligned_seq = str([r.seq for r in alignment if r.id == ref_record.id][0])

        data = {
            "Name": [],
            "MSA Identity %": [],
            "Gapped Identity %": [],
            "Gap-Penalty Identity": [],
            "Jalview Identity %": [],
            "Individual Alignment %": [] if compute_individual_alignments else None
        }

        for record in alignment:
            if record.id == ref_record.id:
                continue
            msa_id = compute_identity(ref_aligned_seq, str(record.seq))
            gapped_id = compute_gapped_identity(ref_aligned_seq, str(record.seq))
            gap_penalty_id = compute_gap_penalty_identity(ref_aligned_seq, str(record.seq))
            jalview_id, _, _ = compute_jalview_identity(ref_aligned_seq, str(record.seq))
            ind_align_id = compute_pairwise_identity(str(ref_record.seq), str(record.seq)) if compute_individual_alignments else None

            data["Name"].append(record.id)
            data["MSA Identity %"].append(msa_id)
            data["Gapped Identity %"].append(gapped_id)
            data["Gap-Penalty Identity"].append(gap_penalty_id)
            data["Jalview Identity %"].append(jalview_id)
            if compute_individual_alignments:
                data["Individual Alignment %"].append(ind_align_id)

        df_results = pd.DataFrame(data)
        st.dataframe(df_results.style.applymap(color_identity, subset=[col for col in df_results.columns if "%" in col]))

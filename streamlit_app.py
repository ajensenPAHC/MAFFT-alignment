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
from matplotlib import pyplot as plt
import base64
import re

st.set_page_config(page_title="Amino Acid Sequence Analyzer", layout="wide")
st.title("🧬 Amino Acid Sequence Analyzer and Classifier")

st.markdown("""
### 🔬 Application Flow
1. **Upload Excel File**: Provide sequence and sample metadata.
2. **Select Name & Sequence Columns**: Choose how each sequence is labeled and where the amino acid data is.
3. **Choose Sequences**: Either by row range or specific row numbers.
4. **Input Amino Acid Positions**: These are alignment-based positions used for comparison against the reference.
5. **Provide Reference Sequence**: Selected from Excel or uploaded as FASTA.
6. **Alignment via Clustal Omega Web API**: Protein sequences aligned with default Clustal Omega settings.
7. **Identity Calculations**:
   - **MSA Identity %**: Jalview-style match rate using ClustalO alignment.
   - **Pairwise Identity Metrics**: Custom scoring comparisons between reference and test sequences.
8. **Amino Acid Comparison at Specific Positions**: Shows differences at positions input by the user.
9. **Downloadable Outputs**: Alignment image, CLUSTAL file, identity comparison CSV.
""")

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
            continue  # skip gap-gap columns
        aligned += 1
        if a == b:
            matches += 1
    return round((matches / aligned) * 100, 2) if aligned > 0 else 0.0, matches, aligned

uploaded_file = st.file_uploader("Upload Excel file", type=[".xlsx"])
ref_fasta = st.file_uploader("(Optional) Upload reference sequence (FASTA format)", type=[".fasta"])

if uploaded_file:
    excel = pd.ExcelFile(uploaded_file)
    sheet_name = st.selectbox("Select sheet", excel.sheet_names)
    df = excel.parse(sheet_name)
    df.index += 2

    st.write("Preview of selected sheet:")
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

    if st.button("Submit Sequences for Alignment"):
        def parse_positions(pos_string):
            pos_set = set()
            for part in pos_string.split(','):
                if '-' in part:
                    start, end = part.split('-')
                    pos_set.update(range(int(start), int(end)+1))
                else:
                    pos_set.add(int(part))
            return sorted(pos_set)

        aa_positions = parse_positions(aa_pos_input)

        records = []
        ref_seq = None

        for idx in selected_rows:
            row = df.loc[idx]
            name = "_".join([str(row[col]) for col in name_cols])
            raw_sequence = str(row[seq_col]).replace("\n", "").strip()
            cleaned_seq = clean_sequence(raw_sequence)
            if not cleaned_seq:
                st.warning(f"Sequence at row {idx} for '{name}' was empty or invalid and will be skipped.")
                continue
            record = SeqRecord(Seq(cleaned_seq), id=name, description="")
            if not use_uploaded_ref and idx == ref_row:
                ref_seq = record
            else:
                records.append(record)

        if use_uploaded_ref and ref_fasta:
            try:
                ref_seq = list(SeqIO.parse(ref_fasta, "fasta"))[0]
                ref_seq.seq = Seq(clean_sequence(str(ref_seq.seq)))
            except Exception as e:
                st.error(f"Error reading uploaded FASTA: {e}")
                st.stop()

        if not ref_seq:
            st.error("Reference sequence not found. Please verify selection or upload.")
            st.stop()

        all_records = [ref_seq] + records

        with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".fasta") as fasta_file:
            SeqIO.write(all_records, fasta_file.name, "fasta")
            fasta_path = fasta_file.name

        with open(fasta_path, 'r') as f:
            st.text_area("Generated FASTA Preview", f.read(), height=200)

        st.info("Submitting alignment job to Clustal Omega Web API...")
        with open(fasta_path, 'r') as f:
            fasta_str = f.read()

        submit_url = "https://www.ebi.ac.uk/Tools/services/rest/clustalo/run"
        result_url = "https://www.ebi.ac.uk/Tools/services/rest/clustalo/result"

        response = requests.post(submit_url, data={
            'sequence': fasta_str,
            'stype': 'protein',
            'email': 'anonymous@example.com'
        })

        if response.status_code != 200:
            st.error("Failed to submit job to Clustal Omega. Please try again later.")
            st.stop()

        job_id = response.text.strip()
        status_url = f"https://www.ebi.ac.uk/Tools/services/rest/clustalo/status/{job_id}"

        with st.spinner("Waiting for Clustal Omega to complete alignment..."):
            for _ in range(40):
                status = requests.get(status_url).text.strip()
                if status == "FINISHED":
                    break
                time.sleep(3)

        alignment = requests.get(f"{result_url}/{job_id}/aln-clustal_num").text
        st.subheader("Alignment Preview")
        st.code(alignment, language="text")

        clustal_io = StringIO(alignment)
        alignments = list(SeqIO.parse(clustal_io, "clustal"))

        ref_aligned_seq = next((str(rec.seq) for rec in alignments if rec.id == ref_seq.id), None)
        if not ref_aligned_seq:
            st.error("Reference sequence not found in alignment.")
            st.stop()

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
    print("Debug: Def line reached")
    try:
        print("Debug: Entered compute_gap_penalty_identity")
        print("Debug: Attempting alignment")
        aligned_score = aligner.align(seq1.replace('-', ''), seq2.replace('-', ''))[0].score
        print(f"Debug: Alignment score: {aligned_score}")
        score_ratio = round(aligned_score / len(seq1.replace('-', '')), 2)
        print(f"Debug: Score ratio: {score_ratio}")
        return score_ratio
    except Exception as e:
        print(f"Alignment error: {e}")
        return 0.0

        alignment_dict = {rec.id: str(rec.seq) for rec in alignments}

        def map_ref_positions(seq):
            mapping = {}
            pos = 0
            for idx, aa in enumerate(seq):
                if aa != '-':
                    pos += 1
                    mapping[pos] = idx
            return mapping

        ref_map = map_ref_positions(ref_aligned_seq)

        full_results = []

        for rec in alignments:
            test_seq = str(rec.seq)
            msa_identity = compute_identity(ref_aligned_seq, test_seq)
            gapped_identity = compute_gapped_identity(ref_aligned_seq, test_seq)
            alignment_score = compute_gap_penalty_identity(ref_aligned_seq, test_seq)
            eff_identity, match_count, eff_len = compute_jalview_identity(ref_aligned_seq, test_seq)

            row = {
                "ID": rec.id,
                "MSA Identity %": msa_identity,
                "Gapped Identity %": gapped_identity,
                "Alignment Score / Len": alignment_score,
                "Jalview Identity %": eff_identity,
                "Effective Matches": match_count,
                "Effective Aligned Positions": eff_len,
            }

            for pos in aa_positions:
                align_idx = ref_map.get(pos)
                label = f"Pos {pos} (Aligned:{ref_aligned_seq[align_idx] if align_idx is not None else '-'}@{(align_idx+1) if align_idx is not None else 'N/A'})"
                row[label] = test_seq[align_idx] if align_idx is not None else '[Invalid]'

            full_results.append(row)

        df_results = pd.DataFrame(full_results)
        df_results = df_results.astype(str)

        sort_option = st.selectbox("Sort results by:", ["ID", "MSA Identity %", "Gapped Identity %", "Jalview Identity %", "Alignment Score / Len"])
        df_results = df_results.sort_values(by=sort_option, ascending=(sort_option == "ID"))
        df_results = pd.concat([df_results[df_results['ID'] == ref_seq.id], df_results[df_results['ID'] != ref_seq.id]])

        def color_identity(val):
    try:
        val = float(val)
        norm_val = val / 100
        rgba = cm.Blues(norm_val)
        bg_color = f"rgba({int(255*rgba[0])},{int(255*rgba[1])},{int(255*rgba[2])}, {rgba[3]})"
        text_color = "#000" if val > 85 else "#FFF"
        return f"background-color: {bg_color}; color: {text_color}"
    except Exception as e:
        print(f"Color mapping error: {e}")
        return ""

        styled_df = df_results.style.map(color_identity, subset=["MSA Identity %", "Gapped Identity %", "Jalview Identity %", "Alignment Score / Len"])
        st.dataframe(styled_df, use_container_width=True)

        csv = df_results.to_csv(index=False).encode()
        st.download_button("Download Comparison CSV", csv, "comparison_results.csv", "text/csv")

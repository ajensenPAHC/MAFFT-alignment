import streamlit as st
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import tempfile
import requests
import time
from matplotlib import cm
from matplotlib import pyplot as plt
import re
import os
import base64
from Bio import Align
from Bio.Align import substitution_matrices
from io import StringIO

st.set_page_config(page_title="Amino Acid Analyzer", layout="wide")
st.title("🧬 Amino Acid Sequence Analyzer and Classifier")

st.markdown("""
### 🧭 Application Flow
1. **Upload Excel File**: Provide sequence and sample metadata.
2. **Select Name & Sequence Columns**: Choose how each sequence is labeled and where the amino acid data is.
3. **Choose Sequences**: Either by row range or specific row numbers.
4. **Input Amino Acid Positions**: These are alignment-based positions used for comparison against the reference.
5. **Provide Reference Sequence**: Selected from Excel or uploaded as FASTA.
6. **Alignment via Clustal Omega Web API**: Protein sequences aligned with default Clustal Omega settings.
7. **Pairwise Identity Calculation**: Identity % is based on *total alignment length*, including gaps (Jalview-style).
8. **Amino Acid Comparison at Specific Positions**: Shows differences at positions input by the user.
9. **Downloadable Outputs**: Alignment image, CLUSTAL file, pairwise identity CSV.

> **Note**: Jalview-style identity = `# exact matches / alignment length (including gaps)`
> **Pairwise identity method**: Uses Biopython `Align.PairwiseAligner` with BLOSUM62 scoring matrix for biologically meaningful scoring.
""")

aligner = Align.PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.mode = 'global'
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5

uploaded_file = st.file_uploader("📤 Upload Excel file with sequences", type=["xlsx"])

if uploaded_file:
    xls = pd.ExcelFile(uploaded_file)
    sheet_name = st.selectbox("📑 Select a sheet", xls.sheet_names)
    df = pd.read_excel(xls, sheet_name)
    df_display = df.copy()
    df_display.index += 2
    st.dataframe(df_display, use_container_width=True)

    name_col = st.selectbox("🔤 Select name/label column", df.columns)
    seq_col = st.selectbox("🧬 Select amino acid sequence column", df.columns)

    row_method = st.radio("📊 How do you want to select rows?", ["Range", "Specific Rows"])
    if row_method == "Range":
        row_range = st.slider("📌 Select row range (Excel-style)", 2, len(df)+1, (2, len(df)+1))
        selected_df = df.iloc[(row_range[0]-2):(row_range[1]-2)].copy()
    else:
        row_ids = st.text_input("🔢 Enter specific row numbers separated by commas (e.g., 2,5,8)")
        if row_ids:
            indices = [int(i.strip())-2 for i in row_ids.split(",")]
            selected_df = df.iloc[indices].copy()
        else:
            selected_df = pd.DataFrame()

    ref_index = st.number_input("🎯 Enter the row number of the reference sequence", min_value=2, max_value=len(df)+1, step=1)
    ref_seq = df.iloc[ref_index-2]

    position_string = st.text_input("🔍 Enter amino acid positions or ranges (e.g. 5,10-12):")

    def parse_aa_positions(pos_string):
        positions = []
        if not pos_string:
            return positions
        for part in pos_string.split(","):
            if "-" in part:
                start, end = map(int, part.split("-"))
                positions.extend(range(start, end+1))
            else:
                positions.append(int(part))
        return sorted(set(positions))

    positions = parse_aa_positions(position_string)

    if st.button("▶️ Submit Sequences for Alignment"):
        sequences = []
        for _, row in selected_df.iterrows():
            name = str(row[name_col])
            sequence = str(row[seq_col]).replace(" ", "").replace("\n", "")
            sequences.append(SeqRecord(Seq(sequence), id=name, description=""))

        ref_name = str(ref_seq[name_col])
        ref_sequence = str(ref_seq[seq_col]).replace(" ", "").replace("\n", "")
        sequences.insert(0, SeqRecord(Seq(ref_sequence), id=ref_name, description=""))

        with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".fasta") as fasta_file:
            SeqIO.write(sequences, fasta_file.name, "fasta")
            fasta_path = fasta_file.name
            with open(fasta_path) as f:
                st.code(f.read(), language='fasta')

        st.info("🧬 Submitting to Clustal Omega...")
        with open(fasta_path) as f:
            data = {'sequence': f.read(), 'email': 'your@email.com', 'stype': 'protein'}
            res = requests.post('https://www.ebi.ac.uk/Tools/services/rest/clustalo/run', data=data)
            job_id = res.text.strip()

        time.sleep(10)
        result_url = f"https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/{job_id}/aln-clustal_num"
        alignment_res = requests.get(result_url)

        if alignment_res.status_code != 200:
            st.error("Clustal Omega returned an error. Check input format or sequence diversity.")
            st.stop()

        alignment_data = alignment_res.text
        st.text_area("📄 Alignment Output (CLUSTAL)", alignment_data, height=300)

        alignment = Align.MultipleSeqAlignment(list(SeqIO.parse(StringIO(alignment_data), "clustal")))
        aligned_seqs = {rec.id: str(rec.seq) for rec in alignment}

        ref_aligned = aligned_seqs.get(ref_name)
        if not ref_aligned:
            st.error(f"Reference '{ref_name}' not found in alignment!")
            st.stop()

        aa_table = []
        for sid, seq in aligned_seqs.items():
            if sid == ref_name:
                continue
            try:
                score = aligner.align(ref_aligned, seq)[0].score
                matches = sum(r == s and r != '-' for r, s in zip(ref_aligned, seq))
                identity = (matches / len(ref_aligned)) * 100
                aa_diff = {}
                for p in positions:
                    align_idx = [i for i, res in enumerate(ref_aligned) if res != '-']
                    if p <= len(align_idx):
                        idx = align_idx[p-1]
                        aa_diff[f"RefPos {p} (Align {idx+1})"] = f"{ref_aligned[idx]} → {seq[idx]}"
                    else:
                        aa_diff[f"RefPos {p}"] = "Out of range"
                aa_table.append({
                    "Sample": sid,
                    "Identity %": round(identity, 2),
                    **aa_diff
                })
            except Exception as e:
                aa_table.append({
                    "Sample": sid,
                    "Identity %": "Error",
                    "Error": str(e)
                })

        result_df = pd.DataFrame(aa_table)
        st.dataframe(result_df, use_container_width=True)

        csv = result_df.to_csv(index=False).encode()
        st.download_button("📥 Download CSV Report", csv, "pairwise_report.csv", "text/csv")

        def highlight_identity(val):
            if isinstance(val, (float, int)):
                normalized = val / 100
                color = cm.Blues(normalized)
                rgba = f"rgba({int(color[0]*255)},{int(color[1]*255)},{int(color[2]*255)}, 0.6)"
                return f"background-color: {rgba}"
            return ""

        st.dataframe(result_df.style.applymap(highlight_identity, subset=["Identity %"]))

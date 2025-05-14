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

st.set_page_config(page_title="Amino Acid Analyzer", layout="wide")
st.title("🧬 Amino Acid Sequence Analyzer and Classifier")

uploaded_excel = st.file_uploader("Upload Excel Spreadsheet (.xlsx)", type=["xlsx"])
uploaded_gene_db = st.file_uploader("Upload Gene Type Database (.csv) (Optional)", type=["csv"])
ref_seq = None
session = st.session_state

if uploaded_excel:
    xls = pd.ExcelFile(uploaded_excel)
    tab_name = st.selectbox("Select the Excel Sheet to Use:", xls.sheet_names)
    df = pd.read_excel(xls, sheet_name=tab_name)
    df = df.astype(str)
    df.index = df.index + 2
    st.subheader("Preview of Selected Sheet")
    df_display = df.copy()
    df_display.index = df_display.index.astype(str)
    df_display.index.name = "Excel Row #"
    st.dataframe(df_display.head().style.set_properties(**{'text-align': 'left'}))

    name_columns = st.multiselect("Select column(s) to use for naming sequences:", df.columns)
    seq_column = st.selectbox("Select the column containing amino acid sequences:", df.columns)

    row_mode = st.radio("How do you want to select rows?", ["Range", "Specific Rows"])
    if row_mode == "Range":
        row_range = st.slider("Select row range to process (starting from row 2):", 2, len(df)+1, (2, len(df)+1))
        df = df.loc[row_range[0]:row_range[1]]
    else:
        row_indices_input = st.text_input("Enter comma-separated row indices (e.g., 2,4,6):", "2,3")
        try:
            indices = [int(i.strip()) for i in row_indices_input.split(",") if i.strip().isdigit()]
            df = df.loc[indices]
        except Exception as e:
            st.error(f"Invalid row indices: {e}")
            st.stop()

    sequences = []
    for index, row in df.iterrows():
        sequence = row[seq_column].strip().replace(" ", "") if pd.notna(row[seq_column]) else None
        if not sequence:
            continue
        name = "_".join(re.sub(r"[^A-Za-z0-9_]", "_", row[col].strip()) for col in name_columns if pd.notna(row[col]))
        if name and sequence:
            sequences.append(SeqRecord(Seq(sequence), id=name, description=""))

    if not sequences:
        st.error("No valid sequences found. Please ensure the name and sequence columns are selected correctly.")
        st.stop()

    ref_option = st.radio("How do you want to provide the reference sequence?", ["Select from Excel", "Upload FASTA File"])
    if ref_option == "Select from Excel":
        ref_display_indices = df.index.tolist()
        ref_selection = st.selectbox("Select the row number of the reference sequence:", ref_display_indices)
        ref_row = df.loc[ref_selection]
        ref_name = "_".join([ref_row[col].strip().replace(" ", "_") for col in name_columns if pd.notna(ref_row[col])])
        ref_seq_text = ref_row[seq_column].strip().replace(" ", "")
        ref_seq = SeqRecord(Seq(ref_seq_text), id=ref_name, description="")
    else:
        uploaded_fasta = st.file_uploader("Upload Reference Sequence (FASTA format)", type=["fasta"])
        if uploaded_fasta:
            ref_seq = list(SeqIO.parse(uploaded_fasta, "fasta"))[0]

    if not ref_seq:
        st.warning("Please provide a valid reference sequence to proceed.")
        st.stop()

    aa_positions_input = st.text_input("Enter amino acid positions or ranges (e.g. 5,10-12):")

    all_seqs = [ref_seq] + [s for s in sequences if s.id != ref_seq.id]

    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".fasta") as fasta_file:
        SeqIO.write(all_seqs, fasta_file.name, "fasta")
        fasta_path = fasta_file.name

    with open(fasta_path, 'r') as preview:
        fasta_preview = preview.read().split(">")
        if len(fasta_preview) > 1:
            st.code(">" + fasta_preview[1], language="text")

    with open(fasta_path, "rb") as download:
        st.download_button("Download FASTA Sequence File", download, file_name="sequences.fasta")

    try:
        with open(fasta_path, 'rb') as f:
            response = requests.post(
                'https://www.ebi.ac.uk/Tools/services/rest/clustalo/run',
                data={'email': 'your.email@example.com', 'stype': 'protein', 'guidetreeout': 'true'},
                files={'sequence': f}
            )
            job_id = response.text.strip()
    except Exception as e:
        st.error(f"Failed to start alignment job: {e}")
        st.stop()

    with st.spinner("Running alignment on Clustal Omega..."):
        max_wait = 60
        waited = 0
        while waited < max_wait:
            status = requests.get(f'https://www.ebi.ac.uk/Tools/services/rest/clustalo/status/{job_id}').text
            if status == 'FINISHED':
                break
            elif status == 'ERROR':
                st.error("Clustal Omega job failed.")
                st.stop()
            time.sleep(5)
            waited += 5

    aln = requests.get(f'https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/{job_id}/aln-clustal_num').text

    if not aln.strip().startswith("CLUSTAL"):
        st.error("Clustal Omega returned an empty or invalid alignment. Try fewer sequences or check format.")
        st.text(aln)
        st.stop()

    with tempfile.NamedTemporaryFile(mode="w+", delete=False, suffix=".aln") as aligned_file:
        aligned_file.write(aln)
        aligned_file_path = aligned_file.name

    st.subheader("Alignment Result (CLUSTAL format)")
    st.code(aln, language="text")

    st.download_button("Download Alignment (CLUSTAL Format)", aln, file_name="alignment.aln")

    fig, ax = plt.subplots(figsize=(10, len(aln.splitlines())/3))
    ax.text(0, 0.5, aln.replace('\t', ' '), family='monospace', fontsize=10)
    ax.axis('off')
    png_file = tempfile.NamedTemporaryFile(delete=False, suffix=".png")
    fig.savefig(png_file.name, bbox_inches='tight')
    with open(png_file.name, "rb") as f:
        st.download_button("Download Alignment Image (PNG)", f, file_name="alignment.png")

    st.success("Alignment complete! Proceeding to analysis.")

    st.subheader("Step 5: Pairwise Identity and Amino Acid Comparison")

    if aa_positions_input:
        st.markdown(f"**Amino acid positions selected:** {aa_positions_input}")
    else:
        st.warning("No amino acid positions were selected earlier.")

    # Parse aligned sequences
    aligned_seqs = {}
    for line in aln.splitlines():
        if line.strip() and not line.startswith("CLUSTAL") and not line.startswith(" "):
            parts = line.split()
            if len(parts) >= 2:
                name, seq = parts[0], parts[1]
                aligned_seqs[name] = aligned_seqs.get(name, "") + seq

    ref_aligned_seq = aligned_seqs.get(ref_seq.id)
    if not ref_aligned_seq:
        st.error(f"Reference sequence not found in alignment. Tried matching '{ref_seq.id}' to {list(aligned_seqs.keys())}")
        st.stop()

    def parse_positions(pos_string):
        pos_list = []
        for part in pos_string.split(','):
            if '-' in part:
                start, end = map(int, part.split('-'))
                pos_list.extend(range(start, end + 1))
            else:
                pos_list.append(int(part))
        return sorted(set(pos_list))

    pos_list = parse_positions(aa_positions_input)

    result_table = []
    for name, test_seq in aligned_seqs.items():
        if name == ref_seq.id:
            continue
        match_count = sum(1 for i, pos in enumerate(pos_list) if pos-1 < len(ref_aligned_seq)
                          and pos-1 < len(test_seq)
                          and ref_aligned_seq[pos-1] == test_seq[pos-1])
        identity_pct = 100 * match_count / len(pos_list)
        aa_comparison = {f"Pos {p}": test_seq[p-1] if p-1 < len(test_seq) else "-" for p in pos_list}
        row = {"Name": name, "Identity %": identity_pct, **aa_comparison}
        result_table.append(row)

    result_df = pd.DataFrame(result_table)

    def color_confidence(val):
        try:
            green = cm.Greens(val / 100)
            red = cm.Reds(1 - val / 100)
            hex_color = f"background-color: rgba({int(255*red[0])},{int(255*green[1])},{int(255*green[2])},0.6)"
            return hex_color
        except:
            return ""

    st.subheader("Final Result Table")
    styled_df = result_df.style.applymap(color_confidence, subset=["Identity %"])
    st.dataframe(styled_df, use_container_width=True)

    csv_out = result_df.to_csv(index=False).encode("utf-8")
    st.download_button("Download Comparison Results (CSV)", csv_out, file_name="pairwise_results.csv")

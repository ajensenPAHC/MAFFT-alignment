import streamlit as st
import pandas as pd
import requests
import time
from io import StringIO

MAFFT_API_URL = "https://mafft.cbrc.jp/alignment/server/"

def send_to_mafft_fasta(fasta_text):
    """Submits sequences to MAFFT's online server and retrieves the aligned result."""
    payload = {
        "sequence": fasta_text,
        "alignment_method": "auto",
        "output": "fasta"
    }
    response = requests.post(MAFFT_API_URL, data=payload)

    if response.status_code == 200 and "Your job is registered" in response.text:
        job_id = response.text.split("job = ")[-1].split()[0]
        result_url = f"{MAFFT_API_URL}{job_id}/{job_id}.fasta"
        for _ in range(30):
            time.sleep(2)
            r = requests.get(result_url)
            if r.status_code == 200 and r.text.startswith(">"):
                return r.text
        st.error("Timed out waiting for MAFFT alignment.")
    else:
        st.error("Failed to submit sequences to MAFFT.")
    return None

def format_fasta(df, id_col, seq_col):
    fasta_lines = []
    for _, row in df.iterrows():
        fasta_lines.append(f">{row[id_col]}")
        fasta_lines.append(row[seq_col])
    return "\n".join(fasta_lines)

st.set_page_config(page_title="Web MAFFT Sequence Aligner", layout="centered")

st.title("🧬 MAFFT Web-Based Amino Acid Aligner")
uploaded_file = st.file_uploader("Upload an Excel file with amino acid sequences", type=["xlsx"])

if uploaded_file:
    df = pd.read_excel(uploaded_file, sheet_name=None)
    sheet_names = list(df.keys())
    selected_sheet = st.selectbox("Select sheet", sheet_names)
    sheet_df = df[selected_sheet]

    st.write("### Preview of Selected Sheet (first few rows):")
    st.dataframe(sheet_df.head().reset_index(drop=True).rename(index=lambda x: x + 2))

    id_col = st.selectbox("Select column for sequence ID", sheet_df.columns)
    seq_col = st.selectbox("Select column for sequence", sheet_df.columns)

    if st.button("Align Sequences via MAFFT Web Server"):
        fasta = format_fasta(sheet_df, id_col, seq_col)
        st.text_area("Generated FASTA", fasta, height=200)
        aligned_result = send_to_mafft_fasta(fasta)

        if aligned_result:
            st.success("Alignment successful!")
            st.text_area("Aligned Sequences (FASTA)", aligned_result, height=300)
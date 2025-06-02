import streamlit as st
import pandas as pd
from Bio import AlignIO, SeqIO, Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import requests, tempfile, time, re
from io import StringIO

st.set_page_config(page_title="Amino Acid Analyzer", layout="wide")
st.title("🧬 Amino Acid Sequence Alignment & Identity Tool")

# === Utilities ===
def clean_sequence(seq):
    return ''.join([aa for aa in seq.upper() if aa in "ACDEFGHIKLMNPQRSTVWY-"])

def sanitize_id(base_id, used):
    base = re.sub(r'[^A-Za-z0-9_]', '_', base_id)[:30]
    while base in used:
        base += "_"
    used.add(base)
    return base

def parse_positions(input_str):
    pos = set()
    for part in input_str.split(','):
        if '-' in part:
            a, b = part.split('-')
            pos.update(range(int(a), int(b)+1))
        else:
            pos.add(int(part.strip()))
    return sorted(pos)

def compute_jalview_identity(seq1, seq2):
    match, total = 0, 0
    for a, b in zip(seq1, seq2):
        if a == '-' and b == '-': continue
        total += 1
        if a == b: match += 1
    return round((match / total) * 100, 2) if total > 0 else 0

def strip_clustal_consensus(text):
    return "\n".join(
        line for line in text.splitlines()
        if line.strip() and not re.match(r'^[\s.:*]+$', line)
    ) + "\n"

def clustalo_api_fasta(fasta_str):
    run_url = "https://www.ebi.ac.uk/Tools/services/rest/clustalo/run"
    payload = {"email": "anon@example.com", "sequence": fasta_str, "stype": "protein", "outfmt": "clustal"}
    job_id = requests.post(run_url, data=payload).text.strip()
    status_url = run_url.replace("/run", f"/status/{job_id}")
    result_url = run_url.replace("/run", f"/result/{job_id}/aln-clustal")
    for _ in range(30):
        if requests.get(status_url).text.strip() == "FINISHED": break
        time.sleep(2)
    return requests.get(result_url).text.strip()

# === Upload ===
uploaded_file = st.file_uploader("Upload Excel with sequences", type="xlsx")
uploaded_fasta = st.file_uploader("Optional: Reference FASTA file", type="fasta")

if uploaded_file:
    df = pd.read_excel(uploaded_file, sheet_name=0).astype(str)
    df.index += 2
    st.dataframe(df.head(10))

    name_cols = st.multiselect("Columns to name sequences", df.columns)
    seq_col = st.selectbox("Column with amino acid sequences", df.columns)
    row_input = st.text_input("Row numbers (e.g. 2-5,7,9)", "2-4")
    ref_row = st.number_input("Reference row (optional)", min_value=2, step=1)
    use_uploaded_ref = st.checkbox("Use uploaded FASTA instead")

    selected_rows = []
    for part in row_input.split(','):
        if '-' in part:
            a, b = map(int, part.split('-'))
            selected_rows.extend(range(a, b+1))
        else:
            selected_rows.append(int(part.strip()))

    sequences, ids, used = [], [], set()
    for row in selected_rows:
        if row in df.index and pd.notna(df.at[row, seq_col]):
            raw_seq = clean_sequence(str(df.at[row, seq_col]))
            label = "_".join([str(df.at[row, col]) for col in name_cols])
            sid = sanitize_id(f"{row}_{label}", used)
            sequences.append(SeqRecord(Seq(raw_seq), id=sid, description=""))
            ids.append((sid, row))

    if not sequences:
        st.error("No valid sequences found.")
        st.stop()

    if use_uploaded_ref and uploaded_fasta:
        ref_record = list(SeqIO.parse(uploaded_fasta, "fasta"))[0]
    elif not use_uploaded_ref and ref_row in df.index:
        ref_seq = clean_sequence(str(df.at[ref_row, seq_col]))
        ref_record = SeqRecord(Seq(ref_seq), id="REF", description="")
    else:
        st.error("Reference sequence not found.")
        st.stop()

    all_records = [ref_record] + sequences
    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as tmp:
        SeqIO.write(all_records, tmp, "fasta")
        fasta_str = tmp.read() if tmp.mode == 'r' else open(tmp.name).read()

    st.code(fasta_str, language="text")

    alignment_text = clustalo_api_fasta(fasta_str)
    st.subheader("🔬 Clustal Omega Alignment")
    st.code(alignment_text)

    cleaned = strip_clustal_consensus(alignment_text)
    try:
        parsed = list(AlignIO.parse(StringIO(cleaned), "clustal"))
        alignment = parsed[0]
    except Exception as e:
        st.error("❌ Failed to parse alignment")
        st.exception(e)
        st.stop()

    ref_seq = next((r.seq for r in alignment if r.id == ref_record.id), None)
    if ref_seq is None:
        st.error(f"Reference ID '{ref_record.id}' not found in alignment.")
        st.stop()

    results = {"Name": [], "Identity %": []}
    for rec in alignment:
        if rec.id == ref_record.id:
            continue
        identity = compute_jalview_identity(str(ref_seq), str(rec.seq))
        label = next((f"{i[0]} (Row {i[1]})" for i in ids if i[0] == rec.id), rec.id)
        results["Name"].append(label)
        results["Identity %"].append(identity)

    df_out = pd.DataFrame(results)
    st.dataframe(df_out.style.background_gradient(cmap="Blues", subset=["Identity %"]))
    st.download_button("📥 Download CSV", df_out.to_csv(index=False), file_name="results.csv")

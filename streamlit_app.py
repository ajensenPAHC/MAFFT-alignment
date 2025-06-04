# Streamlit Amino Acid Sequence Analyzer – Fully Merged

import streamlit as st
import pandas as pd
from Bio import AlignIO, SeqIO, Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import substitution_matrices
import requests, tempfile, time, re, os
from io import StringIO
from matplotlib import cm

st.set_page_config(page_title="Amino Acid Analyzer", layout="wide")
st.title("🧬 Amino Acid Sequence Alignment & Identity Tool")

# === Config ===
pairwise_toggle = st.checkbox("Include individual pairwise alignments (Clustal Omega)")
aa_position_filter = st.text_input("Optional: Filter by AA positions (e.g. 5,10-12)", "")

aligner = Align.PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.mode = 'global'
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5

valid_aa_chars = set("ACDEFGHIKLMNPQRSTVWY-")

# === Utilities ===
def clean_sequence(seq):
    return ''.join([aa for aa in seq.upper() if aa in valid_aa_chars])

def sanitize_id(base_id, used):
    base = re.sub(r'[^A-Za-z0-9_]', '_', base_id)[:30]
    count = 1
    while base in used:
        base = f"{base[:27]}_{count}"
        count += 1
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

def clustalo_api_fasta(fasta_str):
    run_url = "https://www.ebi.ac.uk/Tools/services/rest/clustalo/run"
    payload = {"email": "anon@example.com", "sequence": fasta_str, "stype": "protein", "outfmt": "clustal"}
    try:
        job_id = requests.post(run_url, data=payload).text.strip()
        if not job_id:
            raise ValueError("Failed to receive job ID from ClustalO.")
        status_url = run_url.replace("/run", f"/status/{job_id}")
        result_url = run_url.replace("/run", f"/result/{job_id}/aln-clustal")
        for _ in range(30):
            status = requests.get(status_url).text.strip()
            if status == "FINISHED": break
            time.sleep(2)
        return requests.get(result_url).text.strip()
    except Exception as e:
        st.error(f"❌ Clustal Omega failed: {e}")
        return ""

def strip_clustal_consensus(text):
    return "\n".join(
        line for line in text.splitlines()
        if line.strip() and not re.match(r'^[\s.:*]+$', line)
    ) + "\n"

def color_identity(val):
    try:
        val = float(val)
        norm_val = val / 100
        rgba = cm.Blues(norm_val)
        bg_color = f"rgba({int(255*rgba[0])},{int(255*rgba[1])},{int(255*rgba[2])},{rgba[3]})"
        text_color = "#FFF" if val >= 70 else "#000"
        return f"background-color: {bg_color}; color: {text_color}"
    except:
        return ""

# === File Upload ===
uploaded_file = st.file_uploader("Upload Excel with sequences", type="xlsx")
uploaded_fasta = st.file_uploader("Optional: Reference FASTA file", type="fasta")
submit_button = st.button("🚀 Run Alignment")

if uploaded_file and submit_button:
    try:
        df = pd.read_excel(uploaded_file, sheet_name=0).astype(str)
        df.index += 2
        st.dataframe(df.head(10).astype(str))

        name_cols = st.multiselect("Columns to name sequences", df.columns)
        seq_col = st.selectbox("Column with amino acid sequences", df.columns)

        selection_mode = st.radio("Select rows by", ["Slider", "Manual"])
        if selection_mode == "Slider":
            min_row, max_row = int(df.index.min()), int(df.index.max())
            row_start, row_end = st.slider("Select row range", min_value=min_row, max_value=max_row, value=(min_row, min_row + 2))
            selected_rows = list(range(row_start, row_end + 1))
        else:
            row_input = st.text_input("Row numbers (e.g. 2-5,7,9)", "2-4")
            selected_rows = []
            for part in row_input.split(','):
                if '-' in part:
                    a, b = map(int, part.split('-'))
                    selected_rows.extend(range(a, b+1))
                else:
                    selected_rows.append(int(part.strip()))

        ref_row = st.number_input("Reference row (optional)", min_value=2, step=1)
        use_uploaded_ref = st.checkbox("Use uploaded FASTA instead")

        sequences, ids, used = [], [], set()
        for row in selected_rows:
            if row in df.index and pd.notna(df.at[row, seq_col]):
                raw_seq = clean_sequence(str(df.at[row, seq_col]))
                if len(raw_seq) < 5:
                    continue
                label = "_".join([str(df.at[row, col]) for col in name_cols if col in df.columns and pd.notna(df.at[row, col])])
                sid = sanitize_id(f"{row}_{label}", used)
                sequences.append(SeqRecord(Seq(raw_seq), id=sid, description=""))
                ids.append((sid, row))

        if not sequences:
            st.error("No valid sequences found.")
            st.stop()

        if use_uploaded_ref and uploaded_fasta:
            ref_record = list(SeqIO.parse(uploaded_fasta, "fasta"))[0]
        elif not use_uploaded_ref and ref_row in df.index:
            raw_ref = clean_sequence(str(df.at[ref_row, seq_col]))
            ref_record = SeqRecord(Seq(raw_ref), id="REF", description="")
        else:
            st.error("Reference sequence not found.")
            st.stop()

        all_records = [ref_record] + sequences
        with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".fasta") as tmp:
            SeqIO.write(all_records, tmp, "fasta")
            tmp_path = tmp.name
        with open(tmp_path, 'r') as f:
            fasta_str = f.read()
        os.remove(tmp_path)

        st.code(fasta_str, language="text")
        alignment_text = clustalo_api_fasta(fasta_str)
        if not alignment_text.startswith("CLUSTAL"):
            st.error("❌ Invalid alignment output")
            st.text_area("Raw Output", alignment_text)
            st.stop()

        cleaned = strip_clustal_consensus(alignment_text)
        try:
            parsed = list(AlignIO.parse(StringIO(cleaned), "clustal"))
            alignment = parsed[0]
        except Exception as e:
            st.error("❌ Failed to parse alignment")
            st.exception(e)
            st.text_area("Clustal Raw Output", alignment_text)
            st.stop()

        ref_seq = next((r.seq for r in alignment if r.id == ref_record.id), None)
        if ref_seq is None:
            st.error(f"Reference ID '{ref_record.id}' not found in alignment.")
            st.stop()

        results = {"Name": [], "MSA Identity %": []}
        for rec in alignment:
            if rec.id == ref_record.id:
                continue
            aligned_ref = str(ref_seq)
            aligned_seq = str(rec.seq)
            if aa_position_filter:
                try:
                    positions = parse_positions(aa_position_filter)
                    aligned_ref = ''.join([aligned_ref[pos-1] for pos in positions if pos <= len(aligned_ref)])
                    aligned_seq = ''.join([aligned_seq[pos-1] for pos in positions if pos <= len(aligned_seq)])
                except Exception as e:
                    st.warning(f"⚠️ Invalid AA position filter: {e}")
            identity = compute_jalview_identity(aligned_ref, aligned_seq)
            label = next((f"{i[0]} (Row {i[1]})" for i in ids if i[0] == rec.id), rec.id)
            results["Name"].append(label)
            results["MSA Identity %"].append(identity)

        if pairwise_toggle:
            st.subheader("🔁 Running pairwise alignments...")
            for sid, row in ids:
                if sid == ref_record.id:
                    continue
                rec = next((r for r in sequences if r.id == sid), None)
                if rec:
                    pair_fasta = f">REF\n{str(ref_record.seq)}\n>{rec.id}\n{str(rec.seq)}"
                    raw_pair_align = clustalo_api_fasta(pair_fasta)
                    try:
                        cleaned_pair = strip_clustal_consensus(raw_pair_align)
                        pair_parsed = list(AlignIO.parse(StringIO(cleaned_pair), "clustal"))
                        if pair_parsed:
                            pair_seqs = {r.id: str(r.seq) for r in pair_parsed[0]}
                            ref_aln = pair_seqs.get("REF")
                            test_aln = pair_seqs.get(rec.id)
                            if ref_aln and test_aln:
                                if aa_position_filter:
                                    try:
                                        positions = parse_positions(aa_position_filter)
                                        ref_aln = ''.join([ref_aln[p-1] for p in positions if p <= len(ref_aln)])
                                        test_aln = ''.join([test_aln[p-1] for p in positions if p <= len(test_aln)])
                                    except: pass
                                identity = compute_jalview_identity(ref_aln, test_aln)
                                idx = results["Name"].index(f"{rec.id} (Row {row})")
                                results.setdefault("Pairwise %", ["" for _ in range(len(results["Name"]))])
                                results["Pairwise %"][idx] = identity
                    except Exception as e:
                        st.warning(f"⚠️ Failed pairwise for {rec.id}: {e}")
                    time.sleep(1)

        df_out = pd.DataFrame(results).astype(str)
        st.dataframe(df_out.style.map(color_identity, subset=[col for col in df_out.columns if "%" in col]))
        st.download_button("📥 Download CSV", data=df_out.to_csv(index=False).encode('utf-8'), file_name="results.csv", mime="text/csv")

    except Exception as e:
        st.error("💥 Unexpected error during processing")
        st.exception(e)

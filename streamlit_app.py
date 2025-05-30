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

def parse_positions(pos_string):
    pos_set = set()
    for part in pos_string.split(','):
        part = part.strip()
        if not part:
            continue
        try:
            if '-' in part:
                start, end = part.split('-')
                pos_set.update(range(int(start), int(end)+1))
            else:
                pos_set.add(int(part))
        except (ValueError, TypeError):
            continue
    return sorted(pos_set)

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
        try:
            row_range = st.slider("Select row range (inclusive)", min_value=2, max_value=int(df.index.max()), value=(2, 5))
            selected_rows = list(range(row_range[0], row_range[1] + 1))
        except Exception as e:
            st.error(f"Row range selection failed: {e}")
            st.stop()
    
    else:
        selected_rows_input = st.text_input("Enter specific rows or ranges (e.g., 2,4-6,8)", "2,3,4")
        selected_rows = parse_rows_input(selected_rows_input)

    ref_row = None
    use_uploaded_ref = st.checkbox("Use uploaded FASTA as reference instead of selecting from Excel")
    if not use_uploaded_ref:
        ref_row = st.number_input("Enter the row number of the reference sequence (≥2)", min_value=2, step=1)

    aa_pos_input = st.text_input("Enter amino acid positions or ranges (e.g. 5,10-12)", "5,10-12")
    compute_individual_alignments = st.checkbox("⚠️ Include individual pairwise alignments for each sequence against the reference (slower but more accurate)")

    if st.button("Submit Sequences for Alignment"):
        st.info("Alignment in progress...")

        aa_positions = parse_positions(aa_pos_input)
        names, sequences, invalid_rows = [], [], []
        seen_ids = set()
        row_map = {}
        with st.spinner("Cleaning and preparing sequences..."):
            for i in selected_rows:
                raw_value = df.at[i, seq_col] if pd.notna(df.at[i, seq_col]) else ""
                cleaned_seq = clean_sequence(str(raw_value))
                if cleaned_seq:
                    name = "_".join(str(df.at[i, col]) for col in name_cols)
                    safe_id = sanitize_id(f"{i}_{name}", seen_ids)
                    row_map[safe_id] = i
                    names.append(safe_id)
                    sequences.append(cleaned_seq)
                else:
                    invalid_rows.append(i)

        if invalid_rows:
            st.warning(f"Skipping rows with empty/invalid sequences: {invalid_rows}")

        if not sequences:
            st.error("❌ No valid sequences to process. Please check your input.")
            st.stop()

        fasta_id_map = {}
        records = []
        for name, seq in zip(names, sequences):
            fasta_id = name[:20]  # Truncate to avoid Biopython Clustal parser errors
            fasta_id_map[fasta_id] = name
            records.append(SeqRecord(Seq(seq), id=fasta_id, description=""))

        if use_uploaded_ref:
            ref_record = list(SeqIO.parse(ref_fasta, "fasta"))[0]
        else:
            if ref_row not in selected_rows or ref_row in invalid_rows:
                st.error("❌ Reference row is invalid or missing a sequence.")
                st.stop()
            ref_idx = selected_rows.index(ref_row) - sum(r < ref_row for r in invalid_rows)
            ref_record = records[ref_idx]

        all_records = [ref_record] + [r for r in records if r.id != ref_record.id]
        id_map = {r.id: f"{fasta_id_map[r.id]} (Row {row_map[fasta_id_map[r.id]]})" for r in all_records}

        with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".fasta") as fasta_file:
            SeqIO.write(all_records, fasta_file, "fasta")
            fasta_path = fasta_file.name

        with open(fasta_path, 'r') as preview:
            st.subheader("🗒 Preview of FASTA File Sent to Alignment Server")
            st.code(preview.read(), language="text")

        pairwise_identities = {}

        # Submit MSA job to Clustal Omega
        msa_url = "https://www.ebi.ac.uk/Tools/services/rest/clustalo/run"
        msa_payload = {
            'email': 'anonymous@example.com',
            'sequence': open(fasta_path).read(),
            'stype': 'protein',
            'outfmt': 'clustal'
        }
        msa_response = requests.post(msa_url, data=msa_payload)
        if not msa_response.ok:
            st.error("❌ Failed to submit MSA to Clustal Omega.")
            st.stop()
        job_id = msa_response.text.strip()

        status_url = msa_url.replace("/run", f"/status/{job_id}")
        for _ in range(30):
            status = requests.get(status_url).text.strip()
            if status == "FINISHED":
                break
            time.sleep(2)

        result_url = msa_url.replace("/run", f"/result/{job_id}/aln-clustal")
        aln_text = requests.get(result_url).text.strip()
        if not aln_text.startswith("CLUSTAL"):
            st.error("❌ Clustal Omega result is not a valid Clustal alignment.")
            st.text_area("Raw Output", aln_text, height=300)
            st.stop()

        st.subheader("🔌 Clustal Omega Alignment Preview")
        st.code(aln_text, language="text")

        try:
            alignment = AlignIO.read(StringIO(aln_text), "clustal")
        except AssertionError:
            st.error("❌ Clustal Omega returned a malformed alignment (possible consensus error).")
            st.text_area("Raw Clustal Output", aln_text, height=300)
            st.stop()
        except Exception as e:
            st.error("❌ Failed to parse Clustal Omega alignment.")
            st.exception(e)
            st.text_area("Raw Clustal Output", aln_text, height=300)
            st.stop()
        ref_id = ref_record.id
        matching_seqs = [r.seq for r in alignment if r.id == ref_id]

        if not matching_seqs:
            st.error(f"❌ Reference ID '{ref_id}' not found in Clustal alignment. Available IDs:
{', '.join([r.id for r in alignment])}")
            st.text_area("Raw Clustal Output", aln_text, height=300)
            st.stop()

        ref_aligned_seq = str(matching_seqs[0])
        ref_map = map_ref_positions(ref_aligned_seq)

        # Prepare result structure
        data = {
            "Name": [ref_record.id],
            "MSA Pairwise Identity %": [100.0]
        }
        if compute_individual_alignments:
            data["Individual Alignment %"] = [100.0]

        for pos in aa_positions:
            align_idx = ref_map.get(pos)
            ref_aa = ref_aligned_seq[align_idx] if align_idx is not None else "-"
            data[f"AA @ Pos {pos}\n(MSA:{align_idx+1 if align_idx is not None else 'N/A'})"] = [ref_aa]

        for record in alignment:
            if record.id == ref_trunc:
                continue
            msa_id = compute_gapped_identity(ref_aligned_seq, str(record.seq))
            row = {
                "Name": id_map.get(record.id, record.id),
                "MSA Pairwise Identity %": msa_id
            }
            if compute_individual_alignments:
                row["Individual Alignment %"] = None

            for pos in aa_positions:
                align_idx = ref_map.get(pos)
                test_aa = str(record.seq[align_idx]) if align_idx is not None else "-"
                data[f"AA @ Pos {pos}\n(MSA:{align_idx+1 if align_idx is not None else 'N/A'})"].append(test_aa)

            for key in ["Name", "MSA Pairwise Identity %"] + (["Individual Alignment %"] if compute_individual_alignments else []):
                data[key].append(row[key])

        df_results = pd.DataFrame(data)
        styled_table = df_results.style.map(color_identity, subset=[col for col in df_results.columns if "%" in col])
        placeholder_table = st.empty()
        placeholder_table.dataframe(styled_table)

        csv = df_results.to_csv(index=False)
        st.download_button("📅 Export Results as CSV", data=csv.encode('utf-8'), file_name="alignment_results.csv", mime="text/csv")

        # Optional individual alignment processing
        if compute_individual_alignments:
            import tracemalloc
            tracemalloc.start()
            start_time = time.time()
            progress = st.progress(0)
            for i, record in enumerate(records):
                if record.id == ref_record.id:
                    continue
                try:
                    score = clustalo_pairwise_alignment(str(ref_record.seq), str(record.seq))
                    idx = df_results.index[df_results['Name'] == record.id].tolist()
                    if idx:
                        df_results.at[idx[0], "Individual Alignment %"] = score
                except Exception as e:
                    st.warning(f"⚠️ Alignment failed for {record.id}: {e}")
                progress.progress((i + 1) / len(records))
                styled_table = df_results.style.map(color_identity, subset=[col for col in df_results.columns if "%" in col])
                placeholder_table.dataframe(styled_table)
            current, peak = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            st.success(f"⏱️ Individual alignments completed in {time.time() - start_time:.2f} seconds")
            st.info(f"💾 Memory used: {current / 10**6:.2f} MB (Peak: {peak / 10**6:.2f} MB)")

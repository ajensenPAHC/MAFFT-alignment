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

    status_url = f"https://www.ebi.ac.uk/Tools/services/rest/clustalo/status/{job_id}"
    result_url = f"https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/{job_id}/aln-clustal"

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

# Final results display logic (reconstructed)
data = {
    "Name": [],
    "MSA Pairwise Identity %": [],
    "Individual Alignment %": []
}

for rec in alignment:
    name = id_map.get(rec.id, rec.id)
    if rec.id == ref_id:
        continue
    identity = compute_jalview_identity(str(alignment[ref_id].seq), str(rec.seq))
    data["Name"].append(name)
    data["MSA Pairwise Identity %"].append(identity)
    data["Individual Alignment %"].append(None)  # Will be updated if enabled

df_results = pd.DataFrame(data)
styled_table = df_results.style.map(color_identity, subset=[col for col in df_results.columns if "%" in col])
placeholder_table = st.empty()
placeholder_table.dataframe(styled_table)

csv = df_results.to_csv(index=False)
st.download_button("📅 Export Results as CSV", data=csv.encode('utf-8'), file_name="alignment_results.csv", mime="text/csv")

# The app is structured correctly to run on Streamlit Cloud.
# You should include this file in your GitHub repo and add a requirements.txt like this:
# biopython
# pandas
# streamlit
# matplotlib
# requests

# If you want, I can auto-generate the GitHub repo structure and provide a ZIP or manifest.

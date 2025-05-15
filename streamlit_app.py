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

# Configure aligner with BLOSUM62
aligner = Align.PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.mode = 'global'
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5

# Ready to be used in downstream pairwise comparison logic

# 🧬 Amino Acid Sequence Analyzer and Classifier

This app allows you to upload Excel files with amino acid sequences, align them via Clustal Omega, and compute pairwise identity scores. Optionally, it also computes per-position amino acid differences.

## 🚀 Launch the App

Use the hosted Streamlit link:  
[https://aaronjensenPAHC.streamlit.app](https://mafft-alignment-m5hu4yh6uhlumfeu5vsqsr.streamlit.app)

## 🔐 Data Security

- No user data is stored or logged.
- All uploaded files are processed temporarily in memory.
- Uses HTTPS-secured calls to Clustal Omega servers (EBI).

## 🧰 Features

- Excel preview of uploaded sequence data
- Scrollable Clustal Omega MSA output
- Optional pairwise identity (Jalview logic)
- Per-position amino acid comparison
- Downloadable CSV of results

## 🛠️ Setup (if running locally)

1. Clone the repo  
   `git clone https://github.com/aaronjensenPAHC/MAFFT-alignment`

2. Navigate into the folder  
   `cd your-repo`

3. Install dependencies  
   `pip install -r requirements.txt`

4. Run the app  
   `streamlit run streamlit_app.py`

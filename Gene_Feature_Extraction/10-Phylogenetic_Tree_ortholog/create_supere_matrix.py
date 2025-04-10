#!/usr/bin/env python3
import os
import subprocess
import pandas as pd
from collections import defaultdict

########################################################
#  CONFIGURATION: Absolute paths (WSL format)
########################################################

PROTEINORTHO_TSV = "/home/sathya/proteinortho/myproject.proteinortho.tsv"
FASTA_TAKIN       = "/home/sathya/interproscan_protein/takin_proteins.fasta"
FASTA_BUFFALO     = "/home/sathya/interproscan_protein/water_buffalo_proteins.fasta"
FASTA_YAK         = "/home/sathya/interproscan_protein/wild_yak_proteins.fasta"

GROUP_FASTA_DIR   = "group_fastas"       # For individual orthogroup FASTAs
ALIGNED_DIR       = "aligned_fastas"     # For each group's MAFFT alignment
OUTPUT_SUPERMAT   = "core_orthologs_supermatrix.fasta"

# MAFFT command – full path as installed in your WSL environment.
MAFFT_CMD         = "/home/sathya/miniconda3/bin/mafft"

########################################################
#  HELPER: (Not required in WSL if using Linux paths directly,
#  but we include it for consistency if needed.)
########################################################

def windows_to_wsl_path(win_path):
    """
    Converts a Windows path to a WSL-compatible path.
    Example: "D:\folder\subfolder" -> "/mnt/d/folder/subfolder"
    (Not needed if using absolute Linux paths.)
    """
    drive, rest = os.path.splitdrive(win_path)
    drive_letter = drive[0].lower()
    rest = rest.replace("\\", "/")
    return f"/mnt/{drive_letter}{rest}"

########################################################
#  FUNCTION: Load FASTA into Dictionary
########################################################

def load_fasta_to_dict(fasta_path):
    """
    Loads a FASTA file into a dictionary: {sequenceID -> sequence}.
    """
    seqdict = {}
    with open(fasta_path, "r") as f:
        seq_id = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq_id and seq_lines:
                    seqdict[seq_id] = "".join(seq_lines)
                seq_id = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)
        if seq_id and seq_lines:
            seqdict[seq_id] = "".join(seq_lines)
    return seqdict

########################################################
#  MAIN PIPELINE
########################################################

def main():
    # 1) Load proteomes into memory
    print("[INFO] Loading proteomes into memory...")
    takin_dict   = load_fasta_to_dict(FASTA_TAKIN)
    buffalo_dict = load_fasta_to_dict(FASTA_BUFFALO)
    yak_dict     = load_fasta_to_dict(FASTA_YAK)

    # 2) Create output directories
    os.makedirs(GROUP_FASTA_DIR, exist_ok=True)
    os.makedirs(ALIGNED_DIR, exist_ok=True)

    # 3) Parse the Proteinortho TSV file
    df_rows = []
    with open(PROTEINORTHO_TSV, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("#") or not line:
                continue
            row = line.split("\t")
            df_rows.append(row)
    df = pd.DataFrame(df_rows)
    # Assign column names (assuming 6 initial columns: adjust if needed)
    df.columns = ["GroupID", "Connectivity", "SpeciesCount", "TakinIDs", "BuffaloIDs", "YakIDs"] + df.columns[6:].tolist()

    # 4) Filter groups having exactly one protein ID per species
    selected_rows = []
    for idx, row in df.iterrows():
        takin_id   = str(row["TakinIDs"]).strip()
        buffalo_id = str(row["BuffaloIDs"]).strip()
        yak_id     = str(row["YakIDs"]).strip()
        if "," in takin_id or "," in buffalo_id or "," in yak_id:
            continue
        if takin_id.lower() == "nan" or buffalo_id.lower() == "nan" or yak_id.lower() == "nan":
            continue
        selected_rows.append((row["GroupID"], takin_id, buffalo_id, yak_id))
    print(f"[INFO] Found {len(selected_rows)} orthogroups with exactly 1 protein ID per species.")

    # 5) Create FASTA files for each group and run MAFFT with --anysymbol
    aligned_files = []
    for (grp, tk, bf, yk) in selected_rows:
        fasta_out = os.path.join(GROUP_FASTA_DIR, f"group_{grp}.fasta")
        with open(fasta_out, "w") as outF:
            seq_tk = takin_dict.get(tk)
            seq_bf = buffalo_dict.get(bf)
            seq_yk = yak_dict.get(yk)
            if not seq_tk or not seq_bf or not seq_yk:
                continue
            outF.write(f">{grp}_Takin\n{seq_tk}\n")
            outF.write(f">{grp}_Buffalo\n{seq_bf}\n")
            outF.write(f">{grp}_Yak\n{seq_yk}\n")

        # No conversion needed here since our paths are already in Linux format.
        fasta_wsl = os.path.abspath(fasta_out)
        aligned_out = os.path.join(ALIGNED_DIR, f"group_{grp}_aligned.fasta")
        cmd = [MAFFT_CMD, "--anysymbol", "--auto", fasta_wsl]
        try:
            with open(aligned_out, "w") as alnF:
                subprocess.run(cmd, stdout=alnF, stderr=subprocess.PIPE, check=True, text=True)
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] MAFFT failed for group {grp}.")
            print("Command:", " ".join(cmd))
            print("Error message:", e.stderr)
            continue
        aligned_files.append(aligned_out)

    print(f"[INFO] Aligned {len(aligned_files)} groups. Alignments are in '{ALIGNED_DIR}'.")

    # 6) Concatenate individual alignments into a supermatrix
    supermatrix = defaultdict(str)  # {Species: concatenated sequence}
    aligned_count = 0
    for aln_file in aligned_files:
        seqdict = {}
        with open(aln_file, "r") as f:
            seq_id = None
            seq_buf = []
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if seq_id and seq_buf:
                        seqdict[seq_id] = "".join(seq_buf)
                    seq_id = line[1:]
                    seq_buf = []
                else:
                    seq_buf.append(line)
            if seq_id and seq_buf:
                seqdict[seq_id] = "".join(seq_buf)
        try:
            tk_key = [k for k in seqdict if k.endswith("Takin")][0]
            bf_key = [k for k in seqdict if k.endswith("Buffalo")][0]
            yk_key = [k for k in seqdict if k.endswith("Yak")][0]
        except IndexError:
            continue
        supermatrix["Takin"]   += seqdict[tk_key]
        supermatrix["Buffalo"] += seqdict[bf_key]
        supermatrix["Yak"]     += seqdict[yk_key]
        aligned_count += 1

    print(f"[INFO] Successfully concatenated {aligned_count} alignments into a supermatrix.")

    # 7) Write the supermatrix FASTA file
    with open(OUTPUT_SUPERMAT, "w") as outF:
        for sp in ["Takin", "Buffalo", "Yak"]:
            outF.write(f">{sp}\n{supermatrix[sp]}\n")

    print(f"[INFO] SUPER-MATRIX saved to {OUTPUT_SUPERMAT} with sequences for Takin, Buffalo, and Yak.")

if __name__ == "__main__":
    main()

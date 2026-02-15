# Assignment-2-Maan-Soni  
## Bioinformatics (BBL434)

---

## 📌 Overview

This project implements a **Nucleotide Sequence Alignment Tool** using:

- **Smith–Waterman Local Alignment**
- **Affine Gap Penalty**
- **Biopython for FASTA parsing**

The implementation follows the assignment instruction and is validated against the EBI LALIGN online tool.

---

## 🧬 Algorithm Used

The tool performs:

### ✅ Local Alignment (Smith–Waterman)
- Finds the highest scoring local alignment region.
- Alignment can start and end at any position.

### ✅ Affine Gap Penalty
Gap penalty model:

Gap Penalty = Gap Opening + (k × Gap Extension)


Where:
- Gap Opening = penalty for starting a gap
- Gap Extension = penalty for extending an existing gap
- k = length of gap

Three matrices are used:
- **M** → Match/Mismatch
- **Ix** → Gap in Sequence 1
- **Iy** → Gap in Sequence 2

---

## 📥 Inputs

The program takes:

1. Two FASTA files (command-line arguments)
2. Scoring parameters (user input during execution):
   - Match score
   - Mismatch score
   - Gap opening penalty
   - Gap extension penalty

---

## 📤 Output

The program outputs:

- Best Local Alignment (Sequence 1)
- Best Local Alignment (Sequence 2)
- Alignment Score

---

## 🛠 Installation

### Step 1 — Install Python (if not installed)
Download from:
https://www.python.org/downloads/

### Step 2 — Install Biopython

pip install biopython

---

## ▶️ How to Run

From the project directory:

python affine-alignment.py seq1.fasta seq2.fasta


You will then be prompted to enter:

Match score:
Mismatch score:
Gap Opening penalty:
Gap Extension penalty:


---

## 🧪 Testing & Validation

The implementation was tested using the following EBI tool:

https://www.ebi.ac.uk/jdispatcher/psa/lalign?stype=dna&gapext=0

### Validation Steps:
1. Enter the same sequences in the EBI LALIGN tool.
2. Use identical scoring parameters.
3. Compare:
   - Highest Waterman–Eggert score
   - Alignment region
   - Gap placement

The alignment score and region matched the EBI output.

---

## ⏱ Time and Space Complexity

- **Time Complexity:** O(nm)  
- **Space Complexity:** O(nm)

Where:
- n = length of sequence 1
- m = length of sequence 2

---

## 📂 Project Structure

Assignment-2-Maan-Soni/
│
├── affine-alignment.py
├── seq1.fasta
├── seq2.fasta
└── README.md


---

## 🎓 Course Information

Course: **BBL434 – Bioinformatics**  
Assignment: **Assignment 2**  
Student: **Maan Soni**

---

## 📖 References

- Smith TF, Waterman MS. Identification of common molecular subsequences.
- EMBOSS Needle and LALIGN documentation.
- Biopython Documentation.

---
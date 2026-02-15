import sys
from Bio import SeqIO


def read_fasta(filename):
    """
    Read first sequence from a FASTA file using Biopython.
    """
    record = SeqIO.read(filename, "fasta")
    return str(record.seq)


def affine_local_alignment(seq1, seq2, match, mismatch, gap_open, gap_extend):
    """
    Smith-Waterman Local Alignment with Affine Gap Penalty
    """

    n = len(seq1)
    m = len(seq2)

    NEG_INF = -10**9

    M = [[0] * (m + 1) for _ in range(n + 1)]
    Ix = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    Iy = [[NEG_INF] * (m + 1) for _ in range(n + 1)]

    max_score = 0
    max_i, max_j = 0, 0
    max_state = "M"

    for i in range(1, n + 1):
        for j in range(1, m + 1):

            score = match if seq1[i - 1] == seq2[j - 1] else mismatch

            M[i][j] = max(
                0,
                M[i - 1][j - 1] + score,
                Ix[i - 1][j - 1] + score,
                Iy[i - 1][j - 1] + score
            )

            Ix[i][j] = max(
                M[i][j - 1] - (gap_open + gap_extend),
                Ix[i][j - 1] - gap_extend
            )

            Iy[i][j] = max(
                M[i - 1][j] - (gap_open + gap_extend),
                Iy[i - 1][j] - gap_extend
            )

            current = max(M[i][j], Ix[i][j], Iy[i][j])

            if current > max_score:
                max_score = current
                max_i, max_j = i, j
                if current == M[i][j]:
                    max_state = "M"
                elif current == Ix[i][j]:
                    max_state = "Ix"
                else:
                    max_state = "Iy"

    # Traceback
    align1 = ""
    align2 = ""

    i, j = max_i, max_j
    state = max_state

    while i > 0 and j > 0:

        if state == "M":
            if M[i][j] == 0:
                break

            score = match if seq1[i - 1] == seq2[j - 1] else mismatch

            if M[i][j] == M[i - 1][j - 1] + score:
                state = "M"
            elif M[i][j] == Ix[i - 1][j - 1] + score:
                state = "Ix"
            else:
                state = "Iy"

            align1 = seq1[i - 1] + align1
            align2 = seq2[j - 1] + align2
            i -= 1
            j -= 1

        elif state == "Ix":
            if Ix[i][j] == M[i][j - 1] - (gap_open + gap_extend):
                state = "M"
            else:
                state = "Ix"

            align1 = "-" + align1
            align2 = seq2[j - 1] + align2
            j -= 1

        else:  # Iy
            if Iy[i][j] == M[i - 1][j] - (gap_open + gap_extend):
                state = "M"
            else:
                state = "Iy"

            align1 = seq1[i - 1] + align1
            align2 = "-" + align2
            i -= 1

    return align1, align2, max_score


def main():

    if len(sys.argv) != 3:
        print("\nUsage: python affine-alignment.py <seq1.fasta> <seq2.fasta>\n")
        sys.exit(1)

    file1 = sys.argv[1]
    file2 = sys.argv[2]

    seq1 = read_fasta(file1)
    seq2 = read_fasta(file2)

    print("\n==============================")
    print("   SCORING MATRIX INPUT")
    print("==============================")
    print("Please enter the scoring parameters:\n")

    match = int(input("Match score: "))
    mismatch = int(input("Mismatch score: "))
    gap_open = int(input("Gap Opening penalty: "))
    gap_extend = int(input("Gap Extension penalty: "))

    print("\n==============================")
    print("   COMPUTING ALIGNMENT...")
    print("==============================")

    align1, align2, score = affine_local_alignment(
        seq1, seq2, match, mismatch, gap_open, gap_extend
    )

    print("\n===== BEST LOCAL ALIGNMENT =====\n")
    print(align1)
    print(align2)
    print("\nAlignment Score:", score)


if __name__ == "__main__":
    main()

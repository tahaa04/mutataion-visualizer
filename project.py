import streamlit as st
from Bio.Data import CodonTable
from Bio.Align import substitution_matrices
from rdkit import Chem
from rdkit.Chem import Draw


# ------------------ CONSTANT DATA --------------------

AA_FULL_NAMES = {
    "A": "Alanine", "R": "Arginine", "N": "Asparagine", "D": "Aspartic Acid",
    "C": "Cysteine", "Q": "Glutamine", "E": "Glutamic Acid", "G": "Glycine",
    "H": "Histidine", "I": "Isoleucine", "L": "Leucine", "K": "Lysine",
    "M": "Methionine", "F": "Phenylalanine", "P": "Proline",
    "S": "Serine", "T": "Threonine", "W": "Tryptophan",
    "Y": "Tyrosine", "V": "Valine"
}

AA_SMILES = {
    "A": "C[C@H](N)C(=O)O",
    "R": "N[C@H](CCCNC(=N)N)C(=O)O",
    "N": "N[C@H](CC(=O)N)C(=O)O",
    "D": "N[C@H](CC(=O)O)C(=O)O",
    "C": "N[C@H](CS)C(=O)O",
    "Q": "N[C@H](CCC(=O)N)C(=O)O",
    "E": "N[C@H](CCC(=O)O)C(=O)O",
    "G": "NCC(=O)O",
    "H": "N[C@H](Cc1ncc[nH]1)C(=O)O",
    "I": "CC[C@H](C)[C@H](N)C(=O)O",
    "L": "CC(C)C[C@H](N)C(=O)O",
    "K": "NCCCC[C@H](N)C(=O)O",
    "M": "CSCC[C@H](N)C(=O)O",
    "F": "N[C@H](Cc1ccccc1)C(=O)O",
    "P": "O=C(O)[C@H]1CCCN1",
    "S": "N[C@H](CO)C(=O)O",
    "T": "CC(O)[C@H](N)C(=O)O",
    "W": "N[C@H](Cc1c[nH]c2ccccc12)C(=O)O",
    "Y": "N[C@H](Cc1ccc(O)cc1)C(=O)O",
    "V": "CC(C)[C@H](N)C(=O)O"
}

CANONICAL_AAS = set(AA_FULL_NAMES.keys())
GENETIC_CODE = CodonTable.unambiguous_rna_by_name["Standard"]
BLOSUM62 = substitution_matrices.load("BLOSUM62")


# --------------------functions---------------------------------

def validate_codon(codon):
    """Validate RNA codon."""
    if len(codon) != 3:
        raise ValueError("Codon must be exactly 3 nucleotides long.")
    if not set(codon).issubset({"A", "U", "G", "C"}):
        raise ValueError("Invalid RNA codon (use A, U, G, C only).")
    if codon in GENETIC_CODE.stop_codons:
        raise ValueError("Stop codon does not encode an amino acid.")


def codon_to_aa(codon):
    """Translate RNA codon to amino acid (1-letter)."""
    return GENETIC_CODE.forward_table[codon]


def draw_structure(aa):
    """Return RDKit image for an amino acid."""
    mol = Chem.MolFromSmiles(AA_SMILES[aa])
    return Draw.MolToImage(mol, size=(300, 300))


def main():
    st.set_page_config(page_title="Mutation Visualizer", layout="wide")

    st.markdown(
        """
        <style>
            :root {
                --text-blue: #1f4fd8;
            }

            .stApp {
                background-color: black;
                color: var(--text-blue);
            }

            h1, h2, h3, p, label, div {
                color: var(--text-blue) !important;
            }

            input {
                background-color: #111 !important;
                color: var(--text-blue) !important;
                caret-color: var(--text-blue);
            }

            /* Buttons styled as cards */
            div.stButton > button {
                width: 100%;
                border: 1px solid #2a5d9f;
                background: #0b1220;
                color: var(--text-blue);
                border-radius: 12px;
                padding: 0.6rem;
                text-align: left;
            }

            div.stButton > button:hover {
                background: #0f1a33;
                border-color: #4da3ff;
            }

            /* Remove default Streamlit gray */
            textarea, input[type="text"] {
                color: var(--text-blue) !important;
            }
        </style>
        """,
        unsafe_allow_html=True
    )


    if "selected_aa" not in st.session_state:
        st.session_state.selected_aa = None

    st.markdown("<h1 style='text-align:center;'>Codon → Amino Acid Explorer</h1>", unsafe_allow_html=True)

    _, center, _ = st.columns([1, 2, 1])
    with center:
        codon = st.text_input("Enter RNA codon (e.g. AUG)", max_chars=3).upper()

    if not codon:
        return

    try:
        validate_codon(codon)
        aa = codon_to_aa(codon)
    except ValueError as e:
        st.error(str(e))
        return

    aa_name = AA_FULL_NAMES[aa]

    if st.session_state.selected_aa is None or st.session_state.selected_aa == aa:
        choices = [x for x in BLOSUM62.alphabet if x in CANONICAL_AAS and x != aa]
        st.session_state.selected_aa = max(choices, key=lambda x: BLOSUM62[aa, x])

    spacer_l, left, _, right, spacer_r = st.columns([1, 3, 1, 3, 1])

    with left:
        st.subheader("Original amino acid")
        st.markdown(f"**{aa_name} ({aa})**")
        st.image(draw_structure(aa))

    with right:
        sel = st.session_state.selected_aa
        st.subheader("Selected amino acid")
        st.markdown(f"**{AA_FULL_NAMES[sel]} ({sel})**")
        st.image(draw_structure(sel))

    st.markdown("---")
    st.subheader("BLOSUM62 substitution scores")

    others = sorted(
        [x for x in BLOSUM62.alphabet if x in CANONICAL_AAS and x != aa],
        key=lambda x: BLOSUM62[aa, x],
        reverse=True
    )

    cols_per_row = 7
    for i in range(0, len(others), cols_per_row):
        row = st.columns(cols_per_row)
        for j, other in enumerate(others[i:i + cols_per_row]):
            score = int(BLOSUM62[aa, other])
            label = f"{other} — {AA_FULL_NAMES[other]}\nScore: {score}"
            with row[j]:
                if st.button(label, key=f"{aa}_{other}"):
                    st.session_state.selected_aa = other
                    st.rerun()



if __name__ == "__main__":
    main()

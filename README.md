    # MUTATION VISUALIZER

    #### Video Demo:  https://www.youtube.com/watch?v=uTbsKOjm7rY

    #### Description:

    Mutation visualizer is a python based web app made with the streamlit library which converts the 3 nucleotide RNA codon entered by the user into its corresponding amino acid, and displays substitution probability scores for the other essential amino acids, which when selected displays it's structure alongside the structure of the initial amino acid corresponding to the codon entered.
    The program employs subtitution matrices (to access the BLOSUM62 matrix) and codon table (codon- amino acid) from the biopython library, and uses the RDKit library to draw the chemical structures of the amino acids.
    Other amino acids along with their substitution probability scores (from the BLOSUM62 matrix) are displayed towards the bottom of the web app. When a user selects one of them, the structures of the amino acid corresponding the codon and the selected amino acid are diplayed side-by-side.

    The idea is for the user to observe and study the similarities and differences in the structures of both amino acids, with respect to their substitution likelihood scores. It is not uncommon to observe similarity in structure of amino acids with higher probabilistic score of substitution.
    A point substitution may or may not be deleterious- which maybe reflective of the structural difference it exhibits- which can be observed as a result of the side-by-side design of the program. An amino acid with different structural properties than of the one it replaces may affect the folding and function of the entire protein (as discussed in the demo video).

    Functions such as codon validation, codon-to-amino-acid translation, and molecular structure generation are defined independently of Streamlit. This makes the logic easier to test, reuse, and reason about, and ensures the project meets CS50’s requirement for automated testing using pytest.

    The program with the 'codon validation' function, validates the user input, and deals with any inout that is not an RNA codon.
    Another important design choice was restricting substitution analysis to canonical amino acids only. While BLOSUM matrices may include ambiguous symbols such as B, Z, X, or stop characters, these do not correspond to concrete chemical structures. Therefore, the project explicitly filters substitution scores to the 20 standard amino acids for which names and molecular representations are defined.

    The project consists of the following files:

- `project.py`: Contains the main Streamlit application as well as all core logic functions, including codon validation, codon translation, and molecular structure rendering.
- `test_project.py`: Contains pytest-based unit tests for the core logic functions. These tests verify correct behavior for valid and invalid codons, correct amino acid translation, and successful molecular structure generation.
- `README.md`: This documentation file, which explains the project’s purpose, structure, and design decisions.
- `requirements.txt`: This text file contains the dependencies used to build and run this project.

    This project uses external libraries including Streamlit, Biopython, RDKit, and pytest. These were installed locally using a Conda environment; however, the core logic of the project is independent of the installation method. The automated tests focus only on non-interactive logic and do not require the Streamlit interface to be executed.

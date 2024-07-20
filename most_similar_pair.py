## Molecular fingerprints are molecular descriptors that encode molecular features of a molecule
## in the form of a binary digit(0 or 1). A bit is ON if a certain fragment is found in molecular
## structure. Similarity is based on determining the number of bits that are common to two structure.
## The Tanimoto coefficient is a popular similarity index to compare bitstrings by applying similarity
## indices. 

## Tanimoto coefficient: T(AB) = (c) / (a + b - c) 
## a and b are the number of bits set to 1 in molecule A and B, respectively.
## c is the number of set bits common to both A and B.
## The value of T (similarity) is between zero and one. Zero means maximal dissimalirity.


from rdkit import Chem #to read and write SMILES strings
from rdkit.Chem import rdFingerprintGenerator 
from rdkit import DataStructs # to work with fingerprint and similarity data structures
from itertools import combinations

# Read SMILES file and convert to RDKit molecules
def read_smiles_file(file_path):
    try:
        with open(file_path, 'r') as f:
            smiles_list = f.readlines()
        molecules = [Chem.MolFromSmiles(smiles.strip()) for smiles in smiles_list]
        return molecules
    except FileNotFoundError:
        print(f"File '{file_path}' not found.")
        return None
        
# Function to calculate Tanimoto similarity between two molecules using RDKit's Morgan fingerprints
def calculate_similarity(mol1, mol2):
    generator = rdFingerprintGenerator.GetMorganGenerator(radius = 2)
    fp1 = generator.GetSparseCountFingerprint(mol1)
    fp2 = generator.GetSparseCountFingerprint(mol2)
    similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
    return similarity
    
# Function to find the most similar compounds
def find_most_similar_compounds(molecules):
    if molecules is None:
        return None, None

    max_similarity = 0
    most_similar_pair = None

    for pair in combinations(enumerate(molecules), 2):
        index1, mol1 = pair[0]
        index2, mol2 = pair[1]
        similarity = calculate_similarity(mol1, mol2)
        if similarity > max_similarity:
            max_similarity = similarity
            most_similar_pair = (index1, index2)
    
    return most_similar_pair, max_similarity
    
# Example usage
if __name__ == "__main__":
    file_path = 'test_smiles.smi'  # Adjust the file path here
    molecules = read_smiles_file(file_path)
    most_similar_pair, similarity_score = find_most_similar_compounds(molecules)
    if most_similar_pair is not None:
        index1, index2 = most_similar_pair
        print(f"The two most similar compounds are at indices {index1} and {index2}, with a similarity score of {similarity_score}.")
    else:
        print("Failed to find the most similar compounds.")

## EXAMPLE OUTPUT
# The two most similar compounds are at indices 11 and 14, with a similarity score of 0.7692307692307693. (radius = 0)
# The two most similar compounds are at indices 11 and 14, with a similarity score of 0.6727272727272727. (radius = 1)
# The two most similar compounds are at indices 11 and 14, with a similarity score of 0.5294117647058824. (radius = 2)

## A smaller radius captures local atomic environments, while a larger radius captures more extended structural features.

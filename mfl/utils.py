from rdkit import Chem
from rdkit.Chem import Descriptors


def isValid(s: str) -> bool:
    m = Chem.MolFromSmiles(s)
    return m is not None


def canon(s: str) -> str:
    return Chem.MolToSmiles(Chem.MolFromSmiles(s))

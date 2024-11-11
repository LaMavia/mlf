from rdkit import Chem


def isValid(s: str) -> bool:
    return Chem.MolFromSmiles(s) is not None


def canon(s: str) -> str:
    return Chem.MolToSmiles(Chem.MolFromSmiles(s))

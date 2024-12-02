from rdkit import Chem


def isValid(s: str) -> bool:
    try:
        m = Chem.MolFromSmiles(s)
        if m is None:
            return False
        Chem.SanitizeMol(m, sanitizeOps=Chem.SANITIZE_ALL)
        return True
    except Exception:
        return False


def canon(s: str) -> str:
    return Chem.MolToSmiles(Chem.MolFromSmiles(s))

from rdkit import Chem
from molvs import Validator, validate_smiles
from molvs.validations import (
    VALIDATIONS,
    DichloroethaneValidation,
    FragmentValidation,
    IsNoneValidation,
    IsotopeValidation,
    NeutralValidation,
)

validator = Validator(
    [
        IsNoneValidation,
        DichloroethaneValidation,
        FragmentValidation,
        NeutralValidation,
        IsotopeValidation,
    ],
    stdout=True,
)


def isValid(s: str) -> bool:
    try:
        m = Chem.MolFromSmiles(s)
        if m is None:
            return False
        Chem.SanitizeMol(m, sanitizeOps=Chem.SANITIZE_ALL)
        errors = validator.validate(s)
        print(f"[VALIDATION] errors: {errors}")

        return True
    except Exception:
        return False


print(
    isValid(
        "CCC1CCC2CC3CCCC1CCCC(C2)C3C(C)CN(CCCC(C)(CC)C(C)C)Cc1ccc(CNC(=O)c2csc3ncNc(=O)c23)cc1"
    )
)

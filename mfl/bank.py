from rdkit import Chem
from mfl.lexer import Lexem, lex


class BankEntry:
    def __init__(self, smiles: list[str]) -> None:
        self.mol: list[Chem.Mol] = [Chem.MolFromSmiles(s) for s in smiles]
        self.lexems: list[list[Lexem]] = [lex(s) for s in smiles]
        self.splittable_indices: list[list[int]] = self._calcSplittableIndices()

    def _calcSplittableIndices(self) -> list[list[int]]:
        return [
            [
                i
                for i, lexem in enumerate(lexems)
                if (len(lexem) == 2 and lexem[1] == 0)
                or (len(lexem) == 3 and lexem[2] == 0)
            ]
            + [len(lexems)]
            for lexems in self.lexems
        ]

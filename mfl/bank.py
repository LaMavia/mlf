from rdkit import Chem
from mfl.lexer import MUTABLE_LEXEMS, Lexem, lex


class BankEntry:
    def __init__(self, smiles: list[str]) -> None:
        self.mols: list[Chem.Mol] = [Chem.MolFromSmiles(s) for s in smiles]
        if any(m is None for m in self.mols):
            raise ValueError()

        self.lexems: list[list[Lexem]] = [lex(s) for s in smiles]
        self.splittable_indices: list[list[int]] = self._calcSplittableIndices()
        self.mutable_indices: list[list[int]] = self._calcMutableIndices()

    # @staticmethod
    # def fromSmiles(smiles: list[str]) -> "BankEntry":

    def _calcMutableIndices(self) -> list[list[int]]:
        return [
            [i for i, lexem in enumerate(lexems) if lexem[0] in MUTABLE_LEXEMS]
            for lexems in self.lexems
        ]

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

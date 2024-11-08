from collections.abc import Callable
from typing import Literal
from rdkit import Chem, logging
from rdkit.DataStructs import TanimotoSimilarity
import numpy as np
import numpy.typing as npt
import json
from logging import Logger
import random as rnd

from mfl.bank import BankEntry
from mfl.lexer import lex, serialize

logger = Logger("mlf", logging.DEBUG)

N = 1_000


# class BankEntry:
#     _raws: list[str]
#
#     @staticmethod
#     def setRaws(raws: list[str]):
#         BankEntry._raws = raws
#
#     def __init__(self, smiles: list[str]):
#         assert len(BankEntry._raws) == len(smiles) + 1
#
#         self.smiles: list[str] = smiles
#         self.mols: list[Chem.Mol] = [Chem.MolFromSmiles(s) for s in smiles]
#         self.rings: list[list[tuple[int, int]]] = []
#         self.splittable_indices: list[set[int]] = []
#
#         self._calcSplittableIncides()
#
#     def _calcSplittableIncides(self):
#         for smile in self.smiles:
#             splittable_indices = {0, len(smile)}
#             n_current_rings = 0
#             current_rings: dict[int, int] = {}  # cycle_bond_number => starting_index
#             i = 0
#             while i < len(smile):
#                 logger.debug(f"i={i}")
#                 c = smile[i]
#                 multi_mol = findMultiAtom(smile, i, BankEntry.MULTI_ATOMS)
#                 if c in BankEntry.SKIP_CHARS:
#                     i += 1
#                 elif c == "[":  # skip bracket-enclosed atoms
#                     splittable_indices.add(i)
#                     j = i + 1
#                     depth = 1  # [] nesting depth
#                     while j < len(smile) and depth > 0:
#                         cj = smile[j]
#                         if cj == "[":
#                             depth += 1
#                         elif cj == "]":
#                             depth -= 1
#                         j += 1
#                     i = j
#                 elif c == "%" or c in BankEntry.DIGITS:
#                     # find the last digit of the number
#                     # to get the whole number
#                     if c == "%":
#                         i += 1
#                     j = i + 1
#                     while j < len(smile) and smile[j] in BankEntry.DIGITS:
#                         j += 1
#                     num = int(smile[i:j])
#                     if num not in current_rings:
#                         # open a new ring
#                         n_current_rings += 1
#                         current_rings[num] = i
#                     else:
#                         n_current_rings -= 1
#                         current_rings.pop(num)
#
#                     # skip the parsed number
#                     i = j
#                 elif multi_mol is not None:
#                     splittable_indices.add(i)
#                     i += len(multi_mol)
#                 elif n_current_rings == 0:  # mark the index as splittable
#                     splittable_indices.add(i)
#                     i += 1
#                 else:
#                     i += 1
#                 logger.debug(f"i={i}")
#
#             self.splittable_indices.append(splittable_indices)
#
#     def toCSV(self) -> str:
#         return f"""{json.dumps(self.smiles)} {json.dumps(BankEntry._raws)}"""
#
#     def toSmiles(self) -> str:
#         acc = BankEntry._raws[0]
#         for i in range(len(self.smiles)):
#             acc += self.smiles[i] + BankEntry._raws[i + 1]
#
#         return acc


# def initBank() -> list[BankEntry]:
#     return []
#
#
# def calculateDistances(
#     bank: list[BankEntry],
# ) -> tuple[npt.NDArray[np.float64], np.float64]:
#     """
#     Takes the current bank of molecules, and returns their distances, and the average distance
#     """
#     n = len(bank)
#     fingerprints = np.array(
#         [[Chem.RDKFingerprint(mol) for mol in m.mols] for m in bank]
#     )
#     sum_distances = np.zeros_like(fingerprints)
#     mol_range = range(fingerprints.shape[1])
#
#     for i in range(n):
#         for j in range(i):
#             for k in mol_range:
#                 d = 1 - TanimotoSimilarity(fingerprints[i, k], fingerprints[j, k])
#                 sum_distances[i, k] += d
#                 sum_distances[j, k] += d
#
#     return sum_distances / (n - 1), np.mean(sum_distances) / n
#
#
# BankEntry.setRaws(    raws=["", "", ""])
# m = BankEntry(
#     smiles=[
#         "[Na]ClCCBr[Mg]",
#         "CCCC",
#     ],
# )
#
# mi = max([max(indices) for indices in m.splittable_indices]) + 1
# for j, (indices, s) in enumerate(zip(m.splittable_indices, m.smiles)):
#     for i in sorted(m.splittable_indices[j]):
#         print(f"L={s[:i]}{' ' * (mi+1 - i)}R={s[i:]}")
#     print("")
#
#
# def removeEmptyBranches(smile: str) -> str:
#     for i in range(len(smile)):
#         if smile[i : i + 2] == "()":
#             smile = smile[:i] + smile[i + 2 :]
#             # there can be at most 1 empty branch,
#             # since neither part has empty branches
#             return smile
#     return smile
#
# def balanceBranches(smile: str) -> str:
#     SKIP_CHARS = BankEntry.SKIP_CHARS
#     unmatched_open_branches = 0
#     n = len(smile)
#
#     for i in range(n):
#         c = smile[i]
#         multimol = findMultiAtom(smile, i, BankEntry.MULTI_ATOMS)
#         if c in SKIP_CHARS:
#
#
#
# def crossMolecules(a: BankEntry, b: BankEntry) -> list[BankEntry]:
#     for gi, (lsmile, rsmile) in enumerate(zip(a.smiles, b.smiles, strict=True)):
#         ais = list(a.splittable_indices[gi])
#         bis = list(b.splittable_indices[gi])
#         rnd.shuffle(ais)
#         rnd.shuffle(bis)
#
#         for li, ri in zip(ais, bis):
#             pass
#     rnd.shuffle(list())
#
#     return []


def canon(s: str) -> str:
    return Chem.MolToSmiles(Chem.MolFromSmiles(s))


SMILES = [
    "OCCc1c(C)[n+](cs1)Cc2cnc(C)nc2N",
    "[Cu+2].[O-]S(=O)(=O)[O-]",
    "O=Cc1ccc(O)c(OC)c1",
    "CCc(c1)ccc2[n+]1ccc3c2[nH]c4c3cccc4",
    "CN1CCC[C@H]1c2cccnc2",
    r"CCC[C@@H](O)CC\C=C\C=C\C#CC#C\C=C\CO",
    r"CCC[C@@H](O)CC/C=C/C=C/C#CC#C/C=C/CO",
    r"OC[C@@H](O1)[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)1",
    r"OC[C@@H](O1)[C@@H](O)[C@H](O)[C@@H]2[C@@H]1c3c(O)c(OC)c(O)cc3C(=O)O2",
    r"CC(C)(O1)C[C@@H](O)[C@@]1(O2)[C@@H](C)[C@@H]3CC=C4[C@]3(C2)C(=O)C[C@H]5[C@H]4CC[C@@H](C6)[C@]5(C)Cc(n7)c6nc(C[C@@]89(C))c7C[C@@H]8CC[C@@H]%10[C@@H]9C[C@@H](O)[C@@]%11(C)C%10=C[C@H](O%12)[C@]%11(O)[C@H](C)[C@]%12(O%13)[C@H](O)C[C@@]%13(C)CO",
]

m = BankEntry([SMILES[0]])

for lexems, indices in zip(m.lexems, m.splittable_indices):
    for i in indices:
        print(f"L={serialize(lexems[:i])}\t{serialize(lexems[i:])}")

# print(*lex(SMILES[0]), sep="\n")
#
# for s in SMILES:
#     print(s)
#     tokens = lex(s)
#     s1 = serialize(tokens)
#     assert canon(s) == canon(s1)

# assert canon(s) == canon(s1)

# def iterate(bank: set[str])

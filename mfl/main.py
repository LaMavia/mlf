from typing import Generator
from rdkit import Chem, logging
from rdkit.DataStructs import TanimotoSimilarity
import numpy as np
import numpy.typing as npt
from logging import Logger
import random as rnd
from rdkit.rdBase import DisableLog

from mfl.bank import BankEntry
from mfl.lexer import (
    PAREN_CLOSE,
    PAREN_OPEN,
    Lexem,
    repOfLexem,
    serialize,
    typeOfLexem,
)

logger = Logger("mlf", logging.DEBUG)
DisableLog("rdApp.*")


N = 1_000

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


def initBank() -> list[BankEntry]:
    return [BankEntry([s]) for s in SMILES]


def calculateDistances(
    bank: list[BankEntry],
) -> tuple[npt.NDArray[np.float64], np.float64]:
    """
    Takes the current bank of molecules, and returns their distances, and the average distance
    """
    n = len(bank)
    fingerprints = np.array(
        [[Chem.RDKFingerprint(mol) for mol in m.mols] for m in bank]
    )
    sum_distances = np.zeros_like(fingerprints)
    mol_range = range(fingerprints.shape[1])

    for i in range(n):
        for j in range(i):
            for k in mol_range:
                d = 1 - TanimotoSimilarity(fingerprints[i, k], fingerprints[j, k])
                sum_distances[i, k] += d
                sum_distances[j, k] += d

    return sum_distances / (n - 1), np.mean(sum_distances) / n


def removeEmptyBranches(lexems: list[Lexem]) -> list[Lexem]:
    n = len(lexems)
    for i in range(n - 1):
        left = lexems[i]
        right = lexems[i + 1]
        if repOfLexem(left) == PAREN_OPEN and repOfLexem(right) == PAREN_CLOSE:
            return lexems[:i] + lexems[i + 1 :]

    return lexems


def parenthesisSum(lexems: list[Lexem]) -> int:
    sum = 0

    for lexem in lexems:
        t = typeOfLexem(lexem)
        if t == PAREN_OPEN:
            sum += 1
        elif t == PAREN_CLOSE:
            sum -= 1

    return sum


def invalidClosingParentheses(
    lexems: list[Lexem],
) -> Generator[int, None, None]:
    sum = 0

    for i, lexem in enumerate(lexems):
        t = typeOfLexem(lexem)
        if t == PAREN_OPEN:
            sum += 1
        elif t == PAREN_CLOSE:
            sum -= 1
        if sum < 0:
            yield i
            sum = 0


def invalidOpeningParentheses(
    lexems: list[Lexem],
) -> Generator[int, None, None]:
    sum = 0

    for i in range(len(lexems) - 1, -1, -1):
        t = typeOfLexem(lexems[i])
        if t == PAREN_OPEN:
            sum += 1
        elif t == PAREN_CLOSE:
            sum -= 1
        if sum > 0:
            yield i
            sum = 0


def reduceOpeningBalancing(lexems: list[Lexem], insertions: set[int]) -> list[Lexem]:
    """
    Reduce empty branches for the case `S < 0`,
    and insert insertions.
    """

    acc: list[Lexem] = []
    for i, lexem in enumerate(lexems):
        if i not in insertions:
            acc.append(lexem)
            continue
        if typeOfLexem(lexem) == PAREN_CLOSE:
            continue
        acc.append((PAREN_OPEN, PAREN_OPEN, lexem[2]))
        acc.append(lexem)

    return acc


def reduceClosingBalancing(lexems: list[Lexem], insertions: set[int]) -> list[Lexem]:
    """
    Reduce empty branches for the case `S > 0`,
    and insert insertions.
    """

    acc: list[Lexem] = []
    for i, lexem in enumerate(lexems):
        if i not in insertions:
            acc.append(lexem)
            continue
        if typeOfLexem(lexem) == PAREN_OPEN:
            continue
        acc.append(lexem)
        acc.append((PAREN_CLOSE, PAREN_CLOSE, lexem[2]))

    return acc


def balanceBranches(lexems: list[Lexem]) -> Generator[list[Lexem], None, None]:
    # find unbalanced closing parentheses
    nm1 = len(lexems) - 1
    sum = parenthesisSum(lexems)
    if sum == 0:
        yield lexems
    if sum < 0:
        """
        Balance too many closing parenthesis
        """
        inserted_openings: set[int] = set()
        for pi in invalidClosingParentheses(lexems):
            while True:
                insertion_index = rnd.randint(0, pi)
                if insertion_index not in inserted_openings:
                    inserted_openings.add(insertion_index)
                    break
        yield reduceOpeningBalancing(lexems, inserted_openings)
    else:
        """
        Balance too many closing parenthesis
        """
        inserted_openings: set[int] = set()
        for pi in invalidOpeningParentheses(lexems):
            while True:
                insertion_index = rnd.randint(pi, nm1)
                if insertion_index not in inserted_openings:
                    inserted_openings.add(insertion_index)
                    break
        yield reduceClosingBalancing(lexems, inserted_openings)


def crossMolecules(a: BankEntry, b: BankEntry) -> Generator[BankEntry, None, None]:
    for gi, (left_lexems, right_lexems) in enumerate(
        zip(a.lexems, b.lexems, strict=True)
    ):
        ais = a.splittable_indices[gi]
        bis = b.splittable_indices[gi]
        rnd.shuffle(ais)
        rnd.shuffle(bis)
        smis: set[str] = set()

        for li, ri in zip(ais, bis):
            lexems = left_lexems[:li] + right_lexems[ri:]
            lexems = removeEmptyBranches(lexems)
            for balanced_lexems in balanceBranches(lexems):
                smi = serialize(balanced_lexems)
                if not isValid(smi):
                    continue
                smis.add(canon(smi))

        print(*sorted(smis), sep="\n")

    yield a


def isValid(s: str) -> bool:
    return Chem.MolFromSmiles(s) is not None


def canon(s: str) -> str:
    return Chem.MolToSmiles(Chem.MolFromSmiles(s))


mols = initBank()
for i in range(len(mols)):
    for j in range(i):
        for m in crossMolecules(mols[i], mols[j]):
            ...

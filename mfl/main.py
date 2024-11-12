from rdkit import Chem, logging
from rdkit.DataStructs import TanimotoSimilarity
import numpy as np
import numpy.typing as npt
from logging import Logger
from rdkit.rdBase import DisableLog
from itertools import islice

from mfl.bank import BankEntry
from mfl.mutations import removeAtom, replaceAtom

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


mols = initBank()
for i in range(len(mols)):
    m = mols[i]
    print("original:", *m.smiles, sep="\n", end="\n<<<<<<\n")
    for mutated in islice(removeAtom(m), 5):
        print(*mutated.smiles, sep="\n")
    print("\n")
    # for j in range(i):
    #     for m in crossMolecules(mols[i], mols[j]):
    #         ...

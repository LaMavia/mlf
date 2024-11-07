from collections.abc import Callable
from rdkit import Chem, logging
from rdkit.DataStructs import TanimotoSimilarity
import numpy as np
import numpy.typing as npt
import json
from logging import Logger

logger = Logger("mlf", logging.DEBUG)

N = 1_000


def findBy[A](pred: Callable[[A], bool], xs: list[A]) -> A | None:
    for x in xs:
        if pred(x):
            return x
    return None


class BankEntry:
    SKIP_CHARS = ["-", "=", "#", ".", "$", ":", "/", "\\", "+"]
    DIGITS = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]
    MULTI_ATOMS = ["Cl", "Br"]

    def __init__(self, smiles: list[str], raws: list[str]):
        assert len(raws) == len(smiles) + 1

        self.smiles: list[str] = smiles
        self.mols: list[Chem.Mol] = [Chem.MolFromSmiles(s) for s in smiles]
        self._raws: list[str] = raws
        self.rings: list[list[tuple[int, int]]] = []
        self.splittable_indices: list[set[int]] = []

        self._calcSplittableIncides()

    def _calcSplittableIncides(self):
        for smile in self.smiles:
            splittable_indices = {0, len(smile)}
            n_current_rings = 0
            current_rings: dict[int, int] = {}  # cycle_bond_number => starting_index
            i = 0
            while i < len(smile):
                logger.debug(f"i={i}")
                c = smile[i]
                multi_mol = findBy(
                    lambda m: m == smile[i : i + len(m)], BankEntry.MULTI_ATOMS
                )
                if c in BankEntry.SKIP_CHARS:
                    i += 1
                elif c == "[":  # skip bracket-enclosed atoms
                    splittable_indices.add(i)
                    j = i + 1
                    depth = 1  # [] nesting depth
                    while j < len(smile) and depth > 0:
                        cj = smile[j]
                        if cj == "[":
                            depth += 1
                        elif cj == "]":
                            depth -= 1
                        j += 1
                    i = j
                elif c == "%" or c in BankEntry.DIGITS:
                    # find the last digit of the number
                    # to get the whole number
                    if c == "%":
                        i += 1
                    j = i + 1
                    while j < len(smile) and smile[j] in BankEntry.DIGITS:
                        j += 1
                    num = int(smile[i:j])
                    if num not in current_rings:
                        # open a new ring
                        n_current_rings += 1
                        current_rings[num] = i
                    else:
                        n_current_rings -= 1
                        current_rings.pop(num)

                    # skip the parsed number
                    i = j
                elif multi_mol is not None:
                    splittable_indices.add(i)
                    i += len(multi_mol)
                elif n_current_rings == 0:  # mark the index as splittable
                    splittable_indices.add(i)
                    i += 1
                else:
                    i += 1
                logger.debug(f"i={i}")

            self.splittable_indices.append(splittable_indices)

    def toCSV(self) -> str:
        return f"""{json.dumps(self.smiles)} {json.dumps(self._raws)}"""

    def toSmiles(self) -> str:
        acc = self._raws[0]
        for i in range(len(self.smiles)):
            acc += self.smiles[i] + self._raws[i + 1]

        return acc


def initBank() -> list[BankEntry]:
    return []


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


m = BankEntry(
    smiles=[
        "[Na]ClCCBr[Mg]",
        "CCCC",
    ],
    raws=["", "", ""],
)

mi = max([max(indices) for indices in m.splittable_indices]) + 1
for j, (indices, s) in enumerate(zip(m.splittable_indices, m.smiles)):
    for i in sorted(m.splittable_indices[j]):
        print(f"L={s[:i]}{' ' * (mi+1 - i)}R={s[i:]}")
    print("")

# def crossMolecules()


# def iterate(bank: set[str])

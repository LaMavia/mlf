from itertools import islice
from typing import Generator
from rdkit import Chem
from rdkit.DataStructs import TanimotoSimilarity
from mfl.bank import BankEntry
import numpy as np
import numpy.typing as npt
from pathlib import Path
import pandas as pd
# import random as rnd

from mfl.cross import crossMolecules
from mfl.lexer import serialize
from mfl.mutations import addAtom, removeAtom, replaceAtom
from mfl.utils import canon, isValid

# SMILES = [
#     "OCCc1c(C)[n+](cs1)Cc2cnc(C)nc2N",
#     "[Cu+2].[O-]S(=O)(=O)[O-]",
#     "O=Cc1ccc(O)c(OC)c1",
#     "CCc(c1)ccc2[n+]1ccc3c2[nH]c4c3cccc4",
#     "CN1CCC[C@H]1c2cccnc2",
#     r"CCC[C@@H](O)CC\C=C\C=C\C#CC#C\C=C\CO",
#     r"CCC[C@@H](O)CC/C=C/C=C/C#CC#C/C=C/CO",
#     r"OC[C@@H](O1)[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)1",
#     r"OC[C@@H](O1)[C@@H](O)[C@H](O)[C@@H]2[C@@H]1c3c(O)c(OC)c(O)cc3C(=O)O2",
#     r"CC(C)(O1)C[C@@H](O)[C@@]1(O2)[C@@H](C)[C@@H]3CC=C4[C@]3(C2)C(=O)C[C@H]5[C@H]4CC[C@@H](C6)[C@]5(C)Cc(n7)c6nc(C[C@@]89(C))c7C[C@@H]8CC[C@@H]%10[C@@H]9C[C@@H](O)[C@@]%11(C)C%10=C[C@H](O%12)[C@]%11(O)[C@H](C)[C@]%12(O%13)[C@H](O)C[C@@]%13(C)CO",
# ]
#
#
# def initBank() -> list[BankEntry]:
#     return [BankEntry([s]) for s in SMILES]
#


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


class CSA:
    MAX_ADDED = 200
    MAX_REPLACED = 200
    MAX_REMOVED = 200

    GROUP_DEL = "&"

    def __init__(
        self,
        population_dir: Path,
        round: int,
        last_generation_file_path: Path,
        Rd: float,
        DCut: float | None,
        template: str,
        seed_size: int,
    ) -> None:
        """
        Population file format (csv): SMILES, groups(string[&]), score
        Template - string with `?` used as placeholders
        """
        self.population_dir = population_dir
        self.round = round
        self.Rd = Rd
        self.DCut = DCut
        self.template_raws = template.split("?")
        self.seed_size = seed_size

        self.population = self._calcPopulation(pd.read_csv(last_generation_file_path))

    def _calcPopulation(self, df: pd.DataFrame) -> list[BankEntry]:
        df["score"] = pd.to_numeric(df["score"], errors="coerce")
        # df.drop_duplicates(subset=["groups"], inplace=True)
        df.sort_values(by=["score"], ascending=False, inplace=True)

        return [BankEntry(groups.split(CSA.GROUP_DEL)) for groups in df["groups"]]

    def _generateFromPair(
        self, a: BankEntry, b: BankEntry
    ) -> Generator[list[str], None, None]:
        for crossed in crossMolecules(a, b):
            for added in islice(addAtom(crossed), CSA.MAX_ADDED):
                yield sorted(canon(serialize(_)) for _ in added.lexems)
            for replaced in islice(replaceAtom(crossed), CSA.MAX_REPLACED):
                yield sorted(canon(serialize(_)) for _ in replaced.lexems)
            for removed in islice(removeAtom(crossed), CSA.MAX_REPLACED):
                yield sorted(canon(serialize(_)) for _ in removed.lexems)

    @staticmethod
    def _groupsToSmile(template_raws: list[str], groups: list[str]) -> str:
        acc = template_raws[0]
        for i in range(len(groups)):
            acc += groups[i] + template_raws[i + 1]
        return acc

    def generateNextPopulation(self):
        generated: set[str] = set()
        best = self.population[0]
        seed = self.population[: self.seed_size]
        for entry in seed:
            n = 0
            for groups in self._generateFromPair(best, entry):
                n += 1
                generated.add(CSA.GROUP_DEL.join(groups))
            if n == 0:
                for groups in self._generateFromPair(entry, best):
                    generated.add(CSA.GROUP_DEL.join(groups))

        data: dict[str, list[str]] = {"SMILES": [], "groups": []}
        for group_string in generated:
            groups = group_string.split(CSA.GROUP_DEL)
            smiles = CSA._groupsToSmile(self.template_raws, groups)
            if isValid(smiles):
                data["groups"].append(group_string)
                data["SMILES"].append(smiles)

        gen_num = f"{self.round}".rjust(3, "0")
        df = pd.DataFrame(data)
        df.to_csv(self.population_dir / f"{gen_num}.csv", index=False)

    @staticmethod
    def genInitial(template_raws: list[str], groups: list[list[str]], out: Path):
        data: dict[str, list[str]] = {"SMILES": [], "groups": [], "score": []}
        for gs in groups:
            smiles = CSA._groupsToSmile(template_raws=template_raws, groups=gs)
            if isValid(smiles):
                data["groups"].append(CSA.GROUP_DEL.join(gs))
                data["SMILES"].append(smiles)
                data["score"].append("0")

        df = pd.DataFrame(data)
        df.to_csv(out, index=False)

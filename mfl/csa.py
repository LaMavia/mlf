from itertools import chain, islice
import json
from typing import Generator, Iterable
from rdkit import Chem
from rdkit.DataStructs import SparseBitVect, TanimotoSimilarity
from mfl.bank import BankEntry
import numpy as np
import numpy.typing as npt
from pathlib import Path
import pandas as pd

# import random as rnd
from collections import OrderedDict

from mfl.cross import crossMolecules
from mfl.lexer import ATOMS, BONDS, Lexem, isDigit, serialize
from mfl.mutations import addAtom, changeBond, removeAtom, replaceAtom
from mfl.utils import canon, isValid


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


def calculateAvgDistance(df: pd.DataFrame) -> float:
    fingerprints = [
        Chem.RDKFingerprint(Chem.MolFromSmiles(r["SMILES"])) for _, r in df.iterrows()
    ]
    n = len(fingerprints)

    dist: float = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            dist += 1 - TanimotoSimilarity(fingerprints[i], fingerprints[j])

    return dist / n


def itFst[A](it: Iterable[A]) -> A | None:
    for r in it:
        return r


class CSA:
    GROUP_DEL = "&"
    ROUND_OFFSET = 3

    def __init__(
        self,
        population_dir: Path,
        round: int,
        last_generation_file_path: Path,
        meta_file_path: Path,
        Rd: float,
        template: str,
        seed_size: int,
        max_added: int,
        max_removed: int,
        max_replaced: int,
        max_crossed: int,
        n_bank: int,
        dummy_score: bool,
        double_cross: bool,
    ) -> None:
        """
        Population file format (csv): SMILES, groups(string[&]), score
        Template - string with `?` used as placeholders
        """
        df = pd.read_csv(last_generation_file_path)
        try:
            with open(meta_file_path) as f:
                self.meta = json.load(f)
        except Exception:
            print(f"[Meta] Failed to open/parse «{meta_file_path}», using default meta")
            self.meta = {}

        DCut = self.meta.get("dcut")
        if DCut is None or round < CSA.ROUND_OFFSET:
            DCut = calculateAvgDistance(df) / 2
            self.meta["dcut"] = DCut

        self.meta_file_path = meta_file_path
        self.population_dir = population_dir
        self.round = round
        self.Rd = Rd
        self.DCut = max(DCut * Rd**5, DCut * Rd ** (max(round - CSA.ROUND_OFFSET, 0)))
        self.template_raws = template.split("?")
        self.seed_size = seed_size
        self.max_added = max_added
        self.max_removed = max_removed
        self.max_replaced = max_replaced
        self.max_crossed = max_crossed
        self.n_bank = n_bank
        self.dummy_score = dummy_score
        self.double_cross = double_cross

        self.population = self._calcPopulation(df)

    def _sortPopulation(self, df: pd.DataFrame) -> pd.DataFrame:
        df.sort_values(by=["score"], ascending=False, inplace=True)

        return df

    def _mapPopulation(self, df: pd.DataFrame) -> pd.DataFrame:
        df["score"] = pd.to_numeric(df["score"], errors="coerce")
        df = self._sortPopulation(df)

        return df

    def _preprocessRow(self, row: pd.Series) -> dict:
        return {
            **row,
            "mol": (mol := Chem.MolFromSmiles(row["SMILES"])),
            "fingerprint": Chem.RDKFingerprint(mol),
        }

    def _initBank(self, df: pd.DataFrame) -> tuple[dict[str, dict], pd.DataFrame]:
        df = df.sample(frac=1).reset_index(drop=True)
        initial_rows = df.iloc[: self.n_bank]
        remaining_rows = df.iloc[self.n_bank :]

        return {
            row["SMILES"]: self._preprocessRow(row)
            for _, row in initial_rows.iterrows()
        }, remaining_rows

    def _nearestNeighbour(
        self, bank: dict[str, dict], fingerprint: SparseBitVect
    ) -> tuple[dict, float]:
        nn = itFst(bank.values())
        assert nn is not None
        min_dist = 1 - TanimotoSimilarity(nn["fingerprint"], fingerprint)

        for e in bank.values():
            dist = TanimotoSimilarity(e["fingerprint"], fingerprint)
            if dist <= min_dist:
                nn = e
                min_dist = dist

        return nn, min_dist

    def _getBankMin(self, bank: dict) -> dict:
        mn: dict | None = None
        for r in bank.values():
            if mn is None or r["score"] <= mn["score"]:
                mn = r

        assert mn is not None
        return mn

    def _calcPopulation(self, df: pd.DataFrame) -> list[BankEntry]:
        """
        # min_score, min_mol
        if mol.score < min_score:
            discard mol
            continue

        (nn, dist) = nearest_neighbour(bank, mol)
        if dist < self.DCut:
            if nn.score > mol.score:
                discard mol
                continue
            replace(bank, nn, mol)
        elif dist > self.DCut:
            replace(bank, min_mol, mol)
            min_mol = mol
            min_score = mol.score
        """
        dcut = self.DCut
        assert dcut is not None
        df = self._mapPopulation(df)

        bank, remaining_rows = self._initBank(df)
        min_mol = self._getBankMin(bank)

        for _, r in remaining_rows.iterrows():
            if (s := r["score"]) < min_mol["score"]:
                print(f"Score too low: {s}")
                continue

            r = self._preprocessRow(r)
            nn, dist = self._nearestNeighbour(bank, r["fingerprint"])
            if dist < dcut:
                print(f"dist < cut: {dist} < {dcut}", end=": ")
                if nn["score"] > r["score"]:
                    print("worse score, discarding")
                    continue
                print(f"""better score: {nn['SMILES']} -> {r['SMILES']}""")
                del bank[nn["SMILES"]]
                bank[r["SMILES"]] = r
                if nn["SMILES"] == min_mol["SMILES"]:
                    min_mol = r
            elif dist > dcut:
                print(f"replacing min: {min_mol['SMILES']} -> {r['SMILES']}")
                del bank[min_mol["SMILES"]]
                bank[r["SMILES"]] = r
                min_mol = r

        print("Bank size:", len(bank))
        bank = sorted(bank.values(), key=lambda r: r["score"], reverse=True)
        return [BankEntry(row["groups"].split(CSA.GROUP_DEL)) for row in bank]

    @staticmethod
    def _canonLexems(lexems: list[list[Lexem]]) -> list[str]:
        return sorted(canon(serialize(_)) for _ in lexems)

    def _generateFromPair(
        self, a: BankEntry, b: BankEntry
    ) -> Generator[list[str], None, None]:
        pairs = [(a, b), (b, a)] if self.double_cross else [(a, b)]

        for a, b in pairs:
            for crossed in islice(crossMolecules(a, b), self.max_crossed):
                yield CSA._canonLexems(crossed.lexems)

                for added in islice(addAtom(crossed), self.max_added):
                    yield CSA._canonLexems(added.lexems)
                for replaced in islice(replaceAtom(crossed), self.max_replaced):
                    yield CSA._canonLexems(replaced.lexems)
                for removed in islice(removeAtom(crossed), self.max_removed):
                    yield CSA._canonLexems(removed.lexems)
                for bonded in islice(changeBond(crossed), self.max_replaced):
                    yield CSA._canonLexems(bonded.lexems)

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
            if n == 0 and not self.double_cross:
                for groups in self._generateFromPair(entry, best):
                    generated.add(CSA.GROUP_DEL.join(groups))

        data: dict[str, list[str]] = {"SMILES": [], "groups": [], "score": []}
        for group_string in chain(
            [CSA.GROUP_DEL.join(CSA._canonLexems(best.lexems))], generated
        ):
            groups = group_string.split(CSA.GROUP_DEL)
            smiles = CSA._groupsToSmile(self.template_raws, groups)
            if isValid(smiles):
                data["groups"].append(group_string)
                data["SMILES"].append(smiles)
                data["score"].append(
                    str(CSA._bondScore(smiles)) if self.dummy_score else ""
                )

        gen_num = f"{self.round}".rjust(3, "0")
        df = pd.DataFrame(data)
        df.sort_values(by=["score"], ascending=False, inplace=True)
        df.drop_duplicates(["SMILES"])
        df.to_csv(self.population_dir / f"{gen_num}.csv", index=False)
        self._saveMeta()

    def _saveMeta(self):
        with open(self.meta_file_path, "w") as f:
            json.dump(self.meta, f)

    @staticmethod
    def _carbonScore(s: str) -> float:
        carbons = 0.0
        other = 0.0

        i = 0
        while i < len(s):
            c = s[i]
            if c not in ATOMS:
                i += 1
                continue
            if c == "c":
                carbons += 1
            elif c == "C" and i < len(s) - 1 and s[i + 1] == "l":
                other += 1
            elif c == "C":
                carbons += 1
            else:
                other += 1

            i += 1

        return carbons / (carbons + other)

    @staticmethod
    def _bondScore(s: str) -> float:
        count = 1.0
        sum = 0.0
        map = {"=": 2, "#": 4, "$": 8, ":": 1.5}

        for c in s:
            if c in BONDS:
                sum += map.get(c, 0)
                count += 1
            if c in ATOMS:
                sum += 1
                count += 1

        return sum / count

    @staticmethod
    def genInitial(
        template_raws: list[str], groups: list[list[str]], out: Path, dummy_score: bool
    ):
        data: dict[str, list[str]] = {"SMILES": [], "groups": [], "score": []}
        for gs in groups:
            smiles = CSA._groupsToSmile(template_raws=template_raws, groups=gs)
            if isValid(smiles):
                data["groups"].append(CSA.GROUP_DEL.join(gs))
                data["SMILES"].append(smiles)
                data["score"].append(str(CSA._bondScore(smiles)) if dummy_score else "")
            else:
                raise ValueError(f"Fail to create an entry for {gs}")

        df = pd.DataFrame(data)
        df.to_csv(out, index=False)

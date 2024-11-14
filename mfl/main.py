from os.path import dirname
from pathlib import Path
from rdkit import Chem, logging
from rdkit.DataStructs import TanimotoSimilarity
import numpy as np
import numpy.typing as npt
from logging import Logger
from rdkit.rdBase import DisableLog
from itertools import islice

from mfl.bank import BankEntry
from mfl.csa import CSA
from mfl.mutations import removeAtom, replaceAtom, addAtom

logger = Logger("mlf", logging.DEBUG)
DisableLog("rdApp.*")


# N = 1_000
#
#
# mols = initBank()
# for i in range(len(mols)):
#     m = mols[i]
#     print("original:", *m.smiles, sep="\n", end="\n<<<<<<\n")
#     for mutated in islice(addAtom(m), 5):
#         print(*mutated.smiles, sep="\n")
#     print("\n")
#     # for j in range(i):
#     #     for m in crossMolecules(mols[i], mols[j]):
#     #         ...

root = Path(dirname(__file__)) / ".."
template: str = "?N(?)Cc1ccc(CNC(=O)c2csc3ncNc(=O)c23)cc1"

# CSA.genInitial(
#     out=root / "populations" / "000.csv",
#     template_raws=template.split("?"),
#     groups=[
#         ["NC1CCNCC1", "NC1CCNCC1"],
#         ["C1CCNCC1", "C1CCNCC1"],
#         ["Nc1nc2cc(Br)ccc2s1", "[H-4]"],
#     ],
# )

csa = CSA(
    population_dir=root / "populations",
    round=1,
    last_generation_file_path=root / "populations" / "000.csv",
    Rd=0.98,
    DCut=1,
    template=template,
    seed_size=600,
)

csa.generateNextPopulation()

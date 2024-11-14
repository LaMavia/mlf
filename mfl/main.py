from os.path import dirname
from pathlib import Path
from rdkit import logging
from logging import Logger
from rdkit.rdBase import DisableLog

from mfl.csa import CSA

logger = Logger("mlf", logging.DEBUG)
DisableLog("rdApp.*")


root = Path(dirname(__file__)) / ".."
template: str = "?N(?)Cc1ccc(CNC(=O)c2csc3ncNc(=O)c23)cc1"

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

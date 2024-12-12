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

CSA.genInitial(
    out=root / "populations" / "000.csv",
    template_raws=template.split("?"),
    groups=[
        ["NC1CCN(C)CC1", "NC1CCN(C)CC1"],
        ["c1ccccc1C", "[H]"],
        ["C(CN)OCCO", "[H]"],
        ["CCOCOCNC", "[H]"],
        ["CCCCCCCC", "[H]"],
        ["CCCC", "[H]"],
        ["CCCCC", "[H]"],
        ["CCCCCCCCCC", "[H]"],
        ["CINCCCCCCCN", "[H]"],
        ["CINCCCCCN", "[H]"],
        ["CINCCCCN", "[H]"],
    ],
    dummy_score=True,
)

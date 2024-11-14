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
        ["CNC(=O)c1cc(=O)c(O)cn1O", "[H-3]"],
        ["CNC(=O)CCCCCNC(=O)c1cc(=O)c(O)cn1O", "[H-3]"],
        ["CCCCNC", "[H-3]"],
        ["CCCCCNC", "[H-3]"],
        ["CCOCCOCCNC", "[H-3]"],
        ["CC12CC3CC(CC(C3)C1)C2", "[H-3]"],
    ],
)

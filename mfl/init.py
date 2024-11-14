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
        # ["NC1CCNCC1", "NC1CCNCC1"],
        # ["C1CCNCC1", "C1CCNCC1"],
        # ["Nc1nc2cc(Br)ccc2s1", "[H-3]"],
        # ["NCCCCCCCNCl", "[H-3]"],
        # ["NCCOCCOCC", "[H-3]"],
        ["CNC(=O)CCCCCNC(=O)C1=CC(=O)C(=CN1O)O", "[H-3]"],
        ["CNC(=O)C1=CC(=O)C(=CN1O)O", "[H-3]"],
    ],
)

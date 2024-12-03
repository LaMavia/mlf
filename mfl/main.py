from pathlib import Path
from rdkit import logging
from logging import Logger
from rdkit.rdBase import DisableLog
from argparse import ArgumentParser

from mfl.csa import CSA

logger = Logger("mlf", logging.DEBUG)
DisableLog("rdApp.*")

parser = ArgumentParser()
parser.add_argument("-p", "--population", type=str, required=True)
parser.add_argument("-m", "--meta-path", type=str, required=True)
parser.add_argument("-d", "--out-dir", type=str, required=True)
parser.add_argument("-r", "--round", type=int, required=True)
parser.add_argument("-s", "--seed-size", type=int)
parser.add_argument("--dummy-score", action="store_true")
parser.add_argument("--double-cross", action="store_true")
parser.add_argument("--max-added", type=int, default=100)
parser.add_argument("--max-replaced", type=int, default=100)
parser.add_argument("--max-removed", type=int, default=100)
parser.add_argument("--max-crossed", type=int, default=100)
parser.add_argument("--bank-size", type=int, default=100)

args = parser.parse_args()

template: str = "?N(?)Cc1ccc(CNC(=O)c2csc3ncNc(=O)c23)cc1"

csa = CSA(
    population_dir=Path(args.out_dir),
    round=args.round,
    last_generation_file_path=Path(args.population),
    Rd=0.98,
    meta_file_path=args.meta_path,
    template=template,
    seed_size=args.seed_size,
    max_added=args.max_added,
    max_replaced=args.max_replaced,
    max_removed=args.max_removed,
    max_crossed=args.max_crossed,
    n_bank=args.bank_size,
    dummy_score=args.dummy_score,
    double_cross=args.double_cross,
)

csa.generateNextPopulation()

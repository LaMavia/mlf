from typing import Generator
from rdkit import Chem
import random as rnd

import numpy as np

from mfl.bank import BankEntry
from mfl.cross import removeEmptyBranches
from mfl.lexer import ATOM, PAREN_CLOSE, PAREN_OPEN, Lexem, depthOfLexem, serialize
from mfl.utils import isValid

MAX_MUTATIONS = 30
MAX_MUTATION_TRIES = 30
BRANCHING_CHANCE = 0.2
MUTATION_ATOMS = 500 * (50 * ["C"] + ["N", "P", "O", "S", "Cl"]) + [
    "C%sCC%s",
    "C%sCO%s",
    "C%sCN%s",
    "C%sCCC%s",
    "C%sCCO%s",
    "C%sCCN%s",
    "C%sCCCC%s",
    "c%sccoc%s",
    "c%scc[nH]c%s",
    "c%sccsc%s",
    "c%scnc[nH]%s",
    "c%scn[nH]c%s",
    "C%sCCCCC%s",
    "c%sccccc%s",
    "c%sccncc%s",
    "c%scnccc%s",
    "c%scnccn%s",
    "c%sncncn%s",
    "c%scocn%s",
    "c%snncn%s",
    "c%scscn%s",
]


def _replaceAtomAux(mol: Chem.RWMol) -> Generator[Chem.RWMol, None, None]:
    """
    Modified from https://github.com/duaibeom/MolFinder/blob/e946cacc7fd13b4f457b22b7ee16645b98464360/bin/ModSMI.py#L235-L303
    """

    #                 C   N  P   O  S   F  Cl  I
    normal_replace = [6, 15, 8, 16, 9, 17, 53]
    #               C  N  P   O  S
    arom_replace = [6, 7, 15, 8, 16]
    mutations = 0
    max_len = mol.GetNumAtoms()
    indices = list(range(max_len))
    while mutations < MAX_MUTATIONS and len(indices) > 0:
        atomi = rnd.choice(indices)
        atom = mol.GetAtomWithIdx(atomi)

        valence = atom.GetExplicitValence()
        s = atom.GetSymbol().upper()
        if s == "C" and rnd.random() < 0.99999:
            continue
        if atom.GetIsAromatic():
            if valence == 3:
                new_atomic_num_index = np.random.randint(0, 3)
            elif valence == 2:
                new_atomic_num_index = np.random.randint(1, 5)
            else:
                indices.remove(atomi)
                continue

            new_atom = Chem.Atom(arom_replace[new_atomic_num_index])
            new_atom.SetIsAromatic(True)
            mol.ReplaceAtom(atomi, new_atom)
        else:
            if valence == 4:
                new_atomic_num_index = np.random.randint(0, 1)
            elif valence == 3:
                new_atomic_num_index = np.random.randint(0, 2)
            elif valence == 2:
                new_atomic_num_index = np.random.randint(0, 4)
            elif valence == 1:
                new_atomic_num_index = np.random.randint(0, 7)
            else:
                indices.remove(atomi)
                continue

            mol.ReplaceAtom(atomi, Chem.Atom(normal_replace[new_atomic_num_index]))

        mutations += 1

        try:
            Chem.SanitizeMol(mol)
            yield mol
        except Chem.rdchem.KekulizeException:
            pass


def replaceAtom(entry: BankEntry) -> Generator[BankEntry, None, None]:
    for mutated_mols in zip(*(_replaceAtomAux(Chem.RWMol(mol)) for mol in entry.mols)):
        smis = [Chem.MolToSmiles(mol) for mol in mutated_mols]
        yield BankEntry(smis)


def removeAtom(entry: BankEntry) -> Generator[BankEntry, None, None]:
    while True:
        out_smis: list[str] = []
        for lexems, mutable_indices in zip(entry.lexems, entry.mutable_indices):
            successful = False
            if len(mutable_indices) == 0:
                return
            for _ in range(MAX_MUTATION_TRIES):
                i: int = np.random.choice(mutable_indices)
                smi = serialize(removeEmptyBranches(lexems[:i] + lexems[i + 1 :]))
                if successful := isValid(smi):
                    out_smis.append(smi)
                    break
            if not successful:
                out_smis.append(serialize(lexems))

        yield BankEntry(out_smis)


def addAtom(entry: BankEntry) -> Generator[BankEntry, None, None]:
    while True:
        out_smis: list[str] = []
        for lexems in entry.lexems:
            successful = False
            n = len(lexems)
            for _ in range(MAX_MUTATION_TRIES):
                if n == 0:
                    i = 0
                    depth = 0
                else:
                    i: int = np.random.randint(0, n)
                    depth = depthOfLexem(lexems[i])
                branch = np.random.random() < BRANCHING_CHANCE
                atom = np.random.choice(MUTATION_ATOMS)

                try:
                    d = depth + 1
                    ring_num = f"%{d}" if d >= 10 else f"{d}"
                    atom = atom % (ring_num, ring_num)
                except TypeError:
                    if len(atom) > 2:
                        print(f"Failed to format ring {atom}")

                mutation: list[Lexem] = (
                    [
                        (PAREN_OPEN, PAREN_OPEN, depth),
                        (ATOM, atom, depth),
                        (PAREN_CLOSE, PAREN_CLOSE, depth),
                    ]
                    if branch
                    else [(ATOM, atom, depth)]
                )

                smi = serialize(lexems[:i] + mutation + lexems[i:])
                if successful := isValid(smi):
                    out_smis.append(smi)
                    break

            if not successful:
                out_smis.append(serialize(lexems))

        yield BankEntry(out_smis)


def changeBond(entry: BankEntry) -> Generator[BankEntry, None, None]:
    smis: list[str] = []
    bond_types: list[Chem.BondType] = [
        Chem.BondType.SINGLE,
        Chem.BondType.DOUBLE,
        Chem.BondType.TRIPLE,
        Chem.BondType.QUADRUPLE,
    ]
    for mol in entry.mols:
        smis: list[str] = [Chem.MolToSmiles(mol)]
        indices = list(range(mol.GetNumBonds()))
        if len(indices) == 0:
            return
        rnd.shuffle(indices)
        for i in indices:
            successful = False
            modified_mol = Chem.Mol(mol, quickCopy=True)
            b = modified_mol.GetBondWithIdx(i)
            type = b.GetBondType()
            rnd.shuffle(bond_types)

            for bt in bond_types:
                if bt == type:
                    continue
                try:
                    b.SetBondType(bt)
                    Chem.SanitizeMol(modified_mol, sanitizeOps=Chem.SANITIZE_ALL)
                    smis.append(Chem.MolToSmiles(modified_mol))
                    successful = True
                except Exception:
                    pass
            if successful:
                break

    if len(smis) == len(entry.mols):
        yield BankEntry(smis)

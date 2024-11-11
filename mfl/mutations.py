from typing import Generator
from rdkit import Chem

import numpy as np

from mfl.cross import removeEmptyBranches
from mfl.lexer import ATOM, Lexem, lex, serialize, typeOfLexem
from mfl.utils import isValid

MAX_MUTATIONS = 30


def replaceAtom(lexems: list[Lexem]) -> Generator[list[Lexem], None, None]:
    """
    Modified from https://github.com/duaibeom/MolFinder/blob/e946cacc7fd13b4f457b22b7ee16645b98464360/bin/ModSMI.py#L235-L303
    """
    #                 C  B  N  P   O  S   F  Cl  Br  I
    normal_replace = [6, 5, 7, 15, 8, 16, 9, 17, 35, 53]
    #               C  N  P   O  S
    arom_replace = [6, 7, 15, 8, 16]

    mol = Chem.RWMol(Chem.MolFromSmiles(serialize(lexems)))
    max_len = mol.GetNumAtoms()

    mutations = 0
    while mutations < MAX_MUTATIONS:
        atomi = np.random.randint(0, max_len)
        atom = mol.GetAtomWithIdx(atomi)

        valence = atom.GetExplicitValence()
        if atom.GetIsAromatic():
            if valence == 3:
                new_atomic_num_index = np.random.randint(0, 3)
            elif valence == 2:
                new_atomic_num_index = np.random.randint(1, 5)
            else:
                continue

            new_atom = Chem.Atom(arom_replace[new_atomic_num_index])
            new_atom.SetIsAromatic(True)
            mol.ReplaceAtom(atomi, new_atom)
        else:
            if valence == 4:
                new_atomic_num_index = np.random.randint(0, 1)
            elif valence == 3:
                new_atomic_num_index = np.random.randint(0, 4)
            elif valence == 2:
                new_atomic_num_index = np.random.randint(0, 6)
            elif valence == 1:
                new_atomic_num_index = np.random.randint(0, 10)
            else:
                continue

            mol.ReplaceAtom(atomi, Chem.Atom(normal_replace[new_atomic_num_index]))

        mutations += 1

        try:
            Chem.SanitizeMol(mol)
            yield lex(Chem.MolToSmiles(mol))
        except Chem.rdchem.KekulizeException:
            pass


def removeAtom(lexems: list[Lexem], temp: float) -> Generator[list[Lexem], bool, None]:
    """
    :param lexems List of lexems describing the molecule
    :param temp
    """
    mutation_threshold = max(min(0.5 - 1 / temp, 1), 0)
    for i, lexem in enumerate(lexems):
        if typeOfLexem(lexem) != ATOM:
            continue
        if np.random.random() >= mutation_threshold:
            continue
        mutated_lexems = removeEmptyBranches(lexems[:i] + lexems[i + 1 :])
        if isValid(serialize(mutated_lexems)) and (yield mutated_lexems):
            lexems = mutated_lexems

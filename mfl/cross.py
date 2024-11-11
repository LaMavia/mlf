from typing import Generator
from mfl.bank import BankEntry
from mfl.lexer import PAREN_CLOSE, PAREN_OPEN, Lexem, repOfLexem, serialize, typeOfLexem
import random as rnd

from mfl.utils import canon, isValid


def removeEmptyBranches(lexems: list[Lexem]) -> list[Lexem]:
    n = len(lexems)
    for i in range(n - 1):
        left = lexems[i]
        right = lexems[i + 1]
        if repOfLexem(left) == PAREN_OPEN and repOfLexem(right) == PAREN_CLOSE:
            return lexems[:i] + lexems[i + 1 :]

    return lexems


def parenthesisSum(lexems: list[Lexem]) -> int:
    sum = 0

    for lexem in lexems:
        t = typeOfLexem(lexem)
        if t == PAREN_OPEN:
            sum += 1
        elif t == PAREN_CLOSE:
            sum -= 1

    return sum


def invalidClosingParentheses(
    lexems: list[Lexem],
) -> Generator[int, None, None]:
    sum = 0

    for i, lexem in enumerate(lexems):
        t = typeOfLexem(lexem)
        if t == PAREN_OPEN:
            sum += 1
        elif t == PAREN_CLOSE:
            sum -= 1
        if sum < 0:
            yield i
            sum = 0


def invalidOpeningParentheses(
    lexems: list[Lexem],
) -> Generator[int, None, None]:
    sum = 0

    for i in range(len(lexems) - 1, -1, -1):
        t = typeOfLexem(lexems[i])
        if t == PAREN_OPEN:
            sum += 1
        elif t == PAREN_CLOSE:
            sum -= 1
        if sum > 0:
            yield i
            sum = 0


def reduceOpeningBalancing(lexems: list[Lexem], insertions: set[int]) -> list[Lexem]:
    """
    Reduce empty branches for the case `S < 0`,
    and insert insertions.
    """

    acc: list[Lexem] = []
    for i, lexem in enumerate(lexems):
        if i not in insertions:
            acc.append(lexem)
            continue
        if typeOfLexem(lexem) == PAREN_CLOSE:
            continue
        acc.append(lexem)
        acc.append((PAREN_OPEN, PAREN_OPEN, lexem[2]))

    return acc


def reduceClosingBalancing(lexems: list[Lexem], insertions: set[int]) -> list[Lexem]:
    """
    Reduce empty branches for the case `S > 0`,
    and insert insertions.
    """

    acc: list[Lexem] = []
    for i, lexem in enumerate(lexems):
        if i not in insertions:
            acc.append(lexem)
            continue
        if typeOfLexem(lexem) == PAREN_OPEN:
            continue
        acc.append(lexem)
        acc.append((PAREN_CLOSE, PAREN_CLOSE, lexem[2]))

    return acc


def balanceBranches(lexems: list[Lexem]) -> Generator[list[Lexem], None, None]:
    # find unbalanced closing parentheses
    nm1 = len(lexems) - 1
    sum = parenthesisSum(lexems)
    if sum == 0:
        yield lexems
    if sum < 0:
        """
        Balance too many closing parenthesis
        """
        inserted_openings: set[int] = set()
        for pi in invalidClosingParentheses(lexems):
            while True:
                insertion_index = rnd.randint(0, pi)
                if insertion_index not in inserted_openings:
                    inserted_openings.add(insertion_index)
                    break

        yield reduceOpeningBalancing(lexems, inserted_openings)
    else:
        """
        Balance too many closing parenthesis
        """
        inserted_openings: set[int] = set()
        for pi in invalidOpeningParentheses(lexems):
            while True:
                insertion_index = rnd.randint(pi, nm1)
                if insertion_index not in inserted_openings:
                    inserted_openings.add(insertion_index)
                    break

        yield reduceClosingBalancing(lexems, inserted_openings)


def crossMolecules(a: BankEntry, b: BankEntry) -> Generator[BankEntry, None, None]:
    for gi, (left_lexems, right_lexems) in enumerate(
        zip(a.lexems, b.lexems, strict=True)
    ):
        ais = a.splittable_indices[gi]
        bis = b.splittable_indices[gi]
        rnd.shuffle(ais)
        rnd.shuffle(bis)
        smis: set[str] = set()

        for li, ri in zip(ais, bis):
            lexems = left_lexems[:li] + right_lexems[ri:]
            lexems = removeEmptyBranches(lexems)
            for balanced_lexems in balanceBranches(lexems):
                smi = serialize(balanced_lexems)
                if not isValid(smi):
                    # print("invalid", smi, sep="\n")
                    continue
                smis.add(canon(smi))

                # *sorted(smis)
        # print(f"len={len(smis)}", sep="\n")
        # if len(smis) == 0:
        #     smis.add(a.smiles)

    yield a

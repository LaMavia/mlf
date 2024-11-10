from typing import Literal, Callable

ATOM = 0
BOND = 1
NUM = 2
PAREN_OPEN = "("
PAREN_CLOSE = ")"

UPPER_ATOMS = ["B", "C", "N", "O", "P", "S", "F", "I", "H"]
ATOMS = UPPER_ATOMS + [c.lower() for c in UPPER_ATOMS]
BONDS = ". - = # $ : / \\".split(" ")
DIGITS = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]
MULTI_ATOMS = ["Cl", "Br"]
NUM_PREF = "%"
BRACKET_OPEN = "["
BRACKET_CLOSE = "]"

MUTABLE_LEXEMS = [ATOM]

Lexem = (
    tuple[Literal[0], str, int]  # Atom(rep, rings)
    | tuple[Literal[1], str, int]  # Bond(rep, rings)
    | tuple[Literal[2], str, int]  # Num(rep)
    | tuple[Literal["(", ")"], Literal["(", ")"], int]  # Paren
)


def repOfLexem(lexem: Lexem) -> str:
    return lexem[1]


def typeOfLexem(lexem: Lexem):
    return lexem[0]


def findBy[A](pred: Callable[[A], bool], xs: list[A]) -> A | None:
    for x in xs:
        if pred(x):
            return x
    return None


def findMultiAtom(smile: str, i: int, multi_atoms: list[str]) -> str | None:
    return findBy(lambda m: m == smile[i : i + len(m)], multi_atoms)


def isDigit(c: str) -> bool:
    return c >= "0" and c <= "9"


def lex(smile: str) -> list[Lexem]:
    lexems: list[Lexem] = []
    i = 0
    n = len(smile)
    n_current_rings = 0
    current_rings: set[str] = set()
    while i < n:
        c = smile[i]
        multi = findMultiAtom(smile, i, MULTI_ATOMS)

        if multi is not None:
            """
            Parse multi-character, bracketless atoms
            """
            lexems.append((ATOM, multi, n_current_rings))
            i += len(multi)
        elif c == BRACKET_OPEN:
            """
            Parse bracketted atoms
            """
            depth = 1
            j = i + 1
            acc = c
            while j < n and depth > 0:
                cj = smile[j]
                if cj == BRACKET_OPEN:
                    depth += 1
                elif cj == BRACKET_CLOSE:
                    depth -= 1
                j += 1
                acc += cj

            lexems.append((ATOM, acc, n_current_rings))
            i = j
        elif c in BONDS:
            """
            Parse bonds
            """
            lexems.append((BOND, c, n_current_rings))
            i += 1
        elif c == NUM_PREF:
            """
            Parse a potentially long number, starting with "%"
            """
            j = i + 1
            num_acc = NUM_PREF
            ring_depth: int
            while j < n and isDigit(cj := smile[j]):
                num_acc += cj
                j += 1
            if num_acc in current_rings:
                current_rings.remove(num_acc)
                ring_depth = n_current_rings
                n_current_rings -= 1
            else:
                current_rings.add(num_acc)
                n_current_rings += 1
                ring_depth = n_current_rings

            lexems.append((NUM, num_acc, ring_depth))
            i = j
        elif isDigit(c):
            """
            Parse a single digit
            """
            num = c
            ring_depth: int
            if num in current_rings:
                current_rings.remove(num)
                ring_depth = n_current_rings
                n_current_rings -= 1
            else:
                current_rings.add(num)
                n_current_rings += 1
                ring_depth = n_current_rings

            lexems.append((NUM, num, ring_depth))
            i += 1
        elif c == PAREN_OPEN or c == PAREN_CLOSE:
            """
            Parse parentheses
            """
            lexems.append((c, c, n_current_rings))
            i += 1
        elif c in ATOMS:
            lexems.append((ATOM, c, n_current_rings))
            i += 1
        else:
            raise ValueError(f"Unexpected char: «{c}» at position {i} in «{smile}»")

    return lexems


def serialize(tokens: list[Lexem]) -> str:
    acc = ""
    for lexem in tokens:
        acc += repOfLexem(lexem)

    return acc

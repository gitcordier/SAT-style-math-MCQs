# -*- coding: utf-8 -*-
"""
Ex_1_11.py
==========
Solver for the SAT-style maths MCQ:

    *Find two integers a ≤ b whose product is 255 and whose absolute
    difference is 2.*

Two complementary approaches are implemented:

1. **Brute-force** — exhaustive search over all integer pairs.
2. **Algebraic** — closed-form derivation using the sum/difference
   substitution ``s = a + b``, ``d = |a - b|``.

Only Python standard-library modules are used (``math``, ``pydoc``).

Example
-------
Run directly::

    python Ex_1_11.py

Or generate the HTML reference::

    python -m pydoc -w Ex_1_11.py

"""

import math
import pydoc


# ---------------------------------------------------------------------------
# Core solvers
# ---------------------------------------------------------------------------

def brute_force(product: int, diff: int) -> list[tuple[int, int]]:
    """Return all integer pairs ``(a, b)`` with ``a ≤ b`` satisfying the
    two constraints.

    The search space is ``1 ≤ a ≤ b ≤ product`` (any factor of *product*
    is bounded by *product* itself).

    Parameters
    ----------
    product:
        The required value of ``a * b``.
    diff:
        The required value of ``|a - b|``.

    Returns
    -------
    list[tuple[int, int]]
        Every ``(a, b)`` pair that satisfies both constraints, in
        ascending order of *a*.

    Examples
    --------
    >>> brute_force(255, 2)
    [(15, 17)]
    >>> brute_force(0, 0)
    []
    """
    solutions: list[tuple[int, int]] = []
    for a in range(1, product + 1):
        for b in range(a, product + 1):
            if a * b == product and abs(a - b) == diff:
                solutions.append((a, b))
    return solutions


def algebraic(product: int, diff: int) -> tuple[int, int] | None:
    """Return the unique integer pair ``(a, b)`` with ``a ≤ b`` via the
    closed-form sum/difference substitution, or ``None`` if no integer
    solution exists.

    **Derivation**

    Let ``s = a + b`` and ``d = b - a = diff`` (since ``a ≤ b``).  Then::

        a = (s - d) / 2
        b = (s + d) / 2

    Substituting into ``a * b = product``::

        (s² - d²) / 4 = product
        s² = 4 * product + d²

    So *s* is the positive integer square root of ``4 * product + d²``
    (which must be a perfect square for an integer solution to exist).

    Parameters
    ----------
    product:
        The required value of ``a * b``.
    diff:
        The required value of ``|a - b|`` (must be non-negative).

    Returns
    -------
    tuple[int, int] | None
        The pair ``(a, b)`` if an integer solution exists, otherwise
        ``None``.

    Examples
    --------
    >>> algebraic(255, 2)
    (15, 17)
    >>> algebraic(100, 3)   # no integer solution
    """
    discriminant = 4 * product + diff ** 2   # s² = 4P + d²
    s = math.isqrt(discriminant)
    if s * s != discriminant:                # not a perfect square
        return None
    if (s - diff) % 2 != 0:                 # a must be a whole number
        return None
    a = (s - diff) // 2
    b = (s + diff) // 2
    return (a, b)


# ---------------------------------------------------------------------------
# Verification helpers
# ---------------------------------------------------------------------------

def verify(a: int, b: int, product: int, diff: int) -> bool:
    """Check that the pair ``(a, b)`` satisfies both constraints.

    Parameters
    ----------
    a, b:
        The candidate pair (order does not matter for the check).
    product:
        Expected value of ``a * b``.
    diff:
        Expected value of ``|a - b|``.

    Returns
    -------
    bool
        ``True`` iff ``a ≤ b``, ``a * b == product`` and
        ``|a - b| == diff``.

    Examples
    --------
    >>> verify(15, 17, 255, 2)
    True
    >>> verify(17, 15, 255, 2)   # violates a ≤ b
    False
    >>> verify(15, 25, 255, 2)
    False
    """
    return a <= b and a * b == product and abs(a - b) == diff


def check_mcq_options(
    options: list[tuple[int, int]], product: int, diff: int
) -> dict[tuple[int, int], bool]:
    """Evaluate every MCQ option against the two constraints.

    Parameters
    ----------
    options:
        List of ``(a, b)`` candidate pairs drawn from the MCQ.
    product:
        Required product.
    diff:
        Required absolute difference.

    Returns
    -------
    dict[tuple[int, int], bool]
        Mapping from each option to ``True`` (valid) or ``False``
        (invalid).

    Examples
    --------
    >>> opts = [(15,25),(19,23),(23,25),(13,15),(15,17),(15,15),(17,19),(17,15)]
    >>> check_mcq_options(opts, 255, 2)  # doctest: +NORMALIZE_WHITESPACE
    {(15, 25): False, (19, 23): False, (23, 25): False, (13, 15): False,
     (15, 17): True,  (15, 15): False, (17, 19): False, (17, 15): False}
    """
    return {opt: verify(opt[0], opt[1], product, diff) for opt in options}


# ---------------------------------------------------------------------------
# Pretty-printing
# ---------------------------------------------------------------------------

def _section(title: str) -> None:
    """Print a section separator with *title*."""
    print(f"\n{'=' * 54}")
    print(f"  {title}")
    print(f"{'=' * 54}")


def report(product: int, diff: int) -> None:
    """Print a full pedagogical report for the given ``(product, diff)``
    problem, covering:

    * Brute-force results
    * Algebraic derivation with intermediate values
    * MCQ option-by-option verdict

    Parameters
    ----------
    product:
        The required product ``a * b``.
    diff:
        The required absolute difference ``|a - b|``.
    """
    # ---- brute force -------------------------------------------------------
    _section("1 · Brute-force search")
    found = brute_force(product, diff)
    if found:
        for a, b in found:
            print(f"  (a, b) = ({a}, {b})")
            print(f"  ✓  {a} × {b} = {a * b}    |{a} - {b}| = {abs(a - b)}")
    else:
        print("  No integer solution found.")

    # ---- algebraic ---------------------------------------------------------
    _section("2 · Algebraic derivation")
    d = diff
    discriminant = 4 * product + d ** 2
    s_float = math.sqrt(discriminant)
    s = math.isqrt(discriminant)

    print(f"""
  Constraints :  a × b = {product}
                 b - a = {d}          (a ≤ b)

  Substitution:  s = a + b,   d = b - a = {d}

  Product rule:  a × b = (s² - d²) / 4 = {product}
                 ⟹  s² = 4 × {product} + {d}² = {discriminant}
                 ⟹  s  = √{discriminant} ≈ {s_float:.4f}  →  {s}
    """)

    sol = algebraic(product, diff)
    if sol:
        a, b = sol
        print(f"  a = (s - d) / 2 = ({s} - {d}) / 2 = {a}")
        print(f"  b = (s + d) / 2 = ({s} + {d}) / 2 = {b}")
        print(f"\n  ✓  Solution: (a, b) = ({a}, {b})")
        print(f"     Memo: 255 = 2⁸ - 1 = (2⁴ - 1)(2⁴ + 1) = 15 × 17")
    else:
        print("  No integer solution (discriminant is not a perfect square).")

    # ---- MCQ verdict -------------------------------------------------------
    _section("3 · MCQ option-by-option verdict")
    labels = "ABCDEFGH"
    options = [(15,25),(19,23),(23,25),(13,15),(15,17),(15,15),(17,19),(17,15)]
    verdicts = check_mcq_options(options, product, diff)

    for label, (opt, ok) in zip(labels, verdicts.items()):
        mark = "✓  CORRECT" if ok else "✗  excluded"
        print(f"  {label}. {str(opt):<12} {mark}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    PRODUCT = 255
    DIFF    = 2

    print(__doc__)
    report(PRODUCT, DIFF)

    # Generate the pydoc HTML reference beside the script
    #pydoc.writedoc("Ex_1_11.py")
    print("\nDone!")
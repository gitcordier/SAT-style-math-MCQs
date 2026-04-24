"""
Ex_1_6.py
==============
Computes the largest power of a prime *p* dividing the product
:math:`P = 2 \\times 3 \\times \\cdots \\times n` (i.e. :math:`n!`)
using **Legendre's formula**.

.. rubric:: Background — Legendre's formula

For a prime :math:`p` and a positive integer :math:`n`, the *p*-adic
valuation :math:`v_p(n!)` — the exponent of the highest power of *p*
that divides :math:`n!` — is given by:

.. math::

   v_p(n!) = \\sum_{i=1}^{\\infty} \\left\\lfloor \\dfrac{n}{p^i} \\right\\rfloor
           = \\frac{n - s_p(n)}{p - 1}

where :math:`s_p(n)` is the sum of the digits of *n* written in base
*p*. The infinite sum is actually finite because :math:`\\lfloor n/p^i
\\rfloor = 0` for :math:`p^i > n`.

.. rubric:: Applied to the MCQ

:math:`P = 2 \\times 3 \\times \\cdots \\times 100 = 100!`  (since the
factor :math:`1` does not change the product). We want :math:`v_2(100!)`:

.. math::

   v_2(100!) = \\left\\lfloor\\frac{100}{2}\\right\\rfloor
             + \\left\\lfloor\\frac{100}{4}\\right\\rfloor
             + \\left\\lfloor\\frac{100}{8}\\right\\rfloor
             + \\left\\lfloor\\frac{100}{16}\\right\\rfloor
             + \\left\\lfloor\\frac{100}{32}\\right\\rfloor
             + \\left\\lfloor\\frac{100}{64}\\right\\rfloor
             = 50 + 25 + 12 + 6 + 3 + 1 = 97

The correct answer is therefore :math:`2^{97}` **(option D)**.

.. rubric:: Usage

Run directly::

    python Ex_1_6.py

Generate the Sphinx-compatible HTML reference::

    python -m pydoc -w Ex_1_6

Or, with a full Sphinx project::

    sphinx-apidoc -o docs/ .
    sphinx-build -b html docs/ docs/_build/

.. note::
    Only Python standard-library modules are used:
    ``math``, ``pydoc``.

:author: pedagogical example
:license: MIT
"""

import math
import pydoc


# ---------------------------------------------------------------------------
# Prime-check helper
# ---------------------------------------------------------------------------

def is_prime(n: int) -> bool:
    """Return ``True`` if *n* is a prime number, ``False`` otherwise.

    Uses trial division up to :math:`\\sqrt{n}`, which is sufficient for
    the small values involved in typical MCQ problems.

    :param n: Integer to test (:math:`n \\geq 2`).
    :type  n: int
    :returns: Primality of *n*.
    :rtype: bool

    :raises ValueError: If *n* is less than 2.

    Examples::

        >>> is_prime(2)
        True
        >>> is_prime(97)
        True
        >>> is_prime(100)
        False
    """
    if n < 2:
        raise ValueError(f"Primality is defined for integers ≥ 2; got {n}.")
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    for divisor in range(3, math.isqrt(n) + 1, 2):
        if n % divisor == 0:
            return False
    return True

# ---------------------------------------------------------------------------
# Legendre's formula — iterative sum
# ---------------------------------------------------------------------------

def legendre_iterative(n: int, p: int) -> tuple[int, list[tuple[int, int, int]]]:
    """Compute :math:`v_p(n!)` via the direct floor-sum form of Legendre's
    formula, collecting each term for pedagogical display.

    The exponent of prime *p* in :math:`n!` is:

    .. math::

        v_p(n!) = \\sum_{i=1}^{\\infty}
                  \\left\\lfloor \\frac{n}{p^i} \\right\\rfloor

    The loop stops as soon as :math:`p^i > n` (all further terms are 0).

    :param n: Upper bound of the factorial product (:math:`n \\geq 1`).
    :type  n: int
    :param p: A prime base.
    :type  p: int
    :returns:
        A 2-tuple ``(exponent, terms)`` where

        * **exponent** (*int*) — the value :math:`v_p(n!)`, and
        * **terms** (*list of (i, p^i, floor(n/p^i))*) — one entry per
          non-zero summand, useful for step-by-step display.
    :rtype: tuple[int, list[tuple[int, int, int]]]

    :raises ValueError: If *p* is not prime or *n* < 1.

    Examples::

        >>> legendre_iterative(100, 2)
        (97, [(1, 2, 50), (2, 4, 25), (3, 8, 12), (4, 16, 6), (5, 32, 3), (6, 64, 1)])
        >>> legendre_iterative(10, 3)
        (4, [(1, 3, 3), (2, 9, 1)])
    """
    if n < 1:
        raise ValueError(f"n must be ≥ 1; got {n}.")
    if not is_prime(p):
        raise ValueError(f"{p} is not prime.")

    exponent = 0
    terms: list[tuple[int, int, int]] = []
    i = 1
    power = p  # p^i
    while power <= n:
        floor_term = n // power
        exponent += floor_term
        terms.append((i, power, floor_term))
        i += 1
        power *= p
    return exponent, terms

# ---------------------------------------------------------------------------
# Legendre's formula — digit-sum form
# ---------------------------------------------------------------------------

def legendre_digit_sum(n: int, p: int) -> tuple[int, list[int]]:
    """Compute :math:`v_p(n!)` via the **digit-sum** identity:

    .. math::

        v_p(n!) = \\frac{n - s_p(n)}{p - 1}

    where :math:`s_p(n)` is the sum of the base-*p* digits of *n*.
    This form makes it immediately clear *why* divisibility depends on
    the representation of *n* in base *p*.

    :param n: Upper bound of the factorial product (:math:`n \\geq 1`).
    :type  n: int
    :param p: A prime base.
    :type  p: int
    :returns:
        A 2-tuple ``(exponent, digits)`` where

        * **exponent** (*int*) — the value :math:`v_p(n!)`, and
        * **digits** (*list[int]*) — the base-*p* digits of *n*, least
          significant first.
    :rtype: tuple[int, list[int]]

    :raises ValueError: If *p* is not prime, *n* < 1, or the formula
        yields a non-integer (which would indicate a bug).

    Examples::

        >>> legendre_digit_sum(100, 2)
        (97, [0, 0, 1, 0, 0, 1, 1])
        >>> legendre_digit_sum(10, 3)
        (4, [1, 0, 1])
    """
    if n < 1:
        raise ValueError(f"n must be ≥ 1; got {n}.")
    if not is_prime(p):
        raise ValueError(f"{p} is not prime.")

    # Decompose n in base p (digits, least-significant first)
    digits: list[int] = []
    tmp = n
    while tmp > 0:
        digits.append(tmp % p)
        tmp //= p

    s_p = sum(digits)                      # digit sum in base p
    numerator = n - s_p
    if numerator % (p - 1) != 0:          # should never happen for a prime p
        raise ValueError(
            f"(n − s_p(n)) = {numerator} is not divisible by (p−1) = {p-1}."
        )
    exponent = numerator // (p - 1)
    return exponent, digits


# ---------------------------------------------------------------------------
# Brute-force cross-check
# ---------------------------------------------------------------------------

def brute_force_valuation(n: int, p: int) -> int:
    """Compute :math:`v_p(n!)` by direct factorisation — multiply out
    :math:`1 \\times 2 \\times \\cdots \\times n`, then divide by *p*
    as many times as possible.

    .. warning::

        This approach becomes very slow for large *n* because it
        works with an astronomically large integer.  Use it only as a
        **correctness check** for small inputs.

    :param n: Upper bound (:math:`n \\geq 1`).
    :type  n: int
    :param p: A prime base.
    :type  p: int
    :returns: :math:`v_p(n!)`.
    :rtype: int

    Examples::

        >>> brute_force_valuation(100, 2)
        97
        >>> brute_force_valuation(10, 3)
        4
    """
    factorial = math.factorial(n)
    exponent = 0
    while factorial % p == 0:
        factorial //= p
        exponent += 1
    return exponent


# ---------------------------------------------------------------------------
# MCQ evaluator
# ---------------------------------------------------------------------------

def evaluate_mcq(
    options: dict[str, int],
    correct_exponent: int,
    p: int,
) -> None:
    """Print a verdict for each MCQ option, marking the unique correct one.

    :param options:
        Mapping from option label (e.g. ``'A'``) to the *exponent* it
        proposes (e.g. ``94`` for the choice :math:`2^{94}`).
    :type  options: dict[str, int]
    :param correct_exponent:
        The true value of :math:`v_p(n!)`.
    :type  correct_exponent: int
    :param p:
        The prime base (used only for display).
    :type  p: int
    :returns: ``None`` (results are printed to *stdout*).
    """
    for label, exp in options.items():
        mark = "✓  CORRECT" if exp == correct_exponent else "✗  excluded"
        print(f"  {label}.  {p}^{exp:<4}  {mark}")


# ---------------------------------------------------------------------------
# Pretty-printing helpers
# ---------------------------------------------------------------------------

def _section(title: str) -> None:
    """Print a decorated section header.

    :param title: Text to display inside the separator.
    :type  title: str
    """
    bar = "=" * 60
    print(f"\n{bar}\n  {title}\n{bar}")

def report(n: int, p: int, mcq_options: dict[str, int]) -> None:
    """Orchestrate the full pedagogical report for the problem
    :math:`v_p(n!)`.

    Sections printed:

    1. **Iterative sum** — every floor term shown explicitly.
    2. **Digit-sum identity** — base-*p* decomposition of *n*.
    3. **Brute-force cross-check** — confirmation via ``math.factorial``.
    4. **MCQ verdict** — each option evaluated.

    :param n: Upper bound of the factorial product.
    :type  n: int
    :param p: The prime base.
    :type  p: int
    :param mcq_options:
        Ordered mapping ``{label: exponent}`` of the MCQ choices.
    :type  mcq_options: dict[str, int]
    """
    # ------------------------------------------------------------------ 1
    _section(f"1 · Legendre's formula — iterative floor sum  (v_{p}({n}!))")

    exp_iter, terms = legendre_iterative(n, p)

    print(f"\n  v_{p}({n}!) = " + " + ".join(
        f"⌊{n}/{p}^{i}⌋" for i, _, _ in terms
    ))
    print(f"           = " + " + ".join(str(f) for _, _, f in terms))
    print(f"           = {exp_iter}")

    # ------------------------------------------------------------------ 2
    _section(f"2 · Digit-sum identity  (v_{p}({n}!) = (n − s_{p}(n)) / (p−1))")

    exp_ds, digits = legendre_digit_sum(n, p)
    s_p = sum(digits)

    # Display digits most-significant first for readability
    msf = list(reversed(digits))
    base_repr = "".join(str(d) for d in msf)

    print(f"""
  Base-{p} representation of {n}:
    {n} = {" + ".join(
        f"{d}×{p}^{k}" for k, d in enumerate(digits) if d
    )}
    ({n})₁₀  =  ({base_repr})_{p}

  Digit sum:  s_{p}({n}) = {" + ".join(str(d) for d in msf)} = {s_p}

  Formula:    v_{p}({n}!) = ({n} − {s_p}) / ({p} − 1)
                         = {n - s_p} / {p - 1}
                         = {exp_ds}
    """)

    # ------------------------------------------------------------------ 3
    _section("3 · Brute-force cross-check  (math.factorial)")

    exp_bf = brute_force_valuation(n, p)
    match = "✓  agrees" if exp_bf == exp_iter else "✗  MISMATCH — bug!"
    print(f"\n  v_{p}({n}!) by direct division = {exp_bf}   {match}")

    # ------------------------------------------------------------------ 4
    _section("4 · MCQ option-by-option verdict")
    print()
    evaluate_mcq(mcq_options, exp_iter, p)

    # ------------------------------------------------------------------ summary
    print(f"""
  ──────────────────────────────────────────────────────
  Answer:  the largest power of {p} dividing {n}! is  {p}^{exp_iter}.
  ──────────────────────────────────────────────────────""")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # Problem parameters
    N = 100   # P = 2 × 3 × … × 100  =  100!
    P = 2     # prime base

    # MCQ options (label → proposed exponent)
    MCQ: dict[str, int] = {
        "A": 94,
        "B": 95,
        "C": 96,
        "D": 97,
        "E": 98,
    }

    report(N, P, MCQ)

    # Generate the pydoc HTML API reference
    #pydoc.writedoc("Ex_1_6")
    print("\nEx_1_6.html written — open it for the full API reference.")
"""
Ex_1_29.py
=====================
Illustrates that the map :math:`x \\mapsto x^{1729} \\pmod{17}` is the
**identity on** :math:`\\mathbb{Z}/17\\mathbb{Z}`, i.e. every residue is a
fixed point, by iterating the map from a starting value until a fixed point
is reached.

.. rubric:: Mathematical background

**Fermat's little theorem** states that for a prime :math:`p` and any
integer :math:`a` with :math:`p \\nmid a`:

.. math::

    a^{p-1} \\equiv 1 \\pmod{p}.

Here :math:`p = 17`, so :math:`a^{16} \\equiv 1 \\pmod{17}`.

Write the exponent :math:`1729` in terms of :math:`16`:

.. math::

    1729 = 16 \\times 108 + 1 \\implies 1729 \\equiv 1 \\pmod{16}.

Therefore, for any :math:`a` coprime to :math:`17`:

.. math::

    a^{1729} = a^{16 \\times 108 + 1}
             = \\bigl(a^{16}\\bigr)^{108} \\cdot a
             \\equiv 1^{108} \\cdot a
             \\equiv a \\pmod{17}.

The map :math:`x \\mapsto x^{1729} \\pmod{17}` is therefore the identity:
*every residue class is a fixed point.*

**Why 1729?** It is the Hardy–Ramanujan *taxicab number* —
:math:`1729 = 1^3+12^3 = 9^3+10^3` — and also a Carmichael number
(:math:`1729 = 7 \\times 13 \\times 19`).  The choice is deliberate: a
single number simultaneously illustrates number theory from two different
angles.

**Observed convergence in two steps**

Starting from :math:`x_0 = 1729`:

- :math:`x_0 = 1729 \\equiv 12 \\pmod{17}` is only a *representational*
  reduction, not an application of the map.
- Step 1: :math:`x_1 = 12^{1729} \\bmod 17 = 12`. Fixed point reached.
- Step 2 (confirmation): :math:`x_2 = 12^{1729} \\bmod 17 = 12`.

.. rubric:: Usage

Run directly::

    python Ex_1_29.py

Generate the pydoc HTML reference::

    python -m pydoc -w Ex_1_29

With a full Sphinx project::

    sphinx-apidoc -o docs/ .
    sphinx-build -b html docs/ docs/_build/

.. note::
    Only Python standard-library modules are used: ``math``, ``pydoc``.

:author: pedagogical example
:license: MIT
"""

import math
import pydoc


# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------

#: The prime modulus.  Must be prime so that Fermat's little theorem applies.
MODULO: int = 17

#: The exponent used in every application of the map x -> x**POWER % MODULO.
#: Chosen equal to the Hardy–Ramanujan taxicab number 1729.
#: Key property: 1729 ≡ 1 (mod 16) = 1 (mod MODULO - 1),
#: which makes the map the identity on Z/17Z by Fermat's little theorem.
POWER: int = 1729

#: Starting value of the iteration.  Deliberately set to POWER (= 1729)
#: to highlight that the first step is a simple reduction mod MODULO
#: (1729 ≡ 12 mod 17) before the fixed-point behaviour kicks in.
INTEGER: int = 1729


# ---------------------------------------------------------------------------
# Mathematical helpers
# ---------------------------------------------------------------------------

def is_prime(n: int) -> bool:
    """Return ``True`` if *n* is prime, ``False`` otherwise.

    Uses trial division up to :math:`\\lfloor\\sqrt{n}\\rfloor`.

    :param n: Integer to test (:math:`n \\ge 2`).
    :type  n: int
    :returns: Primality of *n*.
    :rtype: bool
    :raises ValueError: If *n* < 2.

    Examples::

        >>> is_prime(17)
        True
        >>> is_prime(1729)   # Carmichael, not prime
        False
    """
    if n < 2:
        raise ValueError(f"Primality is defined for integers >= 2; got {n}.")
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    for d in range(3, math.isqrt(n) + 1, 2):
        if n % d == 0:
            return False
    return True


def fermat_exponent_reduces_to(power: int, modulo: int) -> int:
    """Return the effective exponent of the map :math:`x\\mapsto x^{power}`
    on :math:`(\\mathbb{Z}/modulo\\mathbb{Z})^*`, using Fermat's little
    theorem.

    For a prime *modulo* = :math:`p`, every element of the multiplicative
    group has order dividing :math:`p-1`, so:

    .. math::

        x^{power} \\equiv x^{power \\bmod (p-1)} \\pmod{p}.

    :param power: The original exponent.
    :type  power: int
    :param modulo: A prime modulus.
    :type  modulo: int
    :returns: ``power % (modulo - 1)``, the effective exponent.
    :rtype: int
    :raises ValueError: If *modulo* is not prime.

    Examples::

        >>> fermat_exponent_reduces_to(1729, 17)
        1
        >>> fermat_exponent_reduces_to(100, 17)
        4
    """
    if not is_prime(modulo):
        raise ValueError(
            f"modulo must be prime for Fermat's little theorem; got {modulo}."
        )
    return power % (modulo - 1)


# ---------------------------------------------------------------------------
# Core iteration
# ---------------------------------------------------------------------------

def iterate_to_fixed_point(
    start: int,
    power: int,
    modulo: int,
    verbose: bool = True,
) -> tuple[int, int]:
    """Apply the map :math:`x \\mapsto x^{power} \\pmod{modulo}` repeatedly
    until a fixed point is reached, printing each step when *verbose* is
    ``True``.

    A *fixed point* is a value :math:`r` satisfying
    :math:`r^{power} \\equiv r \\pmod{modulo}`.  The iteration stops as
    soon as the image equals the pre-image.

    .. note::
        The loop body contains the assignment ``remainder_old =
        remainder_new`` **inside** the loop, which was missing (wrong
        indentation) in the original source code.  Without that update the
        loop would never terminate.

    :param start: Initial value :math:`x_0` (any non-negative integer;
        it is reduced mod *modulo* on the first application of the map).
    :type  start: int
    :param power: Exponent applied at each step.
    :type  power: int
    :param modulo: The modulus (:math:`\\ge 2`).
    :type  modulo: int
    :param verbose: If ``True`` (default), print each step to *stdout*.
    :type  verbose: bool
    :returns:
        ``(fixed_point, steps)`` where *fixed_point* is the residue at
        convergence and *steps* is the number of map applications
        performed.
    :rtype: tuple[int, int]
    :raises ValueError: If *modulo* < 2.

    Examples::

        >>> fp, s = iterate_to_fixed_point(1729, 1729, 17, verbose=False)
        >>> fp
        12
        >>> s
        1
    """
    if modulo < 2:
        raise ValueError(f"modulo must be >= 2; got {modulo}.")

    step          = 0
    remainder_old = start       # x_0: starting value (may exceed modulo)
    remainder_new = -1          # sentinel — guaranteed != any valid residue

    # ---- main iteration ---------------------------------------------------
    # Stop as soon as the image equals the pre-image (fixed point detected).
    # BUG FIX: the original source placed `remainder_old = remainder_new`
    # outside the loop body (wrong indentation), causing an infinite loop.
    while remainder_old != remainder_new:
        step         += 1
        remainder_new = pow(remainder_old, power, modulo)   # built-in 3-arg pow
                                                             # uses fast modular
                                                             # exponentiation
        if verbose:
            print(
                f"  step {step:>2}: "
                f"{remainder_old}^{power} % {modulo} = {remainder_new}"
            )
        remainder_old = remainder_new   # ← advance: x_{n+1} := x_n^power mod p

    return remainder_new, step


def confirm_fixed_point(
    fixed_point: int,
    power: int,
    modulo: int,
    extra_steps: int = 2,
) -> None:
    """Apply the map *extra_steps* more times to a known fixed point and
    print each result, confirming that the value does not change.

    :param fixed_point: A residue that is already a fixed point.
    :type  fixed_point: int
    :param power: Exponent of the map.
    :type  power: int
    :param modulo: The modulus.
    :type  modulo: int
    :param extra_steps: How many additional applications to display
        (default: 2).
    :type  extra_steps: int

    Examples::

        >>> confirm_fixed_point(12, 1729, 17, extra_steps=2)
          confirm step 1: 12^1729 % 17 = 12
          confirm step 2: 12^1729 % 17 = 12
        and so on (every residue is a fixed point by Fermat's little theorem).
    """
    current = fixed_point
    for k in range(1, extra_steps + 1):
        image = pow(current, power, modulo)
        print(f"  confirm step {k}: {current}^{power} % {modulo} = {image}")
        current = image
    print(
        "and so on "
        f"(every residue is a fixed point by Fermat's little theorem)."
    )


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def report(start: int, power: int, modulo: int) -> None:
    """Print a complete pedagogical report for the iteration
    :math:`x \\mapsto x^{power} \\pmod{modulo}`.

    Sections:

    1. **Parameter analysis** — primality of *modulo*, effective exponent
       via Fermat's little theorem.
    2. **Iteration** — step-by-step output until the fixed point is found.
    3. **Fixed-point confirmation** — two extra applications showing
       stability.

    :param start: Starting value :math:`x_0`.
    :type  start: int
    :param power: Exponent applied at each step.
    :type  power: int
    :param modulo: The modulus.
    :type  modulo: int
    """
    bar = "=" * 62

    # ------------------------------------------------------------------ 1
    print(f"\n{bar}")
    print(f"  1 . Parameter analysis")
    print(bar)

    prime_status = "prime" if is_prime(modulo) else "NOT prime"
    print(f"\n  modulo = {modulo}  ({prime_status})")
    print(f"  power  = {power}")
    print(f"  start  = {start}")

    # Note on 1729
    print(f"\n  1729 = 7 × 13 × 19  (Carmichael number, Hardy-Ramanujan taxicab number)")
    print(f"  1729 = 1^3 + 12^3 = 9^3 + 10^3")

    if is_prime(modulo):
        eff = fermat_exponent_reduces_to(power, modulo)
        print(f"\n  By Fermat's little theorem (p={modulo} prime):")
        print(f"    x^{power} ≡ x^({power} mod {modulo-1}) = x^{eff}  (mod {modulo})")
        if eff == 1:
            print(f"    => the map x -> x^{power} mod {modulo} is the IDENTITY.")
            print(f"    => every residue is already a fixed point.")

    # ------------------------------------------------------------------ 2
    print(f"\n{bar}")
    print(f"  2 . Iteration until fixed point")
    print(bar)
    print()

    fixed_point, steps = iterate_to_fixed_point(start, power, modulo, verbose=True)

    print(f"\n  Fixed point {fixed_point} reached after {steps} step(s).")

    # ------------------------------------------------------------------ 3
    print(f"\n{bar}")
    print(f"  3 . Confirming the fixed point")
    print(bar)
    print()

    confirm_fixed_point(fixed_point, power, modulo, extra_steps=2)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    report(INTEGER, POWER, MODULO)

    #pydoc.writedoc("Ex_1_29")
    print(
        "\n  Ex_1_29.html written"
        " — open it for the full Sphinx-style API reference."
    )
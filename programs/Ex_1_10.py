"""
Ex_1_10.py
===============

Find the integer roots of the polynomial

.. math::

    P(X) = X^2 + 6X - 55

by **Newton-Raphson iteration**, then round the numerical zeros to the
nearest integer and verify them exactly.

Context
-------
The problem asks for integer pairs :math:`(a, b)` such that
:math:`a \\cdot b = 55` and :math:`a - b = -6`.

Re-writing the system with :math:`u = a` and :math:`v = -b` gives

.. math::

    \\begin{cases}
        u \\cdot v = -55 \\\\
        u + v     = -6
    \\end{cases}

which means :math:`u` and :math:`v` are the two roots of

.. math::

    P(X) = X^2 - (u+v)X + uv = X^2 + 6X - 55

The integer pairs :math:`(a, b)` are then recovered as
:math:`(u,\\, -v)` for each root :math:`v = -b`.

Algorithm
---------
Newton-Raphson starting from two distinct guesses:

.. math::

    X_{k+1} = X_k - \\frac{P(X_k)}{P'(X_k)}

where :math:`P'(X) = 2X + 6`.  The float result is rounded to the
nearest integer and re-injected into :math:`P` to verify it is an
exact zero.

Usage::

    python roots_newton.py

:author: pedagogical example
:date:   2026
"""


# ═══════════════════════════════════════════════════════════════════
# §0  Constants
# ═══════════════════════════════════════════════════════════════════

TOLERANCE = 1e-12   #: Newton convergence threshold on |P(Xₖ)|
MAX_ITER  = 200     #: safety cap on iterations


# ═══════════════════════════════════════════════════════════════════
# §1  Manual integer rounding (no math import)
# ═══════════════════════════════════════════════════════════════════

def iround(x):
    """
    Round a float to the nearest integer without importing ``math``.

    Uses the standard "round-half-up" rule:
    :math:`\\lfloor x + 0.5 \\rfloor`.

    :param x: value to round
    :type  x: float
    :returns: nearest integer
    :rtype:   int

    >>> iround(4.9999999)
    5
    >>> iround(-10.9999999)
    -11
    """
    if x >= 0:
        return int(x + 0.5)
    else:
        return -int(-x + 0.5)


# ═══════════════════════════════════════════════════════════════════
# §2  Polynomial and its derivative
# ═══════════════════════════════════════════════════════════════════

def P(X):
    """
    Evaluate the polynomial :math:`P(X) = X^2 + 6X - 55`.

    :param X: evaluation point
    :type  X: float
    :returns: :math:`P(X)`
    :rtype:   float

    >>> P(5)
    0
    >>> P(-11)
    0
    """
    return X**2 + 6*X - 55


def dP(X):
    """
    Evaluate the derivative :math:`P'(X) = 2X + 6`.

    :param X: evaluation point
    :type  X: float
    :returns: :math:`P'(X)`
    :rtype:   float

    >>> dP(5)
    16
    >>> dP(-11)
    -16
    """
    return 2*X + 6


# ═══════════════════════════════════════════════════════════════════
# §3  Newton-Raphson root finder
# ═══════════════════════════════════════════════════════════════════

def newton(x0, tol=TOLERANCE, max_iter=MAX_ITER):
    """
    Run Newton-Raphson to find a root of :math:`P`, starting from *x0*.

    .. rubric:: Iteration

    .. math::

        X_{k+1} = X_k - \\frac{P(X_k)}{P'(X_k)}

    Stops when :math:`|P(X_k)| < \\text{tol}`.

    :param x0:       initial guess
    :type  x0:       float
    :param tol:      convergence tolerance on :math:`|P(X_k)|`
    :type  tol:      float
    :param max_iter: maximum iterations allowed
    :type  max_iter: int
    :returns: ``(root_float, n_iterations)``
    :rtype:   tuple[float, int]
    :raises ValueError: if :math:`P'(X_k) = 0` (stationary point hit)
    :raises RuntimeError: if no convergence within *max_iter*

    .. note::
        :math:`P'(X) = 0` at :math:`X = -3` (the vertex of the
        parabola).  Avoid starting there or the iteration is undefined.
    """
    X = float(x0)
    print(f"\n    Starting from X₀ = {x0}")
    print(f"    {'k':>4}  {'Xₖ':>14}  {'P(Xₖ)':>16}  {'P'(Xₖ)':>12}")
    print("    " + "─" * 54)

    for k in range(max_iter):
        Pval  = P(X)
        dPval = dP(X)
        print(f"    {k:>4}  {X:>14.8f}  {Pval:>16.8f}  {dPval:>12.6f}")

        if abs(Pval) < tol:
            print(f"    ✓ Converged in {k} step(s).  Float root ≈ {X:.10f}")
            return X, k

        if abs(dPval) < 1e-14:
            raise ValueError(
                f"P'(X) ≈ 0 at X = {X:.6f} — stationary point hit. "
                "Choose a different starting guess."
            )

        X = X - Pval / dPval        # Newton update

    raise RuntimeError(
        f"Newton did not converge in {max_iter} iterations from X₀={x0}."
    )


# ═══════════════════════════════════════════════════════════════════
# §4  Rounding and exact verification
# ═══════════════════════════════════════════════════════════════════

def round_and_verify(float_root, label):
    """
    Round *float_root* to the nearest integer and verify it is an
    exact zero of :math:`P`.

    An integer :math:`n` is an exact root iff :math:`P(n) = 0`
    (integer arithmetic, no floating-point error).

    :param float_root: numerical root returned by Newton
    :type  float_root: float
    :param label:      human-readable label for display
    :type  label:      str
    :returns: ``(int_root, is_exact)``
    :rtype:   tuple[int, bool]
    """
    n   = iround(float_root)        # round to integer
    val = P(n)                      # exact integer evaluation
    ok  = (val == 0)
    status = "✓  exact zero" if ok else f"✗  P({n}) = {val} ≠ 0"
    print(f"\n    {label}: float ≈ {float_root:.10f}  →  integer {n:+d}")
    print(f"    P({n:+d}) = {n}² + 6·{n} - 55 = {n**2} + {6*n} - 55 = {val}   {status}")
    return n, ok


# ═══════════════════════════════════════════════════════════════════
# §5  Recover (a, b) pairs from roots of P
# ═══════════════════════════════════════════════════════════════════

def recover_pairs(roots):
    """
    Convert roots :math:`\\{u, v\\}` of :math:`P` into integer pairs
    :math:`(a, b)` satisfying :math:`ab = 55` and :math:`a - b = -6`.

    The algebraic reduction gives :math:`u = a` and :math:`v = -b`,
    so :math:`(a, b) = (u, -v)` for each root :math:`v = -b`.

    :param roots: list of integer roots of :math:`P`
    :type  roots: list[int]
    :returns: list of ``(a, b)`` pairs
    :rtype:   list[tuple[int, int]]
    """
    pairs = []
    for u in roots:
        for v in roots:
            if u != v:
                a, b = u, -v
                if a * b == 55 and a - b == -6:
                    pairs.append((a, b))
    return pairs


# ═══════════════════════════════════════════════════════════════════
# §6  Main
# ═══════════════════════════════════════════════════════════════════

def main():
    """
    Full pipeline:

    1. Display :math:`P` and explain the context.
    2. Run Newton from two initial guesses (one per root).
    3. Round each float root to an integer.
    4. Verify each integer is an exact zero of :math:`P`.
    5. Reconstruct the :math:`(a, b)` pairs.
    6. Final check: :math:`ab = 55` and :math:`a - b = -6`.
    """
    sep = "═" * 62

    # ── Banner ────────────────────────────────────────────────────
    print(sep)
    print("  Root finder for  P(X) = X² + 6X - 55  (Newton-Raphson)")
    print(sep)
    print("""
  Context: find integer pairs (a, b) with  ab = 55  and  a - b = -6.
  Setting  u = a,  v = -b  transforms the system into:
      u · v = -55   and   u + v = -6
  so u, v are the roots of  P(X) = X² + 6X - 55  (Vieta's formulas).
""")

    print(f"  P(X)  = X² + 6X - 55")
    print(f"  P'(X) = 2X + 6")
    print(f"  Vertex at X = -3  (avoid as starting guess!)\n")

    # ── Step 1: Newton from two different guesses ─────────────────
    print("  ── Step 1 : Newton iterations ────────────────────────")

    guesses = [
        ( 10.0, "root_1  (positive side)"),
        (-15.0, "root_2  (negative side)"),
    ]

    float_roots = []
    for x0, label in guesses:
        fr, _ = newton(x0)
        float_roots.append((fr, label))

    # ── Step 2: round and verify ──────────────────────────────────
    print("\n  ── Step 2 : Round to integer & verify ────────────────")
    int_roots = []
    all_exact = True
    for fr, label in float_roots:
        n, ok = round_and_verify(fr, label)
        int_roots.append(n)
        if not ok:
            all_exact = False

    print()
    if all_exact:
        print(f"  Both integer roots verified: {int_roots}")
    else:
        print("  WARNING: rounding produced a non-exact root — check guesses.")
        return

    # ── Step 3: recover (a, b) pairs ─────────────────────────────
    print("\n  ── Step 3 : Recover (a, b) pairs ─────────────────────")
    pairs = recover_pairs(int_roots)

    print(f"\n  {'(a, b)':<14}  {'a·b':>6}  {'a-b':>6}  {'ab=55?':>8}  {'a-b=-6?':>10}")
    print("  " + "─" * 52)
    for (a, b) in pairs:
        prod = a * b
        diff = a - b
        ok_prod = "✓" if prod == 55  else "✗"
        ok_diff = "✓" if diff == -6  else "✗"
        print(f"  ({a:+d}, {b:+d}){'':6}  {prod:>6}  {diff:>6}  {ok_prod:>8}  {ok_diff:>10}")

    # ── Step 4: summary ───────────────────────────────────────────
    print(f"\n  ── Answer ────────────────────────────────────────────")
    print(f"  P has two integer roots: {int_roots[0]} and {int_roots[1]}.")
    print(f"  The system  ab = 55,  a - b = -6  has {len(pairs)} solution(s):\n")
    for (a, b) in pairs:
        print(f"    (a, b) = ({a:+d}, {b:+d})")
    print()
    print("  → Answer G  (two solutions).")
    print()


if __name__ == "__main__":
    main()

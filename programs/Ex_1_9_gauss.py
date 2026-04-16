"""
PA_gauss.py
===========

Algebraic solver for a linear system using **Gaussian elimination
with partial pivoting** (no external math library).

The system solved here is:

.. math::

    \\begin{cases}
        a + b + c &= 450 \\\\
        a + b     &= 240 \\\\
        a - b     &= 40
    \\end{cases}

which can be written in matrix form :math:`A \\mathbf{x} = \\mathbf{b}`:

.. math::

    \\begin{pmatrix}
        1 &  1 & 1 \\\\
        1 &  1 & 0 \\\\
        1 & -1 & 0
    \\end{pmatrix}
    \\begin{pmatrix} a \\\\ b \\\\ c \\end{pmatrix}
    =
    \\begin{pmatrix} 450 \\\\ 240 \\\\ 40 \\end{pmatrix}

.. note::
    **Gaussian elimination** transforms the augmented matrix ``[A | b]``
    into an upper-triangular form, then recovers the solution via
    **back-substitution**.  Partial pivoting selects the largest
    available pivot at each step to improve numerical stability.

Usage::

    python PA_gauss.py

:author: pedagogical example
:date:   2026
"""


# ---------------------------------------------------------------------------
# Helper utilities (no import of math libraries)
# ---------------------------------------------------------------------------

def copy_matrix(M):
    """
    Return a deep copy of a 2-D list.

    :param M: original matrix (list of lists of floats)
    :type  M: list[list[float]]
    :returns: independent copy of *M*
    :rtype:   list[list[float]]

    .. rubric:: Why do we need this?

    Python lists are mutable objects.  Passing ``M`` directly and
    modifying it in-place would destroy the original data.  By working
    on a copy we keep the original intact for display purposes.
    """
    return [row[:] for row in M]


def print_matrix(M, label="Matrix"):
    """
    Pretty-print a 2-D list with a label.

    :param M:     matrix to display
    :type  M:     list[list[float]]
    :param label: title printed above the matrix
    :type  label: str
    """
    print(f"\n  {label}")
    for row in M:
        formatted = "  | " + "  ".join(f"{v:10.4f}" for v in row) + "  |"
        print(formatted)
    print()


# ---------------------------------------------------------------------------
# Step 1 - build the augmented matrix [A | b]
# ---------------------------------------------------------------------------

def build_augmented_matrix(A, b):
    """
    Assemble the augmented matrix ``[A | b]``.

    Each row of the result is ``A[i] + [b[i]]``, i.e. the coefficients
    of the equation followed by its right-hand side.

    :param A: coefficient matrix  (n X n)
    :type  A: list[list[float]]
    :param b: right-hand side vector (length n)
    :type  b: list[float]
    :returns: augmented matrix (n X (n+1))
    :rtype:   list[list[float]]

    **Example** - for a 2X2 system::

        A = [[2, 1],    b = [5,
             [1, 3]]        7]

        augmented = [[2, 1, 5],
                     [1, 3, 7]]
    """
    n = len(A)
    augmented = []
    for i in range(n):
        # Concatenate row of A with the corresponding right-hand side value
        augmented.append(A[i][:] + [b[i]])
    return augmented


# ---------------------------------------------------------------------------
# Step 2 - forward elimination (with partial pivoting)
# ---------------------------------------------------------------------------

def forward_elimination(M):
    """
    Transform the augmented matrix *M* to **upper-triangular** form
    using partial pivoting.

    .. rubric:: Algorithm

    For each column ``k`` (the *pivot column*):

    1. **Partial pivoting** - among the rows ``k, k+1, …, n-1``,
       find the one with the largest absolute value in column ``k``
       and swap it with row ``k``.  This avoids dividing by a very
       small (or zero) pivot, which would cause numerical blow-up.

    2. **Elimination** - for every row ``i > k``, subtract a multiple
       of row ``k`` from row ``i`` so that ``M[i][k]`` becomes 0:

       .. math::

           \\text{factor} = \\frac{M[i][k]}{M[k][k]}, \\quad
           \\text{row}_i \\leftarrow \\text{row}_i
                          - \\text{factor} \\times \\text{row}_k

    :param M: augmented matrix (n X (n+1)), modified **in place**
    :type  M: list[list[float]]
    :raises ValueError: if a zero pivot is encountered (singular matrix)
    """
    n = len(M)                       # number of equations

    for k in range(n):               # k = current pivot column
        # --- 2a. Partial pivoting: find the row with the max |value| ---
        max_val = abs(M[k][k])
        max_row = k
        for i in range(k + 1, n):
            if abs(M[i][k]) > max_val:
                max_val = abs(M[i][k])
                max_row = i

        # Swap row k with the max-row (if different)
        if max_row != k:
            M[k], M[max_row] = M[max_row], M[k]
            print(f"    [pivot] swapped row {k} ↔ row {max_row}")

        # Guard against a singular matrix
        if abs(M[k][k]) < 1e-12:
            raise ValueError(
                f"Zero pivot encountered at column {k}: "
                "the system has no unique solution."
            )

        print(f"    Pivot element at column {k}: {M[k][k]:.4f}")

        # --- 2b. Eliminate all entries below the pivot ---
        for i in range(k + 1, n):
            factor = M[i][k] / M[k][k]          # multiplier λ
            print(f"      row[{i}] ← row[{i}] - {factor:.4f} X row[{k}]")
            for j in range(k, n + 1):           # update all columns
                M[i][j] -= factor * M[k][j]


# ---------------------------------------------------------------------------
# Step 3 - back-substitution
# ---------------------------------------------------------------------------

def back_substitution(M):
    """
    Recover the solution vector from an upper-triangular augmented matrix.

    .. rubric:: Algorithm

    Starting from the **last** equation (which has only one unknown)
    and working upward:

    .. math::

        x_i = \\frac{M[i][n] - \\sum_{j=i+1}^{n-1} M[i][j] \\cdot x_j}
                    {M[i][i]}

    :param M: upper-triangular augmented matrix (n X (n+1))
    :type  M: list[list[float]]
    :returns: solution vector ``[x_0, x_1, …, x_{n-1}]``
    :rtype:   list[float]
    """
    n = len(M)
    x = [0.0] * n                        # initialise solution vector

    for i in range(n - 1, -1, -1):      # loop from last row upward
        # Start with the right-hand side of row i
        x[i] = M[i][n]
        # Subtract contributions of already-known unknowns
        for j in range(i + 1, n):
            x[i] -= M[i][j] * x[j]
        # Divide by the diagonal coefficient
        x[i] /= M[i][i]
        print(f"    Back-sub: x[{i}] = {x[i]:.6f}")

    return x


# ---------------------------------------------------------------------------
# Step 4 - verification
# ---------------------------------------------------------------------------

def verify(A, x, b, var_names):
    """
    Check that the computed solution satisfies every equation.

    For each equation ``i`` this function computes
    :math:`\\sum_j A[i][j] \\cdot x[j]` and compares it to ``b[i]``.

    :param A:         coefficient matrix (n X n)
    :type  A:         list[list[float]]
    :param x:         solution vector (length n)
    :type  x:         list[float]
    :param b:         right-hand side vector (length n)
    :type  b:         list[float]
    :param var_names: human-readable names for the unknowns
    :type  var_names: list[str]
    """
    print("\n  ── Verification ──────────────────────────────────")
    n = len(A)
    all_ok = True
    for i in range(n):
        lhs = sum(A[i][j] * x[j] for j in range(n))
        error = abs(lhs - b[i])
        status = "✓" if error < 1e-9 else "✗"
        print(f"  Eq.{i+1}: computed = {lhs:.6f},  expected = {b[i]:.6f},  "
              f"err = {error:.2e}  {status}")
        if error >= 1e-9:
            all_ok = False

    print()
    if all_ok:
        print("  All equations satisfied. Solution is EXACT (up to float precision).")
    else:
        print("  WARNING: one or more equations are NOT satisfied.")


# ---------------------------------------------------------------------------
# Main routine
# ---------------------------------------------------------------------------

def solve_gauss(A, b, var_names=None):
    """
    Solve the linear system :math:`A \\mathbf{x} = \\mathbf{b}` using
    Gaussian elimination with partial pivoting followed by back-substitution.

    :param A:         coefficient matrix (n X n)
    :type  A:         list[list[float]]
    :param b:         right-hand side vector (length n)
    :type  b:         list[float]
    :param var_names: optional labels for the unknowns (default: x0, x1, …)
    :type  var_names: list[str] or None
    :returns: solution vector
    :rtype:   list[float]

    .. rubric:: Summary of steps

    1. Build augmented matrix ``[A | b]``
    2. Forward elimination → upper-triangular form
    3. Back-substitution → solution vector
    4. Verification against the original system
    """
    n = len(A)
    if var_names is None:
        var_names = [f"x{i}" for i in range(n)]

    print("=" * 60)
    print("  PA - Algebraic solver: Gaussian Elimination")
    print("=" * 60)

    # Display initial system
    aug = build_augmented_matrix(A, b)
    print_matrix(aug, label="Initial augmented matrix [A | b]")

    # --- Forward elimination ---
    print("  ── Step 1 : Forward Elimination ──────────────────")
    working = copy_matrix(aug)           # work on a copy, keep aug intact
    forward_elimination(working)
    print_matrix(working, label="Upper-triangular matrix after elimination")

    # --- Back-substitution ---
    print("  ── Step 2 : Back Substitution ────────────────────")
    x = back_substitution(working)

    # --- Display solution ---
    print("\n  ── Solution ──────────────────────────────────────")
    for name, val in zip(var_names, x):
        print(f"    {name} = {val:.6f}")

    # --- Verify ---
    verify(A, x, b, var_names)

    return x


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # -----------------------------------------------------------------------
    # Define the system
    #   a + b + c = 450
    #   a + b     = 240
    #   a - b     = 40
    # -----------------------------------------------------------------------
    A = [
        [1,  1, 1],   # equation 1 : a + b + c = 450
        [1,  1, 0],   # equation 2 : a + b     = 240
        [1, -1, 0],   # equation 3 : a - b     = 40
    ]

    b = [450, 240, 40]

    var_names = ["a", "b", "c"]

    # Solve and display
    solution = solve_gauss(A, b, var_names)

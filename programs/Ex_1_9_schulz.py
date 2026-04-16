"""
PN_schulz.py
============

Numerical solver based on **Schulz's iterative matrix inversion**
(a matrix-level Newton method), applied to the system:

.. math::

    \\begin{cases}
        a + b + c &= 450 \\\\
        a + b     &= 240 \\\\
        a - b     &= 40
    \\end{cases}

Strategy
--------
Instead of solving :math:`A\\mathbf{x} = \\mathbf{b}` directly, we
**invert** :math:`A` iteratively, then recover the solution as
:math:`\\mathbf{x} = A^{-1}\\mathbf{b}`.

We write :math:`A^{-1}` column by column:

.. math::

    A^{-1} = \\begin{bmatrix} \\mathbf{X} & \\mathbf{Y} & \\mathbf{Z} \\end{bmatrix}

where :math:`A\\mathbf{X} = \\mathbf{e}_1`,
:math:`A\\mathbf{Y} = \\mathbf{e}_2`,
:math:`A\\mathbf{Z} = \\mathbf{e}_3`
(standard basis vectors).  All three columns are found **simultaneously**
by the matrix iteration described below.

The Schulz Iteration
--------------------
The scalar Newton method for computing :math:`1/a` is:

.. math::

    x_{k+1} = x_k(2 - a\\,x_k)

which doubles the number of correct digits each step.
**Schulz (1933)** lifted this to matrices: given an initial matrix
approximation :math:`X_0 \\approx A^{-1}`, iterate

.. math::

    X_{k+1} = X_k\\,(2I - A\\,X_k)

Under the convergence condition
:math:`\\|I - A X_0\\|_F < 1`, the sequence
:math:`X_k \\to A^{-1}` **quadratically** (error squares at each step).

Initialization
--------------
A safe initial guess is :math:`X_0 = \\alpha A^{\\top}` with

.. math::

    \\alpha = \\frac{1}{\\|A\\|_F^2}

because one can show :math:`\\|I - \\alpha A A^{\\top}\\|_F < 1`
whenever :math:`\\alpha \\in \\left(0,\\; 2/\\|A\\|_F^2\\right)`.
(See Ben-Israel & Cohen, 1966.)

References
----------
.. [Schulz1933]
   G. Schulz, "Iterative Berechnung der reziproken Matrix,"
   *Zeitschrift für Angewandte Mathematik und Mechanik*, **13** (1), 57–59, 1933.
   https://doi.org/10.1002/zamm.19330130111

.. [Hotelling1943]
   H. Hotelling, "Some new methods in matrix calculation,"
   *Annals of Mathematical Statistics*, **14** (1), 1–34, 1943.
   https://doi.org/10.1214/aoms/1177731489

.. [BenIsrael1966]
   A. Ben-Israel & D. Cohen, "On iterative computation of generalized
   inverses and associated projections,"
   *SIAM Journal on Numerical Analysis*, **3** (3), 410–419, 1966.
   https://doi.org/10.1137/0703035

.. [Higham2008]
   N. J. Higham, *Functions of Matrices: Theory and Computation*,
   SIAM, 2008. §6.3 "Iteration for the Matrix Inverse".
   https://doi.org/10.1137/1.9780898717778

.. [BurdenFaires]
   R. L. Burden & J. D. Faires, *Numerical Analysis*, 10th ed.,
   Brooks/Cole, 2015.  Chapter 6 (Direct Methods for Linear Systems)
   and Chapter 10 (Numerical Linear Algebra).

Usage::

    python PN_schulz.py

:author: pedagogical example
:date:   2026
"""


# ═══════════════════════════════════════════════════════════════════════════
# §0 – Hyper-parameters
# ═══════════════════════════════════════════════════════════════════════════

TOLERANCE = 1e-12   #: convergence threshold on  ‖I − A Xₖ‖_F
MAX_ITER  = 60      #: safety cap on Schulz iterations


# ═══════════════════════════════════════════════════════════════════════════
# §1 – Scalar sqrt without math import
# ═══════════════════════════════════════════════════════════════════════════

def scalar_sqrt(s, n_iter=80):
    """
    Compute :math:`\\sqrt{s}` using Newton's method for :math:`f(x) = x^2 - s`:

    .. math::

        x_{k+1} = \\frac{1}{2}\\left(x_k + \\frac{s}{x_k}\\right)

    This is Heron's method (≈ 60 BC), the oldest known root-finding
    algorithm and a special case of Newton's method.

    :param s:      non-negative radicand
    :type  s:      float
    :param n_iter: number of refinement steps (default 80, far more than needed)
    :type  n_iter: int
    :returns: :math:`\\sqrt{s}`
    :rtype:   float
    """
    if s == 0.0:
        return 0.0
    x = s                           # crude initial guess
    for _ in range(n_iter):
        x = 0.5 * (x + s / x)      # Heron's update
    return x


# ═══════════════════════════════════════════════════════════════════════════
# §2 – Matrix primitives  (all implemented from scratch)
# ═══════════════════════════════════════════════════════════════════════════

def identity(n):
    """
    Return the :math:`n \\times n` identity matrix :math:`I_n`.

    :param n: dimension
    :type  n: int
    :returns: :math:`I_n` as a list of lists
    :rtype:   list[list[float]]
    """
    return [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]


def zeros(n, m):
    """
    Return an :math:`n \\times m` zero matrix.

    :param n: number of rows
    :type  n: int
    :param m: number of columns
    :type  m: int
    :returns: zero matrix
    :rtype:   list[list[float]]
    """
    return [[0.0] * m for _ in range(n)]


def transpose(M):
    """
    Return the transpose :math:`M^{\\top}` of a matrix.

    .. math::

        (M^{\\top})_{ij} = M_{ji}

    :param M: input matrix (n × m)
    :type  M: list[list[float]]
    :returns: transposed matrix (m × n)
    :rtype:   list[list[float]]
    """
    n, m = len(M), len(M[0])
    return [[M[i][j] for i in range(n)] for j in range(m)]


def mat_add(A, B):
    """
    Element-wise matrix addition :math:`A + B`.

    :param A: first matrix (n × m)
    :type  A: list[list[float]]
    :param B: second matrix (n × m)
    :type  B: list[list[float]]
    :returns: :math:`A + B`
    :rtype:   list[list[float]]
    """
    n, m = len(A), len(A[0])
    return [[A[i][j] + B[i][j] for j in range(m)] for i in range(n)]


def mat_sub(A, B):
    """
    Element-wise matrix subtraction :math:`A - B`.

    :param A: first matrix (n × m)
    :type  A: list[list[float]]
    :param B: second matrix (n × m)
    :type  B: list[list[float]]
    :returns: :math:`A - B`
    :rtype:   list[list[float]]
    """
    n, m = len(A), len(A[0])
    return [[A[i][j] - B[i][j] for j in range(m)] for i in range(n)]


def mat_scale(s, M):
    """
    Scalar multiplication :math:`s \\cdot M`.

    :param s: scalar multiplier
    :type  s: float
    :param M: matrix (n × m)
    :type  M: list[list[float]]
    :returns: :math:`s M`
    :rtype:   list[list[float]]
    """
    return [[s * M[i][j] for j in range(len(M[0]))] for i in range(len(M))]


def mat_mul(A, B):
    """
    Standard matrix multiplication :math:`C = AB`.

    .. math::

        C_{ij} = \\sum_{k} A_{ik} \\, B_{kj}

    :param A: left matrix  (n × p)
    :type  A: list[list[float]]
    :param B: right matrix (p × m)
    :type  B: list[list[float]]
    :returns: product matrix (n × m)
    :rtype:   list[list[float]]
    :raises ValueError: if inner dimensions do not match
    """
    n, p, m = len(A), len(A[0]), len(B[0])
    if len(B) != p:
        raise ValueError(
            f"mat_mul: incompatible dimensions ({n}×{p}) × ({len(B)}×{m})"
        )
    C = zeros(n, m)
    for i in range(n):
        for j in range(m):
            C[i][j] = sum(A[i][k] * B[k][j] for k in range(p))
    return C


def mat_vec(A, v):
    """
    Matrix–vector product :math:`A \\mathbf{v}`.

    :param A: matrix (n × m)
    :type  A: list[list[float]]
    :param v: vector (length m)
    :type  v: list[float]
    :returns: vector of length n
    :rtype:   list[float]
    """
    n, m = len(A), len(A[0])
    return [sum(A[i][j] * v[j] for j in range(m)) for i in range(n)]


def frobenius_norm(M):
    """
    Compute the Frobenius norm of a matrix:

    .. math::

        \\|M\\|_F = \\sqrt{\\sum_{i,j} M_{ij}^2}

    This is the matrix analogue of the Euclidean vector norm.

    :param M: input matrix (n × m)
    :type  M: list[list[float]]
    :returns: :math:`\\|M\\|_F`
    :rtype:   float
    """
    s = sum(M[i][j] ** 2 for i in range(len(M)) for j in range(len(M[0])))
    return scalar_sqrt(s)


def col(M, j):
    """
    Extract column *j* of matrix *M* as a vector.

    :param M: matrix (n × m)
    :type  M: list[list[float]]
    :param j: column index (0-based)
    :type  j: int
    :returns: column vector of length n
    :rtype:   list[float]
    """
    return [M[i][j] for i in range(len(M))]


def print_matrix(M, label="", indent=4):
    """
    Pretty-print a matrix with an optional label.

    :param M:      matrix to display
    :type  M:      list[list[float]]
    :param label:  title printed above
    :type  label:  str
    :param indent: number of leading spaces
    :type  indent: int
    """
    sp = " " * indent
    if label:
        print(f"\n{sp}{label}")
    for row in M:
        print(sp + "│ " + "  ".join(f"{v:10.6f}" for v in row) + " │")
    print()


# ═══════════════════════════════════════════════════════════════════════════
# §3 – Schulz iteration  (the heart of the program)
# ═══════════════════════════════════════════════════════════════════════════

def schulz_init(A):
    """
    Build the initial approximation :math:`X_0 = \\alpha A^{\\top}`.

    The scalar :math:`\\alpha` is chosen as:

    .. math::

        \\alpha = \\frac{1}{\\|A\\|_F^2}

    **Why does this work?**  Write the error matrix
    :math:`E_k = I - A X_k`.  For the iteration to converge we need
    :math:`\\|E_0\\|_F < 1`.  With the choice above:

    .. math::

        \\|E_0\\|_F
        = \\|I - \\alpha A A^{\\top}\\|_F
        \\le \\left|1 - \\alpha \\sigma_{\\max}^2\\right| < 1

    since :math:`\\sigma_{\\max}^2 \\le \\|A\\|_F^2 = 1/\\alpha`.
    (Proof: Ben-Israel & Cohen, 1966, Theorem 1.)

    :param A: square matrix (n × n)
    :type  A: list[list[float]]
    :returns: tuple ``(X0, alpha)``
    :rtype:   tuple[list[list[float]], float]
    """
    norm_sq = frobenius_norm(A) ** 2       # ‖A‖_F²
    alpha   = 1.0 / norm_sq               # safe scaling factor
    At      = transpose(A)                # A^T
    X0      = mat_scale(alpha, At)        # X_0 = α A^T
    return X0, alpha


def schulz_iterate(A, tol=TOLERANCE, max_iter=MAX_ITER):
    """
    Compute :math:`A^{-1}` via the **Schulz iteration**:

    .. math::

        X_{k+1} = X_k \\,(2I - A\\,X_k)

    This is the exact matrix generalisation of Newton's scalar formula
    for :math:`1/a`:  :math:`x_{k+1} = x_k(2 - a x_k)`.

    .. rubric:: Error recursion

    Defining :math:`E_k = I - A X_k` (the "identity residual"), one
    shows that:

    .. math::

        E_{k+1} = E_k^2

    So :math:`\\|E_k\\|_F \\le \\|E_0\\|_F^{2^k}` — **quadratic convergence**:
    the number of exact digits *doubles* at each step.

    (See [Schulz1933]_, [Higham2008]_ §6.3.)

    :param A:        invertible square matrix (n × n)
    :type  A:        list[list[float]]
    :param tol:      convergence tolerance on :math:`\\|I - A X_k\\|_F`
    :type  tol:      float
    :param max_iter: maximum number of iterations
    :type  max_iter: int
    :returns: approximate inverse :math:`A^{-1}`
    :rtype:   list[list[float]]
    :raises RuntimeError: if iteration does not converge
    """
    n = len(A)
    I = identity(n)

    # ── Initialise ────────────────────────────────────────────────────
    X, alpha = schulz_init(A)

    E0_norm = frobenius_norm(mat_sub(I, mat_mul(A, X)))
    print(f"\n    α = {alpha:.8f}   (X₀ = α·Aᵀ)")
    print(f"    ‖I − A·X₀‖_F = {E0_norm:.6f}   "
          f"({'< 1 ✓ convergence guaranteed' if E0_norm < 1 else '≥ 1 ✗ divergence risk'})")

    print(f"\n    {'Iter':>4}  {'‖I − AXₖ‖_F':>16}  {'‖AXₖ − I‖²→ next':>22}")
    print("    " + "─" * 50)

    # ── Schulz loop ───────────────────────────────────────────────────
    for k in range(max_iter):
        AX     = mat_mul(A, X)           # A Xₖ
        E      = mat_sub(I, AX)          # Eₖ = I − A Xₖ
        norm_E = frobenius_norm(E)       # ‖Eₖ‖_F

        # Predicted norm after next step (quadratic: ‖E_{k+1}‖ ≈ ‖Eₖ‖²)
        predicted = norm_E ** 2
        print(f"    {k:>4}  {norm_E:>16.4e}  {predicted:>22.4e}")

        if norm_E < tol:
            print(f"\n    ✓ Converged after {k} iteration(s).")
            return X

        # Schulz update:  X_{k+1} = X_k (2I − A X_k) = X_k (I + E_k)
        two_I_minus_AX = mat_add(I, E)       # = 2I − AXₖ = I + Eₖ
        X = mat_mul(X, two_I_minus_AX)

    raise RuntimeError(
        f"Schulz iteration did not converge in {max_iter} steps. "
        "Check that A is invertible and ‖I − AX₀‖_F < 1."
    )


# ═══════════════════════════════════════════════════════════════════════════
# §4 – Verification helpers
# ═══════════════════════════════════════════════════════════════════════════

def verify_columns(A, Ainv, var_names):
    """
    Verify each column of :math:`A^{-1}` by checking
    :math:`A \\mathbf{X} = \\mathbf{e}_1`,
    :math:`A \\mathbf{Y} = \\mathbf{e}_2`,
    :math:`A \\mathbf{Z} = \\mathbf{e}_3`.

    This is the column-by-column interpretation:
    the :math:`j`-th column of :math:`A^{-1}` is the unique vector that,
    when multiplied by :math:`A`, produces the :math:`j`-th standard
    basis vector.

    :param A:         coefficient matrix (n × n)
    :type  A:         list[list[float]]
    :param Ainv:      computed inverse matrix (n × n)
    :type  Ainv:      list[list[float]]
    :param var_names: labels for the columns (e.g. ``["X","Y","Z"]``)
    :type  var_names: list[str]
    """
    n = len(A)
    I = identity(n)
    print("\n    ── Column verification: A · col_j(A⁻¹) = eⱼ ──────────")
    all_ok = True
    for j in range(n):
        c   = col(Ainv, j)               # j-th column of A⁻¹
        Ac  = mat_vec(A, c)              # A · column
        ej  = [I[i][j] for i in range(n)]
        err = frobenius_norm([[Ac[i] - ej[i]] for i in range(n)])
        status = "✓" if err < 1e-8 else "✗"
        vec_str = f"[{', '.join(f'{v:.6f}' for v in Ac)}]"
        print(f"    col {var_names[j]}: A·{var_names[j]} = {vec_str}  "
              f"(expected e_{j+1})  err={err:.2e}  {status}")
        if err >= 1e-8:
            all_ok = False
    print()
    return all_ok


def verify_product(A, Ainv):
    """
    Check that :math:`A \\cdot A^{-1} \\approx I`.

    :param A:    coefficient matrix (n × n)
    :type  A:    list[list[float]]
    :param Ainv: computed inverse (n × n)
    :type  Ainv: list[list[float]]
    """
    n   = len(A)
    I   = identity(n)
    AX  = mat_mul(A, Ainv)
    err = frobenius_norm(mat_sub(AX, I))
    print(f"    ‖A · A⁻¹ − I‖_F = {err:.4e}  "
          f"{'✓ inverse is correct' if err < 1e-8 else '✗ check failed'}")


def verify_solution(A, x, b, eq_labels, unknowns):
    """
    Re-inject the solution :math:`\\mathbf{x}` into the original system
    :math:`A\\mathbf{x} = \\mathbf{b}` and report residuals.

    :param A:          coefficient matrix (n × n)
    :type  A:          list[list[float]]
    :param x:          solution vector (length n)
    :type  x:          list[float]
    :param b:          right-hand side vector (length n)
    :type  b:          list[float]
    :param eq_labels:  human-readable equation strings
    :type  eq_labels:  list[str]
    :param unknowns:   names of the unknowns
    :type  unknowns:   list[str]
    """
    print("\n    ── Solution verification: A x = b ──────────────────")
    all_ok = True
    for i, (row, rhs, eq) in enumerate(zip(A, b, eq_labels)):
        lhs = sum(row[j] * x[j] for j in range(len(x)))
        err = abs(lhs - rhs)
        st  = "✓" if err < 1e-9 else "✗"
        print(f"    Eq.{i+1} [{eq}]:  computed = {lhs:.6f},  "
              f"expected = {rhs:.6f},  err = {err:.2e}  {st}")
        if err >= 1e-9:
            all_ok = False
    print()
    if all_ok:
        print("    All equations satisfied – solution is VERIFIED. ✓")
    else:
        print("    WARNING: solution does NOT fully satisfy the system.")


# ═══════════════════════════════════════════════════════════════════════════
# §5 – Main driver
# ═══════════════════════════════════════════════════════════════════════════

def solve_via_schulz_inversion(A, b, col_names, var_names, eq_labels):
    """
    Full pipeline: invert :math:`A` with Schulz, then solve via
    :math:`\\mathbf{x} = A^{-1} \\mathbf{b}`.

    .. rubric:: Pipeline

    1. Compute :math:`X_0 = \\alpha A^{\\top}` (safe initialisation).
    2. Run the Schulz iteration until :math:`\\|I - AX_k\\|_F < \\varepsilon`.
    3. Decompose :math:`A^{-1} = [\\mathbf{X} \\mid \\mathbf{Y} \\mid \\mathbf{Z}]`
       and verify each column satisfies :math:`A \\cdot \\text{col}_j = e_j`.
    4. Compute :math:`\\mathbf{x} = A^{-1} \\mathbf{b}`.
    5. Verify :math:`A \\mathbf{x} = \\mathbf{b}`.

    :param A:          coefficient matrix (n × n)
    :type  A:          list[list[float]]
    :param b:          right-hand side vector (length n)
    :type  b:          list[float]
    :param col_names:  labels for the columns of :math:`A^{-1}`
                       (e.g. ``["X","Y","Z"]``)
    :type  col_names:  list[str]
    :param var_names:  labels for the unknowns (e.g. ``["a","b","c"]``)
    :type  var_names:  list[str]
    :param eq_labels:  human-readable equation strings for display
    :type  eq_labels:  list[str]
    :returns: solution vector
    :rtype:   list[float]
    """
    sep = "═" * 62

    # ── Banner ─────────────────────────────────────────────────────
    print(sep)
    print("  PN – Schulz Iterative Matrix Inversion  [Schulz 1933]")
    print(sep)
    print("""
  Idea: scalar Newton for 1/a is  x_{k+1} = x_k (2 − a x_k).
  Lift to matrices:  X_{k+1} = X_k (2I − A X_k)
  Converges to A⁻¹ quadratically  (error squares each step).
  Then solve Ax = b  via  x = A⁻¹ b.
""")

    # ── Step 0: display A and b ────────────────────────────────────
    print("  ── Step 0 : Initial data ─────────────────────────────")
    print_matrix(A, label="Coefficient matrix  A", indent=4)
    print(f"    Right-hand side  b = {b}\n")

    # ── Step 1: initialise X_0 ────────────────────────────────────
    print("  ── Step 1 : Initialise  X₀ = α·Aᵀ ──────────────────")

    # ── Step 2: Schulz iterations ─────────────────────────────────
    print("\n  ── Step 2 : Schulz iterations ────────────────────────")
    Ainv = schulz_iterate(A)

    # ── Step 3: display A⁻¹ and its columns ───────────────────────
    print("\n  ── Step 3 : Inspect  A⁻¹ = [ X | Y | Z ] ───────────")
    print_matrix(Ainv, label=f"A⁻¹  (columns = {col_names})", indent=4)
    verify_columns(A, Ainv, col_names)
    print("    Sanity check:")
    verify_product(A, Ainv)

    # ── Step 4: compute x = A⁻¹ b ────────────────────────────────
    print("\n  ── Step 4 : Compute  x = A⁻¹ · b ───────────────────")
    x = mat_vec(Ainv, b)
    print(f"    x = A⁻¹ · b  =  A⁻¹ · {b}")
    for name, val in zip(var_names, x):
        print(f"      {name} = {val:.10f}")

    # ── Step 5: verify ────────────────────────────────────────────
    print("\n  ── Step 5 : Verify  A · x = b ───────────────────────")
    verify_solution(A, x, b, eq_labels, var_names)

    return x


# ═══════════════════════════════════════════════════════════════════════════
# Entry point
# ═══════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    # ── System definition ──────────────────────────────────────────────
    #   a + b + c = 450    ← row 0
    #   a + b     = 240    ← row 1
    #   a - b     = 40     ← row 2
    A = [
        [1,  1, 1],
        [1,  1, 0],
        [1, -1, 0],
    ]
    b = [450, 240, 40]

    # ── Labels ────────────────────────────────────────────────────────
    col_names = ["X", "Y", "Z"]           # columns of A⁻¹
    var_names = ["a", "b", "c"]           # unknowns
    eq_labels = [
        "a + b + c = 450",
        "a + b     = 240",
        "a - b     = 40 ",
    ]

    solution = solve_via_schulz_inversion(
        A, b, col_names, var_names, eq_labels
    )

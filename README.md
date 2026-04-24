# SAT-style math MCQs

> **[Download the PDF](SAT_STYLE_MATH_MCQS.pdf)**

SAT-style math MCQs I wrote when I was teaching.
Each exercise comes with a Python program that computes the solution.

The format is MCQ, several exercises are deliberately designed so that the fastest path to the correct answer is to eliminate
all wrong ones, a strategy explicitly discussed in the solutions. However, all solutions provdide *methods*.

Where possible, multiple proofs are given. The companion Python programs verify the analytically derived result by computation.

---

## Structure

**Five chapters, one annex. Chapters 2–5 are under review

| Chapter | Topic                                                |
| ------- | ---------------------------------------------------- |
| 1       | Arithmetic                                           |
| 2       | Linear problems                                      |
| 3       | Probability                                          |
| 4       | Geometry                                             |
| 5       | Analysis                                             |
| Annex   | Extended proofs (Legendre's formula, Taxicab number) |

---

## Chapter 1 — Arithmetic

| Exercise | Topic                     | What it teaches                                                                       |
| -------- | ------------------------- | ------------------------------------------------------------------------------------- |
| 1        | Combinatorics             | Counting pairs via structured case analysis and arithmetic series                     |
| 2        | Modular arithmetic        | Units digit of 3**100 via congruences; two independent proofs                         |
| 3        | Modular arithmetic        | Units digit of a product of even numbers                                              |
| 4        | Binary expansion          | Decomposing 625 into powers of 2                                                      |
| 5        | Euclid's lemma            | Largest power of 7 dividing P = 3 x 5 x 7 x ... × 99                                 |
| 6        | Euclid's lemma            | Largest power of 2 dividing P = 100!                                                  |
| 7        | Geometric series          | Value of 1 + 11 + 111 + ⋯ + 111…1                                                   |
| 8        | Linear systems            | Recovering (a, b, c) from mean and difference constraints                             |
| 9        | Gaussian elimination      | 3×3 system solved by row reduction, shown explicitly as matrix steps                 |
| 10       | Factorisation             | ab = 55, a − b = −6: two solutions found by two methods (arithmetic and polynomial) |
| 11       | Divisibility, elimination | ab = 255,\|a − b\| = 2: solved by progressive elimination of wrong answers           |
| 12       | Consecutive integers      | Product and mean constraints on three consecutive integers                            |
| 13       | Angles                    | Angle between clock hands at 12h05: relative angular velocity                         |
| 29       | Fermat's Little Theorem   | ((1729**1729 mod 17)**1729 mod 17)… , proved four ways (see below)      |

Exercise 29 (Taxicab) appears in Chapter 1 because it is purely arithmetic in nature.

---

## Showcase: Exercise 29 — Taxicab Number

The question asks for the value of

```
(((1729^1729 mod 17)^1729)^1729 mod 17) …
```

expressed in base 17 (heptadecimal).

**Four proofs:**

1. *Fermat's Little Theorem**
2. **Rapid exponentiation**
3. **Polynomial over Z/17Z** — using 1729 = 12³ + 1³, define
   P(X) = X³ − X + 1, then prove that P(12**17) = 0 modulo 17.
4. **Quadratic residues** — Variant of 3. Factoring P(X) = (X − 12)·Q(X), we establish that Q(0) is not a quadratic residue
   mod 7.

The full working is in the [Annex](SAT_STYLE_MATH_MCQS.pdf).

---

## Typesetting

The document is compiled with **LuaLaTeX** and uses:

- **New Computer Modern** (OTF) for text, math, and monospace
- `unicode-math` with `math-style=upright` (French academic convention)
- `polyglossia` with French as the default language
- French interval notation (`]a, b[` for open intervals) implemented via
  custom `xparse` macros with size options and error handling

The macro design treats LaTeX as code: interval commands validate their
arguments and emit errors on misuse, analogous to typed function signatures.

---

## Repository

```
SAT-style-math-MCQs/
├── SAT_STYLE_MATH_MCQS.pdf        ← compiled document (start here)
├── SAT_STYLE_MATH_MCQS.tex        ← root LaTeX file
├── FR/
│   └── chapter_1/
│       └── chapter_1.tex          ← Chapter 1 source
└── programs/
    ├── 1_1.py                     ← solver for Exercise 1
    ├── 1_2.py                     ← solver for Exercise 2
    └── …
```

---

## License

[CC0](LICENSE) — public domain. No attribution required.

from sympy import *

ANSWER = 8


def p_adic_valuation(p, n):
    """
        Compute the p-adic valuation of ``n``: the largest exponent ``i`` such
        that :math:`p^i` divides ``n``.

        :param p: A prime number greater than 1.
        :type p: int
        :param n: A non-zero integer to evaluate.
        :type n: int
        :returns: The largest ``i ≥ 0`` such that ``p**i`` divides ``n``.
        :rtype: int
        :raises ValueError: If ``p <= 1`` or ``n == 0``.

        **Example:**

        .. code-block:: python

            >>> p_adic_valuation(2, 24)   # 24 = 2^3 * 3
            3
            >>> p_adic_valuation(3, 24)   # 24 = 2^3 * 3^1
            1
            >>> p_adic_valuation(5, 24)   # 5 does not divide 24
            0
    """
    if n < 1:
        raise ValueError("n must be positive.")
    elif not isprime(p):
        raise ValueError("p must be a proven prime.")


    i = 0
    while n % p == 0:
        n //= p
        i += 1
    return i
    
def p_adic_valuation_of_product(p, factors):
    """
        Compute the p-adic valuation of a product of integers, given as an iterable
        of factors.

        Exploits the fundamental additive property of p-adic valuations:

        .. math::

            v_p\\!\\left(\\prod_{k} a_k\\right) = \\sum_{k} v_p(a_k)

        This avoids computing the (potentially enormous) product itself, working
        directly on the individual factors instead.

        :param p: A prime number present in ``PRIMES``.
        :type p: int
        :param factors: An iterable of non-zero integers whose product's valuation
                        is to be computed. Accepts any iterable, e.g. a :class:`list`,
                        :class:`tuple`, or :func:`range`.
        :type factors: iterable of int
        :returns: The p-adic valuation of the product of all elements in ``factors``,
                i.e. the largest ``i ≥ 0`` such that ``p**i`` divides the product.
        :rtype: int
        :raises ValueError: If ``p`` is not in ``PRIMES``, or if any element of
                            ``factors`` is zero (propagated from :func:`p_adic_valuation`).

        **Example:**

        .. code-block:: python

            >>> # v_7 of the product of odd numbers in [3, 99]:
            >>> # odd multiples of 7 in range: 7,21,35,49,63,77,91 → 8 factors of 7
            >>> p_adic_valuation_of_product(7, range(3, 100, 2))
            8
            >>> # v_2 of 10! = 2^8 * ...
            >>> p_adic_valuation_of_product(2, range(1, 11))
            8

        .. note::
            The additivity of :math:`v_p` means this function is exact even when the
            product itself would overflow or be computationally expensive to form.

        .. seealso::
            :func:`p_adic_valuation` — the per-factor valuation this function sums over.
    """
    valuation = 0
    for k in factors: #e.g. range(3, 100, 2):
        valuation += p_adic_valuation(p, k)
    #
    return valuation

#
if __name__ == "__main__":
    # Compute then print the result:
    computed = p_adic_valuation_of_product(7, range(3, 100, 2))
    yes_or_no = "" if computed == ANSWER else " not"

    print("Result: %d" %computed)
    print("It was%s the expected answer." %yes_or_no)


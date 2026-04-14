ANSWER = 97

"""
    Expected 2-adic valuation of the product of even numbers in ``[2, 100]``.

    By Legendre's formula: :math:`v_2(100!) = 97`, and since
    :math:`\\prod_{k=1}^{50}(2k) = 2^{50} \\cdot 50!`, we get
    :math:`v_2 = 50 + 47 = 97`.
"""

def dyadic_valuation(n):
    """
        Special case of p-adic evalution, where p=2. Can be seen as the partial function 
        p_adic_valuation(2, n). However, dyadic_valuation(n) implementation only stands on biwise operators; 
        it does not depend on any prior p_adic_valuation. 
        Compute the 2-adic (dyadic) valuation of a positive integer ``n``:
        the largest exponent ``i`` such that :math:`2^i` divides ``n``.

        Equivalently, this is the number of trailing zero bits in the binary
        representation of ``n``.

        Specialises :func:`p_adic_valuation` to ``p = 2``, allowing the use
        of efficient bitwise operations.

        :param n: A strictly positive integer to evaluate.
        :type n: int
        :returns: The largest ``i ≥ 0`` such that ``2**i`` divides ``n``.
                Returns ``0`` if ``n`` is odd.
        :rtype: int
        :raises ValueError: If ``n < 1``.

        **Bitwise identity:**

        .. math::

            v_2(n) = \\text{bit\\_length}(n \\mathbin{\\&} (-n)) - 1

        where :math:`n \\mathbin{\\&} (-n)` isolates the lowest set bit of ``n``.

        **Example:**

        .. code-block:: python

            >>> dyadic_valuation(24)   # 24 = 2^3 * 3
            3
            >>> dyadic_valuation(7)    # 7 is odd
            0
            >>> dyadic_valuation(1)    # 1 = 2^0 * 1
            0

        .. seealso::
            :func:`p_adic_valuation` — the general case for any prime ``p``.
    """
    if n < 1:
        raise ValueError("n must be positive.")
    i = 0
    while n&1 == 0 : # Pythonic alternative to the loop: (n & -n).bit_length() - 1.
        n >>= 1
        i += 1
    return i
    
def dyadic_valuation_of_product(factors):
    """
        Compute the 2-adic (dyadic) valuation of a product of positive integers,
        given as an iterable of factors.

        Applies the additive property of valuations:

        .. math::

            v_2\\!\\left(\\prod_{k} a_k\\right) = \\sum_{k} v_2(a_k)

        :param factors: An iterable of strictly positive integers. Accepts any
                        iterable, e.g. a :class:`list`, :class:`tuple`, or
                        :func:`range`.
        :type factors: iterable of int
        :returns: The 2-adic valuation of the product of all elements in
                ``factors``, i.e. the largest ``i ≥ 0`` such that ``2**i``
                divides the product.
        :rtype: int
        :raises ValueError: If any element of ``factors`` is less than 1
                            (propagated from :func:`dyadic_valuation`).

        **Example:**

        .. code-block:: python

            >>> # v₂ of the product of even numbers in [2, 100]:
            >>> dyadic_valuation_of_product(range(2, 101, 2))
            97
            >>> # v₂ of 8! = 2^7 * ...
            >>> dyadic_valuation_of_product(range(1, 9))
            7

        .. seealso::
            :func:`dyadic_valuation` — the per-factor valuation this function sums over.

            :func:`p_adic_valuation_of_product` — the general case for any prime ``p``.
    """
    valuation = 0
    for k in factors: #e.g. range(2, 101, 2):
        valuation += dyadic_valuation(k)
    #
    return valuation

#
if __name__ == "__main__":
    # Compute then print the result:
    computed = dyadic_valuation_of_product(range(2, 101, 2))
    yes_or_no = "" if computed == ANSWER else " not"

    print("Result: %d" %computed)
    print("It was%s the expected answer." %yes_or_no)


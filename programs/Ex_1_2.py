def power_N_modulo_M(k, N, M):
   """
        Compute the modular exponentiation :math:`k^N mod M` using binary (square-and-multiply)
        exponentiation.

        This algorithm runs in :math:`O(log N)` time by decomposing the exponent ``N`` in binary
        and accumulating partial results, squaring the base at each step. This is significantly
        faster than naïve repeated multiplication, which would run in :math:`O(N)`.

        :param k: The base of the exponentiation.
        :type k: int
        :param N: The exponent. Must be a non-negative integer.
        :type N: int
        :param M: The modulus. Must be a non-zero positive integer.
        :type M: int
        :returns: The value of :math:`k^N mod M`, an integer in the range ``[0, M-1]``.
        :rtype: int
        :raises ValueError: If ``M`` is zero or if ``N`` is negative.

        **Algorithm outline:**

        At each iteration, if the current exponent ``n`` is odd, the current base is folded into
        the running ``result``. The base is then squared and the exponent halved, repeating until
        only one factor remains.

        **Example:**

        .. code-block:: python

            >>> power_N_modulo_M(3, 100, 10)
            1
            >>> power_N_modulo_M(2, 10, 1000)
            24

        .. note::
            Python's built-in :func:`pow` with three arguments (``pow(k, N, M)``) implements
            the same algorithm natively in C and should be preferred in production code.

        .. seealso::
            `Modular exponentiation <https://en.wikipedia.org/wiki/Modular_exponentiation>`_
    """
    if M == 0:
        raise ValueError("Modulus M must be non-zero.")
    if N < 0:
        raise ValueError("Exponent N must be non-negative.")
    if N == 0:
        return 1

    result = 1
    base = k % M
    n = N

    while n > 1:
        if n & 1:          # n is odd
            result = (result * base) % M
            n -= 1         # make n even
        base = (base * base) % M
        n >>= 1            # divide by 2

    return (base * result) % M


if __name__ == "__main__":
    answer = pow(3, 100, 10) # Python built-in
    computed = power_N_modulo_M(3, 100, 10)
    yes_or_no = "" if computed == answer else " not"
    print("Result: %d" %computed)
    print("It was%s the expected answer." %yes_or_no)


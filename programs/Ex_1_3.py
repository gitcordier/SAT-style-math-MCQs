ANSWER = 6

def P_mod_10():
    """
        Compute the last digit of the product of even numbers in the range ``[2, 98]``,
        excluding multiples of 10.

        Formally, this computes:

        .. math::

            P = \\prod_{\\substack{n=2,4,\\ldots,98 \\\\ 10 \\nmid n}} n \\pmod{10}

        Multiples of 10 are excluded because any such factor would immediately force the
        product to zero modulo 10, making the result trivially 0 regardless of the other
        terms.

        :returns: The last decimal digit of the filtered even product, i.e. an integer in
                the range ``[0, 9]``. The expected result is ``6``.
        :rtype: int

        **Mathematical justification:**

        The qualifying even numbers repeat the pattern ``{2, 4, 6, 8}`` modulo 10 in each
        decade. Their product within one decade is :math:`2 \\times 4 \\times 6 \\times 8 = 384
        \\equiv 4 \\pmod{10}`. There are 10 such decades in ``[2, 98]``, so the total is
        :math:`4^{10} = 1048576 \\equiv 6 \\pmod{10}`.

        **Example:**

        .. code-block:: python

            >>> P_mod_10()
            6
    """
    p = 1
    for n in range(2, 100, 2):
        if n % 10 == 0:
            continue
        p = (p * n) % 10 # Avoid large integers!
    return p


if __name__ == "__main__":
    digit = P_mod_10()
    yes_or_no = "" if digit == ANSWER else " not"
    print("Result: %d" % digit)
    print("It was%s the expected answer." % yes_or_no)


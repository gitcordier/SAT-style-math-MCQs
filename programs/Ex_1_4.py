ANSWER = 5 

def hamming_weight_binary(n):
    """
        Return the Hamming weight of n with respect to its binary expansion.

        The Hamming weight is the number of symbols that are different from the 
        zero-symbol of the alphabet used. In this binary case, it is the number 
        of '1' bits.

        :param n: The integer to analyze.
        :type n: int
        :return: The total count of set bits (1s) in the binary representation.
        :rtype: int
    """
    return bin(n).count('1')


if __name__ == "__main__":
    weight = hamming_weight_binary(625)
    yes_or_no = "" if weight == ANSWER else " not"
    print("Result: %d" % weight)
    print("It was%s the expected answer." % yes_or_no)

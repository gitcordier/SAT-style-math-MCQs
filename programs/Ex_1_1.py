ANSWER_FORMULA = (9999 - 1099 + 1) * 90 + (89 * 90) / 2
ANSWER_RAW     = 8945.5 * 90

def count_pairs():
    """
        Count pairs (N, n) that simultaneously satisfy:

        - 1000 ≤ N ≤ 9999
        - 1000 ≤ n ≤ 9999
        - 10 ≤ N - n ≤ 99

        :return: The number of matching pairs.
        :rtype: int
    """


    count = 0  # running total of matching pairs
    for N in range(1000, 10000):       # 1000 ≤ N ≤ 9999
        for n in range(1000, 10000):   # 1000 ≤ n ≤ 9999
            if 10 <= N - n <= 99:
                count += 1
            #
        #
    #

    return count


if __name__ == "__main__":
    # Compute then print the result:
    count = count_pairs()
    yes_or_no = "" if count == ANSWER_FORMULA == ANSWER_RAW else " not"

    print("Result: %d" %count)
    print("It was%s the expected answer." %yes_or_no)



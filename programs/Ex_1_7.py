import time
import pandas as pd

def I(N):
    """
    Compute a repunit-like integer consisting of N ones in base 10.
    
    This naive approach uses a loop to compute the sum of powers of 10.
    For example, I(3) = 10^0 + 10^1 + 10^2 = 111.

    :param N: The number of terms to sum.
    :type N: int
    :return: The computed integer sum.
    :rtype: int
    """
    x = 0
    # Naive summation loop (O(N) time complexity)
    for k in range(0, N):
        x += pow(10,k)
    return x

def S(N):
    """
    Compute the sum of I(k) for k from 1 to N.
    
    This is the naive, brute-force computation method utilizing a generator expression.
    
    :param N: The upper limit of the sequence.
    :type N: int
    :return: The sum of the sequence up to N.
    :rtype: int
    """
    # Sums the results of I(k) iteratively (O(N^2) overall time complexity)
    return sum(I(k) for k in range(1, N+1))

def F(N):
    """
    Compute the sum of I(k) for k from 1 to N using a closed-form mathematical formula.
    
    This represents an optimized, direct calculation avoiding loops. 
    The mathematical derivation relies on geometric series properties.
    
    :param N: The upper limit of the sequence.
    :type N: int
    :return: The mathematically evaluated sum of the sequence (0 if N ≤ 0).
    :rtype: int
    """
    if N  < 0:
        return 0
    # Calculate the geometric component using integer division to avoid float overflow
    geo = (pow(10, N+1) - 1) // 9
    
    # Calculate the arithmetic component
    ari = N + 1
    
    # Final formula evaluated using integer division
    return (geo - ari) // 9 

def compare_execution_time(f, g, N):
    """
    Compare the execution time of two functions for a given input N.
    
    Runs both functions, verifies that they produce the exact same output, 
    and calculates the performance metrics.

    :param f: The first function to evaluate (usually the optimized formula).
    :type f: callable
    :param g: The second function to evaluate (usually the naive computation).
    :type g: callable
    :param N: The input argument to pass to both functions.
    :type N: int
    :raises AssertionError: If the output of function `f` does not match function `g`.
    :return: A dictionary containing the input N, the validated result, execution times, 
             and the performance quotient.
    :rtype: dict
    """
    # Time the optimized function
    start_time = time.perf_counter()
    result_f = f(N)
    end_time = time.perf_counter()
    time_f = end_time - start_time

    # Time the naive function
    start_time = time.perf_counter()
    result_g = g(N)
    end_time = time.perf_counter()
    time_g = end_time - start_time

    # Verify both methods yield the exact same mathematical result
    try:
        assert result_f == result_g
    except AssertionError as e:
        # Use string formatting to output the mismatching values before raising the error
        print("f=%f" % result_f)
        print("g=%f" % result_g)
        raise e
        
    return dict(
        N=N,
        result=result_f, 
        time_f=time_f, 
        time_g=time_g, 
        quotient=time_f / time_g
    )

def record(L):
    """
    Generate a benchmark record comparing functions F and S over a range of inputs.

    :param L: The upper boundary (exclusive) for the range of N to test.
    :type L: int
    :return: A list of dictionaries containing execution metrics for each N.
    :rtype: list of dict
    """
    # List comprehension to build the comparison dataset
    return [compare_execution_time(F, S, k-12) for k in range(0, L)]

if __name__ == "__main__":
    # Generate data for N=0 up to N=15
    benchmark_data = record(16) 
    print(benchmark_data)

    # Convert the list of dictionaries into a Pandas DataFrame
    df = pd.DataFrame(benchmark_data)
    
    # Clean up display precision for timing columns
    df["time_f"] = df["time_f"].round(8)
    df["time_g"] = df["time_g"].round(8)
    
    # Vectorized multiplication to convert quotient to a percentage
    df["quotient"] = df["quotient"].mul(100)
    df["quotient"] = df["quotient"].round(2)

    # Safely cast the integer results using Pandas Nullable Integer type
    df["result"] = df["result"].astype("Int64")

    # Rename columns for clarity in the final report
    df = df.rename(columns={
        "quotient": "quotient time(F)/time(G), in %", 
        "time_f": "with the formula", 
        "time_g": "naive computation"
    })
    
    # Export to CSV without the DataFrame index
    df.to_csv("Ex_1_7_comparison.csv", index=False)
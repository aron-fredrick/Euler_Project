from math import floor
import math
import numpy as np


def find_factors(n):
    # Calculates all the factors for a given n

    factors = []
    for ii in range(1, floor(np.sqrt(n))+1):
        if n%ii == 0:
            factors.append(ii)
            if ii*ii != n:
                factors.append(int(n/ii))
    factors.sort()
    return factors


def isprime(n):
    # checks to see if a number, n, is prime
    prime = True
    if n >= 2:

        for ii in range(2, int(np.ceil(np.sqrt(n)))):
            if n % ii == 0:
                prime = False
                break
    else:
        prime = False
    return prime


def primes(N_max):
    # returns all primes below N_max
    prime_numbers = [2]
    test = True
    N = 3

    while prime_numbers[-1] < N_max:
        isprime = True
        for prime in prime_numbers:
            if N % prime == 0:
                isprime = False
                break
        if isprime:
            if N > N_max:
                break
            prime_numbers.append(N)
        N = N + 2
    return prime_numbers

def get_primes(n):
    # same as primes, but faster
    m = n+1
    numbers = [True] * m
    for i in range(2, int(n**0.5 + 1)):
        if numbers[i]:
            for j in range(i*i, m, i):
                numbers[j] = False
    primes1 = []
    for i in range(2, m):
        if numbers[i]:
            primes1.append(i)
    return primes1

def prime_factors(n):
    # returns all the unique prime factors of a number n
    prime_numbers = [2]
    prime_factors = []
    test = True
    N = 3
    bignumber = n
    upper = np.ceil(bignumber)
    if n%2 == 0:
        prime_factors.append(2)
    while N<=upper:
        isprime = True
        for prime in prime_numbers:
            if N%prime == 0:
                isprime = False
                break
        if isprime:
            prime_numbers.append(N)
            if bignumber%N==0:
                prime_factors.append(N)
                bignumber = bignumber/N
                upper = np.ceil(bignumber)
                N = N - 2
        N = N + 2
    return prime_factors


def primeFactors(n):
    # returns all the unique prime factors of a number n
    prime_factors = []
    # Print the number of two's that divide n
    while n % 2 == 0:
        if 2 not in prime_factors:
            prime_factors.append(2)
        n = n / 2

    # n must be odd at this point
    # so a skip of 2 ( i = i + 2) can be used
    for i in range(3, int(n ** (1 / 2)) + 1, 2):

        # while i divides n , print i amd divide n
        while n % i == 0:
            if i not in prime_factors:
                prime_factors.append(i)
            n = n / i

    # Condition if n is a prime
    # number greater than 2
    if n > 2:
        prime_factors.append(int(n))

    return prime_factors


def isPalindrome(inString):
    # checks to see if a string is a palindrome
    answer = True

    for ii in range(0, math.ceil(len(inString) / 2)):
        if inString[ii] != inString[len(inString) - ii - 1]:
            answer = False
            break

    return answer


def b10tob2(n):
    # converts base 10 to binary
    answer = ''
    while n > 0:
        r = n % 2
        n = int((n - r) / 2)
        answer = str(r) + answer

    return answer


def next_prime(prime_numbers=[]):
    # returns the next prime in the list, assumes the list that is input is correct and ordered
    if len(prime_numbers) <= 1:
        prime_numbers = [2]
        N = 3
    else:
        N = prime_numbers[-1] + 2
    test = True
    assert N % 2 == 1

    while test:
        isprime = True
        for prime in prime_numbers:
            if N % prime == 0:
                isprime = False
                break
        if isprime:
            prime_numbers.append(N)
            test = False
        N = N + 2
    return prime_numbers

def permutes(num_str):
    # returns all the permutations of a string of characters
    p_list = []
    if len(num_str)==1:
        return [num_str]
    if len(num_str)==2:
        return [num_str , num_str[::-1]]
    for ii in range(0,len(num_str)):
        p_list = p_list+[num_str[ii]+s for s in permuts(num_str[0:ii]+num_str[ii+1:])]
    return p_list


def is_permute(str1, str2):
    # checks to see if 2 strings are permutations of each other
    l1 = sorted(str1)
    l2 = sorted(str2)

    if l1 == l2:
        return True
    else:
        return False

def factorial(N):
    # calculates N!
    if N==0:
        return 1
    if N==1:
        return 1
    else:
        return N*factorial(N-1)


def gcd(a, b):
    # calculate greatest common divisor of a and b
    if a == b:
        return a
    elif a < b:
        tempa = b
        b = a
        a = tempa

    rem = a % b

    if rem == 0:
        return b
    else:
        return gcd(b, rem)


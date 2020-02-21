# Author: Silvio Gregorini (silviogregorini@openforce.it)

from math import ceil, sqrt


def get_n_primes(n):
    if n <= 4:
        return {
            0: [],
            1: [2],
            2: [2, 3],
            3: [2, 3, 5],
            4: [2, 3, 5, 7],
        }.get(n, [])
    primes = [2, 3, 5, 7]
    num = 11
    while len(primes) < n:
        if all(num % div != 0 for div in range(2, ceil(sqrt(num) + 1))):
            primes.append(num)
        num += 1
    return primes


def get_nth_prime(n):
    if n <= 0:
        raise ValueError
    return get_n_primes(n)[-1]

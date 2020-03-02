# Author: Silvio Gregorini (silviogregorini@openforce.it)

from math import ceil, sqrt


ALL_PRIMES = {2}


def factorize(num):
    factors = {}
    primes = get_primes_up_to(ceil(sqrt(num)))
    while primes and num != 1:
        prime = primes.pop(0)
        while num % prime == 0:
            num = num / prime
            if prime not in factors:
                factors[prime] = 0
            factors[prime] += 1
    if num != 1:
        factors[int(num)] = 1
    return factors


def get_n_primes(num):
    if num <= 0:
        return []
    primes = []
    while len(primes) < num:
        primes.append(get_next_prime(primes[-1] if primes else 1))
    return primes


def get_next_prime(num):
    num += 1
    if not is_prime(num):
        return get_next_prime(num)
    return num


def get_primes_up_to(num):
    if num <= 1:
        return []
    primes = [2]
    while primes[-1] < num:
        primes.append(get_next_prime(primes[-1]))
    return primes


def is_prime(num):
    if num <= 1:
        raise ValueError
    elif num in ALL_PRIMES:
        return True
    for x in range(2, ceil(sqrt(num)) + 1):
        if num % x == 0:
            return False
    ALL_PRIMES.update([num])
    return True

import pickle
from functools import reduce
from operator import mul

def unpickl(filename):
    with open(filename, 'rb') as f:
        unpickler = pickle.Unpickler(f)
        content = unpickler.load()
    return content

def pickl(filename, content):
    with open(filename, 'wb') as f:
        pickler = pickle.Pickler(f)
        pickler.dump(content)

def filename(basename, q):
    try:
        return basename % q
    except TypeError:
        return basename

primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]

def factorize(n):
    n_orig = n
    if n < 0:
        return [(-1,1)] + factorize(-n)
    if n == 0:
        raise ValueError('0 cannot be factorized')
    result = []
    for prime in primes:
        if n == 1:
            break
        num = 0
        while n % prime == 0:
            num = num + 1
            n //= prime
        if num > 0:
            result.append((prime,num))
    assert reduce(mul, (p**e for p,e in result), 1) == n_orig
    return result
    

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# _authors_: Vozec
# _date_ : 31/10/2022
# edited by: sindreoyen, oct 2024

from utils.logger import found_x
import gmpy2 as gmpy
from gmpy2 import mpz

def egcd(a, b):
    if a == 0:
        return (b, mpz(0), mpz(1))
    else:
        g, y, x = egcd(gmpy.f_mod(b, a), a)
        return (g, x - (b // a) * y, y)

def modInv(a, m):
    # This function finds the modular inverse using the extended gcd function
    _, x, _ = egcd(mpz(a), mpz(m))
    return gmpy.f_mod(x, m)

def crt(Ni, Mi):
    # Calculate the product of all moduli
    M = mpz(1)
    for m in Mi:
        M *= mpz(m)
    
    # Compute the solution for the system of congruences
    X = mpz(0)
    for m, n in zip(Mi, Ni):
        m = mpz(m)
        n = mpz(n)
        temp_mi = M // m
        X += temp_mi * n * modInv(temp_mi, m)
    
    X = gmpy.f_mod(X, M)
    return X

def run(stop, name, g, h, p, all_factor):
    try:
        x = []

        for factor in all_factor:
            if stop.is_cancelled:
                return None

            # Convert necessary values to mpz for gmpy2 compatibility
            factor_0 = mpz(factor[0])
            factor_1 = factor[1]
            p = mpz(p)
            g = mpz(g)
            h = mpz(h)
            res = {}

            # Precompute modular powers for current factor
            for i in range(int(factor_0)):
                if stop.is_cancelled:
                    return None
                res[gmpy.powmod(g, (p - 1) // factor_0 * i, p)] = i

            c_i = []
            h_ = h

            # Compute the coefficients for the current factor
            for j in range(factor_1):
                if stop.is_cancelled:
                    return None
                tp = gmpy.powmod(h_, (p - 1) // pow(factor_0, j + 1), p)
                c_i.append(res[tp])
                
                # Update h_ based on the inverse calculation
                a = gmpy.powmod(g, res[tp] * pow(factor_0, j), p)
                a = modInv(a, p)
                h_ = gmpy.f_mod(h_ * a, p)

            # Sum up the results with the current coefficients
            x.append(sum(gmpy.mul(pow(factor_0, j), c_i[j]) for j in range(factor_1)))

        # Compute the moduli list for CRT
        factors = [pow(mpz(factor[0]), factor[1]) for factor in all_factor]
        return found_x(stop, name, crt(x, factors))

    except Exception as e:
        stop.cancel()
        print(f"An error occurred: {e}")

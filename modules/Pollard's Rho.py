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
        return b, mpz(0), mpz(1)
    else:
        g, y, x = egcd(gmpy.f_mod(b, a), a)
        return g, x - gmpy.f_div(b, a) * y, y

def Xab(x, a, b, G, H, P, Q):   
    sub = gmpy.f_mod(x, 3)
    if sub == 0:
        return gmpy.f_mod(x * G, P), gmpy.f_mod(a + 1, Q), b
    elif sub == 1:
        return gmpy.f_mod(x * H, P), a, gmpy.f_mod(b + 1, Q)
    elif sub == 2:
        return gmpy.f_mod(x * x, P), gmpy.f_mod(a * 2, Q), gmpy.f_mod(b * 2, Q)
    return x, a, b

def run(stop, name, g, h, p, factors):
    if not gmpy.is_prime(p):
        return None

    q = mpz((p - 1) // 2)
    g, h, p = mpz(g), mpz(h), mpz(p)
    x, a, b = gmpy.f_mod(g * h, p), mpz(1), mpz(1)
    X, A, B = x, a, b

    for i in range(1, int(p)):
        if stop.is_cancelled:
            return None
        x, a, b = Xab(x, a, b, g, h, p, q)
        X, A, B = Xab(*Xab(X, A, B, g, h, p, q), g, h, p, q)  # Triple iteration

        if x == X:
            break

    gcd_result, inverse_b_minus_B, _ = egcd(B - b, q)
    if gcd_result != 1:
        return None

    x = gmpy.f_mod(inverse_b_minus_B * (a - A), q)
    if gmpy.powmod(g, x, p) == h:
        return found_x(stop, name, int(x))
    else:
        return found_x(stop, name, int(x + q))

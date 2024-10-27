#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# _authors_: Vozec
# _date_ : 31/10/2022
# edited by: sindreoyen, oct 2024

from utils.logger import found_x
from sympy import ceiling, isprime, sqrt
import gmpy2 as gmpy
from gmpy2 import mpz

def run(stop, name, g, h, p, _):	
	if(not isprime(p)):
		return None

	n = mpz(ceiling(sqrt(p-1)))
	
	for i in range(n):
		if(stop.is_cancelled):
			return None
		ref = { gmpy.powmod(g, i, p): i }

	c = gmpy.powmod(g, n * (p - 2), p)

	for j in range(n):
		if(stop.is_cancelled):return None

		val = gmpy.f_mod(mpz(gmpy.powmod(c, j, p) * h), p)
		if val in ref:
			x = ref[val] + n * j
			return found_x(stop,name,x)

	return None


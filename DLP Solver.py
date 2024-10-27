#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# _authors_: Vozec
# _date_ : 31/10/2022

import threading
from utils.utils  import *
from utils.logger import *
from gmpy2 import mpz


class thread_stop:
	def __init__(self):
		self.is_cancelled = False
	def cancel(self):
		self.is_cancelled = True

def main(args):
	g, h, p = mpz(args.g), mpz(args.h), mpz(args.p)
	logger('Trying to solve the discrete log problem : %s = %s^x mod %s'%(h,g,p),'',1,1,True)

	all_algo = Load_modules('./modules',globals())

	factors = Get_all_factors(p-1)
	logger('Factors of p-1 : %s\n'%(factors),'',1,1,True)

	all_thread = []
	stopper = thread_stop()

	for name, algo in all_algo.items():		
		logger('[+] Starting %s algorithm '%name,'info',0,0)
		t1 = threading.Thread(target=algo.run, args = (stopper,name,g,h,p,factors))
		all_thread.append(t1)

	[t.start() for t in all_thread]
	[t.join() for t in all_thread]

	stopper.cancel()

if __name__ == '__main__':
	header()
	main(parse_args())

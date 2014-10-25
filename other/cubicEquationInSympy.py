# -*- coding: utf-8 -*-
"""
Created on Wed Jan 08 20:37:40 2014

@author: hakan
"""

import sympy as sp

t,a0,a1,a2,a3,x = sp.symbols('t a0 a1 a2 a3 x')

x = a0+a1*t+a2*t**2+a3*t**3

sp.solve(x,t)
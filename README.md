hobj
====

A Python library for computing matching and independence polynomials


Website: https://github.com/pernici/hobj

Author: Mario Pernici <mario.pernici@gmail.com>
        INFN - Sezione di Milano, Italy

hobj is free software released under the New BSD License (see the
LICENSE file for details)

Dependences
===========

There are no dependences.

One can use gmpy setting GMP = 1 in src/domains.py,
but usually the pure Python version is faster
(dup_permanental_minor_poly(m, QQ) is faster with gmpy).


Usage
=====

Use python from the directory src, or add the path to it using
>>> import sys
>>> sys.path.append('path_to_hobj/src')

>>> from hobj import dup_permanental_minor_poly
>>> from domains import ZZ
>>> m = [[i*j if abs(i-j) < 6 else 0 for i in range(20)] for j in range(20)]
>>> sum(dup_permanental_minor_poly(m, ZZ))
11936810897247956264161397956481650508142206788L
>>> dup_permanental_minor_poly(m, ZZ, 1)
11936810897247956264161397956481650508142206788L

>>> from hobj import dup_matching_generating_poly, dup_independence_poly
>>> from graphs_gen import dict_fuller
>>> d = dict_fuller(60)
>>> dup_matching_generating_poly(d, val=1)
1417036634543488
>>> dup_independence_poly(d, val=1)
217727997152

It can be used in Sage:
sage: from hobj import dup_matching_generating_poly
sage: g = graphs.PetersenGraph()
sage: d = g.to_dictionary()
sage: g.matching_polynomial()
x^10 - 15*x^8 + 75*x^6 - 145*x^4 + 90*x^2 - 6
sage: dup_matching_generating_poly(d)
[6, 90, 145, 75, 15, 1]


Notes
=====

The code is organized in a way to be easily included in SymPy.
The low-level representation of SymPy univariate polynomials is used;
the code in src/densearith.py is taken from SymPy.

Eventually hobj.py will appear in a Pull Request for inclusion in SymPy.


Credits
=======

Thanks to Paolo Butera:
this module has been written in preparation of the article

P. Butera, M. Pernici
``Sums of permanental minors using Grassmann algebra''



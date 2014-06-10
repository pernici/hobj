
GMPY = 0

if GMPY == 1:
    from gmpy import mpz, mpq
    class ZZ(object):
        zero = mpz(0)
        one = mpz(1)
        def __new__(self, a):
            return mpz(a)

    class QQ(object):
        zero = mpz(0)
        one = mpz(1)
        def __new__(self, p, q):
            return mpq(p, q)

else:
    from fractions import Fraction
    class ZZ(object):
        zero = 0
        one = 1
        def __new__(self, a):
            return int(a)

    class QQ(object):
        zero = Fraction(0)
        one = Fraction(1)
        def __new__(self, p, q):
            return Fraction(p, q)

"""
Matching polynomials for a sequence of random regular bipartite graphs.

Sequence of regular bipartite graphs constructed in the following way.
Let $G_0$ be a cycle with $k$ vertices, with $k$ even.
Add another cycle with $k$ vertices;
the odd (even) vertices of this cycle are linked respectively
to the even (odd) vertices of the
previous cycle in a random way, obtaining $G_1$;
continue adding cycles in this way, obtaining the sequence ${G_i}$.
The sequence of the corresponding regular bipartite graphs
is obtained linking the vertices of the first and last cycle.


"""
import sys
sys.path.insert(0,'../src')
from copy import deepcopy
from hobj import Hobj
from densearith import dup_add
from random import randint, seed
from active_nodes import d_from_links
from itertools import permutations, product
from domains import ZZ
from compatibility import itervalues

def sum_values(p, K):
    """
    sum the values in ``p``
    """
    nv = []
    for v in itervalues(p):
        nv = dup_add(nv, v, K)
    nv.reverse()
    return nv



def rand_reg6_gen():
    K = ZZ
    p = {0: [K.one]}
    hb = Hobj()
    # first vertical line
    for i in range(6):
        p = hb.iadd_object(p, 1, [i, (i + 1)%6], [], K)
    n = 6
    # vector of allowed permutations
    v = list(product([[0,2,4], [0,4,2]], list(permutations([1,3,5]))))
    v = [[i0,i3,i1,i4,i2,i5] for (i0,i1,i2),(i3,i4,i5) in v]
    N = 10000
    for ii in range(N):
        # add horizontal lines of a strip
        #sys.stderr.write('ii=%d\n' %ii)
        for i in range(6):
            if ii == 0:
                 p = hb.iadd_object(p, 1, [n - 6 + i, n + i], [], K)
            else:
                p = hb.iadd_object(p, 1, [n - 6 + i, n + i], [n - 6 + i], K)
        # add a random vertical line
        a = v[randint(0,len(v) - 1)]
        #a = list(range(6))

        for i in range(6):
            p = hb.iadd_object(p, 1, [n + a[i], n + a[(i + 1)%6]], [], K)
        n += 6
        n1 = n // 6
        if n1 == 2 or n1 % 2 == 1:
            continue
        # closure
        # it is not necessary to copy `p` since `free` is not empty in the
        # first call of iadd_object
        p2 = p
        hb2 = deepcopy(hb)
        for i in range(6):
            p2 = hb2.iadd_object(p2, 1, [i, n - 6 + i], [i, n - 6 + i], K)
        assert len(p2) == 1
        nv = sum_values(p2, K)
        assert len(hb2.links) == len(set(hb2.links))
        d = d_from_links(hb2.links)
        yield n1, dict(d), nv

if __name__ == '__main__':
    try:
        N = int(sys.argv[1])
    except:
        print('prog N')
        sys.exit()
    print('rand_reg6:')
    for n1, d, nv in rand_reg6_gen():
        if n1 > N:
            break
        print('n1=%d sum(nv)=%s' %(n1, sum(nv)))

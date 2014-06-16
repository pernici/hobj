"""
Compute the matching generating polynomials for a sequence of nanotubes
obtained cutting a C_40 fullerene and inserting ``N`` strips of hexagons.
"""
import sys
sys.path.insert(0,'../src')
from time import time
from hobj import Hobj
from graphs_gen import d_from_links
from domains import ZZ
from copy import deepcopy

def nano_d40_seq():
    K = ZZ
    N = 10000
    # initial graph, North cap
    a0 = [(0,2), (2,8), (6,8), (3,6), (0,3)]
    a0 = [(i,j,[]) for i, j in a0]
    a = a0 + [(0, 1, [0]), (3,4,[3]), (6,7,[6]), (8,9,[8]), (2,5,[2]),
    (4,18,[]), (17,18,[]), (1,17,[]), (1,14,[1]), (16,17,[17]), (18,19,[18]),
    (4,21,[4]), (5,13,[]), (12,13,[]), (13,14,[13]),(5,10,[5]), (20,21,[]),
    (19,20,[]),(21,22,[21]), (7,22,[]), (7,25,[7]), (22,23,[22]),
    (23,24,[]), (24,25,[]), (25,26,[25]), (9,26,[]), (9,29,[9]), (26,27,[26]),
    (10,29,[]), (28,29,[29]), (27,28,[]), (10,11,[10]), (11,12,[]),
    (14,15,[14]), (15,16,[])]

    p = {0: [K.one]}
    hb = Hobj()
    for i, j, free in a:
        p = hb.iadd_object(p, 1, [i, j], free, K)

    n = 30
    b = [(-19,4,[-19]), (-18,7,[-18]), (-15,8,[-15]), (7,8,[]),
    (6,7,[7]), (-11,12,[-11]), (-14,11,[-14]), (11,12,[]), (10,11,[11]),
    (8,9,[8]), (9,10,[]), (12,13,[12]), (-10,15,[-10]), (-7,16,[-7]),
    (15,16,[]), (14,15,[15]), (13,14,[]), (16,17,[16]), (-6,19,[-6]),
    (-3,0,[-3]), (0,19,[]), (18,19,[19]), (17,18,[]), (-2,3,[-2]),
    (3,4,[]), (0,1,[0]), (1,2,[]), (2,3,[3]), (4,5,[4]), (5,6,[])]

    c = [(-19,7,[-19]), (-2,5,[-2]), (5,7,[]), (3,5,[5]), (4,7,[7]),
        (-3,6,[-3]), (-6,6,[-6]), (6,8,[6]), (3,8,[]), (8,9,[8]), (-7,9,[-7]),
        (-10,9,[-10,9]), (2,3,[3]), (-11,2,[-11]), (1,2,[2]), (-14,1,[-14]),
        (1,4,[1]), (0,4,[4]), (-18,0,[-18]), (-15,0,[-15,0])]
    for ii in range(N + 1):
        # strip
        for i, j, free in b:
            i += n
            j += n
            free = [k + n for k in free]
            p = hb.iadd_object(p, 1, [i, j], free, K)

        n += 20
        # closure, South cap
        # it is not vecessary to copy `p` since free = c[0][2] = [-19]
        # see docstring of iadd_object
        p2 =  p
        hb2 = deepcopy(hb)
        # strip
        for i, j, free in c:
            i += n
            j += n
            free = [k + n for k in free]
            p2 = hb2.iadd_object(p2, 1, [i, j], free, K)
        assert len(p2) == 1
        nv = [y for y in p2[0]]
        yield (ii + 1, hb2.links, nv)



def test_seq():
    try:
        N = int(sys.argv[1])
    except:
        print('prog N')
        sys.exit()

    for i, links, nv in nano_d40_seq():
        if i > N:
            break
        print('i=%d sum(nv)=%s' %(i, sum(nv)))
        # put 1 here to print the dict of the graph
        if 0:
            print('  d=', d_from_links(links))

if __name__ == '__main__':
    print('d40nano test_seq:')
    test_seq()

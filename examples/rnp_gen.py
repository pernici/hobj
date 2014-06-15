import sys
sys.path.insert(0,'../src')
from hobj import Hobj
from densearith import dup_add
from compatibility import itervalues

def get_sum(p, K):
    """
    nv from p for rectangle non periodic
    """
    nv = []
    for v in itervalues(p):
        nv = dup_add(nv, v, K)
    return nv

def nv_r_nx_ny_np_rec(ny, K):
    """
    generator of ``(nx, ny, nv)``, where `ny` is the input argument;
    ``nv`` is the matching generating polynomial of the grid ``(nx, ny)``
    with open boundary conditions.
    """
    # initial graph
    N = 10000
    p = {0: [K.one]}
    hb = Hobj()
    for i, j in [(i,i+1) for i in range(ny - 1)]:
        p = hb.iadd_object(p, 1, [i, j], [], K)
    for i in range(ny):
        p = hb.iadd_object(p, 1, [i, i + ny], [i], K)
    for i, j in [(k, k+1) for k in range(ny, 2*ny - 1)]:
        p = hb.iadd_object(p, 1, [i, j], [], K)

    n = 2*ny
    for ii in range(N + 1):
        # strip
        for k in range(n - ny, n):
            p = hb.iadd_object(p, 1, [k, k + ny], [k], K)

        for i, j in [(k, k+1) for k in range(n, n + ny - 1)]:
            p = hb.iadd_object(p, 1, [i, j], [], K)

        n += ny

        # closure
        nv = get_sum(p, K)
        yield n//ny, ny, nv



def test1():
    import sys
    from domains import ZZ as K
    from time import time
    try:
        ny = int(sys.argv[1])
    except:
        print('prog ny')
        sys.exit()

    it = nv_r_nx_ny_np_rec(ny, K)
    for n1, n2, nv in it:
        print('n1=%d n2=%d nv=%s' %(n1, n2, nv))
        if n1 >= ny:
            break


if __name__ == '__main__':
    print('rnp_gen:')
    test1()

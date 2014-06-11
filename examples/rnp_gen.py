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
    # initial graph
    N = 10000
    p = {0: [K.one]}
    hb = Hobj()
    for i, j in [(i,i+1) for i in range(ny - 1)]:
        p = hb.iadd_object(p, 1, [i, j], [], K)
    #print 'DB0 links=', links
    for i in range(ny):
        p = hb.iadd_object(p, 1, [i, i + ny], [i], K)
    for i, j in [(k, k+1) for k in range(ny, 2*ny - 1)]:
        p = hb.iadd_object(p, 1, [i, j], [], K)
    #print 'DB1 links=', links

    n = 2*ny
    for ii in range(N + 1):
        #print 'DB1 ii=', ii
        # strip
        for k in range(n - ny, n):
            p = hb.iadd_object(p, 1, [k, k + ny], [k], K)
            #print 'DB2 links=', links

        for i, j in [(k, k+1) for k in range(n, n + ny - 1)]:
        #for i,j in [(n,n+1),(n+1,n+2),(n+2,n+3),(n+3,n+4),(n+4,n+5)]:
            p = hb.iadd_object(p, 1, [i, j], [], K)
            #p = _prm_mul2(p, p1, [, K) # n-1 already eliminated, so useless
                                        # [n-1], but not wrong

        n += ny

        #d = d_from_links(links)
        #print 'DB3 d=', d
        # closure
        #print 'DB2 len(p)=', len(p)
        nv = get_sum(p, K)
        #t2=time()
        #print '%.3f %.3f' %(t1-t0, t2-t1)
        #print '(%d, %d): %s' %(n//ny, ny, nv)
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
        #sys.stderr.write('n1=%d n2=%d %.2f\n' %(n1, n2, t1 - t0))
        print('n1=%d n2=%d nv=%s' %(n1, n2, nv))
        if n1 >= ny:
            break


if __name__ == '__main__':
    print('rnp_gen:')
    test1()

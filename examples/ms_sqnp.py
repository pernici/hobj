import sys
sys.path.insert(0,'../src')

from time import time
import sys
from graphs_gen import sq_d_np
from hobj import dup_independence_poly

def test1():
    try:
        n1 = int(sys.argv[1])
        n2 = int(sys.argv[2])
    except:
        print 'prof n1 n2'
        sys.exit()
    d = sq_d_np(n1, n2)
    vlist = list(range(n1*n2))
    t0 = time()
    p = dup_independence_poly(d, vlist=vlist, val=1)
    t1 = time()
    print 'p=', p
    print '%.2f' %(t1 - t0)


test1()

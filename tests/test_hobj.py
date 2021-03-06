import sys
sys.path.insert(0,'../src')
from active_nodes import (ordered_links, ip_list_objects_from_vlist,
     ip_ordered_vertices)
from hobj import (dup_permanental_minor_poly, gen_hobj,
    dup_matching_generating_poly, dup_independence_poly, hobj_str)

from densearith import dup_valuate
from domains import ZZ

from graphs_gen import sq_mat, dict_fuller, line_graph

SLOW_TEST = 0

def band_mat1(n, k, modulus=11):
    m = [[0]*n for i in range(n)]
    for i in range(n):
        for j in range(n):
            if abs(i - j) < k:
                m[i][j] = i*j % modulus + 1
    return m

def test_dup_permanental_minor_poly():
    m = band_mat1(10, 10)
    r = dup_permanental_minor_poly(m, ZZ)
    # checked with Sage
    # n = len(m); m = Matrix(m)
    # r1 = [m.permanental_minor(i) for i in range(n, -1, -1)]

    r1 = [18029573537492, 51585910340388, 29716353005684, 6172212694596, \
          590164352416, 29473932148, 821101973, 13071249, 116883, 541, 1]
    assert r == r1
    c1 = dup_permanental_minor_poly(m, ZZ, val=1)
    assert sum(r) == c1
    c2 = dup_permanental_minor_poly(m, ZZ, val=-3)
    assert dup_valuate(r, -3) == c2

    m = [[i+2*j for j in range(7)] for i in range(5)]
    r = dup_permanental_minor_poly(m, ZZ)
    assert r == [44862720, 12227040, 915600, 25550, 280, 1]

def test_matching_poly_bipartite():
    m, d1, d2 = sq_mat(6, 6, 'pp')
    p = dup_permanental_minor_poly(m, ZZ)
    # see Lundow "Computation of matching polynomials and the numbers of"
    #  "1-factors in polygraphs Table 6
    # Department of Mathematics, Umea Univ. Research reports No. 12 (1996)

    p1 = [1, 72, 2340, 45456, 589158, 5386752, 35826516, 176198256, 645204321,
          1758028568, 3538275120, 5185123200, 5409088488, 3885146784,
          1829582496, 524514432, 81145872, 5415552, 90176]
    p1.reverse()
    assert p == p1

def test_matching_generating_poly():
    d = dict_fuller(60)
    nv = dup_matching_generating_poly(d)
    # C60 matching polynomial, see 'H. Hosoya Comp & Maths. with Appls.
    # Vol. 12B (1986) 271, Table 1, Fig.19(i)
    assert nv == [12500, 4202760, 269272620, 6946574300, 94541532165, \
      783047312406, 4310718227685, 16742486291340, 47883826976580, \
      104113567937140, 176345540119296, 237135867688980, 256967614454320, \
      227043126274260, 165074851632300, 99463457244844, 49924889888850, \
      20949286202160, 7362904561730, 2168137517940, 534162544380, \
      109742831260, 18697786680, 2619980460, 298317860, 27130596, \
      1922040, 102120, 3825, 90, 1]
    d = dict_fuller(70)
    nv = dup_matching_generating_poly(d)
    # TODO add reference
    assert nv == [52168, 24949770, 2175411410, 76188054720, 1414031839415, 16060741380337, 121977585926655, 657705883833585, 2628626806075280, 8042439217209190, 19311918160511451, 37116989376321905, 58004771574641390, 74650239997165790, 79943103597641110, 71843167484323524, 54553421377655995, 35194282142522435, 19372757557305170, 9127603794188980, 3688877411254177, 1280236743404535, 381580564799125, 97578536447295, 21362502290720, 3989784929340, 632435470975, 84479881155, 9417493125, 864583005, 64182939, 3753855, 166390, 5250, 105, 1]

def test_dup_independence_poly():
    # the following independence polynomials for fullerenes C20,..,C30 agree with
    # M.B. Ahmadi, H. Alimorad Daskhezr, Commun. Math. Comput. Chem. 71 (2014) 355.
    d = dict_fuller(20)
    ip = dup_independence_poly(d)
    assert ip == [5, 320, 1240, 1912, 1510, 660, 160, 20, 1]
    d = dict_fuller(24)
    ip = dup_independence_poly(d)
    assert ip == [264, 2292, 6756, 9918, 8316, 4212, 1304, 240, 24, 1]
    d = dict_fuller(26)
    ip = dup_independence_poly(d)
    assert ip == [2, 218, 2660, 10523, 20298, 22361, 15120, 6461, 1742, 286, 26, 1]
    d = dict_fuller(28)
    ip = dup_independence_poly(d)
    assert ip == [176, 2860, 14840, 36903, 52200, 45684, 25776, 9506, 2268, 336, 28, 1]
    d = dict_fuller(30)
    ip = dup_independence_poly(d)
    assert ip == [150, 2910, 19263, 60650, 108120, 119670, 86435, 41724, 13515, 2890, 390, 30, 1]

    from graphs_gen import sq_d_np
    n1, n2 = 6, 6
    d = sq_d_np(n1, n2)
    vlist = list(range(n1*n2))
    p = dup_independence_poly(d, vlist=vlist)
    a = [2, 40, 416, 3120, 17834, 76716, 245030, 576052, 991750, 1247116,
            1143638, 763144, 368868, 127960, 31320, 5248, 570, 36, 1]
    assert p == a
    r = dup_independence_poly(d, vlist=vlist, val=1)
    # see https://oeis.org/search?q=A006506
    assert r == sum(a) == 5598861
    pr = 73
    r1 = dup_independence_poly(d, vlist=vlist, val=1, pr=pr)
    assert r1 == r % pr
    n1, n2 = 10, 10
    d = sq_d_np(n1, n2)
    vlist = list(range(n1*n2))
    r = dup_independence_poly(d, vlist=vlist, val=1)
    assert r == 2030049051145980050

    
def test_line_graph():
    d = dict_fuller(20)
    d1 = line_graph(d)
    p = dup_matching_generating_poly(d)
    p1 = dup_independence_poly(d1)
    assert p == p1

def is_independent_set(d, a):
    """
    test if vertices in `a` are independent in graph with dict `d`

    Parameters
    ==========

    d : dictionary of the graph
    a : tuple or list of vertices
    """
    n = len(a)
    for i in range(n):
        for j in range(i):
            k1 = a[i]
            k2 = a[j]
            if k1 in d[k2]:
                return False
    return True

def is_independent_sets(p, d, n):
    """
    test if the `n`-sets in `p` are independent in `d`

    Parameters
    ==========

    p : polynomial generated by gen_hobj
    d : dict of the graph
    n " number of elements in set

    """
    from hobj import _monom
    p = {}
    for k, v in p.items():
        c = count_bits_set(k)
        if c == n:
            p1[k] = v
            a1.append(k)
            am = _monom(k)
            b = is_independent_set(d, am)
            if not b:
                return False
    return True

def test_gen_hobj():
    d={0: [1, 2, 3], 1: [0], 2: [0], 3: [0]}
    links0 = [0, d[0][0]]
    vlist = ip_ordered_vertices(d, *links0)
    objects = ip_list_objects_from_vlist(d, vlist)
    p = gen_hobj(objects, vlist)
    assert hobj_str(p) == '1 + x0 + x1 + x2 + x3 + x1*x2 + x1*x3 + x2*x3 + x1*x2*x3'
    ip = dup_independence_poly(d)
    assert ip == [1, 3, 4, 1]

def test_gen_hobj_C34():
    # the independence polynomial for C34 as in
    # M.B. Ahmadi, H. Alimorad Daskhezr, Commun. Math. Comput. Chem. 71 (2014) 355.
    # differs in a few coefficients; in particular they find for the
    # leading coefficient 94, while we get 95. Here the 95 cases found
    # are enumerated, and it is checked that they are independent sets.
    d = dict_fuller(34)
    ip = dup_independence_poly(d)
    n = len(ip) - 1
    assert n == 14
    assert ip[0] == 95
    vlist = ip_ordered_vertices(d)
    objects = ip_list_objects_from_vlist(d, vlist)
    p = gen_hobj(objects, vlist)
    assert is_independent_sets(p, d, n)


if __name__ == '__main__':
    test_dup_permanental_minor_poly()
    test_matching_poly_bipartite()
    test_matching_generating_poly()
    test_dup_independence_poly()
    test_line_graph()
    test_gen_hobj()

    if SLOW_TEST:
        test_gen_hobj_C34()
    print('test_hobj ok')

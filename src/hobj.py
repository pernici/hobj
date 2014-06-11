""" Algebraic formalism for hard objects

  Hard objects are defined as products of ``eta_i`` where ``i`` is
  the index of a node of the graph on which the hard object is defined.
  The ``eta`` elements are commuting and nilpotent.
  They satisfy
  ``<eta_{i_1} ...eta_{i_k}> = 1``
  if ``i_1,..,i_k`` are all different.

  The generating function for counting hard objects is
  ``Z(t) = <Prod (1 + t*Prod eta_i)>``

"""
from densearith import (dup_lshift, dup_add, dup_mul, dup_mul_ground)
from densearith import dup_degree, dup_strip, dup_lshift
from active_nodes import ip_ordered_vertices, ip_list_objects_from_vlist
from compatibility import iteritems
from domains import ZZ

#                  0  1  10 11, 100 101, 110, 111
_bits_set_table = [0, 1, 1, 2,  1,  2,   2,   3  ]

def count_bits_set(n):
    """
    count the number of set bits in a number
    """
    c = 0
    while n:
        c += _bits_set_table[n & 7]
        n >>= 3
    return c

def _is_free_var(i, k, m):
    """
    i  current row
    k  index of the variable
    m  matrix as a list of lists
    returns true if the variable `k` does not occur from row `i` on
    """
    for j in range(i, len(m)):
        if m[j][k] != 0:
            return False
    return True

def _monom(n):
    """
    monomial in `eta` variables from the number `n` encoding it
    """
    v = []
    i = 0
    while n:
        if n % 2:
            v.append(i)
        n = n >> 1
        i += 1
    return v

def _prm_mul(p1, p2, free_vars_indices, K):
    """
    helper function for dup_permanental_minor_poly
    Return the product of ``p1`` and ``p2``

    Parameters
    ==========

    p1 : polynomial
    p2 : linear polynomial
    free_vars_indices : list of variables not used anymore

    Notes
    =====

    ``p1`` is a polynomial in the ``x_i=t*eta_i`` elements, with polynomials
    in ``t`` as coefficients;
    ``p2`` is a linear polynomial in the ``x_i`` elements, with numbers
    as coefficients.
    If after the product an ``eta-object``
    is not anymore used, it is replaced by the variable ``t``

    Examples
    ========

    ``p1 = (1 + 2*t) + (1+t)*x_0; p2 = 1 + x_0 + x_1``
    ``<p1*p2> = _prm_mul(p1, p2, [0, 1], ZZ)``

    >>> from hobj import _prm_mul
    >>> from domains import ZZ
    >>> p1 = {0:[2, 1], 1:[1, 1]}
    >>> p2 = {0:1, 1:1, 2:1}
    >>> _prm_mul(p1, p2, [0, 1], ZZ)
    {0: [1, 6, 5, 1]}
    """
    p = {}
    mask_free = 0
    for i in free_vars_indices:
        mask_free += 1 << i
    if not p2:
        return p
    get = p.get
    for exp1, v1 in iteritems(p1):
        for exp2, v2 in p2.items():
            if exp1 & exp2:
                continue
            exp = exp1 | exp2
            v =  dup_mul_ground(v1, v2, K)
            if exp & mask_free:
                for i in free_vars_indices:
                    if exp & (1 << i):
                        exp = exp ^ (1 << i)
                        v = dup_lshift(v, 1, K)
            p[exp] = dup_add(get(exp, []), v, K)
    return p

def dup_permanental_minor_poly(m, K):
    """
    return the polynomial of the sum of permanental minors of a matrix ``m``

    Let ``perm(i, m)`` be the sum of the permanents
    of all ``i x i`` minors of ``m``

    The polynomial of the sum of permanental minors is
    ``sum_{i=0}^{len(m)-1} perm(i, m)*x**i``


    Parameters
    ==========

    m : matrix in list form

    Examples
    ========

    >>> from hobj import dup_permanental_minor_poly
    >>> from domains import ZZ
    >>> m = [[2,1,2],[3,0,1],[1,1,2]]
    >>> dup_permanental_minor_poly(m, ZZ)
    [15, 36, 13, 1]
    """

    n = len(m)
    ny = len(m[0])
    p = {0:[K.one]}
    done_vars = set()
    for i in range(n):
        p1 = {0: K.one}
        a = m[i]
        for j in range(len(a)):
            if a[j]:
                p1[1<<j] = a[j]

        free_vars_indices = []
        for j in range(ny):
            if j in done_vars:
                continue
            r = _is_free_var(i+1, j, m)
            if r:
                free_vars_indices.append(j)
                done_vars.add(j)
        p = _prm_mul(p, p1, free_vars_indices, K)

    assert len(p) == 1
    nv = [y for y in p[0]]
    return nv


def _get_num_elements(objects):
    """"
    number of elements in `objects`
    """
    s = set()
    for t in objects:
        for i in t:
            s.add(i)
    nvars = len(s)
    b = list(sorted(s))
    b.sort()
    if b != list(range(nvars)):
         raise ValueError('elements should be labelled in 0,...,nvars-1')
    return nvars


def _poly_str(a):
    """
    string representation of univariate polynomial
    """
    n = dup_degree(a)
    v = []
    for i in range(n + 1):
        c = a[i]
        if not c:
            continue
        if i == n:
            m = '1'
        elif i == n - 1:
            m = 't'
        else:
            m = 't^%d' % (n - i)
        if c > 0:
            if m != '1':
                if c == 1:
                    if i == n:
                        sx = '%s' % m
                    else:
                        sx = '+%s' % m
                else:
                    sx = '+%s*%s' % (c, m)
            else:
                sx = '+%s' % c
        else:
            sx = '%s*m' % (c, m)
        v.append(sx)
    s = ''.join(v)
    return s

def hobj_str(p, noval=True):
    """
    string representation of polynomial for hard objects

    Parameters
    ==========

    p : polynomial for hard objects
    noval : if True the coefficients are put to ``1``
    """
    a = []
    items = list(p.items())
    items.sort(key=lambda k: (count_bits_set(k[0]), _monom(k[0])))
    for exp1, v1 in items:
        if exp1 == 0:
            if noval:
                a.append('1')
            else:
                s = _poly_str(v1)
                a.append(s)
        else:
            sx = '*'.join(['x%d' % ii for ii in _monom(exp1)])
            if noval:
                a.append('%s' % sx)
            else:
                s = _poly_str(v1)
                a.append('%s*(%s)' %(sx, s))
    s = ' + '.join(a)
    return s

def hobj_list(p):
    """
    return the list of lists representing independent sets

    Examples
    ========

    >>> from hobj import gen_hobj, hobj_list
    >>> p = gen_hobj([(0,1),(1,2),(2,3),(3,4),(0,4)])
    >>> hobj_list(p)
    [[], [0], [1], [2], [3], [4], [0, 2], [0, 3], [1, 3], [1, 4], [2, 4]]
    """
    items = list(p.items())
    items.sort(key=lambda k: (count_bits_set(k[0]), _monom(k[0])))
    a = [_monom(expv) for expv, y in items]
    return a


def _hobj_mul_val(p1, p2, free_vars_indices, K):
    p = {}
    mask_free = 0
    for i in free_vars_indices:
        mask_free += 1 << i
    if not p2:
        return p
    get = p.get
    for exp1, v1 in iteritems(p1):
        for exp2, v2 in p2.items():
            if exp1 & exp2:
                continue
            exp = exp1 | exp2
            v =  v1*v2
            if exp & mask_free:
                for i in free_vars_indices:
                    if exp & (1 << i):
                        exp = exp ^ (1 << i)
            p[exp] = get(exp, 0) + v
    return p

def gen_hobj(objects, vlist=None):
    """
    polynomial enumerating all hard object configurations

    Parameters
    ==========

    objects : list of tuple of object element indices
    vlist: list of object indices, a permutation of ``range(len(vlist))``

    Notes
    =====

    TODO explain use of vlist (different naming of indices)

    Examples
    ========
    >>> from hobj import gen_hobj, hobj_str
    >>> hobj_str(gen_hobj([(0,1),(1,2),(2,3),(3,4),(0,4)]))
    '1 + x0 + x1 + x2 + x3 + x4 + x0*x2 + x0*x3 + x1*x3 + x1*x4 + x2*x4'
    """
    from domains import ZZ
    nvars = _get_num_elements(objects)
    masks = []
    done_vars = set()
    one = 1
    p = {0:one}
    n = len(objects)
    if not vlist:
        vlist = list(range(n))
    else:
        assert list(sorted(vlist)) == list(range(n))
    # use the bits up to `n` for the hard object variables;
    # use the subsequent bits for the `eta` elements
    for i in range(n):
        expv = 0
        for j in objects[i]:
            expv = expv + (1 << (j+n))
        masks.append(expv)
    for i in range(n):
        p1 = {0: one}
        expv = masks[i] + (1 << vlist[i])
        p1[expv] = one
        free_vars = set()
        for j in range(n, n + nvars):
            if j in done_vars:
                continue
            hit = False
            jmask = 1 << j
            for mask in masks[i + 1:]:
                if mask & jmask:
                    hit = True
                    break
            if not hit:
                done_vars.add(j)
                free_vars.add(j)
        #p = _prm_mul2(p, p1, free_vars, K)
        p = _hobj_mul_val(p, p1, free_vars, ZZ)

    return p



def obj_free(objects):
    """
    return a list of tuples ``(obj, free)``
    """
    nvars = _get_num_elements(objects)
    v = []
    done = set()
    n = len(objects)
    for i in range(n):
        free = []
        obj = objects[i]
        for j in range(nvars):
            if j in done:
                continue
            hit = 0
            for k in range(i + 1, n):
                if j in objects[k]:
                    hit = 1
                    break
            if not hit:
                free.append(j)
                done.add(j)
        v.append((obj, free))
    return v


class Hobj(object):
    """
    class used with iadd_object, iadd_object_val
    TODO: use if also with _monom, hobj_str; then put these functions as methods
    """
    def __init__(self, pr=None):
        self.links = []
        self.dt = {}
        self.freedt = list(range(1000, -1, -1))
        self.pr = pr


    def hobj_str(hob, p, noval=True):
        """
        string representation of polynomial for hard objects

        Parameters
        ==========

        p : polynomial for hard objects
        noval : if True the coefficients are put to ``1``
        """
        dtinv = {}
        for i, j in iteritems(hob.dt):
            dtinv[j] = i

        a = []
        items = p.items()
        items.sort(key=lambda k: (count_bits_set(k[0]), _monom(k[0])))
        for exp1, v1 in items:
            if exp1 == 0:
                if noval:
                    a.append('1')
                else:
                    s = _poly_str(v1)
                    a.append(s)
            else:
                t = _monom(exp1)
                t = [dtinv[i] for i in t]
                sx = '*'.join(['x%d' % ii for ii in t])
                if noval:
                    a.append('%s' % sx)
                else:
                    s = _poly_str(v1)
                    a.append('%s*(%s)' %(sx, s))
        s = ' + '.join(a)
        return s

    def hobj_list(hob, p):
        """
        return the list of lists representing independent sets

        Examples
        ========

        TODO
        """
        dtinv = {}
        for i, j in iteritems(hob.dt):
            dtinv[j] = i


        items = p.items()
        items.sort(key=lambda k: (count_bits_set(k[0]), _monom(k[0])))
        a = [[dtinv[i] for i in _monom(expv)] for expv, y in items]
        return a

    def iadd_object(hb, p, val, obj, free, K, valued=None, pr=None):
        """
        multiply ``p`` by ``(1 + t*val*eta_i*eta_j)``

        Notes
        =====

        ``p`` is in general changed.

        ``free`` is the list of indices of ``eta`` elements which
        are integrated (that is, put to ``1`` after performing the product).

        Examples
        ========
        
        >>> from domains import ZZ
        >>> from hobj import Hobj, obj_free
        >>> p = {0: [ZZ.one]}
        >>> hb = Hobj()
        >>> a = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]
        >>> a = obj_free(a); a
        [((0, 1), []), ((1, 2), [1]), ((2, 3), [2]), ((3, 4), [3]), ((4, 0), [0, 4])]
        >>> for t, free in a:
        ...   p = hb.iadd_object(p, 1, t, free, ZZ)
        ...
        >>> p[0]
        [5, 5, 1]

        """
        links = hb.links
        dt = hb.dt
        freedt = hb.freedt
        exp2 = 0
        for i in obj:
            if i in dt:
                j = dt[i]
            else:
                j = freedt.pop()
                dt[i] = j
            exp2 += 1 << j
        free0 = free
        free = [dt[i] for i in free]
        t = tuple(sorted(obj))
        if t in links:
            raise ValueError('%s in %s' %(t, links))
        links.append(t)
        if free:
            p1 = p
            p = {}
            mask_free = 0
            for i in free:
                mask_free += 1 << i
            get = p.get
            a = [1, (val, exp2)]
            for exp1, v1 in iteritems(p1):
                for ii in range(2):
                    if not ii:
                        exp = exp1
                        v = v1
                    else:
                        if exp1 & exp2:
                            continue
                        exp = exp1 | exp2
                        if valued:
                            if val != 1:
                                v1 = v1*val
                        else:
                            if val != 1:
                                v1 = dup_mul_ground(v1, val, K)
                            v = dup_lshift(v1, 1, K)
                    if exp & mask_free:
                        for i in free:
                            if exp & (1 << i):
                                exp = exp ^ (1 << i)
                    if valued:
                        c = get(exp, 0) + v
                        if pr:
                            c = c % pr
                    else:
                        p[exp] = dup_add(get(exp, []), v, K)
            for exp in free:
                freedt.append(exp)
            return p

        else:
            a = []
            for exp1, v1 in iteritems(p):
                if exp1 & exp2:
                    continue
                exp = exp1 | exp2
                v1 = dup_mul_ground(v1, val, K)
                v = dup_lshift(v1, 1, K)
                try:
                    p[exp] = dup_add(p[exp], v, K)
                except KeyError:
                    a.append((exp, v))
            for exp, v in a:
                p[exp] = v

        return p

    def iadd_object_val(hb, p, val, obj, free, K, pr=None):
        """
        multiply ``p`` by ``(1 + t*val*eta_i*eta_j)``

        Notes
        =====

        ``p`` is in general changed.

        ``free`` is the list of indices of ``eta`` elements which
        are integrated (that is, put to ``1`` after performing the product).

        Examples
        ========

        """
        links = hb.links
        dt = hb.dt
        freedt = hb.freedt
        exp2 = 0
        for i in obj:
            if i in dt:
                j = dt[i]
            else:
                j = freedt.pop()
                dt[i] = j
            exp2 += 1 << j
        free0 = free
        free = [dt[i] for i in free]
        t = tuple(sorted(obj))
        if t in links:
            raise ValueError('%s in %s' %(t, links))
        links.append(t)
        if free:
            p1 = p
            p = {}
            mask_free = 0
            for i in free:
                mask_free += 1 << i
            get = p.get
            a = [1, (val, exp2)]
            for exp1, v1 in iteritems(p1):
                for ii in range(2):
                    if not ii:
                        exp = exp1
                        v = v1
                    else:
                        if exp1 & exp2:
                            continue
                        exp = exp1 | exp2
                        if val != 1:
                            v = v1*val
                    if exp & mask_free:
                        for i in free:
                            if exp & (1 << i):
                                exp = exp ^ (1 << i)
                    c = get(exp, 0) + v
                    if pr:
                        c = c % pr
                    p[exp] = c

            for exp in free:
                freedt.append(exp)
            return p

        else:
            a = []
            for exp1, v1 in iteritems(p):
                if exp1 & exp2:
                    continue
                exp = exp1 | exp2
                v = v1 * val
                if pr:
                    v = v % pr
                try:
                    p[exp] += v
                except KeyError:
                    a.append((exp, v))
            for exp, v in a:
                p[exp] = v

        return p

def d_relabel(d):
    dt = {}
    keys = list(d.keys())
    dt = dict(zip(keys, range(len(d))))
    d1 = {}
    for k, v in d.items():
        d1[dt[k]] = [dt[i] for i in v]
    return d1, dt


def dup_gen_count_hobj(objects, K, val=None, pr=None):
    """
    Counting polynomial for hard object from a list of edges

    Parameters
    ==========

    objects : list of tuples of element indices
    TODO add docs

    Notes
    =====

    The counting polynomial for hard object is
    ``M(x) = sum_i N(i)*x**i``
    where ``N(i)`` is the number of ways in which ``i`` dimers can be
    put in the graph.
    In the case of dimers it is the matching generating polynomial, which
    is related to the matching polynomial ``mu`` by the relation
    ``mu(x) = x**v*M(-x**2)``
    where ``v`` is the number of vertices of the graph.

    For efficiency reasons the list of edges should be appropriately
    ordered to reduce the number of active nodes.

    Examples
    ========

    >>> from domains import ZZ
    >>> from hobj import dup_gen_count_hobj
    >>> dup_gen_count_hobj([(1, 0), (2, 1), (3, 2), (4, 0), (4, 3)], ZZ)
    [5, 5, 1]
    >>> dup_gen_count_hobj([(1, 0), (2, 1), (3, 2), (4, 0), (4, 3)], ZZ, val=1)
    11
    """
    a = obj_free(objects)
    hb = Hobj(pr=pr)
    if not val:
        if pr:
            raise NotImplementedError
        p = {0: [K.one]}
        for obj, free in a:
            p = hb.iadd_object(p, 1, obj, free, K)
    else:
        p = {0: K.one}
        for obj, free in a:
            p = hb.iadd_object_val(p, val, obj, free, K, pr)

    assert len(p) == 1
    return p[0]


def dup_matching_generating_poly(d, val=None, pr=None, links=None, K=ZZ):
    """
    Return the matching polynomial for the graph defined by ``d``

    Parameters
    ==========

    d : dict for the graph
    val : evaluate the polynomial in ``val``
    pr : evaluate the polynomial modulo the prime ``pr``
    links : list of vertices of the graph forming a path

    Notes
    =====

    A simple greedy algorithm tries to find an efficient ordering of
    vertices to compute the independence polynomial.
    To help it, it can be provided an initial path ``links``;
    for instance in a long rectangle ``(m, n)``, with ``m`` much greater
    than ``n``, the path can be a short side ``n``;
    in this case the greedy algorithm tries to deduce an efficient ordering.

    Examples
    ========

    >>> from hobj import dup_matching_generating_poly
    >>> d = {0:[1,4], 1:[0,2], 2:[1,3], 3:[2,4], 4:[0,3]}
    >>> dup_matching_generating_poly(d)
    [5, 5, 1]

    """
    from active_nodes import ordered_links
    bd = True
    if list(sorted(d.keys())) != list(range(len(d))):
        bd = False
        d, dt = d_relabel(d)
        if not links:
            links = [[dt[k] for k in obj] for obj in links]
    if not links:
        k0 = 0
        links = [k0, d[k0][0]]
        ord_links = ordered_links(d, *links)
    else:
        ord_links = links
    p = dup_gen_count_hobj(ord_links, K, val, pr)
    return p


def dup_independence_poly(d, val=None, pr=None, links=None, vlist=None, K=ZZ):
    """
    Return the independence polynomial for the graph defined by ``d``

    Parameters
    ==========

    d : dict for the graph
    val : evaluate the polynomial in ``val``
    pr : evaluate the polynomial modulo the prime ``pr``
    links : list of vertices of the graph forming a path
    vlist : list of vertices of the graph

    Notes
    =====

    A simple greedy algorithm tries to find an efficient ordering of
    vertices to compute the independence polynomial.
    To help it, it can be provided an initial path ``links``;
    for instance in a long rectangle ``(m, n)``, with ``m`` much greater
    than ``n``, the path can be a short side ``n``;
    in this case the greedy algorithm deduces an efficient ordering
    of vertices, that is it adds short vertical lines of vertices
    till all the rectangle is obtained; the number of active nodes is
    ``nu=n``.
    Alternatively, one can provide an ordering of all vertices in the
    graph in ``vlist``.
    Giving no hint, the algorithm might start with an horizontal line
    of vertices, and continue adding horizontal lines; then the number
    of active nodes would be ``nu=m``. The complexity of the
    computation of the independence polynomial depending on ``2**nu``,
    it would take much longer.

    Examples
    ========

    >>> from hobj import dup_independence_poly
    >>> d = {0:[1,3], 1:[0,2], 2:[1,3], 3:[0,2]}
    >>> dup_independence_poly(d)
    [2, 4, 1]
    """
    bd = True
    if list(sorted(d.keys())) != list(range(len(d))):
        bd = False
        d, dt = d_relabel(d)
    if not vlist:
        if not links:
            k0 = 0
            links = [k0, d[k0][0]]
        vlist = ip_ordered_vertices(d, *links)
    elif not bd:
        vlist = [dt[k] for k in vlist]
    objects = ip_list_objects_from_vlist(d, vlist)
    p = dup_gen_count_hobj(objects, K, val, pr)
    return p

def independent_sets_gen(d):
    """
    Generator for the independent sets

    """
    vlist = ip_ordered_vertices(d)
    objects = ip_list_objects_from_vlist(d, vlist)
    p = gen_hobj(objects, vlist)
    for expv, _ in iteritems(p):
        yield _monom(expv)

def independence_sets(d):
    """
    Return the list of independence sets of the graph with dict ``d``.

    Examples
    ========

    >>> from hobj import independence_sets
    >>> d = {0:[1,3], 1:[0,2], 2:[1,3], 3:[0,2]}
    >>> independence_sets(d)
    [[], [0], [1], [2], [3], [0, 2], [1, 3]]
    """
    vlist = ip_ordered_vertices(d)
    objects = ip_list_objects_from_vlist(d, vlist)
    p = gen_hobj(objects, vlist)
    a = hobj_list(p)
    return a


if __name__ == "__main__":
    import doctest
    import sys
    if sys.version_info < (2, 6):
        print('doctests require Fraction, available from Python2.6')
        sys.exit()
    doctest.testmod()

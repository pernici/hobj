from collections import defaultdict, deque


def _add_links1(links, d1, d2):
    """
    add links without increasing the number of active nodes

    Parameters
    ==========
    links : list of edges defining the part of the graph covered
    d1 : adjacency dict of the graph
    d2 : adjacency dict of the part of the graph covered
    """
    hit = True
    while hit:
        hit = False
        for k, a in d2.items():
            # if there is only one edge missing to complete a node, add it
            if len(a) == len(d1[k]) - 1:
                a1 = d1[k]
                k1 = [x for x in a1 if x not in a][0]
                if k > k1:
                    k, k1 = k1, k
                d2[k].append(k1)
                d2[k1].append(k)
                hit = True
                links.append((k, k1))
        # it two nodes adjacent in d1 are present in d2, and there
        # is no edge between them in d2, add it
        for k1, a in d1.items():
            for k2 in a:
                if k1 < k2 and k1 in d2 and \
                    ((k2 not in d2 and len(d2[k1]) == len(d1[k1]) - 1) or \
                    (k2 in d2 and k2 not in d2[k1])):
                    d2[k1].append(k2)
                    d2[k2].append(k1)
                    hit = True
                    links.append((k1, k2))

    return links


def m_from_d(d):
    n = len(d)
    m = [[0]*n for i in range(n)]
    for k1, v in d.items():
        for k2 in v:
            m[k1][k2] = 1
    m1 = zip(*m)
    return m

def d_from_m_bip(m):
    d = defaultdict(list)
    n1 = len(m)
    for i in range(n1):
        for j in range(len(m[0])):
            if m[i][j]:
                d[i].append(j+n1)
                d[j+n1].append(i)
    return d

def d_from_links(links):
    """
    dict for the graph from a list of its edges
    """
    d = defaultdict(list)
    for i, j in links:
        d[i].append(j)
        d[j].append(i)
    return d


def _append_link(dx, links, k1, k2):
    assert not (k1,k2) in links
    assert not (k2,k1) in links
    dx[k1].append(k2)
    dx[k2].append(k1)
    links.append((k1, k2))

def _short_path_active_nodes(G, s, a):
    """
    Return a short path starting in `s` and ending in a node in `a`
    Return None if there is not such a path
    """
    P, Q = {s: None}, deque([s])
    while Q:
        u = Q.popleft()
        for v in G[u]:
            if v in P:
                continue
            P[v] = u
            if v in a:
                a1 = [v]
                while 1:
                    v = P[v]
                    a1.append(v)
                    if v in a:
                        return a1
                assert 0
            Q.append(v)
    return None

def _short_path_active_nodes_all(d, dx, active):
    ax = [0]*len(d)
    # adjacency list of complement of G and Gx
    for k, v in d.items():
        v1 = dx[k]
        ak = [y for y in v if y not in v1]
        ax[k] = ak
    min_length = 10000
    min_path = None
    for i in active:
        a1 = _short_path_active_nodes(ax, i, active)
        if a1 is not None and len(a1) < min_length:
            min_length = len(a1)
            min_path = a1
    return min_path

def _add_paths(d, dx, links, a1):
    added = []
    for i in range(len(a1) - 1):
        k1, k2 = a1[i], a1[i+1]
        if k2 < k1:
            k2, k1 = k1, k2
        _append_link(dx, links, k1, k2)
        added.append((k1, k2))
    r = _add_links1(links, d, dx)
    added.extend(r)
    return added

def ordered_links(d, k0, k1):
    """
    find ordered links starting from the link (k0, k1)

    Parameters
    ==========

    d : dict for the graph
    k0, k1: adjacents nodes of the graphs
    """
    assert k0 in d
    assert k1 in d[k0]
    dx = defaultdict(list)
    links = []
    _append_link(dx, links, k0, k1)
    r = _add_links1(links, d, dx)
    while 1:
        active = [k for k in dx if 0 < len(dx[k]) < len(d[k])]
        if not active:
            break
        a1 = _short_path_active_nodes_all(d, dx, active)
        if a1 is None:
            break
        a2 = _add_paths(d, dx, links, a1)
    return links

def ordered_links_all(d):
    """
    find ordered links
    """
    max_active = 10000
    for k1x, v in d.items():
        for k2x in v:
            if k1x > k2x:
                continue
            linksx = ordered_links(d, k1x, k2x)
            ma = num_active_nodes(d, linksx)
            if ma < max_active:
                max_active = ma
                k1, k2 = k1x, k2x
                links = linksx
    return k1, k2, links, max_active

def num_active_nodes(d, links):
    """
    number of the active nodes for the graph defined by ``links``

    Parameters
    ==========

    d : dict for the graph
    links : list of edges of the graph
    """
    dx = defaultdict(list)
    max_active = 0
    for edge in links:
        k1, k2 = edge
        dx[k1].append(k2)
        dx[k2].append(k1)
        nactivex = len([k for k in dx if 0 < len(dx[k]) < len(d[k])])
        if max_active < nactivex:
            max_active = nactivex
    return max_active

def ip_ordered_vertices(d, *links0):
    """
    find ordered vertices starting from a path ``links0``
    """
    for k in links0:
        assert k in d
    if not links0:
        links0 = [0, d[0][0]]
    res = list(links0)
    dx = defaultdict(list)
    for i in range(len(links0) - 1):
        k0 = links0[i]
        k1 = links0[i + 1]
        dx[k0].append(k1)
        dx[k1].append(k0)
    while 1:
        active = [k for k in d if 0 < len(dx[k]) < len(d[k])]
        if not active:
            break
        a1 = _short_path_active_nodes_all(d, dx, active)
        if a1 is None:
            break
        res.extend(a1[1:-1])
        for i in range(len(a1) - 1):
            k1, k2 = a1[i], a1[i+1]
            if k2 < k1:
                k1, k2 = k2, k1
            dx[k1].append(k2)
            dx[k2].append(k1)
    if len(res) < len(d):
        for k in d:
            if k not in res:
                res.append(k)
    return res

def ip_get_dn(d):
    """
    return the dictionary with items ``((v1,v2), i)``,
    where ``(v1, v2)`` is an edge, and ``i`` the index given to the vertex
    """
    c = 0
    dn = {}
    for k, v in d.items():
        for k2 in v:
            a = [k, k2]
            a.sort()
            t = tuple(a)
            if t in dn:
                continue
            dn[t] = c
            c += 1
    return dn

def ip_add_object(v, d, dn, ac, vlist):
    """
    v : added object (a vertex)
    d : dict of the graph
    dn : dictionary of edge indices
    ac : set of active elements
    vlist : list of added vertices
    """
    vlist.append(v)
    a = d[v]
    for j in a:
        if j in ac:
            ac.remove(j)
        else:
            ac.add(j)

def ip_list_objects_from_vlist(d, vlist, dn=None):
    """
    return list of objects diven ``d`` and a list of vertices

    Parameters
    ==========

    d : dict of the graph
    vlist : list of vertices of the graph
    dn : dict mapping an edge to its index
    """
    if not dn:
        dn = ip_get_dn(d)
    a = []
    for k in vlist:
        w = []
        if k not in d:
            raise KeyError('k=%s not in d' % k)
        b = d[k]
        for k2 in b:
            t = [k, k2]
            t.sort()
            t = tuple(t)
            n = dn[t]
            w.append(n)
        a.append(tuple(w))
    return a

def ip_num_active_elements(objects):
    """
    Returns the number of active elements correspponding to ``objects``

    Parameters
    ==========

    objects : list of tuple of indices of edges

    Notes
    =====

    In the computation of the independence polynomial an object is
    a tuple of edges starting at a vertex.
    Different vertex ordering correspond to different ``objects``;
    for each ordering there is a number of active elements,
    which is the maximum number ``nu`` of ``eta`` elements appearing in the
    computation of the independence polynomial; its complexity depends
    from ``2**nu``.
 
    """
    ac = set()
    nu = 0
    for obj in objects:
        for i in obj:
            if i in ac:
                ac.remove(i)
            else:
                ac.add(i)
        if nu < len(ac):
            nu = len(ac)
    return nu

if __name__ == "__main__":
    import doctest
    import sys
    if sys.version_info < (2, 6):
        print('doctests require Fraction, available from Python2.6')
        sys.exit()
    doctest.testmod()

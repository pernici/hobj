def iteritems(p):
    try:
        it = p.iteritems()
    except AttributeError:
        it = p.items()
    return it

def itervalues(p):
    try:
        it = p.itervalues()
    except AttributeError:
        it = p.values()
    return it

class Variable(object):
    def __init__(self, f):
        self.f = f

    def __and__(self, othr):
        return Variable(lambda x: self.f(x) & othr.f(x))

    def __or__(self, othr):
        return Variable(lambda x: self.f(x) | othr.f(x))

    def __matmul__(self, othr):
        return Variable(lambda *args, **kw: self.f(othr.f(*args, **kw)))

    def __invert__(self):
        return Variable(lambda *args, **kw: ~self.f(*args, **kw))

    def __call__(self, *args, **kw):
        return self.f(*args, **kw)

    def __lt__(self, othr):
        if isinstance(othr, Variable):
            return Variable(lambda x: self.f(x) < othr.f(x))
        else:
            return Variable(lambda x: self.f(x) < othr)

    def __le__(self, othr):
        if isinstance(othr, Variable):
            return Variable(lambda x: self.f(x) <= othr.f(x))
        else:
            return Variable(lambda x: self.f(x) <= othr)

    def __gt__(self, othr):
        if isinstance(othr, Variable):
            return Variable(lambda x: self.f(x) > othr.f(x))
        else:
            return Variable(lambda x: self.f(x) > othr)

    def __ge__(self, othr):
        if isinstance(othr, Variable):
            return Variable(lambda x: self.f(x) >= othr.f(x))
        else:
            return Variable(lambda x: self.f(x) >= othr)

    def __getattr__(self, idx):
        return Variable(lambda x: self.f(x).__getattr__(idx))

def VAR(f):
    return Variable(f)

def ARGVAR(f):
    return lambda **kw: Variable(lambda x: f(x, **kw))

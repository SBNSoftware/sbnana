import operator

def _extract_othr_value(othr, *args, **kw):
    if callable(othr):
        return othr(*args, *kw)
    else:
        return othr

# convert python operator to the right-side equivalent
def flip_op(op):
    # doesn't change
    if op == operator.__and__:
        return operator.__and__
    elif op == operator.__or__:
        return operator.__or__
    elif op == operator.__eq__:
        return operator.__eq__
    # flip comps
    elif op == operator.__lt__:
        return operator.__ge__
    elif op == operator.__gt__:
        return operator.__le__
    elif op == operator.__ge__:
        return operator.__lt__
    elif op == operator.__le__:
        return operator.__gt__
    # r-operators
    elif op == operator.__sub__:
        return operator.__rsub__
    elif op == operator.__add__:
        return operator.__add__
    elif op == operator.__mul__:
        return operator.__rmul__
    elif op == operator.__truediv__:
        return operator.__rtruediv__
    elif op == operator.__floordiv__:
        return operator.__rfloordiv__


class Variable(object):
    def __init__(self, f):
        self.f = f

    def __call__(self, *args, **kw):
        return self.f(*args, **kw)

    def __getattr__(self, idx):
        return Variable(lambda *args, **kw: self.f(*args, **kw).__getattr__(idx))

    # special case for matmul (which Variable defines as a composition)
    def __matmul__(self, othr):
        return Variable(lambda *args, **kw: self(othr(*args, **kw)))

    def __rmatmul__(self, othr):
        return Variable(lambda *args, **kw: othr(self(*args, **kw)))

    # implement unary operators
    def do_uop(self, op):
        return Variable(lambda *args, **kw: op(self(*args, **kw)))

    def __invert__(self):
        return self.do_uop(operator.__invert__)

    def __neg__(self):
        return self.do_uop(operator.__neg__)

    # implement binary operators
    def do_op(self, othr, op):
        # Special case for SystVars -- defer to those  
        if isinstance(othr, SystVariable) or isinstance(othr, MultiVariable):
            return othr.do_op(self, flip_op(op))
        return Variable(lambda *args, **kw: op(self(*args, **kw), _extract_othr_value(othr, *args, **kw)))

    def __and__(self, othr):
        return self.do_op(othr, operator.__and__)

    def __or__(self, othr):
        return self.do_op(othr, operator.__or__)

    def __eq__(self, othr):
        return self.do_op(othr, operator.eq)

    def __lt__(self, othr):
        return self.do_op(othr, operator.lt)

    def __le__(self, othr):
        return self.do_op(othr, operator.le)

    def __gt__(self, othr):
        return self.do_op(othr, operator.gt)

    def __ge__(self, othr):
        return self.do_op(othr, operator.ge)

    def __sub__(self, othr):
        return self.do_op(othr, operator.sub)

    def __add__(self, othr):
        return self.do_op(othr, operator.add)

    def __mul__(self, othr):
        return self.do_op(othr, operator.mul)

    def __truediv__(self, othr):
        return self.do_op(othr, operator.truediv)

    def __floordiv__(self, othr):
        return self.do_op(othr, operator.floordiv)

    # right side versions
    def __rsub__(self, othr):
        return self.do_op(othr, operator.rsub)

    def __radd__(self, othr):
        return self.do_op(othr, operator.radd)

    def __rmul__(self, othr):
        return self.do_op(othr, operator.rmul)

    def __rtruediv__(self, othr):
        return self.do_op(othr, operator.rtruediv)

    def __rfloordiv__(self, othr):
        return self.do_op(othr, operator.rfloordiv)


class MultiVariable(Variable):
    def __init__(self, *fs):
        self.fs = fs

    def __call__(self, *args, **kw):
        return [f(*args, **kw) for f in self.fs]

    def __getattr__(self, idx):
        return type(self)(*[s.__getattr__(idx) for s in self.fs])

    # override special operators
    def __matmul__(self, othr):
        return type(self)(*[s.__matmul__(othr) for s in self.fs])

    def __rmatmul__(self, othr):
        return type(self)(*[s.__rmatmul__(othr) for s in self.fs])

    # override unary operators
    def do_uop(self, op):
        return type(self)(*[s.do_uop(op) for s in self.fs])

    # override binary operators
    def do_op(self, othr, op):
        if isinstance(othr, type(self)):
            return type(self)(*([s.do_op(o, op) for s,o in zip(self.fs, othr.fs)]))

        return type(self)(*[s.do_op(othr, op) for s in self.fs])


class SystVariable(MultiVariable):
    def __init__(self, CV, *systs):
        self.fs = [CV] + list(systs)
        self.correlated = False

    def correlate(self):
        self.correlated = True
        return self

    def nsysts(self):
        return len(self.fs) - 1

    def systs(self):
        return self.fs[1:]

    def cv(self):
        return self.fs[0]

    def __call__(self, *args, **kw):
        return self.fs[0](*args, **kw)

    # override binary operators
    def do_op(self, othr, op):
        # If othr is also a SystVariable, add together the variations
        if isinstance(othr, SystVariable):
            # are we supposed to correlate these two variables?
            if self.correlated or othr.correlated:
                return SystVariable(*([s.do_op(o, op) for s,o in zip(self.fs, othr.fs)]))

            return SystVariable(*[s.do_op(othr.cv(), op) for s in self.fs], *[self.cv().do_op(o, op) for o in othr.systs()])

        return SystVariable(*[s.do_op(othr, op) for s in self.fs])

def VAR(f):
    return Variable(f)

def ARGVAR(f):
    return lambda **kw: Variable(lambda *args: f(*args, **kw))

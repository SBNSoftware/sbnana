def _extract_othr_value(othr, x):
    if callable(othr):
        return othr(x)
    else:
        return othr

class Variable(object):
    def __init__(self, f):
        self.f = f

    def __and__(self, othr):
        return Variable(lambda x: self(x) & _extract_othr_value(othr, x))

    def __or__(self, othr):
        return Variable(lambda x: self(x) | _extract_othr_value(othr, x))

    def __matmul__(self, othr):
        return Variable(lambda *args, **kw: self(othr(*args, **kw)))

    def __rmatmul__(self, othr):
        return Variable(lambda *args, **kw: othr(self(*args, **kw)))

    def __invert__(self):
        return Variable(lambda *args, **kw: ~self(*args, **kw))

    def __call__(self, *args, **kw):
        return self.f(*args, **kw)

    def __lt__(self, othr):
        return Variable(lambda x: self.f(x) < _extract_othr_value(othr, x))

    def __le__(self, othr):
        return Variable(lambda x: self.f(x) <= _extract_othr_value(othr, x))

    def __gt__(self, othr):
        return Variable(lambda x: self.f(x) > _extract_othr_value(othr, x))

    def __ge__(self, othr):
        return Variable(lambda x: self.f(x) >= _extract_othr_value(othr, x))

    def __eq__(self, othr):
        return Variable(lambda x: self.f(x) == _extract_othr_value(othr, x))

    def __getattr__(self, idx):
        return Variable(lambda x: self.f(x).__getattr__(idx))

    def __sub__(self, othr):
        return Variable(lambda x: self.f(x) - _extract_othr_value(othr, x))
    def __add__(self, othr):
        return Variable(lambda x: self.f(x) + _extract_othr_value(othr, x))
    def __mul__(self, othr):
        return Variable(lambda x: self.f(x) * _extract_othr_value(othr, x))
    def __div__(self, othr):
        return Variable(lambda x: self.f(x) / _extract_othr_value(othr, x))
    def __truediv__(self, othr):
        return Variable(lambda x: self.f(x) / _extract_othr_value(othr, x))
    def __floordiv__(self, othr):
        return Variable(lambda x: self.f(x) // _extract_othr_value(othr, x))

def VAR(f):
    return Variable(f)

def ARGVAR(f):
    return lambda **kw: Variable(lambda x: f(x, **kw))

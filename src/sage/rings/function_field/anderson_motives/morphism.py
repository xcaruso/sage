from sage.categories.map import Map

class DrinfeldToAnderson(Map):
    def __init__(self, parent, phi):
        Map.__init__(self, parent)
        self._phi = phi
        self._motive = parent.codomain()
        self._AK = self._motive.base_combined()

    def _call_(self, f):
        phi = self._phi
        r = phi.rank()
        phiT = phi.gen()
        coords = []
        for _ in range(r):
            coords.append([])
        while f:
            f, rem = f.right_quo_rem(phiT)
            for i in range(r):
                coords[i].append(rem[i])
        coords = [self._AK(c) for c in coords]
        return self._motive(coords)

class AndersonToDrinfeld(Map):
    def __init__(self, parent, phi):
        Map.__init__(self, parent)
        self._phi = phi
        self._Ktau = parent.codomain()

    def _call_(self, x):
        phi = self._phi
        phiT = phi.gen()
        S = self._Ktau
        ans = S.zero()
        for i in range(phi.rank()):
            if x[i].denominator() != 1:
                raise ValueError("not in the Anderson motive")
            xi = x[i].numerator()
            s = S.zero()
            for j in range(xi.degree(), -1, -1):
                s = s*phiT + (S(xi[j]) << i)
            ans += s
        return ans

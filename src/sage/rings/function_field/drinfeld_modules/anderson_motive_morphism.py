# *****************************************************************************
#        Copyright (C) 2024 Xavier Caruso <xavier.caruso@normalesup.org>
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#                   http://www.gnu.org/licenses/
# *****************************************************************************

from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.structure.unique_representation import UniqueRepresentation

from sage.categories.morphism import Morphism
from sage.structure.element import coerce_binop

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

from sage.rings.function_field.drinfeld_modules.anderson_motive import AndersonMotive

class AndersonMotiveMorphism(Morphism):
    def __init__(self, parent, M):
        Morphism.__init__(self, parent)
        self._domain = X = parent.domain()
        self._codomain = Y = parent.codomain()
        if not (isinstance(X, AndersonMotive)
            and isinstance(Y, AndersonMotive)
            and X.category() == Y.category()):
            raise ValueError
        q = X._q
        AK = X._AK
        M = M.change_ring(AK)
        tM = M.parent()([AK([x**q for x in entry.list()])
                         for entry in M.list()])
        if X.matrix() * M != tM * Y.matrix():
            raise ValueError("do not commute with tau")
        self._matrix = M

    def _repr_(self):
        if self.is_endomorphism():
            return f'Endomorphism of {self._domain}'
        else:
            return f'Anderson motive morphism:\n' \
                   f'  From: {self._domain}\n'  \
                   f'  To:   {self._codomain}'

    def __hash__(self):
        return hash((self._domain, self._codomain, self._matrix))

    def matrix(self):
        return self._matrix

    def is_zero(self):
        return self._matrix.is_zero()

    def is_identity(self):
        return self._matrix.is_identity()

    @coerce_binop
    def __add__(self, other):
        return self.parent()(self.matrix() + other.matrix())

    def _composition_(self, other, H):
        return H(other.matrix() * self.matrix())

    def norm(self, ideal=True):
        nu = self._matrix.det()
        # We cast to A
        A = self.domain().function_ring()
        if ideal:
            nu = A([c.in_base() for c in nu.monic().list()])
            return A.ideal(nu)
        elif self.domain() is self.codomain():
            return A([c.in_base() for c in nu.list()])
        else:
            raise ValueError("norm is defined as an actual element only for endomorphisms")

    def trace(self):
        if self.domain() is not self.codomain():
            raise ValueError("trace is only defined for endomorphisms")
        tr = self._matrix.trace()
        # We cast to A
        A = self.domain().function_ring()
        return A([c.in_base() for c in tr.list()])

    def characteristic_polynomial(self, var='X'):
        if self.domain() is not self.codomain():
            raise ValueError("characteristic polynomial is only defined for endomorphisms")
        P = self._matrix().charpoly()
        # We cast to the correct parent
        A = self.domain().function_ring()
        parent = PolynomialRing(A, name=var)
        return parent([A([c.in_base() for c in co.list()]) for co in P.list()])

    def charpoly(self, var='X'):
        return self.characteristic_polynomial(var)

    def __eq__(self, other):
        if not isinstance(other, AndersonMotiveMorphism):
            return False
        return (self.domain() is other.domain()
            and self.codomain() is other.codomain()
            and self.matrix() == other.matrix())

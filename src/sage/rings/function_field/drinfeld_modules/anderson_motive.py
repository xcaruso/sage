r"""
This module provides facilities for manipulation Anderson motives
and computing the associated L-series.

    sage: from sage.rings.function_field.drinfeld_modules.anderson_motive import AndersonMotive


.. RUBRIC:: Construction of Anderson motives

We first need to define an Anderson motive.
This is achieved by passing in the matrix of `\tau` to the
constructor :class:`AndersonMotive`::

    sage: q = 5
    sage: Fq = GF(q)
    sage: K.<theta> = Fq[]
    sage: R.<t> = K[]
    sage: tau = matrix(2, 2, [t, 1, theta, 1])
    sage: M = AndersonMotive(tau)
    sage: M
    Anderson motive defined by the matrix
    [    t     1]
    [theta     1]

As a shortcut, the trivial Anderson motive can be defined by
just passing in the base ring::

    sage: F = AndersonMotive(R)
    sage: F
    Anderson motive defined by the matrix
    [1]

We can create (positive or negative) twists.
For example, the Carlitz module can be defined as follows::

    sage: C = F(1)
    sage: C
    Anderson motive defined by the matrix
    [t + 4*theta]

Another twist::

    sage: M(-3)
    Anderson motive defined by the matrix
    [    t/(t^3 + 2*theta*t^2 + 3*theta^2*t + 4*theta^3)     1/(t^3 + 2*theta*t^2 + 3*theta^2*t + 4*theta^3)]
    [theta/(t^3 + 2*theta*t^2 + 3*theta^2*t + 4*theta^3)     1/(t^3 + 2*theta*t^2 + 3*theta^2*t + 4*theta^3)]


.. RUBRIC:: Tensorial constructions

Basic additive and tensorial constructions are implemented.

Direct sums::

    sage: N = M + F(-2)
    sage: N
    Anderson motive defined by the matrix
    [                            t                             1                             0]
    [                        theta                             1                             0]
    [                            0                             0 1/(t^2 + 3*theta*t + theta^2)]

Tensor products::

    sage: M * M
    Anderson motive defined by the matrix
    [    t^2       t       t       1]
    [theta*t       t   theta       1]
    [theta*t   theta       t       1]
    [theta^2   theta   theta       1]

Duals::

    sage: M.dual()
    Anderson motive defined by the matrix
    [    4/(t + 4*theta)     1/(t + 4*theta)]
    [theta/(t + 4*theta)   4*t/(t + 4*theta)]

Symmetric powers::

    sage: M.symmetric_power(2)
    Anderson motive defined by the matrix
    [      t^2       2*t         1]
    [  theta*t t + theta         1]
    [  theta^2   2*theta         1]

Exterior powers::

    sage: N.exterior_power(2)
    Anderson motive defined by the matrix
    [                      t + 4*theta                                 0                                 0]
    [                                0     t/(t^2 + 3*theta*t + theta^2)     1/(t^2 + 3*theta*t + theta^2)]
    [                                0 theta/(t^2 + 3*theta*t + theta^2)     1/(t^2 + 3*theta*t + theta^2)]

As a shortcut, the method :meth:`determinant` computes the maximal
exterior power::

    sage: N.determinant()
    Anderson motive defined by the matrix
    [1/(t + 4*theta)]


.. RUBRIC:: L-series

The L-series assciated to an Anderson motive is computed
thanks to the method :meth:`Lseries`.
This method takes as input a place of `\FF_q[t]` (encoded
either by ``infinity`` for the place at infinity, an element
of `\FF_q` for a rational place ot an irreducible polynomial
for a general finite place) and a precision::

    sage: F(-3).Lseries(infinity, prec=20)
    (4*u^15 + 2*u^19 + O(u^20))*x + 1 + O(u^20)

    sage: F(-3).Lseries(0, prec=20)
    (3*u^18 + O(u^20))*x^2 + (3*u + u^5 + u^17 + O(u^20))*x + 1 + O(u^20)

    sage: F(-3).Lseries(t^2 + t + 1, prec=10)
    (u^2 + 2*u^3 + 3*u^4 + 4*u^6 + 3*u^8 + 3*u^9 + O(u^10))*x^2 + ((3*a + 4) + (3*a + 4)*u + (3*a + 4)*u^2 + (2*a + 1)*u^3 + (a + 3)*u^6 + (a + 3)*u^7 + (4*a + 2)*u^8 + O(u^10))*x + 1 + O(u^10)

In the output:

- the variable `u` corresponds to a uniformizer of the completion
  of `\mathbb F_q[t]` at the given place: when the place is infinity,
  we have `u = 1/t` whereas, when the place is finite given by
  irreducible polynomial `\mathfrak p(t)`, we have `u = \mathfrak p(t)`,

- the variable `a` is the image of `t` in the residue field,

- the variable `x` is the variable of the L-series.

It is possible to evaluate the `L`-series at a given `x` by just
passing in ``x = value``::

    sage: F(-3).Lseries(infinity, prec=100, x=1)
    1 + 4*u^15 + 2*u^19 + 4*u^23 + 4*u^35 + 2*u^39 + 4*u^43 + 4*u^55 + 2*u^59 + 4*u^63 + 4*u^75 + 2*u^79 + 4*u^83 + u^90 + 3*u^94 + 4*u^95 + u^98 + 2*u^99 + O(u^100)

We check that the L-series of a direct sum is the product of the
L-series of the summands::

    sage: N = M(-1) + F(-3)
    sage: N.Lseries(2, prec=20, x=1)
    1 + 3*u + 3*u^3 + 4*u^5 + 3*u^7 + 3*u^8 + 4*u^9 + 3*u^11 + 3*u^12 + 4*u^13 + 3*u^15 + 3*u^16 + 3*u^18 + 3*u^19 + O(u^20)
    sage: M(-1).Lseries(2, prec=20, x=1) * F(-3).Lseries(2, prec=20, x=1)
    1 + 3*u + 3*u^3 + 4*u^5 + 3*u^7 + 3*u^8 + 4*u^9 + 3*u^11 + 3*u^12 + 4*u^13 + 3*u^15 + 3*u^16 + 3*u^18 + 3*u^19 + O(u^20)


TESTS::

    sage: f = M(-7).Lseries(t^2 + t + 1, prec=100, x=1)
    sage: for prec in [10, 20, 50, 80]:
    ....:     assert(f == M(-7).Lseries(t^2 + t + 1, prec=prec, x=1))

::

    sage: f = F(-3).Lseries(1, prec=50)
    sage: for m in range(1, 10):
    ....:     Fm = AndersonMotive(matrix(R, 1, 1, [(t - theta)^m]))
    ....:     assert(f == Fm(-3-m).Lseries(1, prec=50))

::

    sage: for h in range(50):
    ....:     cond1 = h % (q - 1) == 0
    ....:     cond2 = F(-h).Lseries(0, prec=100, x=1) == 0
    ....:     assert(cond1 == cond2)


AUTHORS:

- Xavier Caruso (2023-08): initial version

"""

# *****************************************************************************
#        Copyright (C) 2024 Xavier Caruso <xavier.caruso@normalesup.org>
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#                   http://www.gnu.org/licenses/
# *****************************************************************************


from sage.misc.lazy_attribute import lazy_attribute

from sage.misc.functional import log
from sage.misc.mrange import mrange
from sage.misc.misc_c import prod
from sage.sets.set import Set

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

from sage.rings.integer_ring import ZZ
from sage.rings.infinity import Infinity
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.fraction_field import FractionField_1poly_field

from sage.matrix.constructor import matrix
from sage.matrix.special import identity_matrix, block_diagonal_matrix



class AndersonMotive(Parent, UniqueRepresentation):
    @staticmethod
    def __classcall_private__(self, category, tau=None, twist=0, check=True, normalize=True):
        K = category.base()
        AK = category.base_combined()

        # We normalize the inputs
        twist = ZZ(twist)

        if tau is None:
            tau = identity_matrix(AK, 1)
        elif tau.nrows() != tau.ncols():
            raise ValueError("tau must be a square matrix")
        else:
            tau = tau.change_ring(AK)

        if normalize:
            divisor = category.divisor()
            exponent = Infinity
            for entry in tau.list():
                if not entry:
                    continue
                e = 0
                while entry.degree() > 0 and e < exponent:
                    entry, R = entry.quo_rem(divisor)
                    if R:
                        break
                    e += 1
                exponent = e
                if exponent == 0:
                    break
            if exponent is not Infinity and exponent > 0:
                denom = divisor ** exponent
                tau = tau.parent()([entry // denom for entry in tau.list()])
                twist += exponent

        Kb = K.backend(force=True)
        if (isinstance(Kb, FractionField_1poly_field)
            and category.constant_coefficient() == K.gen()):
            from .anderson_motive_Fqtheta import AndersonMotive_Fqtheta
            cls = AndersonMotive_Fqtheta
        else:
            cls = AndersonMotive
        obj = cls.__classcall__(cls, category, tau, twist)

        if check:
            obj._dettau

        return obj

    def __init__(self, category, tau, twist):
        self._twist = twist
        self._category = category
        self._A = A = category.function_ring()
        self._t_name = A.variable_name()
        self._Fq = Fq = A.base_ring()
        self._q = Fq.cardinality()
        self._deg = ZZ(log(self._q, Fq.characteristic()))
        self._K = self._base = K = category.base()
        self._theta = category.constant_coefficient()
        self._AK = base = category.base_combined()
        self._t = base.gen()
        self._tau = tau
        Parent.__init__(self, category=category)

    @lazy_attribute
    def _dettau(self):
        det = self._tau.det()
        if det == 0:
            raise ValueError("tau does not define an Anderson motive")
        h = det.degree()
        disc, R = det.quo_rem((self._t - self._theta) ** h)
        if R:
            raise ValueError("tau does not define an Anderson motive")
        return disc[0], det.degree()

    def __repr__(self):
        return "Anderson motive of rank %s over %s" % (self.rank(), self._AK)

    def matrix(self):
        return ((self._t - self._theta) ** self._twist) * self._tau

    def rank(self):
        return self._tau.nrows()

    def twist(self, n):
        n = ZZ(n)
        return AndersonMotive(self._category, self._tau, self._twist + n, check=False, normalize=False)

    def __call__(self, n):
        return self.twist(n)

    def hodge_pink_weights(self):
        S = self._tau.smith_form(transformation=False)
        return [self._twist + S[i,i].degree() for i in range(self.rank())]

    def is_effective(self):
        return self._twist >= 0

    def tensor_product(self, other):
        if not isinstance(other, AndersonMotive):
            raise TypeError("cannot compute the tensor product of an Anderson module with something else")
        if self._category is not other._category:
            raise TypeError("Anderson modules must be in the same category")
        if self.rank() == 0:
            return self
        if other.rank() == 0:
            return other
        tau = self._tau.tensor_product(other._tau)
        tau.subdivide()
        twist = self._twist + other._twist
        return AndersonMotive(self._category, tau, twist, check=False)

    def __mul__(self, other):
        return self.tensor_product(other)

    def direct_sum(self, other):
        if not isinstance(other, AndersonMotive):
            raise TypeError("cannot compute the direct sum of an Anderson module with something else")
        if self._category is not other._category:
            raise TypeError("Anderson modules must be in the same category")
        if self.rank() == 0:
            return other
        if other.rank() == 0:
            return self
        t = self._t
        theta = self._theta
        n = self._twist - other._twist
        if n > 0:
            tau = block_diagonal_matrix([(t - theta)**n * self._tau, other._tau])
            twist = other._twist
        else:
            tau = block_diagonal_matrix([self._tau, (t - theta)**(-n) * other._tau])
            twist = self._twist
        tau.subdivide()
        return AndersonMotive(self._category, tau, twist, check=False, normalize=(n == 0))

    def __add__(self, other):
        return self.direct_sum(other)

    def dual(self):
        disc, deg = self._dettau
        scalar = self._K(~disc)
        tau = (-1)**deg * scalar * self._tau.adjugate()
        twist = -deg - self._twist
        return AndersonMotive(self._category, tau, twist, check=False, normalize=True)

    def symmetric_power(self, n):
        # TODO Optimization:
        # instead of running over mrange([r] * n),
        # run over E and then consider all permutations
        n = ZZ(n)
        if n < 0:
            raise ValueError("exponent must be nonnegative")
        if n == 0:
            return AndersonMotive(self._AK)
        AK = self._AK
        r = self.rank()
        S = Set(range(r+n-1))
        E = [ ]
        Ed = { }
        for s in S.subsets(r-1):
            es = (n + r - 2 - s[-1]) * (0,)
            for i in range(r-2, 0, -1):
                es += (s[i] - s[i-1] - 1) * (r-1-i,)
            es += s[0] * (r-1,)
            Ed[es] = len(E)
            E.append(es)
        rows = [ ]
        for es in E:
            row = [AK.zero() for i in range(len(E))]
            for fs in mrange([r] * n):
                j = Ed[tuple(sorted(fs))]
                row[j] += prod(self._tau[es[i], fs[i]] for i in range(n))
            rows.append(row)
        tau = matrix(rows)
        twist = n * self._twist
        return AndersonMotive(self._category, tau, twist, check=False, normalize=True)

    def exterior_power(self, n):
        n = ZZ(n)
        if n < 0:
            raise ValueError("exponent must be nonnegative")
        if n == 0:
            return AndersonMotive(self._AK)
        r = self.rank()
        if n > r:
            return AndersonMotive(identity_matrix(self._AK, 0))
        I = Set(range(r))
        In = [ sorted(list(J)) for J in I.subsets(n) ]
        rows = [ ]
        for J1 in In:
            row = [ ]
            for J2 in In:
                M = self._tau.matrix_from_rows_and_columns(J1, J2)
                row.append(M.det())
            rows.append(row)
        tau = matrix(rows)
        twist = n * self._twist
        return AndersonMotive(self._category, tau, twist, check=False, normalize=True)

    def determinant(self):
        det, deg = self._dettau
        det *= (self._t - self._theta)**deg
        tau = matrix(self._AK, 1, 1, [det])
        twist = self.rank() * self._twist
        return AndersonMotive(self._category, tau, twist, check=False, normalize=True)

    def _Hom_(self, other, category):
        from sage.rings.function_field.drinfeld_modules.anderson_motive_homset import AndersonMotiveHomset
        return AndersonMotiveHomset(self, other, category)

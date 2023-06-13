r"""
Ring of Laurent Polynomials

If `R` is a commutative ring, then the ring of Laurent polynomials in `n`
variables over `R` is `R[x_1^{\pm 1}, x_2^{\pm 1}, \ldots, x_n^{\pm 1}]`.
We implement it as a quotient ring

.. MATH::

    R[x_1, y_1, x_2, y_2, \ldots, x_n, y_n] / (x_1 y_1 - 1, x_2 y_2 - 1, \ldots, x_n y_n - 1).

TESTS::

    sage: P.<q> = LaurentPolynomialRing(QQ)
    sage: qi = q^(-1)
    sage: qi in P
    True
    sage: P(qi)
    q^-1

    sage: A.<Y> = QQ[]
    sage: R.<X> = LaurentPolynomialRing(A)
    sage: matrix(R,2,2,[X,0,0,1])                                                       # optional - sage.modules
    [X 0]
    [0 1]

AUTHORS:

- David Roe (2008-2-23): created
- David Loeffler (2009-07-10): cleaned up docstrings
"""
# ****************************************************************************
#       Copyright (C) 2008 David Roe <roed@math.harvard.edu>,
#                          William Stein <wstein@gmail.com>,
#                          Mike Hansen <mhansen@gmail.com>
#                          Vincent Delecroix <20100.delecroix@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.lazy_import import LazyImport
from sage.rings.polynomial.laurent_polynomial import LaurentPolynomial, LaurentPolynomial_univariate
from sage.rings.polynomial.laurent_polynomial_ring_base import LaurentPolynomialRing_generic
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.structure.element import parent


def is_LaurentPolynomialRing(R):
    """
    Return ``True`` if and only if R is a Laurent polynomial ring.

    EXAMPLES::

        sage: from sage.rings.polynomial.laurent_polynomial_ring import is_LaurentPolynomialRing
        sage: P = PolynomialRing(QQ, 2, 'x')
        sage: is_LaurentPolynomialRing(P)
        doctest:warning...
        DeprecationWarning: is_LaurentPolynomialRing is deprecated; use isinstance(...,
        sage.rings.polynomial.laurent_polynomial_ring_base.LaurentPolynomialRing_generic) instead
        See https://github.com/sagemath/sage/issues/35229 for details.
        False

        sage: R = LaurentPolynomialRing(QQ,3,'x')                                       # optional - sage.modules
        sage: is_LaurentPolynomialRing(R)                                               # optional - sage.modules
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(35229,
                "is_LaurentPolynomialRing is deprecated; use "
                "isinstance(..., sage.rings.polynomial.laurent_polynomial_ring_base."
                "LaurentPolynomialRing_generic) instead")
    return isinstance(R, LaurentPolynomialRing_generic)


_cache = {}
def LaurentPolynomialRing(base_ring, *args, **kwds):
    r"""
    Return the globally unique univariate or multivariate Laurent polynomial
    ring with given properties and variable name or names.

    There are four ways to call the Laurent polynomial ring constructor:

    1. ``LaurentPolynomialRing(base_ring, name,    sparse=False)``
    2. ``LaurentPolynomialRing(base_ring, names,   order='degrevlex')``
    3. ``LaurentPolynomialRing(base_ring, name, n, order='degrevlex')``
    4. ``LaurentPolynomialRing(base_ring, n, name, order='degrevlex')``

    The optional arguments ``sparse`` and ``order`` *must* be explicitly
    named, and the other arguments must be given positionally.

    INPUT:

    - ``base_ring`` -- a commutative ring
    - ``name`` -- a string
    - ``names`` -- a list or tuple of names, or a comma separated string
    - ``n`` -- a positive integer
    - ``sparse`` -- bool (default: False), whether or not elements are sparse
    - ``order`` -- string or
      :class:`~sage.rings.polynomial.term_order.TermOrder`, e.g.,

      - ``'degrevlex'`` (default) -- degree reverse lexicographic
      - ``'lex'`` -- lexicographic
      - ``'deglex'`` -- degree lexicographic
      - ``TermOrder('deglex',3) + TermOrder('deglex',3)`` -- block ordering

    OUTPUT:

    ``LaurentPolynomialRing(base_ring, name, sparse=False)`` returns a
    univariate Laurent polynomial ring; all other input formats return a
    multivariate Laurent polynomial ring.

    UNIQUENESS and IMMUTABILITY: In Sage there is exactly one
    single-variate Laurent polynomial ring over each base ring in each choice
    of variable and sparseness.  There is also exactly one multivariate
    Laurent polynomial ring over each base ring for each choice of names of
    variables and term order.

    ::

        sage: R.<x,y> = LaurentPolynomialRing(QQ, 2); R                                 # optional - sage.modules
        Multivariate Laurent Polynomial Ring in x, y over Rational Field
        sage: f = x^2 - 2*y^-2                                                          # optional - sage.modules

    You can't just globally change the names of those variables.
    This is because objects all over Sage could have pointers to
    that polynomial ring.

    ::

        sage: R._assign_names(['z','w'])                                                # optional - sage.modules
        Traceback (most recent call last):
        ...
        ValueError: variable names cannot be changed after object creation.


    EXAMPLES:

    1. ``LaurentPolynomialRing(base_ring, name, sparse=False)``

       ::

           sage: LaurentPolynomialRing(QQ, 'w')
           Univariate Laurent Polynomial Ring in w over Rational Field

       Use the diamond brackets notation to make the variable
       ready for use after you define the ring::

           sage: R.<w> = LaurentPolynomialRing(QQ)
           sage: (1 + w)^3
           1 + 3*w + 3*w^2 + w^3

       You must specify a name::

           sage: LaurentPolynomialRing(QQ)
           Traceback (most recent call last):
           ...
           TypeError: you must specify the names of the variables

           sage: R.<abc> = LaurentPolynomialRing(QQ, sparse=True); R
           Univariate Laurent Polynomial Ring in abc over Rational Field

           sage: R.<w> = LaurentPolynomialRing(PolynomialRing(GF(7),'k')); R            # optional - sage.rings.finite_rings
           Univariate Laurent Polynomial Ring in w over
            Univariate Polynomial Ring in k over Finite Field of size 7

       Rings with different variables are different::

           sage: LaurentPolynomialRing(QQ, 'x') == LaurentPolynomialRing(QQ, 'y')
           False

    2. ``LaurentPolynomialRing(base_ring, names,   order='degrevlex')``

       ::

           sage: R = LaurentPolynomialRing(QQ, 'a,b,c'); R                              # optional - sage.modules
           Multivariate Laurent Polynomial Ring in a, b, c over Rational Field

           sage: S = LaurentPolynomialRing(QQ, ['a','b','c']); S                        # optional - sage.modules
           Multivariate Laurent Polynomial Ring in a, b, c over Rational Field

           sage: T = LaurentPolynomialRing(QQ, ('a','b','c')); T                        # optional - sage.modules
           Multivariate Laurent Polynomial Ring in a, b, c over Rational Field

       All three rings are identical.

       ::

           sage: (R is S) and  (S is T)                                                 # optional - sage.modules
           True

       There is a unique Laurent polynomial ring with each term order::

           sage: R = LaurentPolynomialRing(QQ, 'x,y,z', order='degrevlex'); R           # optional - sage.modules
           Multivariate Laurent Polynomial Ring in x, y, z over Rational Field
           sage: S = LaurentPolynomialRing(QQ, 'x,y,z', order='invlex'); S              # optional - sage.modules
           Multivariate Laurent Polynomial Ring in x, y, z over Rational Field
           sage: S is LaurentPolynomialRing(QQ, 'x,y,z', order='invlex')                # optional - sage.modules
           True
           sage: R == S                                                                 # optional - sage.modules
           False


    3. ``LaurentPolynomialRing(base_ring, name, n, order='degrevlex')``

       If you specify a single name as a string and a number of
       variables, then variables labeled with numbers are created.

       ::

           sage: LaurentPolynomialRing(QQ, 'x', 10)                                     # optional - sage.modules
           Multivariate Laurent Polynomial Ring in x0, x1, x2, x3, x4, x5, x6, x7, x8, x9
            over Rational Field

           sage: LaurentPolynomialRing(GF(7), 'y', 5)                                   # optional - sage.modules sage.rings.finite_rings
           Multivariate Laurent Polynomial Ring in y0, y1, y2, y3, y4
            over Finite Field of size 7

           sage: LaurentPolynomialRing(QQ, 'y', 3, sparse=True)                         # optional - sage.modules
           Multivariate Laurent Polynomial Ring in y0, y1, y2 over Rational Field

       By calling the
       :meth:`~sage.structure.category_object.CategoryObject.inject_variables`
       method, all those variable names are available for interactive use::

           sage: R = LaurentPolynomialRing(GF(7), 15, 'w'); R                           # optional - sage.modules sage.rings.finite_rings
           Multivariate Laurent Polynomial Ring in w0, w1, w2, w3, w4, w5, w6, w7,
            w8, w9, w10, w11, w12, w13, w14 over Finite Field of size 7
           sage: R.inject_variables()                                                   # optional - sage.modules sage.rings.finite_rings
           Defining w0, w1, w2, w3, w4, w5, w6, w7, w8, w9, w10, w11, w12, w13, w14
           sage: (w0 + 2*w8 + w13)^2                                                    # optional - sage.modules sage.rings.finite_rings
           w0^2 + 4*w0*w8 + 4*w8^2 + 2*w0*w13 + 4*w8*w13 + w13^2
    """
    from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
    from sage.rings.polynomial.multi_polynomial_ring_base import is_MPolynomialRing

    R = PolynomialRing(base_ring, *args, **kwds)
    if R in _cache:
        return _cache[R]   # put () here to re-enable weakrefs

    if is_PolynomialRing(R):
        # univariate case
        P = LaurentPolynomialRing_univariate(R)
    else:
        assert is_MPolynomialRing(R)
        P = LaurentPolynomialRing_mpair(R)

    _cache[R] = P
    return P

def _split_dict_(D, indices, group_by=None):
    r"""
    Split the dictionary ``D`` by ``indices`` and ``group_by``.

    INPUT:

    - ``D`` -- a dictionary.

    - ``indices`` -- a tuple or list of nonnegative integers.

    - ``group_by`` -- a tuple or list of nonnegative integers.
      If this is ``None`` (default), then no grouping is done.

    OUTPUT:

    A dictionary.

    TESTS::

        sage: from sage.rings.polynomial.laurent_polynomial_ring import _split_dict_
        sage: D = {(0,0,0,0): 'a', (1,0,0,0): 'b',
        ....:      (1,0,0,2): 'c', (1,2,0,3): 'd'}
        sage: _split_dict_(D, [1,0,3])
        {(0, 0, 0): 'a', (0, 1, 0): 'b', (0, 1, 2): 'c', (2, 1, 3): 'd'}
        sage: _split_dict_(D, [2,3], [0,1])
        {(0, 0): {(0, 0): 'a'},
         (1, 0): {(0, 0): 'b', (0, 2): 'c'},
         (1, 2): {(0, 3): 'd'}}
        sage: _split_dict_(D, [3,1], [0])
        {(0,): {(0, 0): 'a'}, (1,): {(0, 0): 'b', (2, 0): 'c', (3, 2): 'd'}}

        sage: _split_dict_(D, [0,None,1,3])
        {(0, 0, 0, 0): 'a', (1, 0, 0, 0): 'b',
         (1, 0, 0, 2): 'c', (1, 0, 2, 3): 'd'}
        sage: _split_dict_(D, [0,1], [None,3,None])
        {(0, 0, 0): {(0, 0): 'a', (1, 0): 'b'},
         (0, 2, 0): {(1, 0): 'c'},
         (0, 3, 0): {(1, 2): 'd'}}
        sage: _split_dict_(D, [None,3,1], [0,None])
        {(0, 0): {(0, 0, 0): 'a'},
         (1, 0): {(0, 0, 0): 'b', (0, 2, 0): 'c',
                     (0, 3, 2): 'd'}}

        sage: _split_dict_(D, [0,1])
        Traceback (most recent call last):
        ...
        SplitDictError: split not possible
        sage: _split_dict_(D, [0], [1])
        Traceback (most recent call last):
        ...
        SplitDictError: split not possible
        sage: _split_dict_({}, [])
        {}
    """
    if not D:
        return {}
    if group_by is None:
        group_by = tuple()

    class SplitDictError(ValueError):
        pass

    def get(T, i):
        return T[i] if i is not None else 0

    def extract(T, indices):
        return tuple(get(T, i) for i in indices)

    remaining = sorted(set(range(len(next(iter(D)))))
                       - set(indices) - set(group_by))
    result = {}
    for K, V in D.items():
        if not all(r == 0 for r in extract(K, remaining)):
            raise SplitDictError('split not possible')
        G = extract(K, group_by)
        I = extract(K, indices)
        result.setdefault(G, dict()).update({I: V})
    if not group_by:
        return result.popitem()[1]
    else:
        return result

def _split_laurent_polynomial_dict_(P, M, d):
    r"""
    Helper function for splitting a multivariate Laurent polynomial
    during conversion.

    INPUT:

    - ``P`` -- the parent to which we want to convert.

    - ``M`` -- the parent from which we want to convert.

    - ``d`` -- a dictionary mapping tuples (representing the exponents)
      to their coefficients. This is the dictionary corresponding to
      an element of ``M``.

    OUTPUT:

    A dictionary corresponding to an element of ``P``.

    TESTS::

        sage: L.<a, b, c, d> = LaurentPolynomialRing(ZZ)                                # optional - sage.modules
        sage: M = LaurentPolynomialRing(ZZ, 'c, d')                                     # optional - sage.modules
        sage: N = LaurentPolynomialRing(M, 'a, b')                                      # optional - sage.modules
        sage: M(c/d + 1/c)  # indirect doctest                                          # optional - sage.modules
        c*d^-1 + c^-1
        sage: N(a + b/c/d + 1/b)  # indirect doctest                                    # optional - sage.modules
        a + (c^-1*d^-1)*b + b^-1
    """
    vars_P = P.variable_names()
    vars_M = M.variable_names()
    if not set(vars_M) & set(vars_P):
        raise TypeError('no common variables')

    def index(T, value):
        try:
            return T.index(value)
        except ValueError:
            return None

    def value(d, R):
        assert d
        if len(d) == 1:
            k, v = next(iter(d.items()))
            if all(i == 0 for i in k):
                return R(v)
        return R(M(d))

    group_by = tuple(index(vars_M, var) for var in vars_P)
    indices = list(range(len(vars_M)))
    for g in group_by:
        if g is not None:
            indices[g] = None
    D = _split_dict_(d, indices, group_by)
    try:
        return {k: value(v, P.base_ring()) for k, v in D.items()}
    except (ValueError, TypeError):
        pass
    return sum(P({k: 1}) * value(v, P) for k, v in D.items()).dict()

def from_fraction_field(L, x):
    r"""
    Helper function to construct a Laurent polynomial from an element of its
    parent's fraction field.

    INPUT:

    - ``L`` -- an instance of :class:`LaurentPolynomialRing_generic`
    - ``x`` -- an element of the fraction field of ``L``

    OUTPUT:

    An instance of the element class of ``L``. If the denominator fails to be
    a unit in ``L`` an error is raised.

    EXAMPLES::

        sage: from sage.rings.polynomial.laurent_polynomial_ring import from_fraction_field
        sage: L.<x, y> = LaurentPolynomialRing(ZZ)                                      # optional - sage.modules
        sage: F = L.fraction_field()                                                    # optional - sage.modules
        sage: xi = F(~x)                                                                # optional - sage.modules
        sage: from_fraction_field(L, xi) == ~x                                          # optional - sage.modules
        True
    """
    d = L(x.denominator())
    if d.is_unit():
        n = L(x.numerator())
        return n * d.inverse_of_unit()
    else:
        raise TypeError("fraction must have unit denominator")


class LaurentPolynomialRing_univariate(LaurentPolynomialRing_generic):
    def __init__(self, R):
        """
        EXAMPLES::

            sage: L = LaurentPolynomialRing(QQ,'x')
            sage: type(L)
            <class 'sage.rings.polynomial.laurent_polynomial_ring.LaurentPolynomialRing_univariate_with_category'>
            sage: L == loads(dumps(L))
            True


        TESTS::

            sage: TestSuite(LaurentPolynomialRing(Zmod(4), 'y')).run()
            sage: TestSuite(LaurentPolynomialRing(ZZ, 'u')).run()
            sage: TestSuite(LaurentPolynomialRing(Zmod(4)['T'], 'u')).run()
        """
        if R.ngens() != 1:
            raise ValueError("must be 1 generator")
        LaurentPolynomialRing_generic.__init__(self, R)

    Element = LaurentPolynomial_univariate

    def _repr_(self):
        """
        TESTS::

            sage: LaurentPolynomialRing(QQ,'x')  # indirect doctest
            Univariate Laurent Polynomial Ring in x over Rational Field
        """
        return "Univariate Laurent Polynomial Ring in %s over %s"%(self._R.variable_name(), self._R.base_ring())

    def _element_constructor_(self, x):
        """
        EXAMPLES::

            sage: L = LaurentPolynomialRing(QQ, 'x')
            sage: L(1/2)
            1/2

            sage: L(x + 3/x)
            3*x^-1 + x

        ::

            sage: L(exp(x))
            Traceback (most recent call last):
            ...
            TypeError: unable to convert e^x to a rational

        ::

            sage: U = LaurentPolynomialRing(QQ, 'a')
            sage: V = LaurentPolynomialRing(QQ, 'c')
            sage: L.<a, b, c, d> = LaurentPolynomialRing(QQ)
            sage: M = LaurentPolynomialRing(QQ, 'c, d')                                 # optional - sage.modules
            sage: Mc, Md = M.gens()                                                     # optional - sage.modules
            sage: N = LaurentPolynomialRing(M, 'a, b')                                  # optional - sage.modules
            sage: Na, Nb = N.gens()                                                     # optional - sage.modules
            sage: U(Na)                                                                 # optional - sage.modules
            a
            sage: V(Mc)                                                                 # optional - sage.modules
            c

            sage: M(L(0))                                                               # optional - sage.modules
            0
            sage: N(L(0))                                                               # optional - sage.modules
            0
            sage: L(M(0))                                                               # optional - sage.modules
            0
            sage: L(N(0))                                                               # optional - sage.modules
            0

        ::

            sage: A.<a> = LaurentPolynomialRing(QQ)
            sage: B.<b> = LaurentPolynomialRing(A)
            sage: B(a)
            a
            sage: C.<c> = LaurentPolynomialRing(B)
            sage: B(C(b))
            b
            sage: D.<d, e> = LaurentPolynomialRing(B)                                   # optional - sage.modules
            sage: B(D(b))                                                               # optional - sage.modules
            b

        TESTS:

        Check that conversion back from fraction field does work (:trac:`26425`)::

            sage: R.<t> = LaurentPolynomialRing(ZZ)
            sage: F = FractionField(R)
            sage: R(F(25/(5*t**2)))
            5*t^-2
            sage: R(F(1/(1+t**2)))
            Traceback (most recent call last):
            ...
            TypeError: fraction must have unit denominator
        """
        from sage.structure.element import Expression
        from sage.rings.fraction_field_element import FractionFieldElement
        from sage.rings.localization import LocalizationElement
        if isinstance(x, Expression):
            return x.laurent_polynomial(ring=self)

        elif isinstance(x, LaurentPolynomial):
            P = x.parent()
            if set(self.variable_names()) & set(P.variable_names()):
                if isinstance(x, LaurentPolynomial_univariate):
                    d = {(k,): v for k, v in x.dict().items()}
                else:
                    d = x.dict()
                x = _split_laurent_polynomial_dict_(self, P, d)
                x = {k[0]: v for k, v in x.items()}
            elif P is self.base_ring():
                x = {0: x}
            elif x.is_constant() and self.has_coerce_map_from(x.parent().base_ring()):
                return self(x.constant_coefficient())
            elif len(self.variable_names()) == len(P.variable_names()):
                x = x.dict()

        elif isinstance(x, FractionFieldElement):
            # since the field of fraction of self is defined corresponding to
            # the polynomial ring of self the conversion of its elements back
            # must be treated separately (:trac:`26425`).
            return from_fraction_field(self, x)

        elif isinstance(x, LocalizationElement):
            # see :trac:`33477`.
            F = self.fraction_field()
            return from_fraction_field(self, F(x))

        return self.element_class(self, x)

    def __reduce__(self):
        """
        Used in pickling.

        EXAMPLES::

            sage: L = LaurentPolynomialRing(QQ, 'x')
            sage: loads(dumps(L)) == L
            True
        """
        return LaurentPolynomialRing_univariate, (self._R,)


class LaurentPolynomialRing_mpair(LaurentPolynomialRing_generic):
    def __init__(self, R):
        """
        EXAMPLES::

            sage: L = LaurentPolynomialRing(QQ,2,'x')                                   # optional - sage.modules
            sage: type(L)                                                               # optional - sage.modules
            <class
            'sage.rings.polynomial.laurent_polynomial_ring.LaurentPolynomialRing_mpair_with_category'>
            sage: L == loads(dumps(L))                                                  # optional - sage.modules
            True
        """
        if R.ngens() <= 0:
            raise ValueError("n must be positive")
        if not R.base_ring().is_integral_domain():
            raise ValueError("base ring must be an integral domain")
        LaurentPolynomialRing_generic.__init__(self, R)

    Element = LazyImport('sage.rings.polynomial.laurent_polynomial_mpair', 'LaurentPolynomial_mpair')

    def _repr_(self):
        """
        TESTS::

            sage: LaurentPolynomialRing(QQ,2,'x').__repr__()                            # optional - sage.modules
            'Multivariate Laurent Polynomial Ring in x0, x1 over Rational Field'
            sage: LaurentPolynomialRing(QQ,1,'x').__repr__()                            # optional - sage.modules
            'Multivariate Laurent Polynomial Ring in x over Rational Field'
        """
        return "Multivariate Laurent Polynomial Ring in %s over %s"%(", ".join(self._R.variable_names()), self._R.base_ring())

    def monomial(self, *args):
        r"""
        Return the monomial whose exponents are given in argument.

        EXAMPLES::

            sage: L = LaurentPolynomialRing(QQ, 'x', 2)                                 # optional - sage.modules
            sage: L.monomial(-3, 5)                                                     # optional - sage.modules
            x0^-3*x1^5
            sage: L.monomial(1, 1)                                                      # optional - sage.modules
            x0*x1
            sage: L.monomial(0, 0)                                                      # optional - sage.modules
            1
            sage: L.monomial(-2, -3)                                                    # optional - sage.modules
            x0^-2*x1^-3

            sage: x0, x1 = L.gens()                                                     # optional - sage.modules
            sage: L.monomial(-1, 2) == x0^-1 * x1^2                                     # optional - sage.modules
            True

            sage: L.monomial(1, 2, 3)                                                   # optional - sage.modules
            Traceback (most recent call last):
            ...
            TypeError: tuple key must have same length as ngens
        """
        if len(args) != self.ngens():
            raise TypeError("tuple key must have same length as ngens")

        from sage.rings.polynomial.polydict import ETuple
        m = ETuple(args, int(self.ngens()))
        return self.element_class(self, self.polynomial_ring().one(), m)

    def _element_constructor_(self, x, mon=None):
        """
        EXAMPLES::

            sage: L = LaurentPolynomialRing(QQ,2,'x')                                   # optional - sage.modules
            sage: L(1/2)                                                                # optional - sage.modules
            1/2

            sage: M = LaurentPolynomialRing(QQ, 'x, y')                                 # optional - sage.modules
            sage: var('x, y')                                                           # optional - sage.modules
            (x, y)
            sage: M(x/y + 3/x)                                                          # optional - sage.modules
            x*y^-1 + 3*x^-1

        ::

            sage: M(exp(x))                                                             # optional - sage.modules
            Traceback (most recent call last):
            ...
            TypeError: unable to convert e^x to a rational

        ::

            sage: L.<a, b, c, d> = LaurentPolynomialRing(QQ)                            # optional - sage.modules
            sage: M = LaurentPolynomialRing(QQ, 'c, d')                                 # optional - sage.modules
            sage: Mc, Md = M.gens()                                                     # optional - sage.modules
            sage: N = LaurentPolynomialRing(M, 'a, b')                                  # optional - sage.modules
            sage: Na, Nb = N.gens()                                                     # optional - sage.modules
            sage: M(c/d)                                                                # optional - sage.modules
            c*d^-1
            sage: N(a*b/c/d)                                                            # optional - sage.modules
            (c^-1*d^-1)*a*b
            sage: N(c/d)                                                                # optional - sage.modules
            c*d^-1
            sage: L(Mc)                                                                 # optional - sage.modules
            c
            sage: L(Nb)                                                                 # optional - sage.modules
            b

            sage: M(L(0))                                                               # optional - sage.modules
            0
            sage: N(L(0))                                                               # optional - sage.modules
            0
            sage: L(M(0))                                                               # optional - sage.modules
            0
            sage: L(N(0))                                                               # optional - sage.modules
            0

            sage: U = LaurentPolynomialRing(QQ, 'a')
            sage: Ua = U.gen()
            sage: V = LaurentPolynomialRing(QQ, 'c')
            sage: Vc = V.gen()
            sage: L(Ua)                                                                 # optional - sage.modules
            a
            sage: L(Vc)                                                                 # optional - sage.modules
            c
            sage: N(Ua)                                                                 # optional - sage.modules
            a
            sage: M(Vc)                                                                 # optional - sage.modules
            c

            sage: P = LaurentPolynomialRing(QQ, 'a, b')                                 # optional - sage.modules
            sage: Q = LaurentPolynomialRing(P, 'c, d')                                  # optional - sage.modules
            sage: Q(P.0)                                                                # optional - sage.modules
            a

        ::

            sage: A.<a> = LaurentPolynomialRing(QQ)
            sage: B.<b> = LaurentPolynomialRing(A)
            sage: C = LaurentPolynomialRing(QQ, 'a, b')                                 # optional - sage.modules
            sage: C(B({1: a}))                                                          # optional - sage.modules
            a*b
            sage: D.<d, e> = LaurentPolynomialRing(B)                                   # optional - sage.modules
            sage: F.<f, g> = LaurentPolynomialRing(D)                                   # optional - sage.modules
            sage: D(F(d*e))                                                             # optional - sage.modules
            d*e

        ::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: R.<x,y,z> = LaurentPolynomialRing(QQ)                                 # optional - sage.modules
            sage: mon = ETuple({}, int(3))                                              # optional - sage.modules
            sage: P = R.polynomial_ring()                                               # optional - sage.modules
            sage: R(sum(P.gens()), mon)                                                 # optional - sage.modules
            x + y + z
            sage: R(sum(P.gens()), (-1,-1,-1))                                          # optional - sage.modules
            y^-1*z^-1 + x^-1*z^-1 + x^-1*y^-1

        ::

            sage: RL = R.localization(x+1)                                              # optional - sage.modules
            sage: xi = RL(~x)                                                           # optional - sage.modules
            sage: R(xi) == ~x     # indirect doctests                                   # optional - sage.modules
            True
        """
        from sage.structure.element import Expression
        from sage.rings.fraction_field_element import FractionFieldElement
        from sage.rings.localization import LocalizationElement

        if mon is not None:
            return self.element_class(self, x, mon)

        P = parent(x)
        if P is self.polynomial_ring():
            from sage.rings.polynomial.polydict import ETuple
            return self.element_class( self, x, mon=ETuple({}, int(self.ngens())) )

        elif isinstance(x, Expression):
            return x.laurent_polynomial(ring=self)

        elif isinstance(x, LaurentPolynomial):
            if self.variable_names() == P.variable_names():
                # No special processing needed here;
                #   handled by LaurentPolynomial_mpair.__init__
                pass
            elif set(self.variable_names()) & set(P.variable_names()):
                if isinstance(x, LaurentPolynomial_univariate):
                    d = {(k,): v for k, v in x.dict().items()}
                else:
                    d = x.dict()
                x = _split_laurent_polynomial_dict_(self, P, d)
            elif P is self.base_ring():
                from sage.rings.polynomial.polydict import ETuple
                mz = ETuple({}, int(self.ngens()))
                return self.element_class(self, {mz: x}, mz)
            elif x.is_constant() and self.has_coerce_map_from(P.base_ring()):
                return self(x.constant_coefficient())
            elif len(self.variable_names()) == len(P.variable_names()):
                x = x.dict()

        elif isinstance(x, FractionFieldElement):
            # since the field of fraction of self is defined corresponding to
            # the polynomial ring of self the conversion of its elements back
            # must be treated separately (:trac:`33477`).
            return from_fraction_field(self, x)

        elif isinstance(x, LocalizationElement):
            # see :trac:`33477`.
            F = self.fraction_field()
            return from_fraction_field(self, F(x))

        return self.element_class(self, x)

    def __reduce__(self):
        """
        Used in pickling.

        EXAMPLES::

            sage: L = LaurentPolynomialRing(QQ, 2, 'x')                                 # optional - sage.modules
            sage: loads(dumps(L)) == L                                                  # optional - sage.modules
            True
        """
        return LaurentPolynomialRing_mpair, (self._R,)

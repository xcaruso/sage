r"""
Localization

Localization is an important ring construction tool. Whenever you have
to extend a given integral domain such that it contains the inverses
of a finite set of elements but should allow non injective homomorphic
images this construction will be needed.

Localization plays a quite important role in algebraic geometry as
it corresponds to taking open subsets of the corresponding schemes.

EXAMPLES::

.. RUBRIC:: Ariki-Koike algebras
    sage: LZ = Localization(ZZ, (5,11))
    sage: m = matrix(LZ, [[5, 7], [0,11]])                                              # optional - sage.modules
    sage: m.parent()                                                                    # optional - sage.modules
    Full MatrixSpace of 2 by 2 dense matrices over Integer Ring localized at (5, 11)
    sage: ~m      # parent of inverse is different: see documentation of m.__invert__   # optional - sage.modules
    [  1/5 -7/55]
    [    0  1/11]
    sage: _.parent()                                                                    # optional - sage.modules
    Full MatrixSpace of 2 by 2 dense matrices over Rational Field
    sage: mi = matrix(LZ, ~m)                                                           # optional - sage.modules
    sage: mi.parent()                                                                   # optional - sage.modules
    Full MatrixSpace of 2 by 2 dense matrices over Integer Ring localized at (5, 11)
    sage: mi == ~m                                                                      # optional - sage.modules
    True

The next example defines the most general ring containing the coefficients
of the irreducible representations of the Ariki-Koike algebra corresponding
to the three colored permutations on three elements::

    sage: R.<u0, u1, u2, q> = ZZ[]
    sage: u = [u0, u1, u2]
    sage: S = Set(u)
    sage: I = S.cartesian_product(S)
    sage: add_units = u + [q, q + 1] + [ui - uj for ui, uj in I if ui != uj]
    sage: add_units += [q*ui - uj for ui, uj in I if ui != uj]
    sage: L = R.localization(tuple(add_units)); L                                       # optional - sage.libs.pari
    Multivariate Polynomial Ring in u0, u1, u2, q over Integer Ring localized at
     (q, q + 1, u2, u1, u1 - u2, u0, u0 - u2, u0 - u1, u2*q - u1, u2*q - u0,
      u1*q - u2, u1*q - u0, u0*q - u2, u0*q - u1)

We define the representation matrices (of one of the three dimensional
irreducible representations)::

    sage: m1 = matrix(L, [[u1, 0, 0], [0, u0, 0], [0, 0, u0]])                          # optional - sage.modules
    sage: m2 = matrix(L, [[(u0*q - u0)/(u0 - u1), (u0*q - u1)/(u0 - u1), 0],            # optional - sage.modules
    ....:                 [(-u1*q + u0)/(u0 - u1), (-u1*q + u1)/(u0 - u1), 0],
    ....:                 [0, 0, -1]])
    sage: m3 = matrix(L, [[-1, 0, 0],                                                   # optional - sage.modules
    ....:                 [0, u0*(1 - q)/(u1*q - u0), q*(u1 - u0)/(u1*q - u0)],
    ....:                 [0, (u1*q^2 - u0)/(u1*q - u0), (u1*q^ 2 - u1*q)/(u1*q - u0)]])
    sage: m1.base_ring() == L                                                           # optional - sage.modules
    True

We check the relations of the Ariki-Koike algebra::

    sage: m1*m2*m1*m2 == m2*m1*m2*m1                                                    # optional - sage.modules
    True
    sage: m2*m3*m2 == m3*m2*m3                                                          # optional - sage.modules
    True
    sage: m1*m3 == m3*m1                                                                # optional - sage.modules
    True
    sage: m1**3 - (u0+u1+u2)*m1**2 + (u0*u1+u0*u2+u1*u2)*m1 - u0*u1*u2 == 0             # optional - sage.modules
    True
    sage: m2**2 - (q-1)*m2 - q == 0                                                     # optional - sage.modules
    True
    sage: m3**2 - (q-1)*m3 - q == 0                                                     # optional - sage.modules
    True
    sage: ~m1 in m1.parent()                                                            # optional - sage.modules
    True
    sage: ~m2 in m2.parent()                                                            # optional - sage.modules
    True
    sage: ~m3 in m3.parent()                                                            # optional - sage.modules
    True

We obtain specializations in positive characteristic::

    sage: Fp = GF(17)                                                                   # optional - sage.rings.finite_rings
    sage: f = L.hom((3,5,7,11), codomain=Fp); f                                         # optional - sage.rings.finite_rings
    Ring morphism:
      From: Multivariate Polynomial Ring in u0, u1, u2, q over Integer Ring localized at
            (q, q + 1, u2, u1, u1 - u2, u0, u0 - u2, u0 - u1, u2*q - u1, u2*q - u0,
             u1*q - u2, u1*q - u0, u0*q - u2, u0*q - u1)
      To:   Finite Field of size 17
      Defn: u0 |--> 3
            u1 |--> 5
            u2 |--> 7
            q |--> 11
    sage: mFp1 = matrix({k: f(v) for k, v in m1.dict().items()}); mFp1                  # optional - sage.modules sage.rings.finite_rings
    [5 0 0]
    [0 3 0]
    [0 0 3]
    sage: mFp1.base_ring()                                                              # optional - sage.modules sage.rings.finite_rings
    Finite Field of size 17
    sage: mFp2 = matrix({k: f(v) for k, v in m2.dict().items()}); mFp2                  # optional - sage.modules sage.rings.finite_rings
    [ 2  3  0]
    [ 9  8  0]
    [ 0  0 16]
    sage: mFp3 = matrix({k: f(v) for k, v in m3.dict().items()}); mFp3                  # optional - sage.modules sage.rings.finite_rings
    [16  0  0]
    [ 0  4  5]
    [ 0  7  6]

We obtain specializations in characteristic zero::

    sage: fQ = L.hom((3,5,7,11), codomain=QQ); fQ                                       # optional - sage.rings.finite_rings
    Ring morphism:
      From: Multivariate Polynomial Ring in u0, u1, u2, q over Integer Ring
            localized at (q, q + 1, u2, u1, u1 - u2, u0, u0 - u2, u0 - u1,
            u2*q - u1, u2*q - u0, u1*q - u2, u1*q - u0, u0*q - u2, u0*q - u1)
      To:   Rational Field
      Defn: u0 |--> 3
            u1 |--> 5
            u2 |--> 7
            q |--> 11
    sage: mQ1 = matrix({k: fQ(v) for k, v in m1.dict().items()}); mQ1                   # optional - sage.modules sage.rings.finite_rings
    [5 0 0]
    [0 3 0]
    [0 0 3]
    sage: mQ1.base_ring()                                                               # optional - sage.modules sage.rings.finite_rings
    Rational Field
    sage: mQ2 = matrix({k: fQ(v) for k, v in m2.dict().items()}); mQ2                   # optional - sage.modules sage.rings.finite_rings
    [-15 -14   0]
    [ 26  25   0]
    [  0   0  -1]
    sage: mQ3 = matrix({k: fQ(v) for k, v in m3.dict().items()}); mQ3                   # optional - sage.modules sage.rings.finite_rings
    [    -1      0      0]
    [     0 -15/26  11/26]
    [     0 301/26 275/26]

    sage: S.<x, y, z, t> = QQ[]
    sage: T = S.quo(x + y + z)
    sage: F = T.fraction_field()                                                        # optional - sage.libs.singular
    sage: fF = L.hom((x, y, z, t), codomain=F); fF                                      # optional - sage.libs.singular
    Ring morphism:
      From: Multivariate Polynomial Ring in u0, u1, u2, q over Integer Ring
            localized at (q, q + 1, u2, u1, u1 - u2, u0, u0 - u2, u0 - u1,
            u2*q - u1, u2*q - u0, u1*q - u2, u1*q - u0, u0*q - u2, u0*q - u1)
      To:   Fraction Field of Quotient of Multivariate Polynomial Ring in x, y, z, t
            over Rational Field by the ideal (x + y + z)
      Defn: u0 |--> -ybar - zbar
            u1 |--> ybar
            u2 |--> zbar
            q |--> tbar
    sage: mF1 = matrix({k: fF(v) for k, v in m1.dict().items()}); mF1                   # optional - sage.libs.singular sage.modules
    [        ybar            0            0]
    [           0 -ybar - zbar            0]
    [           0            0 -ybar - zbar]
    sage: mF1.base_ring() == F                                                          # optional - sage.libs.singular sage.modules
    True

.. RUBRIC:: Localizations at zero divisors and simplifications

Localization at zero divisors (which is sometimes required in algebraic
geometry) is also supported.
Here is a first example over `\ZZ`::

    sage: A = Integers(12)
    sage: L = A.localization(2)
    sage: L
    Ring of integers modulo 12 localized at (2)

Note that in `L`, we can carry out the computation
`3 = \frac{12}{2^2} = \frac{0}{2^2} = 0`.
In SageMath, this simplification is not done by default::

    sage: a = L(3)
    sage: a
    3

However, equality test works well::

    sage: a == 0
    True

It is also possible to simplify the ring `L` using the methods
:meth:`flatten` and :meth:`flattening_morphism`::

    sage: L.flatten()
    Ring of integers modulo 3
    sage: L.flattening_morphism()
    Ring morphism:
      From: Ring of integers modulo 12 localized at (2)
      To:   Ring of integers modulo 3
      Defn: 1 |--> 1

Here is a second example over polynomial rings mixing localizations
and quotients::

    sage: A.<x,y> = QQ[]
    sage: B = A.quotient([x*y^2 - x])
    sage: C = B.localization([x, y-1])
    sage: D = C.quotient([x - 2*y + 1])

Following the interpretation of algebraic geometry:

- ``A`` is the ring of functions on the affine plane,

- ``B`` is the ring of functions on the union of the three lines
  `x = 0`, `y = 1`, `y = -1`,

- ``C`` is obtained from ``B`` by taking the open subset defined
  by `x \neq 0` and `y - 1 \neq 0`; hence it is the ring of functions
  on the line `y = -1`,

- ``D`` corresponds to the intersection of the line `y = -1` with the
  line `x - 2y + 1`, which consists of the unique point `(-3, -1)`;
  consequently `D` is the ring of functions over this point, it is
  then isomorphic to the base field `\QQ` itself.

In this case, the methods :meth:`flatten` and :meth:`flattening_morphism`
are able to do this simplification::

    sage: D.flatten()
    Rational Field
    sage: D.flattening_morphism()
    Ring morphism:
      From: Quotient of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x*y^2 - x) localized at (ybar - 1, xbar) by the ideal (xbar - 2*ybar + 1)
      To:   Rational Field
      Defn: xbarbar |--> -3
            ybarbar |--> -1

TESTS::

    sage: TestSuite(L).run()                                                            # optional - sage.libs.singular sage.modules

AUTHORS:

- Sebastian Oehms (2019-12-09): initial version.
- Sebastian Oehms (2022-03-05): fix some corner cases and add :meth:`factor` (:trac:`33463`)
- Xavier Caruso, Thibaut Verron (2023-03): complete rewrite; allow localization at zero divisors
"""

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import Element
from sage.structure.parent_gens import localvars
from sage.structure.richcmp import richcmp, rich_to_bool
from sage.structure.richcmp import op_EQ, op_NE

from sage.categories.homset import Hom
from sage.categories.integral_domains import IntegralDomains
from sage.categories.rings import Rings

from sage.misc.cachefunc import cached_method
from sage.misc.latex import latex

from sage.rings.integer_ring import ZZ
from sage.rings.homset import RingHomset_generic
from sage.rings.ideal import Ideal_generic
from sage.rings.morphism import RingHomomorphism
from sage.rings.ring import CommutativeRing, CommutativeRing


# Elements
##########

def normalize_extra_units(base_ring, add_units, warning=True):
    """
    Function to normalize input data.

    The given list will be replaced by a list of the involved prime factors
    (if possible).

    INPUT:

    - ``base_ring`` -- an instance of :class:`IntegralDomain`
    - ``add_units`` -- list of elements from base ring
    - ``warning`` -- (optional, default: True) to suppress a warning which is thrown if no normalization was possible

    OUTPUT:

    List of all prime factors of the elements of the given list.

    EXAMPLES::

        sage: from sage.rings.localization import normalize_extra_units
        sage: normalize_extra_units(ZZ, [3, -15, 45, 9, 2, 50])
        [2, 3, 5]
        sage: P.<x,y,z> = ZZ[]
        sage: normalize_extra_units(P, [3*x, z*y**2, 2*z, 18*(x*y*z)**2, x*z, 6*x*z, 5])            # optional - sage.libs.pari
        [2, 3, 5, z, y, x]
        sage: P.<x,y,z> = QQ[]
        sage: normalize_extra_units(P, [3*x, z*y**2, 2*z, 18*(x*y*z)**2, x*z, 6*x*z, 5])            # optional - sage.libs.pari
        [z, y, x]

        sage: R.<x, y> = ZZ[]
        sage: Q.<a, b> = R.quo(x**2 - 5)                                                # optional - sage.libs.singular
        sage: p = b**2 - 5                                                              # optional - sage.libs.singular
        sage: p == (b-a)*(b+a)                                                          # optional - sage.libs.singular
        True
        sage: normalize_extra_units(Q, [p])                                             # optional - sage.libs.pari
        doctest:...: UserWarning: Localization may not be represented uniquely
        [b^2 - 5]
        sage: normalize_extra_units(Q, [p], warning=False)                              # optional - sage.libs.pari
        [b^2 - 5]
    """
    # convert to base ring
    add_units = [base_ring(n) for n in add_units]

    # split down to prime factors if possible
    add_units_result = []
    for n in add_units:
        try:
            if n.is_unit():
                continue
            F = list(n.factor())
            add_units_result += [f[0] for f in F]
        except (NotImplementedError, ArithmeticError, AttributeError):
            # if :meth:`is_unit` or :meth:`factor` are not available we can't do any more.
            if warning:
                from warnings import warn
                warn('Localization may not be represented uniquely')
            add_units_result = add_units
            break

    try:
        add_units_result = list(set(add_units_result))  # remove duplicates
        add_units_result = sorted(add_units_result)     # sort
    except (NotImplementedError, ValueError, TypeError):
        pass
    return add_units_result


class LocalizationElement(Element):
    def __init__(self, parent, *x, simplify=True, check=True):
        r"""
        Initialize this element.

        EXAMPLES:

        Examples with simplification::

            sage: L = ZZ.localization([2,3])
            sage: L(4)
            4
            sage: L(4,2)
            2
            sage: L(1,2^2) + L(3,2^2)
            1

            sage: P.<x,y> = QQ[]
            sage: PP.<xb,yb> = P.quotient(x^2*y)
            sage: L = PP.localization([xb])
            sage: L(xb,xb)
            xb/xb
            sage: L(0,xb^18)
            0

            sage: from sage.rings.localization import LocalizationElement
            sage: P.<x,y,z> = GF(5)[]                                                       # optional - sage.rings.finite_rings
            sage: L = P.localization((x, y*z-x))                                            # optional - sage.rings.finite_rings
            sage: LocalizationElement(L, 4/(y*z-x)**2)                                      # optional - sage.rings.finite_rings
            (-1)/(y^2*z^2 - 2*x*y*z + x^2)
            sage: _.parent()                                                                # optional - sage.rings.finite_rings
            Multivariate Polynomial Ring in x, y, z over Finite Field of size 5 localized at (x, y*z - x)

        """
        Element.__init__(self, parent)
        R = parent.numerator_ring()
        if len(x) == 1:
            x = x[0]
            if x.parent() is parent:
                num = x._num
                denom = x._denom
            else:
                try:
                    num = R(x)
                    denom = R.one()
                    check = False
                except (TypeError, ZeroDivisionError, ArithmeticError):
                    try:
                        num = R(x.numerator())
                        denom = R(x.denominator())
                    except (TypeError, AttributeError):
                        raise ValueError("unable to map %s to %s" % (x, parent))
        elif len(x) == 2:
            num = R(x[0])
            denom = R(x[1])
        else:
            raise TypeError
        if R is ZZ and denom < 0:
            num = -num
            denom = -denom
        if check and not parent.ideal(denom).is_one():
            raise ArithmeticError("denominator is not a unit")

        if num.is_zero():
            denom = R.one()
        elif simplify:
            if parent._has_gcd:
                g = num.gcd(denom)
                num = num // g
                denom = denom // g
            elif parent._has_floordiv:
                try:
                    for u in parent._extra_units:
                        while u.divides(num) and u.divides(denom):
                            num = num // u
                            denom = denom // u
                except NotImplementedError:
                    parent._has_floordiv = False
        self._num = num
        self._denom = denom

    def _repr_(self):
        r"""
        Return a string representation of this element.

        EXAMPLES::

            sage: P.<x,y> = QQ[]
            sage: L = P.localization([x+y, y])
            sage: L(1)
            1

            sage: u = L(x+y).inverse_of_unit()
            sage: u
            1/(x + y)

            sage: v = L(1/y)
            sage: v
            1/y

            sage: (x-y) * u * v^2
            (x - y)/(x*y^2 + y^3)
            sage: x * u * v^2
            x/(x*y^2 + y^3)

        TESTS::

            sage: P.<x,y> = QQ[]
            sage: L.<xl,yl> = P.localization([x+y, y])
            sage: L(x,1)
            xl
            sage: L(x)
            xl
            sage: L(x+y).inverse_of_unit()
            1/(xl + yl)
            sage: x
            x
        """
        facts = [ ]
        R = self.parent()
        S = R.numerator_ring()
        with localvars(S, R.variable_names()):
            if self._denom == 1:
                s = str(self._num)
            else:
                if self._num._is_atomic():
                    s = str(self._num)
                else:
                    s = f"({self._num})"
                s += "/"
                denom_string = str(self._denom)
                if self._denom._is_atomic() and not ('*' in denom_string or '/' in denom_string):
                    s += denom_string
                else:
                    s += f"({denom_string})"
        return s

    def _latex_(self):
        r"""
        Return a LaTeX representation of this element.

        EXAMPLES::

            sage: P.<x, y> = QQ[]
            sage: L = P.localization([x+y, y])
            sage: u = L(x+y).inverse_of_unit()
            sage: latex(u)
            \frac{1}{x + y}

            sage: v = L(y).inverse_of_unit()
            sage: latex(v)
            \frac{1}{y}

            sage: latex((x+y)*u^2*v)
            \frac{1}{x y + y^{2}}
            sage: latex(y*u^2*v)
            \frac{1}{x^{2} + 2 x y + y^{2}}

        TESTS::

            sage: P.<x, y> = QQ[]
            sage: L.<xl, yl> = P.localization([x+y,y])
            sage: latex(L(x+y).inverse_of_unit())
            \frac{1}{\mathit{xl} + \mathit{yl}}
            sage: latex(x)
            x
        """
        R = self.parent()
        S = R.numerator_ring()
        with localvars(S,R.variable_names()):
            if self._denom == 1:
                return latex(self._num)
            else:
                return "\\frac{%s}{%s}" % (latex(self._num), latex(self._denom))

    def __hash__(self):
        r"""
        Return a hash of this element.

        Since deciding equality is not easy, the hash is
        always `0`.

        EXAMPLES::

            sage: P.<x,y> = QQ[]
            sage: PP.<xb,yb> = P.quotient(x^2*y)
            sage: L = PP.localization([xb])
            sage: elt = L(xb*yb,xb^4)
            sage: hash(elt) # random
            0

        """
        return 0 # hash((self._num, self._denom))

    def numerator(self):
        r"""
        Return the numerator of this element.

        Note that simplifications can occur but we always have
        the guarantee that ``x = x.numerator() / x.denominator()``.

        EXAMPLES::

            sage: R = ZZ.localization([2,3])
            sage: R(4).numerator()
            4
            sage: R(4, 2^2*3).numerator()
            1

        .. SEEALSO::

            :meth:`denominator`
        """
        return self._num

    def denominator(self):
        r"""
        Return the denominator of this element.

        Note that simplifications can occur but we always have
        the guarantee that ``x = x.numerator() / x.denominator()``.

        EXAMPLES::

            sage: R = ZZ.localization([2,3])
            sage: R(4).denominator()
            1
            sage: R(4, 2^2*3).denominator()
            3

        .. SEEALSO::

            :meth:`numerator`
        """
        return self._denom

    def _add_(self, other):
        r"""
        Return the sum of this element and ``other``.

        EXAMPLES::

            sage: L = ZZ.localization((5, 11))  # indirect doctest
            sage: L(1/5) + L(2/11)
            21/55
        """
        num = self._num * other._denom + self._denom * other._num
        denom = self._denom * other._denom
        return self.__class__(self.parent(), num, denom, check=False)

    def _sub_(self, other):
        r"""
        Return the subtraction of this element and ``other``.

        EXAMPLES::

            sage: L = ZZ.localization((5, 11))
            sage: L(2) - L(1/121)  # indirect doctest
            241/121
        """
        num = self._num * other._denom - self._denom * other._num
        denom = self._denom * other._denom
        return self.__class__(self.parent(), num, denom, check=False)

    def _neg_(self):
        r"""
        Return the opposite of this element.

        EXAMPLES::

            sage: L = ZZ.localization((5, 11))
            sage: -L(1/121)  # indirect doctest
            -1/121
        """
        return self.__class__(self.parent(), -self._num, self._denom, check=False)

    def _mul_(self,other):
        r"""
        Return the product of this element and ``other``.

        EXAMPLES::

            sage: L = ZZ.localization((5, 11))
            sage: L(1/55) * L(2/25)  # indirect doctest
            2/1375
        """
        num = self._num * other._num
        denom = self._denom * other._denom
        return self.__class__(self.parent(), num, denom, check=False)

    def _lmul_(self,a):
        r"""
        Return the product of this element and ``other``.

        EXAMPLES::

            sage: L = ZZ.localization((5, 11))
            sage: 7 * L(3/5)  # indirect doctest
            21/5
        """
        num = self._num * a
        return self.__class__(self.parent(), num, self._denom, check=False)

    def _rmul_(self,a):
        r"""
        Return the product of this element and ``other``.

        EXAMPLES::

            sage: L = ZZ.localization((5, 11))
            sage: L(3/5) * 7  # indirect doctest
            21/5
            sage: P.<X, Y> = QQ['x, y']
            sage: L = P.localization(X - Y)
            sage: x, y = L.gens()
            sage: p = (x^2 - y^2)/(x-y)^2                                               # optional - sage.libs.singular
            sage: p.factor()                                                            # optional - sage.libs.singular
            (1/(x - y)) * (x + y)
        """
        num = a * self._num
        return self.__class__(self.parent(), num, self._denom, check=False)

    def is_unit(self):
        """
        Return ``True`` if this element is a unit.

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: L = P.localization((x, y*z))                                          # optional - sage.libs.pari
            sage: L(y*z).is_unit()                                                      # optional - sage.libs.pari
            True
            sage: L(z).is_unit()                                                        # optional - sage.libs.pari
            True
            sage: L(x*y*z).is_unit()                                                    # optional - sage.libs.pari
            True

        Other examples where we localize at zero divisors::

            sage: P.<x,y> = QQ[]
            sage: PP.<xx,yy> = P.quotient(x*y)
            sage: L = PP.localization([xx])
            sage: L(x).is_unit()
            True
            sage: L(y).is_unit()
            False

            sage: PP.<xx,yy> = P.quotient(x^2)
            sage: L = PP.localization([xx])
            sage: L(x).is_unit()
            True
            sage: L(y).is_unit()
            True
        """
        return self.parent().ideal(self).is_one()

    def inverse_of_unit(self):
        r"""
        Return the inverse of this element.

        EXAMPLES::

            sage: P.<x,y,z> = ZZ[]
            sage: L = Localization(P, x*y*z)
            sage: L(x*y*z).inverse_of_unit()                                            # optional - sage.libs.singular
            1/(x*y*z)
            sage: L(z).inverse_of_unit()                                                # optional - sage.libs.singular
            1/z

        ::

            sage: L = ZZ.localization([2,3])
            sage: L(1).inverse_of_unit()
            1
            sage: L(-1).inverse_of_unit()
            -1
            sage: L(2).inverse_of_unit()
            1/2
            sage: L(-3).inverse_of_unit()
            -1/3

            sage: L(6).inverse_of_unit()
            1/6
            sage: L(1/6).inverse_of_unit()
            6
            sage: L(5).inverse_of_unit()
            Traceback (most recent call last):
            ...
            ArithmeticError: 5 is not invertible
        """
        if not self.is_unit():
            raise ArithmeticError("%s is not invertible" % self)
        return self.__class__(self.parent(), self._denom, self._num, check=False)

    def is_nilpotent(self):
        r"""
        Return ``True`` if this element is nilpotent.

        EXAMPLES::

            sage: A = Integers(36)
            sage: L = A.localization(2)
            sage: L(2).is_nilpotent()
            False
            sage: L(3).is_nilpotent()
            True

        """
        return self.parent()._quotient(self._num).is_nilpotent()

    def __invert__(self):
        r"""
        Return the inverse of this element in the fraction field.

        EXAMPLES::

            sage: P.<x,y,z> = ZZ[]
            sage: L = P.localization(x*y*z)
            sage: u = 1/x  # indirect doctest
            sage: u
            1/x
            sage: u.parent()
            Fraction Field of Multivariate Polynomial Ring in x, y, z over Integer Ring
        """
        K = self.parent().fraction_field()
        return K(self._denom, self._num)

    def factor(self, proof=False):
        r"""
        Return the factorization of this polynomial.

        INPUT:

        - ``proof`` -- a boolean (default: ``False``); passed to the
          corresponding method of the numerator of ``self``

        EXAMPLES::

            sage: P.<x, y> = QQ[]
            sage: L = P.localization(x - y)
            sage: p = L((x^2 - y^2)/(x-y)^2)
            sage: p.factor()
            (1/(x - y)) * (x + y)
        """
        num = self._num
        den = self._denom
        F = self._num.factor(proof=proof)
        P = self.parent()
        fac = [(P(f), e) for (f, e) in F]
        unit = P(self._denom).inverse_of_unit() * F.unit()
        from sage.structure.factorization import Factorization
        return Factorization(fac, unit)

    def _richcmp_(self, other, op):
        """
        EXAMPLES::

            sage: P.<x,y,z> = GF(7)[]
            sage: L = P.localization((x, y, z))
            sage: ~L(y*z/x) == L(x/(y*z))
            True

        ::

            sage: A = Integers(40)
            sage: L = A.localization(2)
            sage: L(5/2) == 0
            True
        """
        if op != op_EQ and op != op_NE:
            return NotImplemented
        if not self.parent()._has_equality_test:
            return NotImplemented
        a = self._num * other._denom
        b = self._denom * other._num
        if (a - b) in self.parent()._ideal:
            return rich_to_bool(op, 0)
        else:
            return rich_to_bool(op, -1)

    def _rational_(self):
        r"""
        Convert this element to a rational.

        This is only possible if its base ring is the ring of integers.

        TESTS::

            sage: L = ZZ.localization(5)
            sage: cp3 = cyclotomic_polynomial(3).change_ring(L)                         # optional - sage.libs.pari
            sage: cp3.splitting_field('t')      #   indirect doctest                    # optional - sage.libs.pari sage.rings.number_field
            Number Field in t with defining polynomial x^2 + x + 1
        """
        # Only if base ring is ZZ
        if not self.parent()._numring is ZZ:
            raise ValueError("not a rational")
        else:
            return self.numerator()/self.denominator()


# Parents
#########

class Localization(CommutativeRing, UniqueRepresentation):
    r"""

    The localization generalizes the construction of the field of
    fractions of an integral domain to an arbitrary ring. Given a
    ring `R` and a subset `S` of `R`, there exists a ring `R[S^{-1}]`
    together with the ring homomorphism `R \longrightarrow R[S^{-1}]`
    that "inverts" `S`; that is, the homomorphism maps elements in `S`
    to unit elements in `R[S^{-1}]` and, moreover, any ring homomorphism
    from `R` that "inverts" `S` uniquely factors through `R[S^{-1}]`.

    The ring `R[S^{-1}]` is called the *localization* of `R` with
    respect to `S`. For example, if `f` an element in `R`, then the
    localization consists of elements of the form `r/f^n, r\in R,
    n \geq 0` (to be precise, `R[f^{-1}] = R[t]/(ft-1)`.

    .. NOTE::

        The base ring should support saturation of ideals in order
        to guarantee an successful application.
        When :meth:`sage.structure.element.CommutativeRingElement.divides`
        and the exact division operator `//`
        (:meth:`sage.structure.element.Element.__floordiv__`) are
        implemented, fractions are simplified automatically.

    INPUT:

    - ``base_ring`` -- a ring

    - ``extra_units`` -- tuple of elements of ``base_ring`` which
      should be turned into units

    - ``names`` -- a list of string or ``None`` (default: ``None``),
      the names of the variables in the localization; if not given,
      the variable names of the base are reused

    - ``normalize`` -- a boolean (default: ``True``),
      passed to :class:`CommutativeRing`

    - ``category`` -- a category or ``None`` (default: ``None``),
      passed to :class:`CommutativeRing`

    - ``warning`` -- a boolean (default: ``True``),
      to suppress a warning which is thrown if self cannot be
      represented uniquely

    REFERENCES:

    - :wikipedia:`Ring_(mathematics)#Localization`

    EXAMPLES::

        sage: L = ZZ.localization((3,5))
        sage: 1/45 in L
        True
        sage: 1/43 in L
        False

        sage: L.localization((7, 11))
        Integer Ring localized at (3, 5) localized at (7, 11)
        sage: _.is_subring(QQ)
        True

        sage: L(~7)
        Traceback (most recent call last):
        ...
        ArithmeticError: denominator is not a unit

        sage: Zp(7).localization((3, 5))
        7-adic Ring with capped relative precision 20

        sage: R.<x> = ZZ[]
        sage: L = R.localization(x**2 + 1)                                              # optional - sage.libs.pari
        sage: s = (x+5)/(x**2+1)
        sage: s in L                                                                    # optional - sage.libs.pari
        True
        sage: t = (x+5)/(x**2+2)
        sage: t in L                                                                    # optional - sage.libs.pari
        False
        sage: L(t)                                                                      # optional - sage.libs.pari
        Traceback (most recent call last):
        ...
        ArithmeticError: denominator is not a unit
        sage: L(s) in R
        False
        sage: y = L(x)                                                                  # optional - sage.libs.pari
        sage: g = L(s)                                                                  # optional - sage.libs.pari
        sage: g.parent()                                                                # optional - sage.libs.pari
        Univariate Polynomial Ring in x over Integer Ring localized at (x^2 + 1,)
        sage: f = (y+5)/(y**2+1); f                                                     # optional - sage.libs.pari
        (x + 5)/(x^2 + 1)
        sage: f == g                                                                    # optional - sage.libs.pari
        True
        sage: (y+5)/(y**2+2)                                                            # optional - sage.libs.pari
        Traceback (most recent call last):
        ...
        ValueError: factor x^2 + 2 of denominator is not a unit

        sage: Lau.<u, v> = LaurentPolynomialRing(ZZ)
        sage: LauL = Lau.localization(u + 1)
        sage: LauL(~u).parent()
        Multivariate Laurent Polynomial Ring in u, v over Integer Ring localized at (u + 1)

    Localization by a divisor of 0::

        sage: P.<x,y> = QQ[]
        sage: PP.<xb,yb> = P.quotient(x^2-y^2)
        sage: L = PP.localization([xb-yb]); L
        Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2 - y^2) localized at (xb - yb)

    TESTS:

    Check that :trac:`33463` is fixed::

        sage: R = ZZ.localization(5)
        sage: R.localization(~5)
        Integer Ring localized at (5)

    """
    Element = LocalizationElement

    def __classcall_private__(cls, base, units, names=None, normalize=True, category=None):
        r"""
        TESTS::

            sage: L1 = ZZ.localization((2, 5))
            sage: L2 = ZZ.localization([2, 5])
            sage: L1 is L2
            True
        """
        from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing_generic
        if not isinstance(units, (list, tuple)):
            units = [units]
        units = tuple(base(u) for u in units)
        if isinstance(names, list):
            names = tuple(names)
        #if isinstance(base, Localization):
        #    # I really dislike this automatic simplification;
        #    # however, I implement it for backward compatibility
        #    units = base._extra_units + tuple(u.numerator() for u in units)
        #    base = base.numerator_ring()
        #if isinstance(base, LaurentPolynomialRing_generic):
        #    # Same remark
        #    units = base.gens() + units
        #    base = base.polynomial_ring()
        if base in IntegralDomains():
            if category is None:
                category = IntegralDomains()
            else:
                category = category.join([IntegralDomains()])
        return UniqueRepresentation.__classcall__(cls, base, units, names, normalize, category)

    def __init__(self, base, units, names=None, normalize=True, category=None):
        """
        Python constructor of Localization.

        TESTS::

            sage: L = Localization(ZZ, (3, 5))
            sage: TestSuite(L).run()

            sage: R.<x> = ZZ[]
            sage: L = R.localization(x**2 + 1)                                          # optional - sage.libs.pari
            sage: TestSuite(L).run()                                                    # optional - sage.libs.pari
        """
        self._numring = base
        self._given_units = units
        self._extra_units = tuple(normalize_extra_units(base, units, warning=False))

        J = base.ideal(0)
        self._is_ideal_domain = False
        self._has_equality_test = False
        if base in IntegralDomains():
            self._has_equality_test = True
            self._is_integral_domain = True
            if 0 in self._extra_units:
                J = base.ideal([1])
        else:
            try:
                for u in self._extra_units:
                    J, _ = J.saturation(u)
                self._has_equality_test = True
            except (AttributeError, NotImplementedError):
                pass

        self._ideal = J
        self._quotient = base.quotient(J)
        self._has_floordiv = hasattr(base.zero(), "_floordiv_")
        self._has_gcd = self._has_floordiv and hasattr(base.zero(), "gcd")

        if names is None:
            self._names = base.variable_names()
        else:
            self._assign_names(names)
        CommutativeRing.__init__(self, base.base_ring(), names=names, normalize=normalize, category=category)
        self.register_coercion(base)
        self._gens = tuple(self(v) for v in base.gens())
        # We create the fraction field if it's easy;
        # this sets at the same time a coercion map to the fraction field
        if self in IntegralDomains():
            try:
                _ = self.fraction_field()
            except (NotImplementedError, TypeError):
                pass

    def _element_constructor_(self, *args, **kwds):
        """
        Return an element of this parent constructed from the given arguments.

        EXAMPLES::

            sage: L = ZZ.localization((5, 2))
            sage: L(1/25)   # indirect doctest
            1/25
            sage: L(1, 20)  # indirect doctest
            1/20
            sage: L(2/10)  # indirect doctest
            1/5

        ::

            sage: L = Integers(12).localization(3)
            sage: L(1/3)
            1/3
            sage: L(1/2)
            Traceback (most recent call last):
            ...
            ArithmeticError: denominator is not a unit
        """
        kwds['check'] = True
        return self.element_class(self, *args, **kwds)

    def _repr_(self):
        r"""
        Return a string representation of this ring.

        EXAMPLES::

            sage: P.<x,y> = QQ[]
            sage: PP.<xb,yb> = P.quotient(x^2-y^2)
            sage: L = PP.localization([xb-yb]); L
            Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2 - y^2) localized at (xb - yb)
            sage: R.<a> = GF(3)[]                                                       # optional - sage.rings.finite_rings
            sage: Localization(R, a**2 - 1)                                             # optional - sage.rings.finite_rings
            Univariate Polynomial Ring in a over Finite Field of size 3
             localized at (a + 1, a + 2)
        """
        if len(self._extra_units) == 1:
            return f"{self._numring} localized at ({self._extra_units[0]})"
        elif len(self._extra_units) > 1:
            return f"{self._numring} localized at {self._extra_units}"
        else:
            return str(self._numring)

    def _latex_(self):
        r"""
        Return a LaTeX representation of this ring.

        EXAMPLES::

            sage: P.<x,y> = QQ[]
            sage: PP.<xb,yb> = P.quotient(x^2-y^2)
            sage: L = PP.localization([xb-yb])
            sage: latex(L)
            \Bold{Q}[x, y]/\left(x^{2} - y^{2}\right)\Bold{Q}[x, y] \left[ \frac{1}{\mathit{xb} - \mathit{yb}}\right]
        """
        s = latex(self._numring) + r"\left["
        invs = [f"\\frac{{1}}{{{latex(u)}}}" for u in self._extra_units]
        s += ", ".join(invs) + r"\right]"
        return s

    def ngens(self):
        """
        Return the number of generators of this ring,
        according to the same method for the base ring.

        EXAMPLES::

            sage: R.<x, y> = ZZ[]
            sage: R.localization((x^2+1, y-1)).ngens()
            2
            sage: ZZ.localization(2).ngens()
            1
        """
        return self._numring.ngens()

    def gens(self):
        """
        Return a tuple whose entries are the generators for this
        object, in order.

        EXAMPLES::

            sage: R.<x, y> = ZZ[]
            sage: Localization(R, (x**2 + 1, y - 1)).gens()                             # optional - sage.libs.pari
            (x, y)
            sage: ZZ.localization(2).gens()
            (1,)
        """
        return self._gens

    def gen(self,i):
        """
        Return the ``i``-th generator of this ring,
        which is the image in this ring of the ``i``-th generator
        of the base ring.

        EXAMPLES::

            sage: R.<x, y> = ZZ[]
            sage: R.localization((x^2+1, y-1)).gen(0)
            x
            sage: ZZ.localization(2).gen(0)
            1
        """
        return self._gens[i]

    def _ideal_class_(self, num_gens):
        r"""
        Use a specialized class for ideals in localizations.

        EXAMPLES::

            sage: R.<x, y> = ZZ[]
            sage: L = R.localization(x + 1)
            sage: I = L.ideal(x)
            sage: I
            Ideal (x) of Multivariate Polynomial Ring in x, y over Integer Ring localized at (x + 1)
            sage: type(I)  # indirect doctest
            <class 'sage.rings.localization.LocalizationIdeal_generic'>
        """
        return LocalizationIdeal_generic

    def units(self, normalize=False):
        r"""
        Return the elements which were inverted to create
        this localizated ring.

        INPUT:

        - ``normalize`` -- a boolean (default: ``False``); if ``False``,
          return the elements as provided by the user; if ``True``,
          normalize the elements (by factoring and sorting them if
          possible)

        EXAMPLES::

            sage: R.<x, y> = ZZ[]
            sage: L = R.localization(x*y)
            sage: L.units()
            (x*y,)
            sage: L.units(normalize=True)
            (y, x)

        Note that the output consists of elements of the base ring::

            sage: u = L.units()[0]
            sage: u.parent()
            Multivariate Polynomial Ring in x, y over Integer Ring

        .. SEEALSO::

            :meth:`inverse_of_units`
        """
        if normalize:
            return self._extra_units
        else:
            return self._given_units

    def inverse_of_units(self, normalize=False):
        r"""
        Return the inverses of the elements which were inverted
        to create this localizated ring.

        INPUT:

        - ``normalize`` -- a boolean (default: ``False``); if ``False``,
          return the inverses of the elements provided by the user;
          if ``True``, normalize first the elements (by factoring and
          sorting them if possible)

        EXAMPLES::

            sage: R.<x, y> = ZZ[]
            sage: L = R.localization(x*y)
            sage: L.inverse_of_units()
            (1/(x*y),)
            sage: L.inverse_of_units(normalize=True)
            (1/y, 1/x)
        """
        if normalize:
            return tuple(self(u).inverse_of_unit() for u in self._extra_units)
        else:
            return tuple(self(u).inverse_of_unit() for u in self._given_units)

    def numerator_kernel(self):
        r"""
        If this ring is defined as `R[S^{-1}]`, return the kernel of the
        defining morphism `R \to R[S^{-1}]`.

        Note that this kernel is nonzero if and only if `S` contains a
        zero divisor.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: L = R.localization(x+1)
            sage: L.numerator_kernel()
            Principal ideal (0) of Univariate Polynomial Ring in x over Rational Field

        ::

            sage: A = Integers(100)
            sage: L = A.localization(5)
            sage: L.numerator_kernel()
            Principal ideal (4) of Ring of integers modulo 100
        """
        return self._ideal

    def numerator_ring(self):
        r"""
        If this ring is defined as `R[S^{-1}]`, return `R`.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: L = R.localization(x+1)
            sage: L
            Univariate Polynomial Ring in x over Rational Field localized at (x + 1)
            sage: L.numerator_ring()
            Univariate Polynomial Ring in x over Rational Field
        """
        return self._numring

    def is_integral_domain(self):
        r"""
        Return ``True`` if this ring is an integral domain.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: L = R.localization(x+1)
            sage: L.is_integral_domain()
            True

        ::

            sage: A = Integers(12)
            sage: L2 = A.localization(2)  # isomorphic to Z/3Z
            sage: L2.is_integral_domain()
            True

            sage: L3 = A.localization(3)  # isomorphic to Z/4Z
            sage: L3.is_integral_domain()
            False
        """
        if self in IntegralDomains():
            return True
        # maybe we should have an option to avoid the next call
        # which can be costly
        return self._ideal.is_prime()

    @cached_method
    def fraction_field(self):
        """
        Return the fraction field of this ring.

        EXAMPLES::

            sage: R.<a> = GF(5)[]
            sage: L = R.localization((a^2-3, a))
            sage: K = L.fraction_field()
            sage: K
            Fraction Field of Univariate Polynomial Ring in a over Finite Field of size 5

            sage: L.is_subring(K)
            True
        """
        try:
            if not self.is_integral_domain():
                raise TypeError("%s is not an integral domain" % self)
        except NotImplementedError:
            pass
        field = self._quotient.fraction_field()
        field.register_coercion(self)
        return field

    def is_field(self, proof=False):
        """
        Return ``True`` if this ring is a field.

        EXAMPLES::

            sage: R.<x, y> = QQ[]
            sage: L = R.localization(x*y)
            sage: L.is_field()
            False

        .. warning::

            This function has known bugs. For example::

                sage: A = Zp(7)
                sage: L = A.localization(7)  # L is Q7
                sage: L.is_field()
                False

        """
        if self.numerator_ring().is_field(proof=proof):
            return True
        if proof is False:
            # It's never a field if the base ring is not semi-local
            return False
        else:
            raise NotImplementedError

    def is_finite(self):
        r"""
        Return ``True`` if this ring is finite.

        EXAMPLES::

            sage: L = ZZ.localization(2)
            sage: L.is_finite()
            False

            sage: A = Integers(12)
            sage: L = A.localization(2)
            sage: L.is_finite()
            True

        """
        return self._quotient.is_finite()

    def characteristic(self):
        """
        Return the characteristic of this ring.

        EXAMPLES::

            sage: R.<a> = GF(5)[]                                                       # optional - sage.rings.finite_rings
            sage: L = R.localization((a**2 - 3, a))                                     # optional - sage.rings.finite_rings
            sage: L.characteristic()                                                    # optional - sage.rings.finite_rings
            5
        """
        return self.base_ring().characteristic()

    def _flattening_function(self):
        r"""
        Return a function from this ring to its flattened version.

        This is an helper function; do not call it directly.

        TESTS::

            sage: A.<x> = QQ[]
            sage: B = A.localization(0)
            sage: B.flatten()  # indirect doctest
            The zero ring

        """
        from sage.rings.quotient_ring import quotient_with_simplification
        R = self._numring
        try:
            Rs, f = R._flattening_function()   # f : R -> Rs
        except (AttributeError, NotImplementedError):
            Rs = R
            f = lambda x: x
        I = Rs.ideal([ f(x) for x in self._ideal.gens() ])
        RsI, to_RsI = quotient_with_simplification(Rs, I)
        try:
            S, to_S = RsI._flattening_function()
            g = lambda x: to_S(to_RsI(f(x)))   # g : R -> S
        except (AttributeError, NotImplementedError):
            S = RsI
            g = lambda x: to_RsI(f(x))         # g : R -> S
        if isinstance(S, Localization):
            base = S._numring
            units = S._extra_units + tuple(g(R(u)).numerator() for u in self._extra_units)
            ring = base.localization(units)

            def isom(x):  # self -> ring
                num = g(x.numerator())      # in S
                denom = g(x.denominator())  # in S
                num2 = num.numerator() * denom.denominator()
                denom2 = num.denominator() * denom.numerator()
                return ring(num2) * ring(denom2).inverse_of_unit()

        else:
            base = S
            units = [ g(u) for u in self._extra_units ]
            ring = base.localization(units)

            def isom(x):  # self -> ring
                num = g(x.numerator())      # in S
                denom = g(x.denominator())  # in S
                return ring(num) * ring(denom).inverse_of_unit()

        return ring, isom

    def flattening_morphism(self):
        r"""
        Return the flattening morphism from this ring to its
        flattened version (see :meth:`flatten`).

        EXAMPLES::

            sage: L1 = ZZ.localization(2)
            sage: L2 = L1.localization(3)
            sage: L2.flattening_morphism()
            Ring morphism:
              From: Integer Ring localized at (2) localized at (3)
              To:   Integer Ring localized at (2, 3)
              Defn: 1 |--> 1

        ::

            sage: A = Integers(12)
            sage: B = A.localization(2)
            sage: B.flattening_morphism()
            Ring morphism:
              From: Ring of integers modulo 12 localized at (2)
              To:   Ring of integers modulo 3
              Defn: 1 |--> 1

        ::

            sage: A.<x, y> = QQ[]
            sage: B = A.localization(x)
            sage: C = B.quotient(x^2*y)
            sage: C.flattening_morphism()
            Ring morphism:
              From: Quotient of Multivariate Polynomial Ring in x, y over Rational Field localized at (x) by the ideal (x^2*y)
              To:   Univariate Polynomial Ring in x over Rational Field localized at (x)
              Defn: xbar |--> x
                    ybar |--> 0

        .. SEEALSO::

            :meth:`flatten`

        """
        ring, isom = self._flattening_function()
        im_gens = [ isom(x) for x in self.gens() ]
        return self._Hom_(ring, category=Rings())(im_gens)

    def flatten(self):
        r"""
        Return the flattened (and possibly simplified) version of this ring.

        EXAMPLES:

        Localizing twice is equivalent to localizing just once::

            sage: L1 = ZZ.localization(2)
            sage: L2 = L1.localization(3)
            sage: L2
            Integer Ring localized at (2) localized at (3)
            sage: L2.flatten()
            Integer Ring localized at (2, 3)

        In the flattened version, quotients always appear before
        localizations::

            sage: A.<x,y> = QQ[]
            sage: L = A.localization(x)
            sage: B = L.quotient(x + y + 1)
            sage: B
            Quotient of Multivariate Polynomial Ring in x, y over Rational Field localized at (x) by the ideal (x + y + 1)
            sage: B.flatten()
            Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x + y + 1) localized at (-ybar - 1)

        A series of quotients are localizations are reduced to a unique
        quotient followed by a localization::

            sage: A.<x, y, z> = QQ[]
            sage: B = A.localization(x)
            sage: C = B.quotient(x^2 + y^2 + z^2 - 1)
            sage: D = C.localization(z - 1)
            sage: E = D.quotient(x + y + 1)
            sage: E
            Quotient of Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field localized at (x) by the ideal (x^2 + y^2 + z^2 - 1) localized at (zbar - 1) by the ideal (xbar + ybar + 1)
            sage: E.flatten()
            Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by the ideal (x + y + 1, 2*y^2 + z^2 + 2*y, x^2 + y^2 + z^2 - 1) localized at (zbar - 1, -ybar - 1)

        Simplifications are made when possible::

            sage: A.<x, y> = QQ[]
            sage: B = A.localization(x)
            sage: C = B.quotient(x^2*y)
            sage: C
            Quotient of Multivariate Polynomial Ring in x, y over Rational Field localized at (x) by the ideal (x^2*y)
            sage: C.flatten()
            Univariate Polynomial Ring in x over Rational Field localized at (x)

        ::

            sage: A = Integers(12)
            sage: B = A.localization(2)
            sage: B
            Ring of integers modulo 12 localized at (2)
            sage: B.flatten()
            Ring of integers modulo 3

            sage: C = B.localization(3)
            sage: C
            Ring of integers modulo 12 localized at (2) localized at (3)
            sage: C.flatten()
            The zero ring

        .. SEEALSO::

            :meth:`flattening_morphism`, :meth:`sage.rings.quotient_ring.QuotientRing_nc.flatten`

        """
        ring, _ = self._flattening_function()
        return ring

    def _Hom_(self, other, category):
        r"""
        Return the homset from this ring to ``other`` in the
        category ``category``.

        EXAMPLES::

            sage: A = ZZ.localization(2)
            sage: Hom(A, QQ)
            Set of Homomorphisms from Integer Ring localized at (2) to Rational Field

        """
        return RingHomset_localized(self, other, category=category)


class LocalizationIdeal_generic(Ideal_generic):
    r"""
    Generic class for ideal in localized rings.

    See :func:`sage.rings.ideal.Ideal()`.
    """
    @cached_method
    def numerator_ideal(self):
        r"""
        If the underlying ring is `R[S^{-1}]`, return the
        intersection of the ideal with `R`; it is an ideal of `R`.

        EXAMPLES::

            sage: P.<x,y> = QQ[]
            sage: PP.<xb, yb> = P.quotient(x^2*y)
            sage: L = PP.localization([x])
            sage: I = L.ideal(y/x^2)
            sage: I
            Ideal (yb/xb^2) of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2*y) localized at (xb)
            sage: I.numerator_ideal()
            Ideal (yb) of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2*y)
        """
        R = self.ring()
        S = R.numerator_ring()
        gens_S = [g.numerator() for g in self.gens()]
        I_S = S.ideal(gens_S)
        for u in R._extra_units:
            I_S, _ = I_S.saturation(u)
        return I_S

    def _contains_(self, elt):
        r"""
        Return ``True`` is ``elt`` belongs to this ideal.

        EXAMPLES::

            sage: P.<x,y> = QQ[]
            sage: PP.<xb, yb> = P.quotient(x^2*y)
            sage: L = PP.localization([x])
            sage: I = L.ideal(y/x^2)
            sage: y in I
            True
            sage: 1 in I
            False

        """
        return elt.numerator() in self.numerator_ideal()

    def is_zero(self):
        r"""
        Return ``True`` if this ideal is zero.

        EXAMPLES::

            sage: A = Integers(12)
            sage: L = A.localization(2)
            sage: I = L.ideal(2)
            sage: I.is_zero()
            False

            sage: J = L.ideal(3)
            sage: J.is_zero()
            True

        """
        return all([ x == 0 for x in self.gens() ])

    def is_one(self):
        r"""
        Return ``True`` if this ideal is the full ideal,
        that is if it contains `1`.

        EXAMPLES::

            sage: A = Integers(12)
            sage: L = A.localization(2)
            sage: I = L.ideal(2)
            sage: I.is_one()
            True

            sage: J = L.ideal(3)
            sage: J.is_one()
            False
        """
        return self.numerator_ideal().is_one()

    def saturation(self, other):
        r"""
        Return the saturation and an upper bound on the saturation exponent
        of this ideal with respect to the ideal ``other``.

        INPUT:

        - ``other`` -- another ideal or an element of the base ring
          (in which case, the saturating ideal is the principal ideal
          generated by the given element)

        OUTPUT:

        - a pair (ideal, integer)

        EXAMPLES::

            sage: R.<x, y> = QQ[]
            sage: L = R.localization(x)
            sage: I = L.ideal(x*y^3)
            sage: I.saturation(y)
            (Ideal (1) of Multivariate Polynomial Ring in x, y over Rational Field localized at (x),
             3)
        """
        R = self.ring()
        S = R.numerator_ring()
        I = self.numerator_ideal()
        if isinstance(other, Ideal_generic):
            if other.ring() is not R :
                raise TypeError
            K = other.numerator_ideal()
        else:
            if not isinstance(other, (list,tuple)):
                other = [other]
            other = [R(o).numerator() for o in other]
            K = S.ideal(other)
        J, n = I.saturation(K)
        return R.ideal([R(g) for g in J.gens()]), n

    def is_prime(self):
        r"""
        Return ``True`` is this ideal is a prime ideal.

        EXAMPLES::

            sage: R.<x, y> = QQ[]
            sage: L = R.localization(x*y)
            sage: I = L.ideal(x*y + x)
            sage: I.is_prime()
            True
        """
        return self.numerator_ideal().is_prime()


# Morphisms
###########

class LocalizationMorphism(RingHomomorphism):
    r"""
    A class for ring homomorphisms whose domain is a localized ring.

    TESTS::

        sage: A.<x,y> = QQ[]
        sage: L.<xl, yl> = A.localization((x, y))
        sage: L.hom([yl, xl])
        Ring endomorphism of Multivariate Polynomial Ring in x, y over Rational Field localized at (y, x)
          Defn: xl |--> yl
                yl |--> xl

        sage: L.hom([y, x])
        Traceback (most recent call last):
        ...
        ValueError: inverted element y is not mapped to a unit
    """
    def __init__(self, parent, numerator_morphism, base_map=None, check=True):
        r"""
        Initialize this morphism.

        INPUT:

        - ``parent`` -- the homset in which this morphism lives

        - ``numerator_morphism`` -- a datum describing the morphism on the
          numerator ring

        - ``base_map`` -- a morphism or ``None`` (default: ``None``), the
          action of the morphism of the base ring of the numerator ring

        - ``check`` -- a boolean (default: ``True``), whether to check if
          the inverted elements are mapped to units in the codomain

        TESTS::

            sage: A.<x,y> = QQ[]
            sage: L = A.localization(x*y)
            sage: f = L.hom([y, x], codomain=L)  # indirect doctest
            sage: TestSuite(f).run()
        """
        domain = parent.domain()
        if not isinstance(domain, Localization):
            raise TypeError
        RingHomomorphism.__init__(self, parent)
        try:
            f = parent.numerator_homset()(numerator_morphism, base_map=base_map, check=check)
        except (ValueError, TypeError):
            # sometimes works better (e.g. the class SplittingAlgebras redefines the method hom)
            f = domain.numerator_ring().hom(numerator_morphism, codomain=parent.codomain(), base_map=base_map, check=check)
        self._numerator_morphism = f
        if check:
            for u in domain._extra_units:
                if not f(u).is_unit():
                    raise ValueError("inverted element %s is not mapped to a unit" % u)

    def numerator_morphism(self):
        r"""
        If the codomain is `R[S^{-1}]`, return the restriction of this
        morphism to `R`.

        EXAMPLES::

            sage: R.<x, y> = QQ[]
            sage: L.<xl, yl> = R.localization(x*y)
            sage: f = L.hom([yl, xl])
            sage: f
            Ring endomorphism of Multivariate Polynomial Ring in x, y over Rational Field localized at (y, x)
              Defn: xl |--> yl
                    yl |--> xl

            sage: f.numerator_morphism()
            Ring morphism:
              From: Multivariate Polynomial Ring in x, y over Rational Field
              To:   Multivariate Polynomial Ring in x, y over Rational Field localized at (y, x)
              Defn: x |--> yl
                    y |--> xl

        """
        return self._numerator_morphism

    @cached_method
    def im_gens(self):
        r"""
        Return the list of images of the generators of the domain
        under this morphism.

        EXAMPLES::

            sage: R.<x, y> = QQ[]
            sage: L.<xl, yl> = R.localization(x*y)
            sage: f = L.hom([y, x], codomain=L)
            sage: f.im_gens()
            [yl, xl]

        """
        return [ self(x) for x in self.domain().gens() ]

    def _call_(self, x):
        r"""
        Return the image of `x` under this morphism.

        EXAMPLES::

            sage: R.<x, y> = QQ[]
            sage: L.<xl, yl> = R.localization(x*y)
            sage: f = L.hom([yl, xl])
            sage: f(xl/yl)  # indirect doctest
            yl/xl

        """
        f = self._numerator_morphism
        num = x.numerator()
        denom = x.denominator()
        return f(num) * f(denom).inverse_of_unit()

    def _repr_defn(self):
        r"""
        Return a string with the definition of this morphism.

        Used in constructing string representation.

        EXAMPLES::

            sage: R.<x, y> = QQ[]
            sage: L.<xl, yl> = R.localization(x*y)
            sage: f = L.hom([yl, xl])
            sage: print(f._repr_defn())
            xl |--> yl
            yl |--> xl

        """
        R = self.domain()
        S = R.numerator_ring()
        with localvars(S, R.variable_names()):
            return self._numerator_morphism._repr_defn()


class RingHomset_localized(RingHomset_generic):
    r"""
    A class for homset when the domain is a localization.

    TESTS::

        sage: R.<x, y> = QQ[]
        sage: L.<xl, yl> = R.localization(x*y)
        sage: H = Hom(L, R)  # indirect doctest
        sage: H
        Set of Homomorphisms from Multivariate Polynomial Ring in x, y over Rational Field localized at (y, x) to Multivariate Polynomial Ring in x, y over Rational Field

    """
    Element = LocalizationMorphism

    @cached_method
    def numerator_homset(self):
        r"""
        If the domain is `R[S^{-1}]` and the codomain is `C`,
        return `\text{Hom}{R, C}`.

        TESTS::

            sage: R.<x, y> = QQ[]
            sage: L = R.localization(x*y)
            sage: H = End(L)
            sage: H
            Set of Homomorphisms from Multivariate Polynomial Ring in x, y over Rational Field localized at (y, x) to Multivariate Polynomial Ring in x, y over Rational Field localized at (y, x)

            sage: H.numerator_homset()
            Set of Homomorphisms from Multivariate Polynomial Ring in x, y over Rational Field to Multivariate Polynomial Ring in x, y over Rational Field localized at (y, x)

        """
        return Hom(self.domain().numerator_ring(), self.codomain())

    def _element_constructor_(self, *args, **kwds):
        r"""
        Construct a morphism in this homset.

        EXAMPLES::

            sage: R.<x, y> = QQ[]
            sage: L.<xl, yl> = R.localization(x*y)
            sage: H = End(L)

        We can pass in the images of the generators::

            sage: H([y, x])
            Ring endomorphism of Multivariate Polynomial Ring in x, y over Rational Field localized at (y, x)
              Defn: xl |--> yl
                    yl |--> xl

        or a morphism at the level of numerators::

            sage: g = R.hom([y, x])
            sage: g
            Ring endomorphism of Multivariate Polynomial Ring in x, y over Rational Field
              Defn: x |--> y
                    y |--> x
            sage: H(g)
            Ring endomorphism of Multivariate Polynomial Ring in x, y over Rational Field localized at (y, x)
              Defn: xl |--> yl
                    yl |--> xl

        """
        return self.element_class(self, *args, **kwds)

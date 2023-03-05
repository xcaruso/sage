r"""
Localization

Localization is an important ring construction tool. Whenever you have
to extend a given integral domain such that it contains the inverses
of a finite set of elements but should allow non injective homomorphic
images this construction will be needed. See the example on
Ariki-Koike algebras below for such an application.

EXAMPLES::

    sage: LZ = ZZ.localization((5,11))
    sage: m = matrix(LZ, [[5, 7], [0,11]])
    sage: m.parent()
    Full MatrixSpace of 2 by 2 dense matrices over Integer Ring localized at (5, 11)
    sage: ~m      # parent of inverse is different: see documentation of m.__invert__
    [  1/5 -7/55]
    [    0  1/11]
    sage: _.parent()
    Full MatrixSpace of 2 by 2 dense matrices over Rational Field
    sage: mi = matrix(LZ, ~m)
    sage: mi.parent()
    Full MatrixSpace of 2 by 2 dense matrices over Integer Ring localized at (5, 11)
    sage: mi == ~m
    True

The next example defines the most general ring containing the coefficients of the irreducible
representations of the Ariki-Koike algebra corresponding to the three colored permutations on
three elements::

    sage: R.<u0, u1, u2, q> = ZZ[]
    sage: u = [u0, u1, u2]
    sage: S = Set(u)
    sage: I = S.cartesian_product(S)
    sage: add_units = u + [q, q + 1] + [ui - uj for ui, uj in I if ui != uj]
    sage: add_units += [q*ui - uj for ui, uj in I if ui != uj]
    sage: L = R.localization(tuple(add_units)); L
    Multivariate Polynomial Ring in u0, u1, u2, q over Integer Ring localized at
    (q, q + 1, u2, u1, u1 - u2, u0, u0 - u2, u0 - u1, u2*q - u1, u2*q - u0,
    u1*q - u2, u1*q - u0, u0*q - u2, u0*q - u1)

Define the representation matrices (of one of the three dimensional irreducible representations)::

    sage: m1 = matrix(L, [[u1, 0, 0],[0, u0, 0],[0, 0, u0]])
    sage: m2 = matrix(L, [[(u0*q - u0)/(u0 - u1), (u0*q - u1)/(u0 - u1), 0],
    ....:                 [(-u1*q + u0)/(u0 - u1), (-u1*q + u1)/(u0 - u1), 0],
    ....:                 [0, 0, -1]])
    sage: m3 = matrix(L, [[-1, 0, 0],
    ....:                 [0, u0*(1 - q)/(u1*q - u0), q*(u1 - u0)/(u1*q - u0)],
    ....:                 [0, (u1*q^2 - u0)/(u1*q - u0), (u1*q^ 2 - u1*q)/(u1*q - u0)]])
    sage: m1.base_ring() == L
    True

Check relations of the Ariki-Koike algebra::

    sage: m1*m2*m1*m2 == m2*m1*m2*m1
    True
    sage: m2*m3*m2 == m3*m2*m3
    True
    sage: m1*m3 == m3*m1
    True
    sage: m1^3 -(u0+u1+u2)*m1^2 +(u0*u1+u0*u2+u1*u2)*m1 - u0*u1*u2 == 0
    True
    sage: m2^2 -(q-1)*m2 - q == 0
    True
    sage: m3^2 -(q-1)*m3 - q == 0
    True
    sage: ~m1 in m1.parent()
    True
    sage: ~m2 in m2.parent()
    True
    sage: ~m3 in m3.parent()
    True

Obtain specializations in positive characteristic::

    sage: Fp = GF(17)
    sage: f = L.hom((3,5,7,11), codomain=Fp); f
    Ring morphism:
      From: Multivariate Polynomial Ring in u0, u1, u2, q over Integer Ring localized at
      (q, q + 1, u2, u1, u1 - u2, u0, u0 - u2, u0 - u1, u2*q - u1, u2*q - u0,
      u1*q - u2, u1*q - u0, u0*q - u2, u0*q - u1)
      To:   Finite Field of size 17
      Defn: u0 |--> 3
            u1 |--> 5
            u2 |--> 7
            q |--> 11
    sage: mFp1 = matrix({k:f(v) for k, v in m1.dict().items()}); mFp1
    [5 0 0]
    [0 3 0]
    [0 0 3]
    sage: mFp1.base_ring()
    Finite Field of size 17
    sage: mFp2 = matrix({k:f(v) for k, v in m2.dict().items()}); mFp2
    [ 2  3  0]
    [ 9  8  0]
    [ 0  0 16]
    sage: mFp3 = matrix({k:f(v) for k, v in m3.dict().items()}); mFp3
    [16  0  0]
    [ 0  4  5]
    [ 0  7  6]

Obtain specializations in characteristic 0::

    sage: fQ = L.hom((3,5,7,11), codomain=QQ); fQ
    Ring morphism:
      From: Multivariate Polynomial Ring in u0, u1, u2, q over Integer Ring localized at
            (q, q + 1, u2, u1, u1 - u2, u0, u0 - u2, u0 - u1, u2*q - u1, u2*q - u0,
            u1*q - u2, u1*q - u0, u0*q - u2, u0*q - u1)
      To:   Rational Field
      Defn: u0 |--> 3
            u1 |--> 5
            u2 |--> 7
            q |--> 11
    sage: mQ1 = matrix({k:fQ(v) for k, v in m1.dict().items()}); mQ1
    [5 0 0]
    [0 3 0]
    [0 0 3]
    sage: mQ1.base_ring()
    Rational Field
    sage: mQ2 = matrix({k:fQ(v) for k, v in m2.dict().items()}); mQ2
    [-15 -14   0]
    [ 26  25   0]
    [  0   0  -1]
    sage: mQ3 = matrix({k:fQ(v) for k, v in m3.dict().items()}); mQ3
    [    -1      0      0]
    [     0 -15/26  11/26]
    [     0 301/26 275/26]

    sage: S.<x, y, z, t> = QQ[]
    sage: T = S.quo(x+y+z)
    sage: F = T.fraction_field()
    sage: fF = L.hom((x, y, z, t), codomain=F); fF
    Ring morphism:
      From: Multivariate Polynomial Ring in u0, u1, u2, q over Integer Ring localized at
            (q, q + 1, u2, u1, u1 - u2, u0, u0 - u2, u0 - u1, u2*q - u1, u2*q - u0,
            u1*q - u2, u1*q - u0, u0*q - u2, u0*q - u1)
      To:   Fraction Field of Quotient of Multivariate Polynomial Ring in x, y, z, t over
            Rational Field by the ideal (x + y + z)
      Defn: u0 |--> -ybar - zbar
            u1 |--> ybar
            u2 |--> zbar
            q |--> tbar
    sage: mF1 = matrix({k:fF(v) for k, v in m1.dict().items()}); mF1
    [        ybar            0            0]
    [           0 -ybar - zbar            0]
    [           0            0 -ybar - zbar]
    sage: mF1.base_ring() == F
    True

TESTS::

    sage: TestSuite(L).run()

AUTHORS:

- Sebastian Oehms (2019-12-09): initial version.
- Sebastian Oehms (2022-03-05): fix some corner cases and add :meth:`factor` (:trac:`33463`)
- Xavier Caruso, Thibaut Verron (2023-03): complete rewrite; allow localization at zero divisors
"""


from sage.categories.homset import Hom
from sage.categories.integral_domains import IntegralDomains
from sage.categories.rings import Rings

from sage.rings.homset import RingHomset_generic
from sage.rings.ideal import Ideal_generic
from sage.rings.integer_ring import ZZ
from sage.rings.morphism import RingHomomorphism
from sage.rings.ring import CommutativeRing, CommutativeRing

from sage.misc.cachefunc import cached_method
from sage.misc.latex import latex

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import Element
from sage.structure.parent_gens import localvars

from sage.structure.richcmp import richcmp, rich_to_bool
from sage.structure.richcmp import op_EQ, op_NE


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
        sage: normalize_extra_units(P, [3*x, z*y^2, 2*z, 18*(x*y*z)^2, x*z, 6*x*z, 5])
        [2, 3, 5, z, y, x]
        sage: P.<x,y,z> = QQ[]
        sage: normalize_extra_units(P, [3*x, z*y^2, 2*z, 18*(x*y*z)^2, x*z, 6*x*z, 5])
        [z, y, x]

        sage: R.<x, y> = ZZ[]
        sage: Q.<a, b> = R.quo(x^2-5)
        sage: p = b^2-5
        sage: p == (b-a)*(b+a)
        True
        sage: normalize_extra_units(Q, [p])
        doctest:...: UserWarning: Localization may not be represented uniquely
        [b^2 - 5]
        sage: normalize_extra_units(Q, [p], warning=False)
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
        except (NotImplementedError, AttributeError):
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

        EXAMPLES:

        Example with simplification::

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
                    units = parent._extra_units
                    for u in units:
                        while u.divides(num) and u.divides(denom):
                            num = num // u
                            denom = denom // u
                except NotImplementedError:
                    parent._has_floordiv = False
        self._num = num
        self._denom = denom

    def _repr_(self):
        r"""
        How this element should be represented

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
        with localvars(S,R.variable_names()):
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

            sage: P.<x,y> = QQ[]
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

        ..TESTS::

            sage: P.<x,y> = QQ[]
            sage: L.<xl,yl> = P.localization([x+y,y])
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

        EXAMPLES::

            sage: P.<x,y> = QQ[]
            sage: PP.<xb,yb> = P.quotient(x^2*y)
            sage: L = PP.localization([xb])
            sage: elt = L(xb*yb,xb^4)
            sage: hash(elt) # random
            211140964612250627
            sage: dict = {elt} # indirect_doctest
            sage: elt in dict
            True

        """
        return hash((self._num, self._denom))

    def numerator(self):
        r"""

        EXAMPLES::

            sage: R = ZZ.localization([2,3])
            sage: R(4).numerator()
            4
            sage: R(4,2^2*3).numerator()
            1

        """
        return self._num

    def denominator(self):
        r"""

        EXAMPLES::

            sage: R = ZZ.localization([2,3])
            sage: R(4).denominator()
            1
            sage: R(4,2^2*3).denominator()
            3
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

    def _sub_(self,other):
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

    def is_unit(self):
        """
        Return ``True`` if this element is a unit.

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: L = P.localization((x, y*z))
            sage: L(y*z).is_unit()
            True
            sage: L(z).is_unit()
            True
            sage: L(x*y*z).is_unit()
            True

        ::

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
            sage: L = P.localization(x*y*z)
            sage: L(x*y*z).inverse_of_unit()
            1/(x*y*z)
            sage: L(z).inverse_of_unit()
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

    def factor(self, proof=None):
        r"""
        Return the factorization of this polynomial.

        INPUT:

        - ``proof`` -- (optional) if given it is passed to the
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
        if proof is not None:
            F = self._num.factor(proof=proof)
        else:
            F = self._num.factor()
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
           sage: L(1/x) < L(3/(x*y*z)^3)
           False
           sage: ~L(y*z/x) == L(x/(y*z))
           True
        """
        if not self.parent()._has_equality_test:
            return NotImplemented
        a = self._num * other._denom
        b = self._denom * other._num
        eq = (a - b) in self.parent()._ideal
        if eq or op == op_EQ or op == op_NE:
            return rich_to_bool(op, 0 if eq else -1)
        else:
            return richcmp(a, b, op)

    def _rational_(self):
        r"""
        Convert this element to a rational.

        This is only possible if its base ring is the ring of integers.

        OUTPUT:

        A rational.

        TESTS::

            sage: L = ZZ.localization(5)
            sage: cp3 = cyclotomic_polynomial(3).change_ring(L)
            sage: cp3.splitting_field('t')      #   indirect doctest
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
    The localization generalizes the construction of the field of fractions of an integral domain to
    an arbitrary ring. Given a (not necessarily commutative) ring `R` and a subset `S` of `R`,
    there exists a ring `R[S^{-1}]` together with the ring homomorphism `R \longrightarrow R[S^{-1}]`
    that "inverts" `S`; that is, the homomorphism maps elements in `S` to unit elements in `R[S^{-1}]`
    and, moreover, any ring homomorphism from `R` that "inverts" `S` uniquely factors through `R[S^{-1}]`.

    The ring `R[S^{-1}]` is called the *localization* of `R` with respect to `S`. For example, if `R` is
    a commutative ring and `f` an element in `R`, then the localization consists of elements of the form
    `r/f, r\in R, n \geq 0` (to be precise, `R[f^{-1}] = R[t]/(ft-1)`.

    The above text is taken from `Wikipedia`. The construction here used for this class relies on the
    construction of the field of fraction and is therefore restricted to integral domains.

    Accordingly, this class is inherited from :class:`IntegralDomain` and can only be used in that context.
    Furthermore, the base ring should support :meth:`sage.structure.element.CommutativeRingElement.divides` and
    the exact division operator `//` (:meth:`sage.structure.element.Element.__floordiv__`) in order to guarantee
    an successful application.

    INPUT:

    - ``base_ring`` -- an instance of :class:`Ring` allowing the construction of :meth:`fraction_field` (that is an integral domain)

    - ``extra_units`` -- tuple of elements of ``base_ring`` which should be turned into units

    - ``names`` -- passed to :class:`IntegralDomain`

    - ``normalize`` -- (optional, default: True) passed to :class:`IntegralDomain`

    - ``category`` -- (optional, default: None) passed to :class:`IntegralDomain`

    - ``warning`` -- (optional, default: True) to suppress a warning which is thrown if self cannot be represented uniquely

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
        sage: L = R.localization(x^2+1)
        sage: s = (x+5)/(x^2+1)
        sage: s in L
        True
        sage: t = (x+5)/(x^2+2)
        sage: t in L
        False
        sage: L(t)
        Traceback (most recent call last):
        ...
        ArithmeticError: denominator is not a unit
        sage: L(s) in R
        False
        sage: y = L(x)
        sage: g = L(s)
        sage: g.parent()
        Univariate Polynomial Ring in x over Integer Ring localized at (x^2 + 1,)
        sage: f = (y+5)/(y^2+1); f
        (x + 5)/(x^2 + 1)
        sage: f == g
        True

        sage: Lau.<u, v> = LaurentPolynomialRing(ZZ)
        sage: LauL = Lau.localization(u+1)
        sage: LauL(~u).parent()
        Multivariate Laurent Polynomial Ring in u, v over Integer Ring localized at (u + 1,)

    Localization by a divisor of 0::

        sage: P.<x,y> = QQ[]
        sage: PP.<xb,yb> = P.quotient(x^2-y^2)
        sage: L = PP.localization([xb-yb]); L
        Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2 - y^2) localized at (xb - yb,)

    Localization by a non-divisor of 0::

        sage: L = PP.localization([xb-yb^2]); L
        Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2 - y^2) localized at (-yb^2 + xb,)

    TESTS:

    Check that :trac:`33463` is fixed::

        sage: R = ZZ.localization(5)
        sage: R.localization(~5)
        Integer Ring localized at (5,)

    """
    Element = LocalizationElement

    def __classcall_private__(cls, base, units, names=None, normalize=True, category=None):
        if not isinstance(units, (list, tuple)):
            units = [units]
        units = tuple(base(u) for u in units)
        if isinstance(names, list):
            names = tuple(names)
        if base in IntegralDomains():
            if category is None:
                category = IntegralDomains()
            else:
                category = category.join([IntegralDomains()])
        return UniqueRepresentation.__classcall__(cls, base, units, names, normalize, category)

    def __init__(self, base, units, names, normalize, category):
        """
        Python constructor of Localization.

        TESTS::

            sage: L = ZZ.localization((3,5))
            sage: TestSuite(L).run()

            sage: R.<x> = ZZ[]
            sage: L = R.localization(x^2+1)
            sage: TestSuite(L).run()
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
            _ = self.fraction_field()

    def _element_constructor_(self, *args, **kwds):
        """
        Return an element of this parent constructed from the given arguments.

        EXAMPLES::

            sage: L = ZZ.localization((5, 2))
            sage: L(1/25)  # indirect doctest
            1/25
            sage: L(1/20)  # indirect doctest
            1/20
            sage: L(1/10)  # indirect doctest
            1/10
        """
        kwds['check'] = True
        return self.element_class(self, *args, **kwds)

    def _repr_(self):
        r"""
        String representation of this ring.

        EXAMPLES::

            sage: P.<x,y> = QQ[]
            sage: PP.<xb,yb> = P.quotient(x^2-y^2)
            sage: L = PP.localization([xb-yb]); L
            Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2 - y^2) localized at (xb - yb,)

        """
        if self._extra_units:
            return f"{self._numring} localized at {self._extra_units}"
        else:
            return str(self._numring)

    def _latex_(self):
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
            sage: R.localization((x^2+1, y-1)).gens()
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
        return LocalizationIdeal_generic

    def units(self, normalize=False):
        if normalize:
            return self._extra_units
        else:
            return self._given_units

    def numerator_kernel(self):
        return self._ideal

    def numerator_ring(self):
        return self._numring

    def is_integral_domain(self):
        if self in IntegralDomains():
            return True
        return self._ideal.is_prime()

    def inverse_of_units(self, normalize=False):
        if normalize:
            return tuple(self(u).inverse_of_unit() for u in self._extra_units)
        else:
            return tuple(self(u).inverse_of_unit() for u in self._given_units)

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
        if not self.is_integral_domain():
            raise TypeError
        field = self._numring.quotient(self._ideal).fraction_field()
        field.register_coercion(self)
        return field

    def is_field(self, proof=False):
        if self.base_ring().is_field(proof=proof):
            return True
        if proof is False:
            # It's never a field if the base ring is not semi-local
            return False
        else:
            raise NotImplementedError

    def characteristic(self):
        """
        Return the characteristic of this ring.

        EXAMPLES::

            sage: R.<a> = GF(5)[]
            sage: L = R.localization((a^2-3, a))
            sage: L.characteristic()
            5
        """
        return self.base_ring().characteristic()

    def _flattening_function(self):
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
            units = S._extra_units + [ g(R(u)).numerator() for u in self._extra_units ]
            ring = localization_with_simplification(base, units)
            def isom(x):  # self -> ring
                num = g(x.numerator())      # in S
                denom = g(x.denominator())  # in S
                num2 = num.numerator() * denom.denominator()
                denom2 = num.denominator() * denom.numerator()
                return ring(num2) * ring(denom2).inverse_of_unit()
        else:
            base = S
            units = [ g(u) for u in self._extra_units ]
            ring = localization_with_simplification(base, units)
            def isom(x):  # self -> ring
                num = g(x.numerator())      # in S
                denom = g(x.denominator())  # in S
                return ring(num) * ring(denom).inverse_of_unit()
        return ring, isom

    def flattening_morphism(self):
        ring, isom = self._flattening_function()
        im_gens = [ isom(x) for x in self.gens() ]
        return self._Hom_(ring, category=Rings())(im_gens)

    def flatten(self):
        ring, _ = self._flattening_function()
        return ring


    def _Hom_(self, other, category):
        return RingHomset_localized(self, other, category=category)


def localization_with_simplification(base, units):
    simplified_units = [ ]
    for u in units:
        try:
            if not u.is_unit():
                simplified_units.append(u)
        except NotImplementedError:
            simplified_units.append(u)
    if simplified_units:
        return Localization(base, simplified_units)
    else:
        return base


class LocalizationIdeal_generic(Ideal_generic):
    @cached_method
    def numerator_ideal(self):
        r"""

        EXAMPLES::

            sage: P.<x,y> = QQ[]
            sage: PP.<xb,yb> = P.quotient(x^2*y)
            sage: L = PP.localization([xb])
            sage: I = L.ideal(L(yb,xb^2))
            sage: I
            Ideal (yb/xb^2) of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2*y) localized at (xb,)
            sage: I.numerator_ideal()
            Ideal (yb) of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2*y)

        """
        R = self.ring()
        S = R.numerator_ring()
        gens_S = [g.numerator() for g in self.gens()]
        I_S = S.ideal(gens_S)
        for u in R._extra_units:
            I_S,_ = I_S.saturation(u)
        return I_S

    def _contains_(self,elt):
        r"""

        EXAMPLES::

            sage: P.<x,y> = QQ[]
            sage: PP.<xb,yb> = P.quotient(x^2*y)
            sage: L = PP.localization([xb])
            sage: I = L.ideal(L(yb,xb^2))
            sage: L(xb*yb,xb^4) in I
            True

        """
        return elt.numerator() in self.numerator_ideal()

    def is_zero(self):
        return all([ x == 0 for x in self.gens() ])

    def is_one(self):
        return self.numerator_ideal().is_one()

    def saturation(self,other):
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
        return self.numerator_ideal().is_prime()


# Morphisms
###########

class LocalizationMorphism(RingHomomorphism):
    def __init__(self, parent, numerator_morphism, base_map=None, check=True):
        domain = parent.domain()
        if not isinstance(domain, Localization):
            raise TypeError
        RingHomomorphism.__init__(self, parent)
        self._numerator_morphism = f = parent.numerator_homset()(numerator_morphism, base_map=base_map, check=check)
        if check:
            for u in domain._extra_units:
                if not f(u).is_unit():
                    raise ValueError("localized element %s is not mapped to a unit" % u)

    def numerator_morphism(self):
        return self._numerator_morphism

    @cached_method
    def im_gens(self):
        return [ self(x) for x in self.domain().gens() ]

    def _call_(self, x):
        f = self._numerator_morphism
        num = x.numerator()
        denom = x.denominator()
        return f(num) * f(denom).inverse_of_unit()

    def _repr_defn(self):
        return self._numerator_morphism._repr_defn()


class RingHomset_localized(RingHomset_generic):
    Element = LocalizationMorphism

    @cached_method
    def numerator_homset(self):
        return Hom(self.domain().numerator_ring(), self.codomain())

    def _element_constructor_(self, *args, **kwds):
        return self.element_class(self, *args, **kwds)

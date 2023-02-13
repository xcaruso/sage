# TODO: Take documentation from localization.py and adapt the tests
# (Main missing: possible normalization of units, and hom)

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

from sage.structure.element import Element, coerce_binop
from sage.structure.parent_gens import localvars

# Elements
##########

def par_str(a):
    if a._is_atomic():
        return str(a)
    else:
        return f"({a})"

def normalize_extra_units(base_ring, add_units, warning=True):
    # Minimal at the moment, we should see what we can do more (TODO)
    # At the very least we should remove duplicates and opposites
    # We should probably also factor when possible?
    # The original class did sort, we cannot do that without breaking inverse_of_units
    return tuple([ base_ring(u) for u in add_units ])

class LocalizedRingElement(Element):
    def __init__(self, parent, *x, simplify=True):
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
            if isinstance(x, LocalizedRingElement):
                num = R(x._num)
                denom = R(x._denom)
            else:
                try:
                    num = R(x)
                    denom = R.one()
                except (TypeError, ZeroDivisionError):
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
        if not (denom.is_one() or parent.ideal(denom).is_one()):
            raise ArithmeticError("denominator is a unit")

        if num.is_zero():
            denom = R.one()
        elif simplify:
            if parent._has_gcd:
                g = num.gcd(denom)
                num = num // g
                denom = denom // g
            elif parent._has_floordiv:
                units = parent._units
                for u in units:
                    while u.divides(num) and u.divides(denom):
                        num = num // u
                        denom = denom // u
        self._num = num
        self._denom = denom

    def _repr_(self):
        r"""
        How this element should be represented

        EXAMPLES::

            sage: P.<x,y> = QQ[]
            sage: L = P.localization([x+y,y])
            sage: u,v = L.inverse_of_units()
            sage: L(1)
            1
            sage: u
            1/(x + y)
            sage: v
            1/y
            sage: L(x-y)*u*v^2
            (x - y)/(x*y^2 + y^3)
            sage: L(x)*u*v^2
            x/(x*y^2 + y^3)

        ..TESTS::

            sage: P.<x,y> = QQ[]
            sage: L.<xl,yl> = P.localization([x+y,y])
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
                return str(self._num)
            else:
                return "%s/%s" % (par_str(self._num), par_str(self._denom))

    def _latex_(self):
        r"""
        How this element should be represented in latex.

        EXAMPLES::

            sage: P.<x,y> = QQ[]
            sage: L = P.localization([x+y,y])
            sage: u,v = L.inverse_of_units()
            sage: latex(L(1))
            1
            sage: latex(u)
            \frac{1}{x + y}
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

        EXAMPLES::

            sage: R = ZZ.localization([2,3])
            sage: R(4,1) + R(6,1)
            10
            sage: R(4,3) + R(6,2)
            13/3
            sage: R(4,2^2*3^3) + R(6,2^3*3^2)
            13/108

        """
        num = self._num * other._denom + self._denom * other._num
        denom = self._denom * other._denom
        return self.__class__(self.parent(), num, denom)

    def _sub_(self,other):
        r"""

        EXAMPLES::

            sage: R = ZZ.localization([2,3])
            sage: R(4,1) - R(6,1)
            -2
            sage: R(4,3) - R(6,2)
            -5/3
            sage: R(4,2^2*3^3) - R(6,2^3*3^2)
            -5/108

        """
        num = self._num * other._denom - self._denom * other._num
        denom = self._denom * other._denom
        return self.__class__(self.parent(), num, denom)

    def _mul_(self,other):
        r"""

        EXAMPLES::

            sage: R = ZZ.localization([2,3])
            sage: R(4,1) * R(6,1)
            24
            sage: R(4,3) * R(6,2)
            4
            sage: R(4,2^2*3^3) * R(6,2^3*3^2)
            1/324
        
        """
        num = self._num * other._num
        denom = self._denom * other._denom
        return self.__class__(self.parent(), num, denom)

    def _lmul_(self,a):
        num = self._num * a
        return self.__class__(self.parent(), num, self._denom)

    def is_unit(self):
        r"""

        EXAMPLES::

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
        EXAMPLES::

            sage: L = ZZ.localization([2,3])
            sage: L(1).inverse_of_unit()
            1
            sage: L(-1).inverse_of_unit() # optional - bug (outputs 1/-1)
            -1
            sage: L(2).inverse_of_unit()
            1/2
            sage: L(-3).inverse_of_unit() # optional - bug (outputs 1/-3)
            -1/3
            sage: a, b = L.inverse_of_units()
            sage: a.inverse_of_unit()
            2
            sage: L(6).inverse_of_unit()
            1/6
            sage: L(5).inverse_of_unit()
            Traceback (most recent call last):
            ...
            ArithmeticError: 5 is not invertible

        """
        if not self.is_unit():
            raise ArithmeticError("%s is not invertible" % self)
        return self.__class__(self.parent(), self._denom, self._num)

    def __invert__(self):
        K = self.parent().fraction_field()
        return K(self._denom, self._num)

    @coerce_binop
    def __eq__(self, other):
        r"""

        EXAMPLES::

            sage: R = ZZ.localization([2,3])
            sage: R(2) == R(4,2)
            True
            sage: R(2) == R(6,2)
            False

            sage: A.<x,y> = QQ[] # 2 variables because no saturate for univariate ideals
            sage: AA.<xx,yy> = A.quotient(x^2)
            sage: L = AA.localization([xx])
            sage: L(1) == L(0)
            True
            sage: L(x) == L(0)
            True
            sage: L(y) == L(0)
            True
        """
        if not self.parent()._has_equality_test:
            raise TypeError("no equality test in this localization")
        diff = self._num * other._denom - self._denom * other._num
        return diff in self.parent()._ideal

    def __ne__(self, other):
        return not self.__eq__(other)

    def _rational_(self):
        # Only if base ring is ZZ
        if not self.parent()._numring is ZZ:
            raise ValueError("not a rational")
        else:
            return self.numerator()/self.denominator()


# Parents
#########

class LocalizedRing(CommutativeRing):
    r"""


    EXAMPLES:

    Integral case::
    
        sage: R = ZZ.localization([2,3]); R
        Integer Ring localized at (2, 3)
        sage: R._numring
        Integer Ring
        sage: R._ideal
        Principal ideal (0) of Integer Ring
        sage: R._units
        (2, 3)

    Localization by a divisor of 0::

        sage: P.<x,y> = QQ[]
        sage: PP.<xb,yb> = P.quotient(x^2-y^2)
        sage: L = PP.localization([xb-yb]); L
        Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2 - y^2) localized at (xb - yb,)
        sage: L._ideal
        Ideal (xb + yb) of Quotient of  Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2 - y^2)
        sage: L._units
        (xb - yb,)

    Localization by a non-divisor of 0::

        sage: L = PP.localization([xb-yb^2]); L
        Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2 - y^2) localized at (-yb^2 + xb,)
        sage: L._ideal 
        Ideal (0) of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2 - y^2)
        sage: L._units
        (-yb^2 + xb,)

    Inverses of elements can be recovered::

        sage: L = PP.localization([xb])
        sage: xinv, = L.inverse_of_units()
        sage: xinv
        1/xb
        sage: xb*xinv
        xb/xb

    """

    Element = LocalizedRingElement

    def __init__(self, base, units, names=None, normalize=True, category=None):
        self._numring = base
        self._units = units = normalize_extra_units(base,units)

        J = base.ideal(0)
        self._is_ideal_domain = False
        try:
            if base in IntegralDomains():
                self._has_equality_test = True
                self._is_integral_domain = True
                if 0 in units:
                    J = base.ideal([1])
                if category is None:
                    category = IntegralDomains()
                else:
                    category = category.join([IntegralDomains()])
            else:
                _ = J.saturation
                self._has_equality_test = True
                for u in units:
                    J, _ = J.saturation(u)
        except (AttributeError, NotImplementedError) as e:
            self._has_equality_test = False

        self._ideal = J
        self._has_floordiv = hasattr(base.zero(), "_floordiv_")
        self._has_gcd = self._has_floordiv and hasattr(base.zero(), "gcd")

        if names is None:
            names = base.variable_names()
        self._assign_names(names)
        CommutativeRing.__init__(self, base.base_ring(), names=names, normalize=normalize, category=category)
        self.register_coercion(base)

        self._gens = tuple(self(v) for v in base.gens())

    def _repr_(self):
        r"""
        How this localized ring should be represented.

        EXAMPLES::

            sage: P.<x,y> = QQ[]
            sage: PP.<xb,yb> = P.quotient(x^2-y^2)
            sage: L = PP.localization([xb-yb]); L
            Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2 - y^2) localized at (xb - yb,)


        """
        return f"{self._numring} localized at {self._units}"

    def _latex_(self):
        s = latex(self._numring) + r"\left["
        invs = [f"\\frac{{1}}{{{latex(u)}}}" for u in self._units]
        s += ", ".join(invs) + r"\right]"
        return s

    def ngens(self):
        return self._numring.ngens()

    def gens(self):
        return self._gens

    def gen(self,i):
        return self._gens[i]

    def _ideal_class_(self, num_gens):
        return LocalizedRingIdeal_generic

    def units(self):
        return self._units

    def numerator_kernel(self):
        return self._ideal

    def numerator_ring(self):
        return self._numring

    def is_integral_domain(self):
        if self in IntegralDomains():
            return True
        return self._ideal.is_prime()

    def inverse_of_units(self):
        return [self(u).inverse_of_unit() for u in self._units]

    @cached_method
    def fraction_field(self):
        if not self.is_integral_domain():
            raise TypeError
        field = self._numring.quotient(self._ideal).fraction_field()
        field.register_coercion(self)
        return field

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
        if isinstance(S, LocalizedRing):
            base = S._numring
            units = S._units + [ g(R(u)).numerator() for u in self._units ]
            ring = localization_with_simplification(base, units)
            def isom(x):  # self -> ring
                num = g(x.numerator())      # in S
                denom = g(x.denominator())  # in S
                num2 = num.numerator() * denom.denominator()
                denom2 = num.denominator() * denom.numerator()
                return ring(num2) * ring(denom2).inverse_of_unit()
        else:
            base = S
            units = [ g(u) for u in self._units ]
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
        return LocalizedRing(base, simplified_units)
    else:
        return base


class LocalizedRingIdeal_generic(Ideal_generic):
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
        for u in R._units:
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

class LocalizedRingMorphism(RingHomomorphism):
    def __init__(self, parent, numerator_morphism, base_map=None, check=True):
        domain = parent.domain()
        if not isinstance(domain, LocalizedRing):
            raise TypeError
        RingHomomorphism.__init__(self, parent)
        self._numerator_morphism = f = parent.numerator_homset()(numerator_morphism, base_map=base_map, check=check)
        if check:
            for u in domain._units:
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
    Element = LocalizedRingMorphism

    @cached_method
    def numerator_homset(self):
        return Hom(self.domain().numerator_ring(), self.codomain())

    def _element_constructor_(self, *args, **kwds):
        return self.element_class(self, *args, **kwds)

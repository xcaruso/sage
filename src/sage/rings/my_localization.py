from sage.rings.ring import CommutativeRing, CommutativeRing
from sage.rings.ring_extension import RingExtension
from sage.rings.ideal import Ideal_generic
from sage.categories.integral_domains import IntegralDomains

from sage.structure.element import Element, coerce_binop


def pow_str(a,b):
    r"""
    Representation of a^b
    """
    if b > 1:
        return f"({a})^{b}"
    else:
        return f"({a})"

class MyLocalizationElement(Element):
    def __init__(self, parent, *x, simplify=True):
        r"""

        EXAMPLES:

        Example with simplification::

            sage: L = MyLocalization(ZZ,[2,3])
            sage: L(4)
            4
            sage: L(4,[1,0])
            2
            sage: L(1,[2,0]) + L(3,[2,0])
            1

            sage: P.<x,y> = QQ[]
            sage: PP.<xb,yb> = P.quotient(x^2*y)
            sage: L = MyLocalization(PP,[xb])
            sage: L._has_simplify
            False
            sage: L(xb,[1])
            xb/((xb))
            sage: L(0,[18])
            0

        """
        Element.__init__(self, parent)
        if len(x)==1:
            x = x[0]
            if isinstance(x, MyLocalizationElement):
                num = x._num
                powers = x._powers[:]
            elif x in self.parent().numerator_ring():
                num = x
                powers = [0] * len(self.parent()._units)
            else:
                raise ValueError    
        else:
            num, powers = x

        if num.is_zero():
            powers = [0] * len(self.parent()._units)
        elif simplify and self.parent()._has_simplify:
            for i in range(len(powers)):
                p = powers[i]
                u = self.parent()._units[i]
                while p > 0 and u.divides(num):
                    num = num // u
                    p -= 1
                powers[i] = p
        self._num = num
        self._powers = powers

    @cached_method
    def _repr_(self):
        r"""


        EXAMPLES::

            sage: R = MyLocalization(ZZ,[2,3])
            sage: R(4,[0,0])
            4
            sage: R(4,[1,2])
            4/((2)*(3)^2)
            sage: R(4,[3,0])
            4/((2)^3)
        
        """
        nunits = len(self._powers)
        pos = [i for i in range(nunits) if self._powers[i] != 0]
        if pos:
            den = pow_str(self.parent()._units[pos[0]], self._powers[pos[0]])
            for i in pos[1:] :
                den += "*" + pow_str(self.parent()._units[pos[i]], self._powers[pos[i]])
            res = f"{self._num}/({den})"
        else:
            res = f"{self._num}"
        import re

        # In case the localization has different names than the original ring,
        # we want to use the new names
        R = self.parent()
        S = R.numerator_ring()
        for i in range(S.ngens()):
            res = re.sub(f"\\b{S.variable_names()[i]}\\b", R.variable_names()[i], res)
        return res

    def _get_common_den(self,other):
        return [max(self._powers[i],other._powers[i]) for i in range(len(self._powers))]

    def __hash__(self):
        r"""

        EXAMPLES::

            sage: P.<x,y> = QQ[]
            sage: PP.<xb,yb> = P.quotient(x^2*y)
            sage: L = MyLocalization(PP,[xb])
            sage: elt = L(xb*yb,[4]) 
            sage: hash(elt) # random
            211140964612250627

        """
        return hash((self._num,tuple(self._powers)))

    def numerator(self):
        r"""

        EXAMPLES::

            sage: R = MyLocalization(ZZ,[2,3])
            sage: R(4).numerator()
            4
            sage: R(4,[2,3]).numerator()
            4

        """
        return self._num

    def denominator(self):
        r"""

        EXAMPLES::

            sage: R = MyLocalization(ZZ,[2,3])
            sage: R(4).denominator()
            1
            sage: R(4,[2,1]).denominator()
            12
        
        """
        return prod(self.parent()._units[i]**self._powers[i]
                    for i in range(len(self.parent()._units)))
    
    def _add_(self,other):
        r"""

        EXAMPLES::

            sage: R = MyLocalization(ZZ,[2,3])
            sage: R(4,[0,0]) + R(6,[0,0])
            10
            sage: R(4,[0,1]) + R(6,[1,0])
            26/((2)*(3))
            sage: R(4,[2,3]) + R(6,[3,2])
            26/((2)^3*(3)^3)

        """
        res_powers = self._get_common_den(other)
        units = self.parent()._units
        res_num = (self._num*prod(units[i]**(res_powers[i]-self._powers[i])
                                  for i in range(len(units)))
                   + other._num*prod(units[i]**(res_powers[i]-other._powers[i])
                                     for i in range(len(units))))
        return self.parent()(res_num,res_powers)

    def _sub_(self,other):
        r"""

        EXAMPLES::

            sage: R = MyLocalization(ZZ,[2,3])
            sage: R(4,[0,0]) - R(6,[0,0])
            -2
            sage: R(4,[0,1]) - R(6,[1,0])
            -10/((2)*(3))
            sage: R(4,[2,3]) - R(6,[3,2])
            -10/((2)^3*(3)^3)
        
        """
        res_powers = self._get_common_den(other)
        units = self.parent()._units
        res_num = (self._num*prod(units[i]**(res_powers[i]-self._powers[i])
                                  for i in range(len(units)))
                   - other._num*prod(units[i]**(res_powers[i]-other._powers[i])
                                     for i in range(len(units))))
        return self.parent()(res_num, res_powers)

    def _mul_(self,other):
        r"""

        EXAMPLES::

            sage: R = MyLocalization(ZZ,[2,3])
            sage: R(4,[0,0]) * R(6,[0,0])
            24
            sage: R(4,[0,1]) * R(6,[1,0])
            24/((2)*(3))
            sage: R(4,[2,3]) * R(6,[3,2])
            24/((2)^5*(3)^5)
       
        
        """
        res_num = self._num * other._num
        res_powers = [self._powers[i] + other._powers[i] for i in range(len(self._powers))]
        return self.parent()(res_num, res_powers)
    
    def _lmul_(self,a):
        res_num = self._num * a
        return self.parent()(res_num, self._powers[:])

    def is_unit(self):
        r"""

        EXAMPLES::

            sage: P.<x,y> = QQ[]
            sage: PP.<xx,yy> = P.quotient(x*y)
            sage: L = MyLocalization(PP,[xx])
            sage: L(x).is_unit()
            True
            sage: L(y).is_unit()
            False
        
            sage: PP.<xx,yy> = P.quotient(x^2)
            sage: L = MyLocalization(PP,[xx])
            sage: L(x).is_unit()
            True
            sage: L(y).is_unit()
            True

        """
        I = self.parent().ideal(self)
        I_numring = I.numerator_ideal()
        return I_numring.is_one()

    def inverse_of_unit(self):
        r"""

        EXAMPLES::

            sage: L.<a,b> = MyLocalization(ZZ,[2,3])
            sage: L(1).inverse_of_unit()
            1
            sage: L(-1).inverse_of_unit()
            1
            sage: L(2).inverse_of_unit()
            1/((2))
            sage: L(-3).inverse_of_unit()
            -1/((3))
            sage: a.inverse_of_unit()
            2
            sage: L(5).inverse_of_unit()
            Traceback (most recent call last):
            ...
            NotImplementedError:

        It is not currently possible to compute the inverse for all units.
        ::
        
            sage: L(6).inverse_of_unit()
            Traceback (most recent call last):
            ...
            NotImplementedError: 

        """
        units = self.parent()._units
        for i in [1,-1]:
            num = i*self._num
            if num == 1:
                res_num = num
                res_powers = [0]*len(units)
                break
            elif num in units:
                res_num = i
                res_powers = [1 if units[i] == num else 0 for i in range(len(units))]
                break
        else:
            raise NotImplementedError
        for i in range(len(units)):
            res_num *= units[i]**self._powers[i]
        return self.parent()(res_num, res_powers)
            
    @coerce_binop
    def __eq__(self,other):
        r"""

        EXAMPLES::

            sage: R = MyLocalization(ZZ,[2,3])
            sage: R(2) == R(4,[1,0])
            True
            sage: R(2) == R(6,[1,0])
            False

            sage: A.<x,y> = QQ[] # 2 variables because no saturate for univariate ideals
            sage: AA.<xx,yy> = A.quotient(x^2)
            sage: L = MyLocalization(AA,[xx])
            sage: L(1) == L(0)
            True
            sage: L(x) == L(0)
            True
            sage: L(y) == L(0)
            True
        
        """
        if not self.parent()._has_equality_test:
            raise TypeError("no equality test in this localization")
        diff = self.numerator()*other.denominator() - self.denominator()*other.numerator()
        return diff in self.parent()._ideal
        

class MyLocalization(CommutativeRing):
    r"""


    EXAMPLES:

    Integral case::
    
        sage: R = MyLocalization(ZZ,[2,3])
        sage: R._numring
        Integer Ring
        sage: R._ideal
        Principal ideal (0) of Integer Ring
        sage: R._units
        [2, 3]

    Localization by a divisor of 0 (not working until saturation in ideals of quotients is implemented)::
    
        sage: P.<x,y> = QQ[]
        sage: PP.<xb,yb> = P.quotient(x^2-y^2)
        sage: L = MyLocalization(PP,[xb-yb])
        sage: L._ideal # incorrect
        Ideal (x + y) of Multivariate Polynomial Ring in x, y over Rational Field
        sage: L._units
        [x - y]

    Localization by a non-divisor of 0::

        sage: L = MyLocalization(PP,[xb-yb^2])
        sage: L._ideal # incorrect
        Ideal (x^2 - y^2) of Multivariate Polynomial Ring in x, y over Rational Field
        sage: L._units
        [-y^2 + x]

    Elements can be given a name for input purposes::

        sage: L.<xinv> = MyLocalization(PP,[xb])
        sage: xinv
        1/((xb))
        sage: xb*xinv
        xb/((xb))

    """

    Element = MyLocalizationElement

    def __init__(self, base, units, names=None, inverse_names=None):
        A = base
        J = base.ideal(0)

        try:
            if A in IntegralDomains():
                self._has_equality_test = True
            else:
                _ = J.saturation
                self._has_equality_test = True
                for u in units:
                    J, _ = J.saturation(u)
        except (AttributeError, NotImplementedError) as e:
            self._has_equality_test = False

        try:
            _ = A.zero()._floordiv_
            self._has_simplify = True
        except AttributeError:
            self._has_simplify = False
    
        self._ideal = J
        self._numring = base
        self._units = units
        CommutativeRing.__init__(self,base.base_ring())
        self.register_coercion(base)

        if bool(inverse_names) and bool(names):
            raise ValueError("both names and inverse_names are provided")

        if inverse_names is None:
            inverse_names = tuple(f"u{randint(1000,9999)}" for _ in range(len(units)))
        if names is None:
            try:
                orig_names = base.variable_names()
            except AttributeError:
                orig_names = base.variable_names()
            names = orig_names + inverse_names
        self._assign_names(names)

        self._gens = (tuple(self(v,[0]*len(units)) for v in base.gens())
                      + tuple(self(base(1),
                                   [1 if j==i else 0 for j in range(len(units))])
                              for i in range(len(units))))

    def _repr_(self):
        return f"Localization of {self._numring} at {self._units}"

    def ngens(self):
        return len(self._gens)

    def gens(self):
        return self._gens

    def gen(self,i):
        return self._gens[i]

    def _ideal_class_(self, num_gens):
        return MyLocalizationIdeal_generic

    def units(self):
        return self._units

    def numerator_ideal(self):
        return self._ideal

    def numerator_ring(self):
        return self._numring
        

class MyLocalizationIdeal_generic(Ideal_generic):
    r"""
    """

    @cached_method
    def numerator_ideal(self):
        r"""

        EXAMPLES::

            sage: P.<x,y> = QQ[]
            sage: PP.<xb,yb> = P.quotient(x^2*y)
            sage: L = MyLocalization(PP,[xb])
            sage: I = L.ideal(L(yb,[2]))
            sage: I
            Ideal (yb/((xb)^2)) of Localization of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2*y) at [xb]
            sage: I.numerator_ideal()
            (Ideal (yb) of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2*y),
             0)

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
            sage: L = MyLocalization(PP,[xb])
            sage: I = L.ideal(L(yb,[2]))
            sage: L(xb*yb,[4]) in I
            True

        """
        return elt.numerator() in self.numerator_ideal()

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
    
    
# class MyLocalizationExtension(CommutativeRing):
#     def __init__(self,base,units):
#         if isinstance(base, QuotientRing_generic):
#             R = base.cover_ring()
#             I = base.defining_ideal()
#             units_lift = [u.lift() for u in units]
#             J = I
#             for u in units_lift:
#                 J = 
#         elif isinstance(base, MyLocalization):
#             R = base.numerator_ring()
#             U = base.units()
#             backend = MyLocalization(R,U+units)
#         else:
#             backend = MyLocalization(base,units)


# P = QQ[("x","y")]
# P.inject_variables()
# PP = P.quotient(x*y, names=("xb","yb"))
# PP.inject_variables()
# P2 = MyLocalization(PP,[xb], names=("xl","yl","u"))
# P2.inject_variables()

###

Z6 = Integers(6)
L = MyLocalization(Z6,[ZZ(3)])
LL = MyLocalization(L,[ZZ(2)])

#print(L.is_zero())

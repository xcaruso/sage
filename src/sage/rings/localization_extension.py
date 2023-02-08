from sage.rings.ring import CommutativeRing, CommutativeRing
from sage.rings.ideal import Ideal_generic


def pow_str(a,b):
    r"""
    Representation of a^b
    """
    if b > 1:
        return f"({a})^{b}"
    else:
        return f"({a})"

class MyLocalizationElement(CommutativeRingElement):
    def __init__(self,parent,*x):
        r"""

        EXAMPLES::

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
        CommutativeRingElement.__init__(self,parent)
        if len(x)==1:
            x = x[0]
            if isinstance(x, MyLocalizationElement):
                num = x._num
                powers = x._powers
            elif x in self.parent().base():
                num = x
                powers = [0] * len(self.parent()._units)
            else:
                raise ValueError    
        else:
            num, powers = x

        if num == 0:
            powers = [0] * len(self.parent()._units)
        elif self.parent()._has_simplify:
            for i in range(len(powers)):
                p = powers[i]
                u = self.parent()._units[i]
                while p > 0 and u.divides(num):
                    num //= u
                    p -= 1
                powers[i] = p
        self._num = num
        self._powers = powers
        
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
            return f"{self._num}/({den})"
        else:
            return f"{self._num}"

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
                    for i in range(self.parent().ngens()))
    
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
        return self.parent()(res_num, self._powers)

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
        sage: R._base
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

    """

    Element = MyLocalizationElement

    def __init__(self, base, units, names="x"):
        A = base
        J = A.ideal(0)
        try:
            _ = J.saturation
            self._has_equality_test = True
            for u in units:
                J, _ = J.saturation(u)
        except AttributeError:
            self._has_equality_test = False

        try:
            _ = A.zero()._floordiv_
            self._has_simplify = True
        except AttributeError:
            self._has_simplify = False
    
        self._ideal = J
        self._base = base
        self._units = units
        CommutativeRing.__init__(self,base)
        self.register_coercion(base)
        self._assign_names('x') # Necessary for over() but why?

        

    def _repr_(self):
        return f"Localization of {self._base} at {self._units}"

    def ngens(self):
        return len(self._units)

    def gen(self,i):
        return self(self._base(1),[1 if j==i else 0 for j in range(self.ngens())])

    def _ideal_class_(self, num_gens):
        return MyLocalizationIdeal_generic

    def base_ring(self):
        return self._base

    def units(self):
        return self._units

    def cover_ideal(self):
        return self._ideal
        

class MyLocalizationIdeal_generic(Ideal_generic):
    r"""
    """

    @cached_method
    def base_ideal(self):
        r"""

        EXAMPLES::

            sage: P.<x,y> = QQ[]
            sage: PP.<xb,yb> = P.quotient(x^2*y)
            sage: L = MyLocalization(PP,[xb])
            sage: I = L.ideal(L(yb,[2]))
            sage: I
            Ideal (yb/((xb)^2)) of Localization of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2*y) at [xb]
            sage: I.base_ideal()
            (Ideal (yb) of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2*y),
             0)

        """
        R = self.ring()
        S = R.base()
        gens_S = [S(g.numerator()) for g in self.gens()]
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
        return elt.numerator() in self.base_ideal()

    
    
# class MyLocalizationExtension(RingExtension):
#     def __init__(self,base,units):
#         if isinstance(base, QuotientRing_generic):
#             R = base.cover_ring()
#             I = base.defining_ideal()
#             # todo...
#         elif isinstance(base, MyLocalization):
#             R = base.ambient_ring()
#             U = base.units()
#             # todo...
#         else:
#             backend = MyLocalization(base,units)

from sage.rings.ring import CommutativeRing, CommutativeRing


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

        try:
            self._num = self.parent()._ambient(num)
        except TypeError:
            self._num = self.parent()._ambient(num.lift())
        self._num = self.parent()._ideal.reduce(self._num)
        self._powers = powers
        
    def _repr_(self):
        r"""


        EXAMPLES:

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

    def numerator(self):
        r"""

        EXAMPLES:

            sage: R = MyLocalization(ZZ,[2,3])
            sage: R(4).numerator()
            4
            sage: R(4,[2,3]).numerator()
            4

        """
        return self._num

    def denominator(self):
        r"""

        EXAMPLES:

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

        EXAMPLES:

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

        EXAMPLES:

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

        EXAMPLES:

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

        EXAMPLES:

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
        diff = self.numerator()*other.denominator() - self.denominator()*other.numerator()
        return diff in self.parent()._ideal
        

class MyLocalization(CommutativeRing):
    r"""


    EXAMPLES:

    Integral case
    
        sage: R = MyLocalization(ZZ,[2,3])
        sage: R._ambient
        Integer Ring
        sage: R._ideal
        Principal ideal (0) of Integer Ring
        sage: R._units
        [2, 3]

    Localization by a divisor of 0
    
        sage: P.<x,y> = QQ[]
        sage: PP.<xb,yb> = P.quotient(x^2-y^2)
        sage: L = MyLocalization(PP,[xb-yb])
        sage: L._ambient
        Multivariate Polynomial Ring in x, y over Rational Field
        sage: L._ideal
        Ideal (x + y) of Multivariate Polynomial Ring in x, y over Rational Field
        sage: L._units
        [x - y]

    Localization by a non-divisor of 0

        sage: L = MyLocalization(PP,[xb-yb^2])
        sage: L._ambient
        Multivariate Polynomial Ring in x, y over Rational Field
        sage: L._ideal
        Ideal (x^2 - y^2) of Multivariate Polynomial Ring in x, y over Rational Field
        sage: L._units
        [-y^2 + x]

    """

    Element = MyLocalizationElement

    def __init__(self, base, units):
        from sage.rings.quotient_ring import QuotientRing_generic
        if base in IntegralDomains():
            A = base
            J = A.ideal(0)
        elif isinstance(base, QuotientRing_generic):
            I = base.defining_ideal()
            A = base.ambient()
            J = I
            units = [u.lift() for u in units]
            for u in units:
                J, _ = J.saturation(u)
        else:
            raise NotImplementedError("..")
        self._ambient = A
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
        return self(self._ambient(1),[1 if j==i else 0 for j in range(self.ngens())])


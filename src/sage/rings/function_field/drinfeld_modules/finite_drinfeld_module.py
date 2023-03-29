r"""
Finite Drinfeld modules

This module provides the class
:class:`sage.rings.function_fields.drinfeld_module.finite_drinfeld_module.FiniteDrinfeldModule`,
which inherits
:class:`sage.rings.function_fields.drinfeld_module.drinfeld_module.DrinfeldModule`.

AUTHORS:

- Antoine Leudière (2022-04)
"""

# *****************************************************************************
#        Copyright (C) 2022 Antoine Leudière <antoine.leudiere@inria.fr>
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#                   http://www.gnu.org/licenses/
# *****************************************************************************

from sage.matrix.constructor import Matrix
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.special import companion_matrix
from sage.modules.free_module_element import vector
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.function_field.drinfeld_modules.drinfeld_module import DrinfeldModule
from sage.functions.other import ceil, sqrt
from sage.misc.misc_c import prod


class FiniteDrinfeldModule(DrinfeldModule):
    r"""
    This class implements finite Drinfeld `\mathbb{F}_q[T]`-modules.

    A *finite Drinfeld module* is a Drinfeld module whose base field is
    finite. In this case, the function field characteristic is a prime
    ideal.

    For general definitions and help on Drinfeld modules, see class
    :class:`sage.rings.function_fields.drinfeld_module.drinfeld_module.DrinfeldModule`.

    .. RUBRIC:: Construction:

    The user does not ever need to directly call
    ``FiniteDrinfeldModule`` --- the metaclass ``DrinfeldModule`` is
    responsible for instantiating ``DrinfeldModule`` or
    ``FiniteDrinfeldModule`` depending on the input::

        sage: Fq = GF(343)
        sage: A.<T> = Fq[]
        sage: K.<z6> = Fq.extension(2)
        sage: phi = DrinfeldModule(A, [z6, 0, 5])
        sage: phi
        Drinfeld module defined by T |--> 5*t^2 + z6

    ::

        sage: isinstance(phi, DrinfeldModule)
        True
        sage: from sage.rings.function_field.drinfeld_modules.finite_drinfeld_module import FiniteDrinfeldModule
        sage: isinstance(phi, FiniteDrinfeldModule)
        True

    The user should never use ``FiniteDrinfeldModule`` to test if a
    Drinfeld module is finite, but rather the ``is_finite`` method::

        sage: phi.is_finite()
        True

    .. RUBRIC:: Complex multiplication of rank two finite Drinfeld modules

    We can handle some aspects of the theory of complex multiplication
    of finite Drinfeld modules. Apart from the method
    ``frobenius_endomorphism``, we only handle rank two Drinfeld
    modules.

    First of all, it is easy to create the Frobenius endomorphism::

        sage: frobenius_endomorphism = phi.frobenius_endomorphism()
        sage: frobenius_endomorphism
        Endomorphism of Drinfeld module defined by T |--> 5*t^2 + z6
          Defn: t^2

    Its characteristic polynomial can be computed::

        sage: chi = phi.frobenius_charpoly()
        sage: chi
        X^2 + (T + 2*z3^2 + 2*z3 + 1)*X + 2*T^2 + (z3^2 + z3 + 4)*T + 2*z3
        sage: frob_pol = frobenius_endomorphism.ore_polynomial()
        sage: chi(frob_pol, phi(T))
        0

    as well as its trace and norm::

        sage: phi.frobenius_trace()
        6*T + 5*z3^2 + 5*z3 + 6
        sage: phi.frobenius_trace() == -chi[1]
        True
        sage: phi.frobenius_norm()
        2*T^2 + (z3^2 + z3 + 4)*T + 2*z3

    We can decide if a Drinfeld module is ordinary or supersingular::

        sage: phi.is_ordinary()
        True
        sage: phi.is_supersingular()
        False

    .. RUBRIC:: Inverting the Drinfeld module

    The morphism that defines a Drinfeld module is injective (see
    [Gos1998]_, cor. 4.5.2). If the Drinfeld module is finite, one can
    retrieve preimages:

        sage: a = A.random_element()
        sage: phi.invert(phi(a)) == a
        True
    """

    def __init__(self, gen, category):
        """
        Initialize `self`.

        Validity of the input is checked in `__classcall_private__`. The
        `__init__` just saves attributes.

        INPUT:

        - ``function_ring`` -- a univariate polynomial ring whose base
          is a finite field

        - ``gen`` -- the generator of the Drinfeld module as a list of
          coefficients or an Ore polynomial

        - ``name`` (default: `'t'`) -- the name of the Ore polynomial
          ring gen

        TESTS::

            sage: Fq = GF(25)
            sage: A.<T> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: gen = [p_root, z12^3, z12^5]
            sage: phi = DrinfeldModule(A, gen)
            sage: ore_polring = phi.ore_polring()
            sage: phi._gen == ore_polring(gen)
            True
        """
        # NOTE: There used to be no __init__ here (which was fine). I
        # added one to ensure that FiniteDrinfeldModule would always
        # have _frobenius_norm and _frobenius_trace attributes.
        super().__init__(gen, category)
        self._frobenius_charpoly = None

    def frobenius_endomorphism(self):
        r"""
        Return the Frobenius endomorphism of the Drinfeld module as a
        morphism object.

        Let `q` be the order of the base field of the function ring. The
        *Frobenius endomorphism* is defined as the endomorphism whose
        defining Ore polynomial is `t^q`.

        EXAMPLES::

            sage: Fq = GF(343)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(2)
            sage: phi = DrinfeldModule(A, [1, 0, z6])
            sage: phi.frobenius_endomorphism()
            Endomorphism of Drinfeld module defined by T |--> z6*t^2 + 1
              Defn: t^2

       TESTS::

            sage: from sage.rings.function_field.drinfeld_modules.morphism import DrinfeldModuleMorphism
            sage: isinstance(phi.frobenius_endomorphism(), DrinfeldModuleMorphism)
            True
        """
        t = self.ore_polring().gen()
        Fq = self._function_ring.base()
        deg = self._base.over(Fq).degree(Fq)
        return self._Hom_(self, category=self.category())(t**deg)

    def frobenius_charpoly(self, var='X', algorithm='crystalline'):
        r"""
        Return the characteristic polynomial of the Frobenius
        endomorphism.

        Let `\mathbb{F}_q` be the base field of the function ring. The
        *characteristic polynomial* `\chi` *of the Frobenius endomorphism*
        is defined in [Gek1991]_. An important feature of this
        polynomial is that it is monic, univariate, and has coefficients
        in the function ring. As in our case the function
        ring is a univariate polynomial ring, it is customary to see the
        characteristic polynomial of the Frobenius endomorphism as a
        bivariate polynomial.

        Let `\chi = X^r + \sum_{i=1}^{r} A_{i}(T)X^{r-i}` be the
        characteristic polynomial of the Frobenius endomorphism, and
        let `t^n` be the Ore polynomial that defines the Frobenius
        endomorphism of `\phi`; by definition, `n` is the degree of `K`
        over the base field `\mathbb{F}_q`. Then we have

        .. MATH::

            \chi(t^n)(\phi(T))
            = t^{nr} + \sum_{i=1}^{r} \phi_{A_{i}}t^{n(r-i)}
            = 0,

        with `\deg(A_i) \leq \frac{ir}{n}`.

        Note that the *Frobenius trace* is defined as `A_1(T)` and the
        *Frobenius norm* is defined as `A_r(T)`.

        INPUT:

        - ``var`` (default: ``'X'``) -- the name of the second variable
        - ``algorithm`` (default: ``'crystalline'``) -- the algorithm
          used to compute the characteristic polynomial

        EXAMPLES::

            sage: Fq = GF(25)
            sage: A.<T> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(A, [p_root, z12^3, z12^5])
            sage: phi.frobenius_charpoly()
            X^2 + ((4*z2 + 4)*T^3 + (z2 + 3)*T^2 + 3*T + 2*z2 + 3)*X + 3*z2*T^6 + (4*z2 + 3)*T^5 + (4*z2 + 4)*T^4 + 2*T^3 + (3*z2 + 3)*T^2 + (z2 + 2)*T + 4*z2

            sage: Fq = GF(343)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(2)
            sage: phi = DrinfeldModule(A, [1, 0, z6])
            sage: chi = phi.frobenius_charpoly()
            sage: chi
            X^2 + ((3*z3^2 + z3 + 4)*T + 4*z3^2 + 6*z3 + 3)*X + (5*z3^2 + 2*z3)*T^2 + (4*z3^2 + 3*z3)*T + 5*z3^2 + 2*z3

        ::

            sage: frob_pol = phi.frobenius_endomorphism().ore_polynomial()
            sage: chi(frob_pol, phi(T))
            0

        ::

            sage: phi.frobenius_charpoly(algorithm="NotImplemented")
            Traceback (most recent call last):
            NotImplementedError: Algorithm "NotImplemented" not implemented

        ALGORITHM:

        By default, this method uses the so-called *crystalline*
        algorithm which computes the characteristic polynomial of the
        Frobenius acting on the crystalline cohomology of the Drinfeld
        module. For further details, see [Ang1997]_. Currently, the only
        alternative is to use the *gekeler* approach based on solving
        the linear system given by `t^{nr} + \sum_{i=1}^{r}
        \phi_{A_{i}}t^{n(r-i)} = 0`. For more details, see [Gek2008]_.
        """
        methods = [method for method in dir(self)
                   if method.startswith('_frobenius_charpoly_')]
        method_name = f'_frobenius_charpoly_{algorithm}'
        if method_name in methods:
            if self._frobenius_charpoly is not None:
                return self._frobenius_charpoly
            self._frobenius_charpoly = getattr(self, method_name)(var)
            return self._frobenius_charpoly
        raise NotImplementedError(f'Algorithm \"{algorithm}\" not implemented')

    def _frobenius_charpoly_crystalline(self, var='X'):
        r"""
        Return the characteristic polynomial of the Frobenius
        endomorphism using Crystalline cohomology.

        The algorithm works for Drinfeld `\mathbb{F}_q[T]`-modules of
        any rank.

        This method is private and should not be directly called.
        Instead, use :meth:`frobenius_charpoly`.

        INPUT:

        - ``var`` (default: ``'X'``) -- the name of the second variable

        EXAMPLES::

            sage: Fq = GF(25)
            sage: A.<T> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(A, [p_root, z12^3, z12^5])
            sage: phi._frobenius_charpoly_crystalline()
            X^2 + ((4*z2 + 4)*T^3 + (z2 + 3)*T^2 + 3*T + 2*z2 + 3)*X + 3*z2*T^6 + (4*z2 + 3)*T^5 + (4*z2 + 4)*T^4 + 2*T^3 + (3*z2 + 3)*T^2 + (z2 + 2)*T + 4*z2

        ::

            sage: Fq = GF(25)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(8)
            sage: phi = DrinfeldModule(A, [z, 4, 1, z, z+1, 2, z+2, 1, 1, 3, 1])
            sage: phi._frobenius_charpoly_crystalline()
            X^10 + X^9 + (3*T + z2 + 1)*X^8 + (4*T^2 + z2*T + 2*z2 + 1)*X^7 + ... + (4*z2 + 4)*T^4 + 4*z2*T^2 + (z2 + 2)*T + z2

        ::

            sage: Fq = GF(27)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(10)
            sage: phi = DrinfeldModule(A, [z, z^2 + z, 2, 1, z, z+1, 2, z+2, 0, 1, 1, z^2 + z])
            sage: phi._frobenius_charpoly_crystalline()
            X^11 + (z3^2 + 2*z3)*X^10 + ((z3 + 1)*T + z3)*X^9 + ((2*z3^2 + z3 + 2)*T^2 + ... + (2*z3^2 + 2*z3 + 2)*T + z3^2

        ALGORITHM:

        Compute the characteristic polynomial of the Frobenius endomorphism
        acting on the crystalline cohomology of a Drinfeld module, which
        is equal to that of the Frobenius endomorphism on the Drinfeld
        module. A recurrence on elements of the cohomology allows us to
        compute a matrix representation of the Frobenius endomorphism
        efficiently using a companion matrix method.
        """
        Fq = self._Fq
        K = self.base_over_constants_field()
        A = self.function_ring()
        q, r, n = Fq.cardinality(), self.rank(), K.degree(Fq)
        nstar = ceil(sqrt(n))
        nquo, nrem = divmod(n, nstar)
        drin_coeff = self.coefficients(sparse=False)
        poly_K = PolynomialRing(K, name=str(A.gen()))
        matrix_poly_K = MatrixSpace(poly_K, r, r)
        mu_coeffs = ((poly_K.gen() - drin_coeff[0])**(n+1)).coefficients(
                                                            sparse=False)

        def companion(order):
            # + [1] is required to satisfy formatting for companion matrix
            M = matrix_poly_K(companion_matrix(
                              [(drin_coeff[i]/drin_coeff[r])**(q**order)
                               for i in range(r)] + [1], format='top'))
            M[0, r-1] += poly_K.gen() / drin_coeff[r]**(q**order)
            return M

        companion_initial = prod([companion(i)
                                 for i in range(nrem, 0, -1)])
        companion_step = prod([companion(i)
                              for i in range(nstar + nrem, nrem, -1)])
        reduced_companions = []
        for k in range(nquo - 1, 0, -1):
            M = Matrix(poly_K, r, r)
            for i, row in enumerate(companion_step):
                for j, entry in enumerate(row):
                    reduction = entry % poly_K([c**(q**(-k*nstar % n))
                                                for c in mu_coeffs])
                    M[i, j] = poly_K([c**(q**(k*nstar))
                                     for c in reduction.coefficients(
                                              sparse=False)])
            reduced_companions.append(M)
        charpoly_coeffs_K = (prod(reduced_companions) * companion_step *
                             companion_initial).charpoly(var).coefficients(
                             sparse=False)
        # The above line obtains a char poly with coefficients in L[T]
        # This maps them into A
        return PolynomialRing(A, name=var)(
                list(map(lambda coeff:
                     A(list(map(lambda x: K(x).vector()[0], coeff))),
                       charpoly_coeffs_K)))

    def _frobenius_charpoly_gekeler(self, var='X'):
        r"""
        Return the characteristic polynomial of the Frobenius
        endomorphism using Gekeler's algorithm.

        The algorithm works for Drinfeld `\mathbb{F}_q[T]`-modules of
        any rank, provided that the constant coefficient is a generator
        of the base field.

        This method is private and should not be directly called.
        Instead, use :meth:`frobenius_charpoly`.

        .. WARNING:

            This algorithm only works in the generic case
            when the corresponding linear system is invertible.
            Notable cases where this fails include Drinfeld
            modules whose minimal polynomial is not equal to
            the characteristic polynomial, and rank 2 Drinfeld
            modules where the degree 1 coefficient of `\phi_T` is 0.

        INPUT:

        - ``var`` (default: ``'X'``) -- the name of the second variable

        OUTPUT: a univariate polynomial with coefficients in the
                function ring

        EXAMPLES::

            sage: Fq = GF(25)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(6)
            sage: phi = DrinfeldModule(A, [z, 4, 1, z])
            sage: phi._frobenius_charpoly_gekeler()
            X^3 + ((z2 + 2)*T^2 + (z2 + 2)*T + 4*z2 + 4)*X^2 + ... + (3*z2 + 2)*T^2 + (3*z2 + 3)*T + 4

            sage: Fq = GF(125)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(2)
            sage: phi = DrinfeldModule(A, [z, 0, z])
            sage: phi._frobenius_charpoly_gekeler()
            Traceback (most recent call last):
            NotImplementedError: 'Gekeler' algorithm failed

        ALGORITHM:

        Construct a linear system based on the requirement that the
        Frobenius satisfies a degree r polynomial with coefficients in
        the function ring. This generalizes the procedure from
        [Gek2008]_ for the rank 2 case.
        """
        Fq = self._Fq
        K = self.base_over_constants_field()
        A = self.function_ring()
        r, n = self.rank(), K.degree(Fq)
        # Compute constants that determine the block structure of the
        # linear system. The system is prepared such that the solution
        # vector has the form [a_0, a_1, ... a_{r-1}]^T with each a_i
        # corresponding to a block of length (n*(r - i))//r + 1
        shifts = [(n*(r - i))//r + 1 for i in range(r)]
        rows, cols = n*r + 1, sum(shifts)
        block_shifts = [0]
        for i in range(r-1):
            block_shifts.append(block_shifts[-1] + shifts[i])
        # Compute the images \phi_T^i for i = 0 .. n.
        gen_powers = [self(A.gen()**i).coefficients(sparse=False)
                      for i in range(0, n + 1)]
        sys, vec = Matrix(K, rows, cols), vector(K, rows)
        vec[rows - 1] = -1
        for j in range(r):
            for k in range(shifts[j]):
                for i in range(len(gen_powers[k])):
                    sys[i + n*j, block_shifts[j] + k] = gen_powers[k][i]
        if sys.right_nullity() != 0:
            raise NotImplementedError("'Gekeler' algorithm failed")
        sol = list(sys.solve_right(vec))
        # The system is solved over L, but the coefficients should all be in Fq
        # We project back into Fq here.
        sol_Fq = list(map(lambda x: K(x).vector()[0], sol))
        char_poly = []
        for i in range(r):
            char_poly.append([sol_Fq[block_shifts[i] + j]
                              for j in range(shifts[i])])
        return PolynomialRing(self._function_ring, name=var)(char_poly + [1])

    def frobenius_norm(self):
        r"""
        Return the Frobenius norm of the Drinfeld module.

        Let `C(X) = \sum_{i=0}^r a_iX^{r-i}` denote the characteristic
        polynomial of the Frobenius endomorphism. The Frobenius norm
        is `(-1)^r a_{r}`. This is an element of the regular function ring
        and if `n` is the degree of the base field over `\mathbb{F}_q`,
        then the Frobenius norm has degree `n`.

        EXAMPLES::

            sage: Fq = GF(343)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(2)
            sage: phi = DrinfeldModule(A, [1, 0, z6])
            sage: B = phi.frobenius_norm()
            sage: B
            (5*z3^2 + 2*z3)*T^2 + (4*z3^2 + 3*z3)*T + 5*z3^2 + 2*z3

        ::

            sage: n = 2  # Degree of the base field over Fq
            sage: B.degree() == n
            True

        ::

            sage: B == phi.frobenius_charpoly()[0]
            True

        ALGORITHM:

            The Frobenius norm is computed using the formula, by
            Gekeler, given in [MS2019]_, Section 4, Proposition 3.
        """
        return ((-1)**(self.rank() % 2)) * \
               self.frobenius_charpoly().coefficients(sparse=False)[0]

    def frobenius_trace(self):
        r"""
        Return Frobenius trace of the Drinfeld module, if the rank is
        two; raise a NotImplementedError otherwise.

        Let `C(X) = \sum_{i=0}^r a_iX^{r-i}` denote the characteristic
        polynomial of the Frobenius endomorphism. The Frobenius trace
        is `-a_{1}`. This is an element of the regular function ring
        and if `n` is the degree of the base field over `\mathbb{F}_q`, 
        then the Frobenius trace has degree at most `\frac{n}{r}`.

        EXAMPLES::

            sage: Fq = GF(343)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(2)
            sage: phi = DrinfeldModule(A, [1, 0, z6])
            sage: A = phi.frobenius_trace()
            sage: A
            (4*z3^2 + 6*z3 + 3)*T + 3*z3^2 + z3 + 4

        ::

            sage: n = 2  # Degree over Fq of the base codomain
            sage: A.degree() <= n/2
            True

        ::

            sage: A == -phi.frobenius_charpoly()[1]
            True

        ALGORITHM:

            We extract the coefficient of `X^{r-1}` from the
            characteristic polynomial.
        """
        return -self.frobenius_charpoly().coefficients(sparse=False)[-2]

    def invert(self, ore_pol):
        r"""
        Return the preimage of the input under the Drinfeld module, if it
        exists.

        INPUT:

        - ``ore_pol`` -- the Ore polynomial whose preimage we want to
          compute

        EXAMPLES::

            sage: Fq = GF(25)
            sage: A.<T> = Fq[]
            sage: K.<z12> = Fq.extension(6)
            sage: p_root = 2*z12^11 + 2*z12^10 + z12^9 + 3*z12^8 + z12^7 + 2*z12^5 + 2*z12^4 + 3*z12^3 + z12^2 + 2*z12
            sage: phi = DrinfeldModule(A, [p_root, z12^3, z12^5])
            sage: a = A.random_element()
            sage: phi.invert(phi(a)) == a
            True
            sage: phi.invert(phi(T)) == T
            True
            sage: phi.invert(phi(Fq.gen())) == Fq.gen()
            True

        When the input is not in the image of the Drinfeld module, an
        exception is raised::

            sage: t = phi.ore_polring().gen()
            sage: phi.invert(t + 1)
            Traceback (most recent call last):
            ...
            ValueError: input must be in the image of the Drinfeld module

            sage: phi.invert(t^3 + t^2 + 1)
            Traceback (most recent call last):
            ...
            ValueError: input must be in the image of the Drinfeld module

        ALGORITHM:

            The algorithm relies on the inversion of a linear algebra
            system. See [MS2019]_, 3.2.5 for details.

        TESTS::

            sage: a = A.random_element()
            sage: cat = phi.category()
            sage: phi_r1 = cat.random_object(1)
            sage: phi_r1.invert(phi_r1(a)) == a
            True
            sage: phi_r2 = cat.random_object(2)
            sage: phi_r2.invert(phi_r2(a)) == a
            True
            sage: phi_r3 = cat.random_object(3)
            sage: phi_r3.invert(phi_r3(a)) == a
            True
            sage: phi_r4 = cat.random_object(4)
            sage: phi_r4.invert(phi_r4(a)) == a
            True
            sage: phi_r5 = cat.random_object(5)
            sage: phi_r5.invert(phi_r5(a)) == a
            True
        """
        deg = ore_pol.degree()
        r = self.rank()
        base_over_Fq = self.base_over_constants_field()
        if ore_pol not in self._ore_polring:
            raise TypeError('input must be an Ore polynomial')
        if ore_pol.degree() == 0:
            coord = base_over_Fq(self._base(ore_pol)).vector()
            if coord.nonzero_positions == [0]:
                return self._Fq(coord[0])
        if ore_pol == 0:
            return self._Fq.zero()
        if deg % r != 0:
            raise ValueError('input must be in the image of the Drinfeld '
                             'module')

        k = deg // r
        T = self._function_ring.gen()
        mat_lines = [[0 for _ in range(k+1)] for _ in range(k+1)]
        for i in range(k+1):
            phi_T_i = self(T**i)
            for j in range(i+1):
                mat_lines[j][i] = phi_T_i[r*j]
        mat = Matrix(mat_lines)
        vec = vector([list(ore_pol)[r*j] for j in range(k+1)])
        coeffs_K = list((mat**(-1)) * vec)
        coeffs_Fq = list(map(lambda x: base_over_Fq(x).vector()[0], coeffs_K))
        pre_image = self._function_ring(coeffs_Fq)

        if self(pre_image) == ore_pol:
            return pre_image
        else:
            raise ValueError('input must be in the image of the Drinfeld '
                             'module')

    def is_isogenous(self, psi):
        r"""
        Return ``True`` when ``self`` is isogenous to the other
        Drinfeld module.

        If the Drinfeld modules do not belong to the same category, an
        exception is raised.

        EXAMPLES::

            sage: Fq = GF(2)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(3)
            sage: psi = DrinfeldModule(A, [z, z + 1, z^2 + z + 1])
            sage: phi = DrinfeldModule(A, [z, z^2 + z + 1, z^2 + z])
            sage: phi.is_isogenous(psi)
            True

        ::

            sage: chi = DrinfeldModule(A, [z, z + 1, z^2 + z])
            sage: phi.is_isogenous(chi)
            False

        ::

            sage: mu = DrinfeldModule(A, [z + 1, z^2 + z + 1, z^2 + z])
            sage: phi.is_isogenous(mu)
            Traceback (most recent call last):
            TypeError: Drinfeld modules are not in the same category

        ::

            sage: mu = 1
            sage: phi.is_isogenous(mu)
            Traceback (most recent call last):
            TypeError: Input must be a Drinfeld module

        ALGORITHM:

        Two Drinfeld A-modules of equal characteristic are isogenous
        if and only if:

        - they have the same rank
        - the characteristic polynomial of the Frobenius endomorphism
          for both Drinfeld modules are equal.
        """
        if not isinstance(psi, DrinfeldModule):
            raise TypeError("Input must be a Drinfeld module")
        if self.category() != psi.category():
            raise TypeError("Drinfeld modules are not in the same category")
        return self.rank() == psi.rank() and \
               self.frobenius_charpoly() == psi.frobenius_charpoly()

    def is_ordinary(self):
        r"""
        Return ``True`` whether the Drinfeld module is ordinary; raise a
        NotImplementedError if the rank is not two.

        A rank two finite Drinfeld module is *ordinary* if and only if
        the function ring-characteristic does not divide the Frobenius
        trace. A *supersingular* rank two finite Drinfeld module is a
        Drinfeld module that is not ordinary.

        A rank two Drinfeld module is *ordinary* if and only if it is
        not supersingular; see :meth:`is_supersingular`.

        EXAMPLES::

            sage: Fq = GF(343)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(2)
            sage: phi = DrinfeldModule(A, [1, 0, z6])
            sage: phi.is_ordinary()
            False
            sage: phi_p = phi(phi.characteristic())
            sage: phi_p  # Purely inseparable
            z6*t^2

        ALGORITHM:

            Compute the Frobenius trace and test if the
            `\mathbb{F}_q[T]` characteristic divides it.

            We could also test if the image of the
            `\mathbb{F}_q[T]`-characteristic under the Drinfeld module
            is purely inseparable; see [Gek1991]_, Proposition 4.1.
        """
        return not self.is_supersingular()

    def is_supersingular(self):
        r"""
        Return ``True`` whether the Drinfeld module is supersingular; raise a
        NotImplementedError if the rank is not two.

        A rank two finite Drinfeld module is *supersingular* if and only
        if the function ring-characteristic divides the Frobenius
        trace. An *ordinary* rank two finite Drinfeld module is a
        Drinfeld module that is not supersingular.

        EXAMPLES::

            sage: Fq = GF(343)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(2)
            sage: phi = DrinfeldModule(A, [1, 0, z6])
            sage: phi.is_supersingular()
            True
            sage: phi(phi.characteristic()) # Purely inseparable
            z6*t^2

        ALGORITHM:

            Compute the Frobenius trace and test if the function
            ring-characteristic divides it.

            We could also test if the image of the function
            ring-characteristic under the Drinfeld module is purely
            inseparable; see [Gek1991]_, Proposition 4.1.
        """
        return self.characteristic().divides(self.frobenius_trace())

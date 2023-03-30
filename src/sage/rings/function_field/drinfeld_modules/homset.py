r"""
Set of morphisms between two Drinfeld modules

This module provides the class
:class:`sage.rings.function_field.drinfeld_module.homset.DrinfeldModuleHomset`.

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

from sage.categories.drinfeld_modules import DrinfeldModules
from sage.categories.homset import Homset
from sage.misc.latex import latex
from sage.rings.function_field.drinfeld_modules.morphism import DrinfeldModuleMorphism
from sage.structure.parent import Parent
from sage.functions.log import logb
from sage.matrix.constructor import Matrix
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing


class DrinfeldModuleHomset(Homset):
    r"""
    This class implements the set of morphisms between two Drinfeld
    `\mathbb{F}_q[T]`-modules.

    INPUT:

    - ``X`` -- the domain

    - ``Y`` -- the codomain

    EXAMPLES::

        sage: Fq = GF(27)
        sage: A.<T> = Fq[]
        sage: K.<z6> = Fq.extension(2)
        sage: phi = DrinfeldModule(A, [z6, z6, 2])
        sage: psi = DrinfeldModule(A, [z6, 2*z6^5 + 2*z6^4 + 2*z6 + 1, 2])
        sage: H = Hom(phi, psi)
        sage: H
        Set of Drinfeld module morphisms from (gen) 2*t^2 + z6*t + z6 to (gen) 2*t^2 + (2*z6^5 + 2*z6^4 + 2*z6 + 1)*t + z6

    ::

        sage: from sage.rings.function_field.drinfeld_modules.homset import DrinfeldModuleHomset
        sage: isinstance(H, DrinfeldModuleHomset)
        True

    There is a simpler syntax for endomorphisms sets::

        sage: E = End(phi)
        sage: E
        Set of Drinfeld module morphisms from (gen) 2*t^2 + z6*t + z6 to (gen) 2*t^2 + z6*t + z6
        sage: E is Hom(phi, phi)
        True

    The domain and codomain must have the same Drinfeld modules
    category::

        sage: rho = DrinfeldModule(A, [Frac(A)(T), 1])
        sage: Hom(phi, rho)
        Traceback (most recent call last):
        ...
        ValueError: Drinfeld modules must be in the same category

    ::

        sage: sigma = DrinfeldModule(A, [1, z6, 2])
        sage: Hom(phi, sigma)
        Traceback (most recent call last):
        ...
        ValueError: Drinfeld modules must be in the same category

    One can create morphism objects by calling the homset::

        sage: identity_morphism = E(1)
        sage: identity_morphism
        Identity morphism of Drinfeld module defined by T |--> 2*t^2 + z6*t + z6

    ::

        sage: t = phi.ore_polring().gen()
        sage: frobenius_endomorphism = E(t^6)
        sage: frobenius_endomorphism
        Endomorphism of Drinfeld module defined by T |--> 2*t^2 + z6*t + z6
          Defn: t^6

    ::

        sage: isogeny = H(t + 1)
        sage: isogeny
        Drinfeld Module morphism:
          From: Drinfeld module defined by T |--> 2*t^2 + z6*t + z6
          To:   Drinfeld module defined by T |--> 2*t^2 + (2*z6^5 + 2*z6^4 + 2*z6 + 1)*t + z6
          Defn: t + 1

    And one can test if an Ore polynomial defines a morphism using the
    ``in`` syntax::

        sage: 1 in H
        False
        sage: t^6 in H
        False
        sage: t + 1 in H
        True
        sage: 1 in E
        True
        sage: t^6 in E
        True
        sage: t + 1 in E
        False

    This also works if the candidate is a morphism object::

        sage: isogeny in H
        True
        sage: E(0) in E
        True
        sage: identity_morphism in H
        False
        sage: frobenius_endomorphism in H
        False
    """

    Element = DrinfeldModuleMorphism
    __contains__ = Parent.__contains__

    def __init__(self, X, Y, category=None, check=True):
        """
        Initialize ``self``.

        INPUT:

        - ``X`` -- the domain of the homset

        - ``Y`` -- the codomain of the homset

        - ``category`` (default: ``None``) -- the Drinfeld modules category of
          the domain and codomain

        - ``check`` (default: ``True``) -- check the validity of the category

        TESTS::

            sage: Fq = GF(27)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(2)
            sage: phi = DrinfeldModule(A, [z6, z6, 2])
            sage: psi = DrinfeldModule(A, [z6, 2*z6^5 + 2*z6^4 + 2*z6 + 1, 2])
            sage: H = Hom(phi, psi)
            sage: H.domain() is phi
            True
            sage: H.codomain() is psi
            True
        """
        if category is None:
            category = X.category()
        if check:
            if X.category() != Y.category() \
                    or not isinstance(X.category(), DrinfeldModules):
                raise ValueError('Drinfeld modules must be in the same '
                                 'category')
            if category != X.category():
                raise ValueError('category should be DrinfeldModules')
        base = category.base()
        super().__init__(X, Y, category=category, base=base, check=check)

    def _latex_(self):
        r"""
        Return a LaTeX representation of the homset.

        EXAMPLES::

            sage: Fq = GF(27)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(2)
            sage: phi = DrinfeldModule(A, [z6, z6, 2])
            sage: psi = DrinfeldModule(A, [z6, 2*z6^5 + 2*z6^4 + 2*z6 + 1, 2])
            sage: H = Hom(phi, psi)
            sage: latex(H)
            \text{Set{ }of{ }Drinfeld{ }module{ }morphisms{ }from{ }(gen){ }}2 t^{2} + z_{6} t + z_{6}\text{{ }to{ }(gen){ }}2 t^{2} + \left(2 z_{6}^{5} + 2 z_{6}^{4} + 2 z_{6} + 1\right) t + z_{6}
        """
        return f'\\text{{Set{{ }}of{{ }}Drinfeld{{ }}module{{ }}morphisms' \
               f'{{ }}from{{ }}(gen){{ }}}}{latex(self.domain().gen())}' \
               f'\\text{{{{ }}to{{ }}(gen){{ }}}}'\
               f'{latex(self.codomain().gen())}'

    def _repr_(self):
        r"""
        Return a string representation of the homset.

        EXAMPLES::

            sage: Fq = GF(27)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(2)
            sage: phi = DrinfeldModule(A, [z6, z6, 2])
            sage: psi = DrinfeldModule(A, [z6, 2*z6^5 + 2*z6^4 + 2*z6 + 1, 2])
            sage: H = Hom(phi, psi)
            sage: H
            Set of Drinfeld module morphisms from (gen) 2*t^2 + z6*t + z6 to (gen) 2*t^2 + (2*z6^5 + 2*z6^4 + 2*z6 + 1)*t + z6
        """
        return f'Set of Drinfeld module morphisms from (gen) '\
               f'{self.domain().gen()} to (gen) {self.codomain().gen()}'

    def __contains__(self, x):
        r"""
        Return ``True`` if the input defines a morphism in the homset.

        INPUT:

        - ``x`` -- an Ore polynomial or a Drinfeld module morphism

        EXAMPLES:

        In the next examples, the input is an Ore polynomial::

            sage: Fq = GF(27)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(2)
            sage: phi = DrinfeldModule(A, [z6, z6, 2])
            sage: psi = DrinfeldModule(A, [z6, 2*z6^5 + 2*z6^4 + 2*z6 + 1, 2])
            sage: H = Hom(phi, psi)
            sage: E = End(phi)
            sage: t = phi.ore_polring().gen()
            sage: 1 in H
            False
            sage: t^6 in H
            False
            sage: t + 1 in H
            True
            sage: 1 in E
            True
            sage: t^6 in E
            True
            sage: t + 1 in E
            False

        Whereas the input is now a Drinfeld module morphism::

            sage: isogeny = H(t + 1)
            sage: isogeny in H
            True
            sage: E(0) in E
            True
            sage: identity_morphism = E(1)
            sage: identity_morphism in H
            False
            sage: frobenius_endomorphism = phi.frobenius_endomorphism()
            sage: frobenius_endomorphism in H
            False
        """
        try:
            x = self(x)
            return True
        except (AttributeError, ValueError, TypeError):
            return False

    def _element_constructor_(self, *args, **kwds):
        r"""
        Return the Drinfeld module morphism defined by the given Ore
        polynomial.

        INPUT: an Ore polynomial

        EXAMPLES::

            sage: Fq = GF(27)
            sage: A.<T> = Fq[]
            sage: K.<z6> = Fq.extension(2)
            sage: phi = DrinfeldModule(A, [z6, z6, 2])
            sage: psi = DrinfeldModule(A, [z6, 2*z6^5 + 2*z6^4 + 2*z6 + 1, 2])
            sage: H = Hom(phi, psi)
            sage: E = End(phi)
            sage: t = phi.ore_polring().gen()
            sage: identity_morphism = E(1)
            sage: identity_morphism
            Identity morphism of Drinfeld module defined by T |--> 2*t^2 + z6*t + z6

        ::

            sage: frobenius_endomorphism = E(t^6)
            sage: frobenius_endomorphism
            Endomorphism of Drinfeld module defined by T |--> 2*t^2 + z6*t + z6
              Defn: t^6

        ::

            sage: isogeny = H(t + 1)
            sage: isogeny
            Drinfeld Module morphism:
              From: Drinfeld module defined by T |--> 2*t^2 + z6*t + z6
              To:   Drinfeld module defined by T |--> 2*t^2 + (2*z6^5 + 2*z6^4 + 2*z6 + 1)*t + z6
              Defn: t + 1
        """
        # NOTE: This used to be self.element_class(self, ...), but this
        # would call __init__ instead of __classcall_private__. This
        # seems to work, but I don't know what I'm doing.
        return DrinfeldModuleMorphism(self, *args, **kwds)

    def basis(self, degree):
        r"""
        Return a basis for the `\mathbb{F}_q`-space of morphisms from self to
        a Drinfeld module `\psi` of degree at most `degree`. A morphism
        `\iota: \phi \to psi` is an element `\iota \in K\{\tau\}` such that
        `iota \phi_T = \psi_T \iota`. The degree of a morphism is the
        `\tau`-degree of `\iota`.

        INPUT:

        - ``'\psi'`` -- the morphism's codomain, a Drinfeld module

        - ``degree`` (default: ``'1'``) -- the maximum degree of the morphism

        OUTPUT: a list of univariate ore polynomials with coefficients in `K`

        EXAMPLES::

            sage: Fq = GF(2)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(3)
            sage: psi = DrinfeldModule(A, [z, z + 1, z^2 + z + 1])
            sage: phi = DrinfeldModule(A, [z, z^2 + z + 1, z^2 + z])
            sage: homset = Hom(phi, psi)
            sage: basis = homset.basis(3)
            sage: iso = basis[0]
            sage: iso*phi.gen() - psi.gen()*iso
            0

        ALGORITHM:

            We return the basis of the kernel of a matrix derived from the
            constraint that `\iota \phi_T = \psi_T \iota`. See [Wes2022]_ for
            details on this algorithm.
        """
        domain, codomain = self.domain(), self.codomain()
        Fq = domain._Fq
        K = domain.base_over_constants_field()
        q = Fq.cardinality()
        char = Fq.characteristic()
        r = domain.rank()
        n = K.degree(Fq)
        # shorten name for readability
        d = degree
        qorder = logb(q, char)
        K_basis = K.basis_over(Fq)
        # This weird casting is done to allow extraction
        # of coefficients in non-sparse mode from field elements
        # There really should be a cleaner way to do this
        pol_var = K_basis[0].polynomial().parent().gen()
        pol_ring = PolynomialRing(Fq, str(pol_var))
        dom_coeffs = domain.coefficients(sparse=False)
        cod_coeffs = codomain.coefficients(sparse=False)
        frob_matrices = []

        # Compute matrices for the Frobenius operator for
        # K over Fq
        for order in range(d + r):
            frob = Matrix(Fq, n, n)
            for i, elem in enumerate(K_basis):
                col = pol_ring(elem.frobenius(qorder*order).polynomial()) \
                      .coefficients(sparse=False)
                col += [0 for _ in range(n - len(col))]
                for j in range(n):
                    frob[j, i] = col[j]
            frob_matrices.append(frob)

        sys = Matrix(Fq, (d + r + 1)*n, (d + 1)*n)
        for k in range(0, d + r + 1):
            for i in range(max(0, k - r), min(k, d) + 1):
                # We represent multiplication and Frobenius
                # as operators acting on K as a vector space
                # over Fq
                # Require matrices act on the right, so we
                # take a transpose of operators here
                oper = K(dom_coeffs[k-i]
                       .frobenius(qorder*i)).matrix().transpose() \
                       - K(cod_coeffs[k-i]).matrix().transpose() \
                       * frob_matrices[k - i]
                for j in range(n):
                    for l in range(n):
                        sys[k*n + j, i*n + l] = oper[j, l]
        sol = sys.right_kernel().basis()
        # Reconstruct the Ore polynomial from the coefficients
        basis = []
        tau = domain.ore_polring().gen()
        for basis_elem in sol:
            basis.append(sum([sum([K_basis[j]*basis_elem[i*n + j]
                               for j in range(n)])*(tau**i)
                               for i in range(d + 1)]))
        return basis

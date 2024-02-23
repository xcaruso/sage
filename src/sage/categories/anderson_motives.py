r"""
Anderson motives

This module provides the class
:class:`sage.category.drinfeld_modules.AndersonMotives`.

AUTHORS:

- Xavier Caruso (2024-03)
"""

# *************************************************************************
#   Copyright (C) 2024      Xavier Caruso <xavier.caruso@normalesup.org>
#
#   Distributed under the terms of the GNU General Public License (GPL)
#                   http://www.gnu.org/licenses/
# *************************************************************************

from sage.misc.latex import latex

from sage.categories.objects import Objects
from sage.categories.category_types import Category_over_base_ring
from sage.categories.homsets import Homsets
from sage.categories.drinfeld_modules import DrinfeldModules

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing


class AndersonMotives(Category_over_base_ring):
    def __init__(self, category):
        if not isinstance(category, DrinfeldModules):
            category = DrinfeldModules(category)
        self._base_field = K = category.base()
        self._base_morphism = category.base_morphism()
        self._function_ring = A = category.function_ring()
        self._base_over_constants_field = Kp = category.base_over_constants_field()
        self._characteristic = category._characteristic
        self._base_combined = AK = PolynomialRing(Kp, A.variable_name())
        self._constant_coefficient = category.constant_coefficient()
        self._divisor = AK.gen() - Kp(self._constant_coefficient)
        super().__init__(base=K)

    def _latex_(self):
        return f'\\text{{Category{{ }}of{{ }}Anderson{{ }}motives{{ }}' \
               f'over{{ }}{latex(self._base_field)}'

    def _repr_(self):
        return f'Category of Anderson motives over {self._base_field}'

    def Homsets(self):
        return Homsets()

    def Endsets(self):
        return Homsets().Endsets()

    def base_morphism(self):
        return self._base_morphism

    def base_over_constants_field(self):
        return self._base_over_constants_field

    def base_combined(self):
        return self._base_combined

    def divisor(self):
        return self._divisor

    def characteristic(self):
        if self._characteristic is None:
            raise NotImplementedError('function ring characteristic not '
                                      'implemented in this case')
        return self._characteristic

    def constant_coefficient(self):
        return self._constant_coefficient

    def function_ring(self):
        return self._function_ring

    def object(self, gen):
        raise NotImplementedError

    def super_categories(self):
        """
        EXAMPLES::

            sage: Fq = GF(11)
            sage: A.<T> = Fq[]
            sage: K.<z> = Fq.extension(4)
            sage: p_root = z^3 + 7*z^2 + 6*z + 10
            sage: phi = DrinfeldModule(A, [p_root, 0, 0, 1])
            sage: C = phi.category()
            sage: C.super_categories()
            [Category of objects]
        """
        return [Objects()]

    class ParentMethods:

        def base(self):
            return self.category().base()

        def base_morphism(self):
            return self.category().base_morphism()

        def base_combined(self):
            return self.category().base_combined()

        def base_over_constants_field(self):
            return self.category().base_over_constants_field()

        def characteristic(self):
            return self.category().characteristic()

        def function_ring(self):
            return self.category().function_ring()

        def constant_coefficient(self):
            return self.category().constant_coefficient()

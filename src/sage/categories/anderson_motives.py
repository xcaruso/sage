# *****************************************************************************
#        Copyright (C) 2024 Xavier Caruso <xavier.caruso@normalesup.org>
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#                   http://www.gnu.org/licenses/
# *****************************************************************************


from sage.misc.latex import latex

from sage.categories.objects import Objects
from sage.categories.category_types import Category_over_base_ring
from sage.categories.homsets import Homsets
from sage.categories.drinfeld_modules import DrinfeldModules

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing


class AndersonMotives(Category_over_base_ring):
    @staticmethod
    def __classcall_private__(cls, category, dispatch=True):
        if isinstance(category, AndersonMotives):
            return category
        if isinstance(category, DrinfeldModules):
            category = DrinfeldModules(category.base())
        else:
            category = DrinfeldModules(category)
        return AndersonMotives.__classcall__(cls, category)

    def __init__(self, category):
        self._base_field = K = category.base()
        self._base_morphism = category.base_morphism()
        self._function_ring = A = category.function_ring()
        self._base_over_constants_field = category.base_over_constants_field()
        self._characteristic = category._characteristic
        self._base_combined = AK = PolynomialRing(K, A.variable_name())  # TODO: find a better name
        self._constant_coefficient = category.constant_coefficient()
        self._divisor = AK.gen() - self._constant_coefficient
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

# sage.doctest: optional - sage.rings.finite_rings
r"""
Set of morphisms between two Drinfeld modules

This module provides the class
:class:`sage.rings.function_field.drinfeld_module.homset.DrinfeldModuleHomset`.

AUTHORS:

- Antoine Leudière (2022-04)
"""

# *****************************************************************************
#        Copyright (C) 2024 Xavier Caruso <xavier.caruso@normalesup.org>
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#                   http://www.gnu.org/licenses/
# *****************************************************************************

from sage.categories.anderson_motives import AndersonMotives
from sage.categories.homset import Homset
from sage.rings.function_field.drinfeld_modules.anderson_motive_morphism import AndersonMotiveMorphism


class AndersonMotiveHomset(Homset):
    Element = AndersonMotiveMorphism

    def __init__(self, X, Y, category=None, check=True):
        if category is None:
            category = X.category()
        if check:
            if X.category() != Y.category() \
                    or not isinstance(X.category(), AndersonMotives):
                raise ValueError('Anderson motives must be in the same category')
            if category != X.category():
                raise ValueError('category should be AndersonMotives')
        base = category.base()
        super().__init__(X, Y, category=category, base=base, check=check)

    def _element_constructor_(self, *args, **kwds):
        return AndersonMotiveMorphism(self, *args, **kwds)

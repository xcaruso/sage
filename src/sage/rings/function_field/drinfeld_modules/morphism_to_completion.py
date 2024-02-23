# *****************************************************************************
#        Copyright (C) 2024 Xavier Caruso <xavier.caruso@normalesup.org>
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#                   http://www.gnu.org/licenses/
# *****************************************************************************

from sage.categories.homset import Hom

from sage.rings.morphism import RingHomomorphism
from sage.rings.integer_ring import ZZ
from sage.rings.infinity import Infinity
from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
from sage.rings.laurent_series_ring import LaurentSeriesRing
from sage.rings.function_field.function_field import FunctionField


def normalize_place(A, place, infty=True):
    if infty and place is Infinity:
        return place
    if place in A.base_ring():
        return A.gen() - place
    elif place in A:
        place = A(place)
        if place.degree() == 0:
            return A.gen() - place
        if place.is_irreducible():
            return place.monic()
    if infty:
        raise ValueError("place must be Infinity or an irreducible polynomial")
    else:
        raise ValueError("place must an irreducible polynomial")


class MorphismToCompletion(RingHomomorphism):
    def __init__(self, domain, place, default_prec, residue_name, uniformizer_name):
        if not isinstance(domain, PolynomialRing_general):
            raise NotImplementedError
        k = base = domain.base_ring()
        place = self._place = normalize_place(domain, place)
        if place is Infinity:
            codomain = LaurentSeriesRing(base, name=uniformizer_name, default_prec=default_prec)
            image = codomain.one() >> 1
        elif place.degree() == 1:
            codomain = LaurentSeriesRing(base, name=uniformizer_name, default_prec=default_prec)
            image = codomain.gen() - place[0]
        else:
            k = base.extension(place, name=residue_name)
            codomain = LaurentSeriesRing(k, name=uniformizer_name, default_prec=default_prec)
            image = codomain(k.gen()).add_bigoh(1)
        parent = Hom(domain, codomain)
        RingHomomorphism.__init__(self, parent)
        self._image = image
        self._k = k
        self._q = k.cardinality()

    def degree(self):
        if self._place is Infinity:
            return ZZ(1)
        else:
            return self._place.degree()

    def residue_field(self):
        return self._k

    def image_generator(self, prec):
        S = self.codomain()
        image = self._image
        current_prec = image.precision_absolute()
        if current_prec < prec:
            # We update the image
            u = S.gen()
            place = self._place
            der = place.derivative()
            while current_prec < prec:
                current_prec *= 2
                image = image.lift_to_precision(current_prec)
                image -= (place(image) - u) / der(image)
            self._image = image
        return image.add_bigoh(prec)

    def _repr_type(self):
        return "Completion"

    def place(self):
        return self._place

    def _call_with_args(self, P, args, kwds):
        if 'prec' in kwds:
            prec = kwds['prec']
        else:
            prec = self.codomain().default_prec()
        image = self.image_generator(prec)
        return P(image)

    def _call_(self, P):
        image = self._image
        if image.precision_absolute() is not Infinity:
            prec = self.codomain().default_prec()
            image = self.image_generator(prec)
        return P(image)

    def teichmuller_lift(self, x, prec):
        place = self._place
        q = self._q
        modulus = place
        Q = self._k(x).polynomial()
        if Q == 0:
            return Q
        Q = self.domain()(Q)
        current_prec = 1
        while current_prec < prec:
            modulus = modulus ** 2
            current_prec *= 2
            R = Q ** (q-2)
            G, D, _ = R.xgcd(modulus)
            if G != 1:
                raise RuntimeError
            Q = 2*Q - D
        return Q % (place**prec)

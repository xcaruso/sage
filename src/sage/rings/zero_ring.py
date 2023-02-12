from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.ring import CommutativeRing
from sage.structure.element import CommutativeRingElement
from sage.rings.integer_ring import ZZ

from sage.structure.richcmp import rich_to_bool


class ZeroRingElement(CommutativeRingElement):
    def __init__(self, parent, x=0):
        CommutativeRingElement.__init__(self, parent)

    def _repr_(self):
        return "0"

    def _latex_(self):
        return "0"

    def is_zero(self):
        return True

    def is_one(self):
        return True

    def _richcmp_(self, other, op):
        return rich_to_bool(op, 0)

    def _add_(self, other):
        return self

    def _sub_(self, other):
        return self

    def _mul_(self, other):
        return self

    def __invert__(self):
        return self

    def is_unit(self):
        return True

    def inverse_of_unit(self):
        return self


class ZeroRing(UniqueRepresentation, CommutativeRing):
    Element = ZeroRingElement

    def __init__(self, base=ZZ, category=None):
        CommutativeRing.__init__(self, base, category=category)

    def _repr_(self):
        return "The zero ring"

    def _latex_(self):
        return "0"

    def _coerce_map_from_(self, R):
        return True

    def is_zero(self):
        return True

    def zero(self):
        return ZeroRingElement(self)

    def one(self):
        return ZeroRingElement(self)

    def ngens(self):
        return ZZ(0)

    def gens(self):
        return tuple()

    def gen(self, i=0):
        raise IndexError

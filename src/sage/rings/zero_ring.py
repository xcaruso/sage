r"""
The zero ring.

AUTHOR:

- Xavier Caruso (2023-03)

"""


from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.ring import CommutativeRing
from sage.structure.element import CommutativeRingElement
from sage.rings.integer_ring import ZZ

from sage.structure.richcmp import rich_to_bool


class ZeroRingElement(CommutativeRingElement):
    r"""
    The unique element of the zero ring.

    TESTS::

        sage: from sage.rings.zero_ring import ZeroRing
        sage: Z = ZeroRing()
        sage: z = Z()
        sage: TestSuite(z).run()

    """
    def __init__(self, parent, x=0):
        r"""
        Initialize this element.

        TESTS::

            sage: from sage.rings.zero_ring import ZeroRing
            sage: Z = ZeroRing()
            sage: Z()
            0
            sage: Z(0)
            0
            sage: Z(1)
            0

        """
        CommutativeRingElement.__init__(self, parent)

    def _repr_(self):
        r"""
        Return a string representation of this element.

        EXAMPLES::

            sage: from sage.rings.zero_ring import ZeroRing
            sage: z = ZeroRing()()
            sage: str(z)
            '0'

        """
        return "0"

    def _latex_(self):
        r"""
        Return a LaTeX representation of this element.

        EXAMPLES::

            sage: from sage.rings.zero_ring import ZeroRing
            sage: z = ZeroRing()()
            sage: latex(z)
            0

        """
        return "0"

    def is_zero(self):
        r"""
        Return ``True`` if this element is zero,
        which is always the case in the zero ring.

        EXAMPLES::

            sage: from sage.rings.zero_ring import ZeroRing
            sage: z = ZeroRing()()
            sage: z.is_zero()
            True

        """
        return True

    def is_one(self):
        r"""
        Return ``True`` if this element is one,
        which is always the case in the zero ring.

        EXAMPLES::

            sage: from sage.rings.zero_ring import ZeroRing
            sage: z = ZeroRing()()
            sage: z.is_one()
            True

        """
        return True

    def _richcmp_(self, other, op):
        r"""
        Compare this element with ``other``.

        EXAMPLES::

            sage: from sage.rings.zero_ring import ZeroRing
            sage: Z = ZeroRing()
            sage: Z(0) == Z(1)
            True

        """
        return rich_to_bool(op, 0)

    def _add_(self, other):
        r"""
        Return the sum of this element and ``other``.

        EXAMPLES::

            sage: from sage.rings.zero_ring import ZeroRing
            sage: Z = ZeroRing()
            sage: Z(0) + Z(1)
            0

        """
        return self

    def _sub_(self, other):
        r"""
        Return the subtraction of this element and ``other``.

        EXAMPLES::

            sage: from sage.rings.zero_ring import ZeroRing
            sage: Z = ZeroRing()
            sage: Z(0) - Z(1)
            0

        """
        return self

    def _mul_(self, other):
        r"""
        Return the product of this element and ``other``.

        EXAMPLES::

            sage: from sage.rings.zero_ring import ZeroRing
            sage: Z = ZeroRing()
            sage: Z(0) * Z(1)
            0

        """
        return self

    def __invert__(self):
        r"""
        Return the inverse of this element.

        EXAMPLES::

            sage: from sage.rings.zero_ring import ZeroRing
            sage: Z = ZeroRing()
            sage: ~Z(0)
            0

        """
        return self

    def is_unit(self):
        r"""
        Return ``True`` if this element is a unit,
        which is always the case in the zero ring.

        EXAMPLES::

            sage: from sage.rings.zero_ring import ZeroRing
            sage: z = ZeroRing()()
            sage: z.is_unit()
            True

        """
        return True

    def inverse_of_unit(self):
        r"""
        Return the inverse of this element.

        EXAMPLES::

            sage: from sage.rings.zero_ring import ZeroRing
            sage: z = ZeroRing()()
            sage: z.inverse_of_unit()
            0

        """
        return self


class ZeroRing(UniqueRepresentation, CommutativeRing):
    r"""
    A class for the zero ring.

    TESTS::

        sage: from sage.rings.zero_ring import ZeroRing
        sage: Z = ZeroRing()
        sage: TestSuite(Z).run()

    """
    Element = ZeroRingElement

    def __init__(self, base=ZZ, category=None):
        r"""
        Initialize this ring.

        TESTS::

            sage: from sage.rings.zero_ring import ZeroRing
            sage: Z = ZeroRing()
            sage: Z
            The zero ring

        """
        CommutativeRing.__init__(self, base, category=category)
        self._names = tuple()

    def _repr_(self):
        r"""
        Return a string representation of this ring

        EXAMPLES::

            sage: from sage.rings.zero_ring import ZeroRing
            sage: Z = ZeroRing()
            sage: str(Z)  # indirect doctest
            'The zero ring'

        """
        return "The zero ring"

    def _latex_(self):
        r"""
        Return a LaTeX representation of this ring

        EXAMPLES::

            sage: from sage.rings.zero_ring import ZeroRing
            sage: Z = ZeroRing()
            sage: latex(Z)  # indirect doctest
            0

        """
        return "0"

    def _coerce_map_from_(self, R):
        r"""
        Return ``True`` if there is a coercion map from
        `R` to this ring.

        This method always returns true since any ring
        canonically maps to the zero ring.

        EXAMPLES::

            sage: from sage.rings.zero_ring import ZeroRing
            sage: Z = ZeroRing()
            sage: Z.has_coerce_map_from(Z)
            True
            sage: Z.has_coerce_map_from(ZZ)
            True
            sage: Z.has_coerce_map_from(QQ)
            True

        """
        return True

    def is_zero(self):
        r"""
        Return ``True`` if this ring is zero.

        EXAMPLES::

            sage: from sage.rings.zero_ring import ZeroRing
            sage: Z = ZeroRing()
            sage: Z.is_zero()
            True

        """
        return True

    def zero(self):
        r"""
        Return the zero element in this ring.

        EXAMPLES::

            sage: from sage.rings.zero_ring import ZeroRing
            sage: Z = ZeroRing()
            sage: Z.zero()
            0

        """
        return ZeroRingElement(self)

    def one(self):
        r"""
        Return the one element in this ring.

        EXAMPLES::

            sage: from sage.rings.zero_ring import ZeroRing
            sage: Z = ZeroRing()
            sage: Z.one()
            0

        """
        return ZeroRingElement(self)

    def ngens(self):
        r"""
        Return the number of generators of this ring.

        EXAMPLES::

            sage: from sage.rings.zero_ring import ZeroRing
            sage: Z = ZeroRing()
            sage: Z.ngens()
            0

        """
        return ZZ(0)

    def gens(self):
        r"""
        Return the generators of this ring.

        EXAMPLES::

            sage: from sage.rings.zero_ring import ZeroRing
            sage: Z = ZeroRing()
            sage: Z.gens()
            ()

        """
        return tuple()

    def gen(self, i=0):
        r"""
        Return the ``i``-th generator of this ring.

        EXAMPLES::

            sage: from sage.rings.zero_ring import ZeroRing
            sage: Z = ZeroRing()
            sage: Z.gen(0)
            Traceback (most recent call last):
            ...
            IndexError: generator not defined

        """
        raise IndexError("generator not defined")

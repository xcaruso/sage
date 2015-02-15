r"""
Arbitrary precision real intervals using Arb

AUTHORS:

- Clemens Heuberger (2014-10-21): Initial version.

This is a rudimentary binding to the optional `Arb library
<http://fredrikj.net/arb/>`_; it may be useful to refer to its
documentation for more details.

You may have to run ``sage -i arb`` to use the arb library.

Parts of the documentation for this module are copied or adapted from
Arb's own documentation, licenced under the GNU General Public License
version 2, or later.

Comparison
==========

.. WARNING::

    Identical :class:`RealBall` objects are understood to give
    permission for algebraic simplification. This assumption is made
    to improve performance.  For example, setting ``z = x*x`` sets `z`
    to a ball enclosing the set `\{t^2 : t \in x\}` and not the
    (generally larger) set `\{tu : t \in x, u \in x\}`.

Two elements are equal if and only if they are the same object
or if both are exact and equal::

    sage: from sage.rings.real_arb import RealBallField # optional - arb
    sage: RBF = RealBallField() # optional - arb
    sage: a = RBF(1) # optional - arb
    sage: b = RBF(1) # optional - arb
    sage: a is b # optional - arb
    False
    sage: a == b # optional - arb
    True
    sage: a = RBF(1/3) # optional - arb
    sage: b = RBF(1/3) # optional - arb
    sage: a.is_exact() # optional - arb
    False
    sage: b.is_exact() # optional - arb
    False
    sage: a is b # optional - arb
    False
    sage: a == b # optional - arb
    False

A ball is non-zero if and only if it does not contain zero. ::

    sage: a = RBF(RIF(-0.5, 0.5)) # optional - arb
    sage: bool(a) # optional - arb
    False
    sage: a != 0 # optional - arb
    False
    sage: b = RBF(1/3) # optional - arb
    sage: bool(b) # optional - arb
    True
    sage: b != 0 # optional - arb
    True

A ball ``left`` is less than a ball ``right`` if all elements of
``left`` are less than all elements of ``right``. ::

    sage: a = RBF(RIF(1, 2)) # optional - arb
    sage: b = RBF(RIF(3, 4)) # optional - arb
    sage: a < b # optional - arb
    True
    sage: a <= b # optional - arb
    True
    sage: a > b # optional - arb
    False
    sage: a >= b # optional - arb
    False
    sage: a = RBF(RIF(1, 3)) # optional - arb
    sage: b = RBF(RIF(2, 4)) # optional - arb
    sage: a < b # optional - arb
    False
    sage: a <= b # optional - arb
    False
    sage: a > b # optional - arb
    False
    sage: a >= b # optional - arb
    False

TESTS::

    sage: from sage.rings.real_arb import RBF
    sage: (RBF(pi) * identity_matrix(QQ, 3)).parent()
    Full MatrixSpace of 3 by 3 dense matrices over Real ball field with 53 bits precision

Classes and Methods
===================
"""
#*****************************************************************************
# Copyright (C) 2014 Clemens Heuberger <clemens.heuberger@aau.at>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                http://www.gnu.org/licenses/
#*****************************************************************************

include 'sage/ext/interrupt.pxi'
include "sage/ext/python.pxi"
include "sage/ext/stdsage.pxi"

import operator

import sage.symbolic.constants

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.real_mpfi import RealIntervalField, RealIntervalField_class
from sage.structure.unique_representation import UniqueRepresentation

cimport sage.rings.integer
cimport sage.rings.rational
cimport sage.structure.element

from sage.libs.arb.arb cimport *
from sage.libs.arb.arf cimport arf_t, arf_init, arf_get_mpfr, arf_set_mpfr, arf_clear, arf_set_mag
from sage.libs.arb.mag cimport mag_t, mag_init, mag_clear, mag_add, mag_set_d, MAG_BITS
from sage.libs.flint.flint cimport flint_free
from sage.libs.flint.fmpz cimport fmpz_t, fmpz_init, fmpz_set_mpz, fmpz_clear
from sage.libs.flint.fmpq cimport fmpq_t, fmpq_init, fmpq_set_mpq, fmpq_clear
from sage.libs.mpfi cimport mpfi_get_left, mpfi_get_right, mpfi_interv_fr
from sage.libs.mpfr cimport mpfr_t, mpfr_init2, mpfr_clear, GMP_RNDN, MPFR_PREC_MIN
from sage.rings.real_double cimport RealDoubleElement
from sage.rings.real_mpfr cimport RealField_class, RealField, RealNumber
from sage.structure.element cimport Element, ModuleElement, RingElement

cdef void mpfi_to_arb(arb_t target, const mpfi_t source, const long precision):
    """
    Convert an MPFI interval to an Arb ball.

    INPUT:

    - ``target`` -- an ``arb_t``.

    - ``source`` -- an ``mpfi_t``.

    - ``precision`` -- an integer `\ge 2`.
    """
    cdef mpfr_t left
    cdef mpfr_t right

    if _do_sig(precision): sig_on()

    mpfr_init2(left, precision)
    mpfr_init2(right, precision)

    mpfi_get_left(left, source)
    mpfi_get_right(right, source)

    arb_set_interval_mpfr(target,
                          left,
                          right,
                          precision)

    mpfr_clear(left)
    mpfr_clear(right)

cdef int arb_to_mpfi(mpfi_t target, arb_t source, const long precision) except -1:
    """
    Convert an Arb ball to an MPFI interval.

    INPUT:

    - ``target`` -- an ``mpfi_t``.

    - ``source`` -- an ``arb_t``.

    - ``precision`` -- an integer `\ge 2`.

    EXAMPLES::

        sage: cython("\n".join([ # optional - arb
        ....:     '#cinclude $SAGE_ROOT/local/include/flint',
        ....:     '#clib arb',
        ....:     'from sage.rings.real_mpfi cimport RealIntervalFieldElement',
        ....:     'from sage.libs.arb.arb cimport *',
        ....:     'from sage.rings.real_arb cimport arb_to_mpfi',
        ....:     'from sage.rings.real_mpfi import RIF',
        ....:     '',
        ....:     'cdef extern from "arb.h":',
        ....:     '    void arb_pow_ui(arb_t y, const arb_t b, unsigned long e, long prec)',
        ....:     '',
        ....:     'cdef RealIntervalFieldElement result',
        ....:     'cdef arb_t arb',
        ....:     'arb_init(arb)',
        ....:     'result = RIF(0)',
        ....:     'arb_set_ui(arb, 65536)',
        ....:     'arb_pow_ui(arb, arb, 65536**3 * 65535, 53)',
        ....:     'arb_to_mpfi(result.value, arb, 53)',
        ....:     'arb_clear(arb)'
        ....: ]))
        Traceback (most recent call last):
        ...
        ArithmeticError: Error converting arb to mpfi. Overflow?
    """
    cdef mpfr_t left
    cdef mpfr_t right

    mpfr_init2(left, precision)
    mpfr_init2(right, precision)

    try:
        sig_on()
        arb_get_interval_mpfr(left, right, source)
        mpfi_interv_fr(target, left, right)
        sig_off()
    except RuntimeError:
        raise ArithmeticError("Error converting arb to mpfi. Overflow?")
    finally:
        mpfr_clear(left)
        mpfr_clear(right)

class RealBallField(UniqueRepresentation, Parent):
    r"""
    An approximation of the field of real numbers using mid-rad intervals, also
    known as balls.

    INPUT:

    - ``precision`` -- an integer `\ge 2`.

    EXAMPLES::

        sage: from sage.rings.real_arb import RealBallField # optional - arb
        sage: RBF = RealBallField() # optional - arb; indirect doctest
        sage: RBF(1) # optional - arb
        1.000000000000000

    ::

        sage: from sage.rings.real_arb import RBF
        sage: (1/2*RBF(1)) + AA(sqrt(2)) - 1 + polygen(QQ, x)
        x + [0.914213562373095 +/- 4.10e-16]
    """
    Element = RealBall

    @staticmethod
    def __classcall__(cls, long precision=53, category=None):
        r"""
        Normalize the arguments for caching.

        TESTS::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: RealBallField(53) is RealBallField() # optional - arb
            True
        """
        return super(RealBallField, cls).__classcall__(cls, precision, category)

    def __init__(self, precision, category):
        r"""
        Initialize the real ball field.

        INPUT:

        - ``precision`` -- an integer `\ge 2`.

        EXAMPLES::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: RBF = RealBallField() # optional - arb
            sage: RBF(1) # optional - arb
            1.000000000000000
            sage: RealBallField(0) # optional - arb
            Traceback (most recent call last):
            ...
            ValueError: Precision must be at least 2.
            sage: RealBallField(1) # optional - arb
            Traceback (most recent call last):
            ...
            ValueError: Precision must be at least 2.
        """
        if precision < 2:
            raise ValueError("Precision must be at least 2.")
        super(RealBallField, self).__init__(
                #category=category or sage.categories.magmas_and_additive_magmas.MagmasAndAdditiveMagmas().Infinite(),
                # FIXME: RBF is not even associative, but CompletionFunctor only works with rings.
                category=category or sage.categories.rings.Rings().Infinite(),
                element_constructor=RealBall)
        self._prec = precision

    def _repr_(self):
        r"""
        String representation of ``self``.

        EXAMPLES::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: RealBallField() # optional - arb
            Real ball field with 53 bits precision
            sage: RealBallField(106) # optional - arb
            Real ball field with 106 bits precision
        """
        return "Real ball field with {} bits precision".format(self._prec)

    def _coerce_map_from_(self, other):
        r"""
        Parents that canonically coerce into real ball fields include:

        - some exact or lazy parents representing subsets of the reals, such as
          `ZZ`, `QQ`, `AA`, and `RLF`;

        - real ball fields with a larger precision.

        TESTS::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: RealBallField().has_coerce_map_from(RealBallField(54))
            True
            sage: RealBallField().has_coerce_map_from(RealBallField(52))
            False
            sage: RealBallField().has_coerce_map_from(RIF) # optional - arb
            False
            sage: RealBallField().has_coerce_map_from(SR) # optional - arb
            False
            sage: RealBallField().has_coerce_map_from(RR) # optional - arb
            False
        """
        from sage.rings.qqbar import AA
        from sage.rings.real_lazy import RLF
        if isinstance(other, RealBallField):
            return (other._prec >= self._prec)
        elif (other is ZZ) or (other is QQ) or (other is AA) or (other is RLF):
            return True
        else:
            return False

    def gens(self):
        r"""
        EXAMPLE::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: RealBallField().gens()
            (1.000000000000000,)
        """
        return (self.one(),)

    def construction(self):
        """
        Return the construction of a real ball field as a completion of the
        rationals.

        EXAMPLES::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: RBF = RealBallField(42) # optional - arb
            sage: functor, base = RBF.construction()
            sage: functor, base
            (Completion[+Infinity], Rational Field)
            sage: functor(base) is RBF
            True
        """
        from sage.categories.pushout import CompletionFunctor
        functor = CompletionFunctor(sage.rings.infinity.Infinity,
                                    self._prec,
                                    {'type': 'Ball'})
        return functor, QQ

    def precision(self):
        """
        Return the bit precision used for operations on elements of this field.

        EXAMPLES::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: RealBallField().precision() # optional - arb
            53
        """
        return self._prec

    def is_exact(self):
        """
        Real ball fields are not exact.

        EXAMPLES::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: RealBallField().is_exact() # optional - arb
            False
        """
        return False

    def characteristic(self):
        """
        Real ball fields have characteristic zero.

        EXAMPLES::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: RealBallField().characteristic() # optional - arb
            0
        """
        return 0


cdef inline bint _do_sig(long prec):
    """
    Whether signal handlers should be installed for calls to arb.

    TESTS::

        sage: from sage.rings.real_arb import RealBallField # optional - arb
        sage: _ = RealBallField()(1).psi() # optional - arb; indirect doctest
        sage: _ = RealBallField(1500)(1).psi() # optional - arb
    """
    return (prec > 1000)

cdef inline long prec(RealBall ball):
    return ball._parent._prec

cdef class RealBall(RingElement):
    """
    Hold one ``arb_t`` of the `Arb library
    <http://fredrikj.net/arb/>`_

    EXAMPLES::

        sage: from sage.rings.real_arb import RealBallField # optional - arb
        sage: a = RealBallField()(RIF(1))                     # optional - arb; indirect doctest
        sage: b = a.psi()                         # optional - arb
        sage: b                                   # optional - arb
        [-0.577215664901533 +/- 3.85e-16]
        sage: b._interval()        # optional - arb
        -0.577215664901533?
    """

    def __cinit__(self):
        """
        Allocate memory for the encapsulated value.

        EXAMPLES::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: RealBallField()(RIF(1)) # optional - arb; indirect doctest
            1.000000000000000
        """
        arb_init(self.value)

    def __dealloc__(self):
        """
        Deallocate memory of the encapsulated value.

        EXAMPLES::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: a = RealBallField()(RIF(1)) # optional - arb; indirect doctest
            sage: del a # optional - arb
        """
        arb_clear(self.value)

    def __init__(self, parent, mid=None, rad=None):
        """
        Initialize the :class:`RealBall` using ``mid``.

        INPUT:

        - ``parent`` -- a :class:`RealBallField`.

        - ``mid`` (optional) --  ball midpoint, see examples below. If omitted,
          initialize the ball to zero, ignoring the ``rad`` argument.

        - ``rad`` (optional) -- a :class:`RealDoubleElement`, ball radius. If
          the midpoint is not exactly representable in floating-point, the
          radius is adjusted to account for the roundoff error.

        EXAMPLES::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: RBF = RealBallField()
            sage: RBF()
            0

        One can create exact real balls using elements of various exact parents,
        or using floating-point numbers::

            sage: RBF(3)
            3.000000000000000
            sage: RBF(3r)
            3.000000000000000
            sage: RBF(1/3)
            [0.3333333333333333 +/- 7.04e-17]
            sage: RBF(3.14)
            [3.140000000000000 +/- 1.25e-16]

        ::

            sage: RBF(3, 0.125r)
            [3e+0 +/- 0.126]
            sage: RBF(pi, 0.125r)
            [3e+0 +/- 0.267]

        Note that integers and floating-point numbers are ''not'' rounded to
        the parent's precision::

            sage: b = RBF(11111111111111111111111111111111111111111111111); b
            [1.111111111111111e+46 +/- 1.12e+30]
            sage: b.mid().exact_rational()
            11111111111111111111111111111111111111111111111

        Similarly, converting a real ball from one real ball field to another
        (with a different precision) only changes the way it is displayed and
        the precision of operations involving it, not the actual representation
        of its center::

            sage: RBF100 = RealBallField(100)
            sage: b100 = RBF100(1/3); b100
            [0.333333333333333333333333333333 +/- 4.65e-31]
            sage: b53 = RBF(b100); b53
            [0.3333333333333333 +/- 3.34e-17]
            sage: RBF100(b53)
            [0.333333333333333333333333333333 +/- 4.65e-31]

        Special values are supported::

            sage: RBF(oo).mid(), RBF(-oo).mid(), RBF(unsigned_infinity).mid()
            (+infinity, -infinity, 0.000000000000000)
            sage: RBF(NaN)
            nan

        Anything that is convertible to a real interval can also be used to
        construct a real ball::

            sage: RBF(RIF(0, 1))                  # optional - arb; indirect doctest
            [+/- 1.01]
            sage: RBF(1)                          # optional - arb
            1.000000000000000
            sage: RBF(x)                          # optional - arb
            Traceback (most recent call last):
            ...
            TypeError: unable to convert x to a RealIntervalFieldElement

        """
        # TODO: init from mid+rad
        # TODO: perhaps dissociate basic constructors (here) from conversions
        #       (→ element_constructor)?
        cdef fmpz_t tmpz
        cdef fmpq_t tmpq
        cdef arf_t  tmpr
        cdef mag_t  tmpm

        super(RealBall, self).__init__(parent)

        if mid is None:
            return

        elif isinstance(mid, RealBall):
            arb_set(self.value, (<RealBall> mid).value) # no rounding!
        elif isinstance(mid, int):
            arb_set_si(self.value, PyInt_AS_LONG(mid)) # no rounding!
        elif isinstance(mid, sage.rings.integer.Integer):
            if _do_sig(prec(self)): sig_on()
            fmpz_init(tmpz)
            fmpz_set_mpz(tmpz, (<sage.rings.integer.Integer> mid).value) # no rounding!
            arb_set_fmpz(self.value, tmpz)
            fmpz_clear(tmpz)
            if _do_sig(prec(self)): sig_off()
        elif isinstance(mid, sage.rings.rational.Rational):
            if _do_sig(prec(self)): sig_on()
            fmpq_init(tmpq)
            fmpq_set_mpq(tmpq, (<sage.rings.rational.Rational> mid).value)
            arb_set_fmpq(self.value, tmpq, prec(self))
            fmpq_clear(tmpq)
            if _do_sig(prec(self)): sig_off()
        elif isinstance(mid, RealNumber):
            if _do_sig(prec(self)): sig_on()
            arf_init(tmpr)
            arf_set_mpfr(tmpr, (<RealNumber> mid).value)
            arb_set_arf(self.value, tmpr) # no rounding!
            arf_clear(tmpr)
            if _do_sig(prec(self)): sig_off()
        elif isinstance(mid, sage.rings.infinity.AnInfinity):
            if isinstance(mid, sage.rings.infinity.PlusInfinity):
                arb_pos_inf(self.value)
            elif isinstance(mid, sage.rings.infinity.MinusInfinity):
                arb_neg_inf(self.value)
            else:
                arb_zero_pm_inf(self.value)

        else:
            if _do_sig(prec(self)): sig_on()
            if isinstance(mid, sage.symbolic.constants.Constant):
                if isinstance(mid, sage.symbolic.constants.NotANumber):
                    arb_indeterminate(self.value)
                elif isinstance(mid, sage.symbolic.constants.Pi):
                    arb_const_pi(self.value, prec(self))
                elif isinstance(mid, sage.symbolic.constants.E):
                    arb_const_e(self.value, prec(self))
                elif isinstance(mid, sage.symbolic.constants.Log2):
                    arb_const_log2(self.value, prec(self))
                elif isinstance(mid, sage.symbolic.constants.Catalan):
                    arb_const_catalan(self.value, prec(self))
                elif isinstance(mid, sage.symbolic.constants.Khinchin):
                    arb_const_khinchin(self.value, prec(self))
                elif isinstance(mid, sage.symbolic.constants.Glaisher):
                    arb_const_glaisher(self.value, prec(self))
            if _do_sig(prec(self)): sig_off()
            if not isinstance(mid, RealIntervalFieldElement):
                try:
                    mid = RealIntervalField(prec(self))(mid)
                except TypeError:
                    raise TypeError("unable to convert {} to a RealIntervalFieldElement".format(mid))
            mpfi_to_arb(self.value,
                        (<RealIntervalFieldElement> mid).value,
                        prec(self))

        if rad is not None:
            if isinstance(rad, float):
                mag_init(tmpm)
                mag_set_d(tmpm, PyFloat_AS_DOUBLE(rad))
                mag_add(arb_radref(self.value), arb_radref(self.value), tmpm)
                mag_clear(tmpm)
            else:
                raise TypeError("rad should be a Python float")


    cdef RealBall _new(self):
        """
        Return a new real ball element with the same parent as ``self``.

        TESTS::

            sage: from sage.rings.real_arb import RealBallField
            sage: RealBallField()(2)**2 # indirect doctest
            4.000000000000000

        """
        cdef RealBall x
        x = RealBall.__new__(RealBall)
        x._parent = self._parent
        return x

    def _repr_(self):
        """
        Return a string representation of ``self``.

        OUTPUT:

        A string.

        EXAMPLES::

           sage: from sage.rings.real_arb import RealBallField # optional - arb
           sage: RealBallField()(RIF(1.9, 2))  # optional - arb
           [2e+0 +/- 0.101]
        """
        cdef char* c_result
        cdef bytes py_string

        c_result = arb_get_str(self.value, (prec(self) * 31) // 100, 0)
        try:
            py_string = c_result
        finally:
            flint_free(c_result)

        return py_string

    cpdef RealIntervalFieldElement _interval(self):
        """
        Return a :mod:`real interval <sage.rings.real_mpfr>` containing this ball.

        OUTPUT:

        A :class:`~sage.rings.real_mpfi.RealIntervalFieldElement`.

        EXAMPLES::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: a = RealBallField()(RIF(2))                     # optional - arb
            sage: a._interval()        # optional - arb
            2
        """
        cdef RealIntervalFieldElement result
        result = RealIntervalField(prec(self))(0)
        arb_to_mpfi(result.value, self.value, prec(self))
        return result

    # Comparisons and predicates

    def is_zero(self):
        """
        Return ``True`` iff the midpoint and radius of this ball are both zero.

        EXAMPLES::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: RBF = RealBallField() # optional - arb
            sage: RBF(0).is_zero() # optional - arb
            True
            sage: RBF(RIF(-0.5, 0.5)).is_zero() # optional - arb
            False
        """
        return arb_is_zero(self.value)

    def __nonzero__(self):
        """
        Return ``True`` iff zero is not contained in the interval represented
        by this ball.

        EXAMPLES::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: RBF = RealBallField() # optional - arb
            sage: bool(RBF(pi)) # optional - arb
            True
            sage: bool(RBF(RIF(-0.5, 0.5))) # optional - arb
            False
        """
        return arb_is_nonzero(self.value)

    def is_exact(self):
        """
        Return ``True`` iff the radius of this ball is zero.

        EXAMPLES::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: RBF = RealBallField() # optional - arb
            sage: RBF(1).is_exact() # optional - arb
            True
            sage: RBF(RIF(0.1, 0.2)).is_exact() # optional - arb
            False
        """
        return arb_is_exact(self.value)

    def __richcmp__(left, right, int op):
        """
        Compare ``left`` and ``right``.

        For more information, see :mod:`sage.rings.real_arb`.

        EXAMPLES::

                sage: from sage.rings.real_arb import RealBallField # optional - arb
                sage: RBF = RealBallField() # optional - arb
                sage: a = RBF(1) # optional - arb
                sage: b = RBF(1) # optional - arb
                sage: a is b # optional - arb
                False
                sage: a == b # optional - arb
                True
                sage: a = RBF(1/3) # optional - arb
                sage: a.is_exact() # optional - arb
                False
                sage: b = RBF(1/3) # optional - arb
                sage: b.is_exact() # optional - arb
                False
                sage: a == b # optional - arb
                False
        """
        return (<Element>left)._richcmp(right, op)

    cdef _richcmp_c_impl(left, Element right, int op):
        """
        Compare ``left`` and ``right``.

        For more information, see :mod:`sage.rings.real_arb`.

        EXAMPLES::

                sage: from sage.rings.real_arb import RealBallField # optional - arb
                sage: RBF = RealBallField() # optional - arb
                sage: a = RBF(1) # optional - arb
                sage: b = RBF(1) # optional - arb
                sage: a is b # optional - arb
                False
                sage: a == b # optional - arb
                True

        TESTS:

            Balls whose intersection consists of one point::

                sage: a = RBF(RIF(1, 2)) # optional - arb
                sage: b = RBF(RIF(2, 4)) # optional - arb
                sage: a < b # optional - arb
                False
                sage: a > b # optional - arb
                False
                sage: a <= b # optional - arb
                False
                sage: a >= b # optional - arb
                False
                sage: a == b # optional - arb
                False
                sage: a != b # optional - arb
                False

            Balls with non-trivial intersection::

                sage: a = RBF(RIF(1, 4)) # optional - arb
                sage: a = RBF(RIF(2, 5)) # optional - arb
                sage: a < b # optional - arb
                False
                sage: a <= b # optional - arb
                False
                sage: a > b # optional - arb
                False
                sage: a >= b # optional - arb
                False
                sage: a == b # optional - arb
                False
                sage: a != b # optional - arb
                False

            One ball contained in another::

                sage: a = RBF(RIF(1, 4)) # optional - arb
                sage: b = RBF(RIF(2, 3)) # optional - arb
                sage: a < b # optional - arb
                False
                sage: a <= b # optional - arb
                False
                sage: a > b # optional - arb
                False
                sage: a >= b # optional - arb
                False
                sage: a == b # optional - arb
                False
                sage: a != b # optional - arb
                False

            Disjoint balls::

                sage: a = RBF(1/3) # optional - arb
                sage: b = RBF(1/2) # optional - arb
                sage: a < b # optional - arb
                True
                sage: a <= b # optional - arb
                True
                sage: a > b # optional - arb
                False
                sage: a >= b # optional - arb
                False
                sage: a == b # optional - arb
                False
                sage: a != b # optional - arb
                True

            Exact elements::

                sage: a = RBF(2) # optional - arb
                sage: b = RBF(2) # optional - arb
                sage: a.is_exact() # optional - arb
                True
                sage: b.is_exact() # optional - arb
                True
                sage: a < b # optional - arb
                False
                sage: a <= b # optional - arb
                True
                sage: a > b # optional - arb
                False
                sage: a >= b # optional - arb
                True
                sage: a == b # optional - arb
                True
                sage: a != b # optional - arb
                False
        """
        cdef RealBall lt, rt
        cdef arb_t difference

        lt = left
        rt = right

        if op == Py_EQ:
            return (lt is rt) or (
                arb_is_exact(lt.value) and arb_is_exact(rt.value)
                and arb_equal(lt.value, rt.value))

        arb_init(difference)
        arb_sub(difference, lt.value, rt.value, prec(lt))

        if op == Py_NE:
            result = arb_is_nonzero(difference)
        elif op == Py_GT:
            result = arb_is_positive(difference)
        elif op == Py_GE:
            result = arb_is_nonnegative(difference)
        elif op == Py_LT:
            result = arb_is_negative(difference)
        elif op == Py_LE:
            result = arb_is_nonpositive(difference)

        arb_clear(difference)
        return result

    # Center and radius

    def mid(self):
        """
        Return the center of this ball.

        EXAMPLES::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: RealBallField(16)(1/3).mid()
            0.3333
            sage: RealBallField(16)(1/3).mid().parent()
            Real Field with 15 bits of precision
        """
        cdef long mid_prec = arb_bits(self.value) or prec(self)
        if mid_prec < MPFR_PREC_MIN:
            mid_prec = MPFR_PREC_MIN
        cdef RealField_class mid_field = RealField(mid_prec)
        cdef RealNumber mid = RealNumber(mid_field, None)
        if _do_sig(mid_prec): sig_on()
        cdef int rnd = arf_get_mpfr(mid.value, arb_midref(self.value), GMP_RNDN)
        if _do_sig(mid_prec): sig_off()
        if rnd != 0:
            raise OverflowError("Unable to represent the center of this ball within the exponent range of RealNumbers")
        return mid

    def rad(self):
        """
        Return the radius of this ball.

        EXAMPLES::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: RealBallField()(1/3).rad()
            5.5511151e-17
            sage: RealBallField()(1/3).rad().parent()
            Real Field with 30 bits of precision
        """
        cdef RealField_class rad_field = RealField(MAG_BITS)
        cdef RealNumber rad = RealNumber(rad_field, None)
        cdef arf_t tmp
        arf_init(tmp)
        arf_set_mag(tmp, arb_radref(self.value))
        cdef int rnd = arf_get_mpfr(rad.value, tmp, GMP_RNDN)
        arf_clear(tmp)
        if rnd != 0:
            raise OverflowError("Unable to represent the radius of this ball within the exponent range of RealNumbers")
        return rad

    # Precision

    def round(self):
        """
        Return a copy of this ball with center rounded to the precision of the
        parent.

        .. SEEALSO:: :meth:`trim`

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: b = RBF(pi.n(100))
            sage: b.mid()
            3.141592653589793238462643383
            sage: b.round().mid()
            3.1415926535898
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_set_round(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def trim(self):
        """
        Return a trimmed copy of this ball.

        Round `self` to a number of bits equal to the accuracy of `self` (as
        indicated by its radius), plus a few guard bits. The resulting ball is
        guaranteed to contain `self`, but is more economical if `self` has less
        than full accuracy.

        .. SEEALSO:: :meth:`round`

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: b = RBF(RIF(3.1415,3.1416))
            sage: b.mid()
            3.14155000000000
            sage: b.trim().mid()
            3.14155000
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_trim(res.value, self.value)
        if _do_sig(prec(self)): sig_off()
        return res

    # Arithmetic

    def __neg__(self):
        """
        Return the opposite of this ball.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: -RBF(1/3)
            [-0.3333333333333333 +/- 7.04e-17]
        """
        cdef RealBall res = self._new()
        arb_neg(res.value, self.value)
        return res

    def __invert__(self):
        """
        Return the inverse of this ball.

        The result is guaranteed to contain the inverse of any point of the
        input ball.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: ~RBF(5)
            [0.2000000000000000 +/- 4.45e-17]
            sage: ~RBF(0)
            [+/- inf]
            sage: RBF(RIF(-0.1,0.1))
            [+/- 0.101]

        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_inv(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    cpdef ModuleElement _add_(self, ModuleElement other):
        """
        Return the sum of two balls, rounded to the ambient field's precision.

        The resulting ball is guaranteed to contain the sums of any two points
        of the respective input balls.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(1) + RBF(1/3)
            [1.333333333333333 +/- 5.37e-16]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_add(res.value, self.value, (<RealBall> other).value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    cpdef ModuleElement _sub_(self, ModuleElement other):
        """
        Return the difference of two balls, rounded to the ambient field's
        precision.

        The resulting ball is guaranteed to contain the differences of any two
        points of the respective input balls.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(1) - RBF(1/3)
            [0.666666666666667 +/- 5.37e-16]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_sub(res.value, self.value, (<RealBall> other).value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    cpdef RingElement _mul_(self, RingElement other):
        """
        Return the product of two balls, rounded to the ambient field's
        precision.

        The resulting ball is guaranteed to contain the products of any two
        points of the respective input balls.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(-2) * RBF(1/3)
            [-0.666666666666667 +/- 4.82e-16]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_mul(res.value, self.value, (<RealBall> other).value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    cpdef RingElement _div_(self, RingElement other):
        """
        Return the quotient of two balls, rounded to the ambient field's
        precision.

        The resulting ball is guaranteed to contain the quotients of any two
        points of the respective input balls.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(pi)/RBF(e)
            [1.15572734979092 +/- 2.49e-15]
            sage: RBF(2)/RBF(0)
            [+/- inf]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_div(res.value, self.value, (<RealBall> other).value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def __pow__(base, expo, _):
        """
        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(e)^17
            [24154952.753575 +/- 3.79e-7]
            sage: RBF(e)^(-1)
            [0.367879441171442 +/- 4.80e-16]
            sage: RBF(e)^(1/2)
            [1.648721270700128 +/- 5.63e-16]
            sage: RBF(e)^RBF(pi)
            [23.140692632779 +/- 3.08e-13]

        ::

            sage: RBF(-1)^(1/3)
            nan
            sage: RBF(0)^(-1)
            [+/- inf]
            sage: RBF(-e)**RBF(pi)
            nan

        TESTS::

            sage: RBF(e)**(2r)
            [7.38905609893065 +/- 5.88e-15]
            sage: RBF(e)**(-1r)
            [0.367879441171442 +/- 4.80e-16]
        """
        cdef fmpz_t tmpz
        if not isinstance(base, RealBall):
            return sage.structure.element.bin_op(base, expo, operator.pow)
        cdef RealBall self = base
        cdef RealBall res = self._new()
        if isinstance(expo, int) and expo > 0:
            if _do_sig(prec(self)): sig_on()
            arb_pow_ui(res.value, self.value, PyInt_AS_LONG(expo), prec(self))
            if _do_sig(prec(self)): sig_off()
        elif isinstance(expo, sage.rings.integer.Integer):
            if _do_sig(prec(self)): sig_on()
            fmpz_init(tmpz)
            fmpz_set_mpz(tmpz, (<sage.rings.integer.Integer> expo).value)
            arb_pow_fmpz(res.value, self.value, tmpz, prec(self))
            fmpz_clear(tmpz)
            if _do_sig(prec(self)): sig_off()
        elif isinstance(expo, RealBall):
            if _do_sig(prec(self)): sig_on()
            arb_pow(res.value, self.value, (<RealBall> expo).value, prec(self))
            if _do_sig(prec(self)): sig_off()
        else:
            return sage.structure.element.bin_op(base, expo, operator.pow)
        return res

    # Special functions

    cpdef RealBall psi(self):
        """
        Compute the digamma function with argument self.

        OUTPUT:

        A :class:`RealBall`.

        EXAMPLES::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: a = RealBallField()(RIF(1))                     # optional - arb
            sage: a.psi()  # optional - arb
            [-0.577215664901533 +/- 3.85e-16]
        """

        cdef RealBall result = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_digamma(result.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

RBF = RealBallField()

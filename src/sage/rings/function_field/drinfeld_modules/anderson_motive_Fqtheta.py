# *****************************************************************************
#        Copyright (C) 2024 Xavier Caruso <xavier.caruso@normalesup.org>
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#                   http://www.gnu.org/licenses/
# *****************************************************************************

from sage.misc.timing import walltime

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute

from sage.functions.other import ceil
from sage.functions.other import binomial
from sage.matrix.constructor import matrix

from sage.rings.infinity import Infinity
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.fraction_field import FractionField_1poly_field
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.laurent_series_ring import LaurentSeriesRing

from sage.categories.anderson_motives import AndersonMotives
from sage.rings.function_field.drinfeld_modules.anderson_motive import AndersonMotive
from sage.rings.function_field.drinfeld_modules.morphism_to_completion import MorphismToCompletion, normalize_place


class AndersonMotive_Fqtheta(AndersonMotive):
    def __init__(self, category, tau, twist):
        AndersonMotive.__init__(self, category, tau, twist)
        K = self._K.backend(force=True)
        self._K_int = Kint = K.ring()
        self._AK_int = PolynomialRing(Kint, name=self._t_name)
        self._theta_int = Kint.gen()
        self._theta_name = Kint.variable_name()

    @lazy_attribute
    def _tau_int(self):
        q = self._q
        K = self._K.backend(force=True)
        Kint = self._K_int
        AK = PolynomialRing(K, name=self._t_name)
        denom = Kint.one()
        tau = self._tau.change_ring(AK)
        for entry in tau.list():
            for coeff in entry.list():
                denom = denom.lcm(coeff.denominator())
        scalar = Kint.one()
        for F, m in denom.factor():
            e = (-m) // (q-1)
            scalar *= F ** (-e*(q-1))
        return (scalar * tau).change_ring(self._AK_int)

    @lazy_attribute
    def _dettau_int(self):
        det = self._tau_int.det()
        return det.leading_coefficient(), det.degree()

    @lazy_attribute
    def _bad_places(self):
        disc, _ = self._dettau_int
        return [F for F, _ in disc.factor()]

    @cached_method
    def _local_maximal_model(self, place):
        q = self._q
        Fq = self._Fq
        r = self.rank()

        F = Fq.extension(place, name='w')
        if place.degree() == 1:
            a = -place[0]
        else:
            a = F.gen()
        A = PolynomialRing(F, self._t_name)
        K = PolynomialRing(F, 'v')
        v = K.gen()
        S = PolynomialRing(K, self._t_name)
        base_map = F.frobenius_endomorphism(self._deg * (place.degree() - 1))
        psiA = A.hom([A.gen()], base_map=base_map)

        tau = matrix(S, r, r,
                     [ S([c(v+a) for c in entry.list()]) for entry in self._tau_int.list() ])
        disc, _ = self._dettau_int
        val = val0 = disc.valuation(place)

        # 1st sequence: the N_i
        while val >= q - 1:
            rows = [ [ A([f[e] for f in c.list()])
                       for c in taurow.list() for e in range(q) ]
                     for taurow in tau.rows() ]
            M = matrix(A, rows)
            basis = M.minimal_kernel_basis()
            dim = basis.nrows()
            if dim == 0:
                break
            M = basis.stack(basis.basis_completion())
            N = M.inverse_of_unit()
            N = N.parent()([psiA(x) for x in N.list()])
            tau = M * tau * N
            for i in range(dim):
                for j in range(dim):
                    tau[i,j] = S([x >> q-1 for x in tau[i,j].list()])
                for j in range(dim, r):
                    tau[i,j] = S([x >> q for x in tau[i,j].list()])
                    tau[j,i] = S([x << 1 for x in tau[j,i].list()])
            val -= (q-1) * dim

        # 2nd sequence: the L_i
        if val >= q - 1:
            # The following can be improved
            while True:
                rows = [ [ A([f[e] for f in c.list()])
                           for c in taurow.list() for e in range(q-1) ]
                         for taurow in tau.rows() ]
                M = matrix(A, rows)
                if M.is_zero():
                    break
                basis = M.minimal_kernel_basis()
                dim = basis.nrows()
                M = basis.stack(basis.basis_completion())
                N = M.inverse_of_unit()
                N = N.parent()([psiA(x) for x in N.list()])
                tau = M * tau * N
                for i in range(dim, r):
                    for j in range(dim, r):
                        tau[i,j] = S([x << q-1 for x in tau[i,j].list()])
                    for j in range(dim):
                        tau[i,j] = S([x << q for x in tau[i,j].list()])
                        tau[j,i] = S([x >> 1 for x in tau[j,i].list()])
                val += (q-1) * (r - dim)
            for i in range(r):
                for j in range(r):
                    tau[i,j] = S([x >> q-1 for x in tau[i,j].list()])
            val -= (q-1) * r

        if val < val0:
            tau0 = matrix(A, r, r,
                          [ A([c(a) for c in entry.list()]) for entry in self._tau_int.list() ])
            tau1 = matrix(A, r, r,
                          [ A([c[0] for c in entry.list()]) for entry in tau.list() ])
            return val, (tau0, tau1)
        else:
            return val, None

    def discriminant(self):
        disc = self._K_int.one()
        for place in self._bad_places:
            val = self._local_maximal_model(place)[0]
            disc *= place ** val
        return disc

    def bad_reduction_places(self):
        return [place for place in self._bad_places if self._local_maximal_model(place)[0]]

    def has_good_reduction(self):
        return not bool(self.bad_reduction_places())

    def has_good_reduction_at(self, place):
        place = normalize_place(self._K_int, place, infty=False)
        return not bool(self._local_maximal_model(place)[0])

    def reduce(self, place, a=None):
        place = normalize_place(self._K_int, place, infty=False)
        val, taus = self._local_maximal_model(place)
        if val > 0:
            raise ValueError("bad reduction")

        if a is None:
            a = 'a'
        Fq = self._Fq
        F = Fq.extension(place, name=a)
        if place.degree() == 1:
            a = -place[0]
        else:
            a = F.gen()
        A = self._A
        F = F.over(A.hom([a]))

        B = PolynomialRing(F, self._t_name)
        r = self.rank()
        if taus is None:
            tau = matrix(B, r, r,
                         [ B([c(a) for c in entry.list()]) for entry in self._tau_int.list() ])
        else:
            _, tau = taus
            tau = matrix(B, r, r,
                         [ B([F(c) for c in entry.list()]) for entry in tau.list() ])

        category = AndersonMotives(F)
        return AndersonMotive(category, tau, self._twist, check=False)


    def Lseries(self, place, prec, x=None, a=None, u=None):
        n = self.rank()
        h = -self._twist
        q = self._q

        tme = walltime()

        # Figure out which precision is needed
        val = 0
        if place is Infinity:
            for entry in self._tau_int.list():
                val = max(val, entry.degree())  # t-degree
            prectau = prec + n*val - min(0, h)
        else:
            prectau = prec

        # Construct the completion
        if a is None:
            a = 'a'
        if u is None:
            u = 'u'
        completion = MorphismToCompletion(self._A, place, prectau, a, u)
        place = completion.place()
        C = completion.codomain()
        k = completion.residue_field()
        t = completion(self._A.gen())   # t in C
        Ctheta = PolynomialRing(C, self._theta_name)
        theta = Ctheta.gen()

        S = PolynomialRing(C, name='x')

        #print(" [%.5f] rings created" % walltime(tme))

        # Correction at bad places
        corr_num = S.one()
        corr_denom = S.one()
        for pl in self._bad_places:
            if place is not Infinity and pl == place(self._theta):
                continue
            _, taus = self._local_maximal_model(pl)
            if taus is None:
                continue
            tau0, tau1 = taus
            d = pl.degree()
            A = tau0.base_ring()
            F = A.base_ring()
            phiA = A.hom([A.gen()], base_map = F.frobenius_endomorphism(self._deg))
            T0 = tau0; T1 = tau1
            for _ in range(1, d):
                tau0 = tau0.parent()([phiA(y) for y in tau0.list()])
                T0 = tau0 * T0
                tau1 = tau1.parent()([phiA(y) for y in tau1.list()])
                T1 = tau1 * T1
            if d > 1:
                ell = PolynomialRing(k, name='w').quo(pl)
                Aell = A.change_ring(ell)
                T0 = T0.change_ring(Aell)
                T1 = T1.change_ring(Aell)
                Cell = PowerSeriesRing(ell, name='x')
            else:
                Cell = C
            v = Cell(t).add_bigoh(prectau)
            T0 = matrix(n, n, [f(v).add_bigoh(prectau) for f in T0.list()])
            T1 = matrix(n, n, [f(v).add_bigoh(prectau) for f in T1.list()])
            scalar = pl(t) ** self._twist
            chi0 = (scalar*T0).charpoly()
            chi1 = (scalar*T1).charpoly()
            if d > 1:
                coeffs = []
                for z in chi0.list():
                    coeff = C([y.lift()[0] for y in z.list()], z.valuation(), z.precision_absolute())
                    coeffs += [coeff] + (d-1)*[C.zero()]
                chi0 = S(coeffs)
                coeffs = []
                for z in chi1.list():
                    coeff = C([y.lift()[0] for y in z.list()], z.valuation(), z.precision_absolute())
                    coeffs += [coeff] + (d-1)*[C.zero()]
                chi1 = S(coeffs)
            corr_num *= chi0
            corr_denom *= chi1

        #print(" [%.5f] bad places handled" % walltime(tme))

        # Computation of rho
        hi = h
        current_prec = 1
        if place is Infinity:
            v = Ctheta(C.gen().add_bigoh(prectau - h))
            rho = v ** h
            while current_prec < prectau - h:
                hj = ceil(hi / q)
                e = q*hj - hi
                vt = 1 - v*theta
                rho *= vt ** e
                hi = hj
                v = v ** q
                current_prec *= q
        else:
            d = completion.degree()
            v = t.add_bigoh(prec)
            a = v[0]
            m = 0
            ev = -h
            rho = Ctheta(1)
            while current_prec < prectau:
                hj = ceil(hi / q)
                e = q*hj - hi
                rho *= (v - theta) ** e
                ev -= e * q**m
                hi = hj
                v = v ** q
                current_prec *= q
                m = (m + 1) % d
            ev %= q**d - 1
            for _ in range(d):
                ev, e = ev.quo_rem(q)
                rho *= (a - theta) ** (e + q - 1)
                a = a ** q

        #print(" [%.5f] rho computed" % walltime(tme))

        # Computation of tau*
        u = Ctheta.gen()
        B = [ ]
        vals = [ ]
        kmax = 1
        for entry in self._tau_int.list():
            e = Ctheta.zero()
            for i in range(entry.degree() + 1):
                e += entry[i](theta) * t**i
            e *= rho
            B.append(e)
            kmax = max(kmax, e.degree() // (q - 1))
        #print(" [%.5f] matrix B computed" % walltime(tme))
        rows = [ ]
        for i in range(n):
            for k in range(kmax):
                row = [ ]
                for j in range(n):
                    s = B[i + n*j]
                    row += [s[q*kp - k + q - 1] for kp in range(kmax)]
                rows.append(row)
        taudual = matrix(rows)

        #print(" [%.5f] tau* computed (size=%s)" % (walltime(tme), len(rows)))

        # Computation of the L-series
        chi = taudual.charpoly()
        #print(" [%.5f] characteristic polynomial computed" % walltime(tme))
        chi *= corr_num
        chi //= corr_denom
        L = S([l.add_bigoh(prec) for l in chi]).reverse()

        #print(" [%.5f] L-series computed" % walltime(tme))

        # Format and return the final result
        if x is None:
            return L
        elif isinstance(x, str):
            return L.change_variable_name(x)
        else:
            return L(completion(x))

    def special_value(self, place, prec, u=None, a=None):
        L = self.Lseries(place, prec, u, a)
        x = L.parent().gen()
        order = 0
        value = L(1)
        while value == 0:
            L = L // (x - 1)
            value = L(1)
            order += 1
        return value, order

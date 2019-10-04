#ifndef CIRCUIT_H
#define CIRCUIT_H

#include <initializer_list>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using std::string;
using std::to_string;

#include "FR.hpp"
#include "constants.hpp"

#define one FieldR::one()

typedef string variable;
typedef std::vector<std::pair<FieldR, variable>> LinearCombination;
typedef LinearCombination LC;
typedef std::tuple<LC, LC, LC> Constraint;

std::pair<FieldR, variable> operator*(const FieldR &l, const variable &r) {
    return std::make_pair(l, r);
}

// very inefficient if the compiler doesn't optimize it
LinearCombination operator+(LinearCombination lc,
                            std::pair<FieldR, variable> x) {
    lc.push_back(x);
    return lc;
}

struct Circuit {
    const variable ONE = "ONE";
    std::set<variable> pub = {ONE};
    std::set<variable> aux;
    std::vector<Constraint> constraints;

    variable get(variable v) {
        aux.insert(v);
        return v;
    }

    variable get_pub(variable v) {
        pub.insert(v);
        return v;
    }

    void addC(LC a, LC b, LC c) { constraints.push_back(make_tuple(a, b, c)); }

    // TODO: deal with replacements of zeros
    void replace_variable(variable v, FieldR constant, bool is_aux) {

        for (Constraint &constraint : constraints) {
            LC lc[3];
            tie(lc[0], lc[1], lc[2]) = constraint;

            for (int i = 0; i < 3; ++i)
                for (auto &p : lc[i]) {
                    if (p.second == v) {
                        p.second = ONE;
                        p.first *= constant;
                    }
                }

            constraint = make_tuple(lc[0], lc[1], lc[2]);
        }

        if (is_aux)
            aux.erase(aux.find(v));
        else
            pub.erase(pub.find(v));
    }
};

void test() {
    LC l;
    variable s;
    l = l + one * s;

    LC() + one *s;
}

void CUBE3(Circuit &c, string sc, variable output, variable x, variable y,
           FieldR z) {

    variable tmp = c.get(sc + ".tmp");

    auto lc1 = LC() + one * x + one * y + z * c.ONE;
    auto lc2 = LC() + one * tmp;
    auto lc3 = LC() + one * output;

    c.addC(lc1, lc1, lc2);
    c.addC(lc1, lc2, lc3);
}

void MiMC(Circuit &c, string sc, variable output, variable x, variable k,
          const FieldR MiMC_constants[162]) {

    variable tmp1 = c.get(sc + ".tmp1");
    variable tmp2 = c.get(sc + ".tmp2");
    variable x1 = c.get(sc + ".x1");
    variable x161 = c.get(sc + ".x161");

    auto lc1 = LC() + 1 * x + 1 * k;
    auto lc2 = LC() + 1 * tmp1;

    c.addC(lc1, lc1, lc2);

    for (int i = 1; i < 161; ++i) {
        variable xold = c.get(sc + to_string(i));
        variable xnew = c.get(sc + to_string(i + 1));
        CUBE3(c, sc + ".cubespace" + to_string(i), xnew, xold, k,
              MiMC_constants[i]);
    }

    auto lc3 = LC() + one * x161 + 1 * k + MiMC_constants[161] * c.ONE;
    auto lc4 = LC() + one * tmp2;
    auto lc5 = LC() + one * output + (-one) * k;

    c.addC(lc3, lc3, lc4);
    c.addC(lc3, lc4, lc5);
}

void CRH(Circuit &c, string sc, variable output, variable m1, variable m2) {
    variable h0 = c.get(sc + ".h0");
    variable h1 = c.get(sc + ".h1");

    MiMC(c, sc + "mimc1", h1, h0, m1, constants::hash);
    LC x, y, z;

    // currently h1 = E_m1(h0)
    // need h1 = E_m1(h0) + h1
    tie(x, y, z) = c.constraints.back();
    z = z + (-one) * h0;
    c.constraints.pop_back();
    c.addC(x, y, z);

    MiMC(c, sc + "mimc2", output, h1, m2, constants::hash);
    // currently output = E_m2(h1)
    // need output = E_m2(h1) + h1
    tie(x, y, z) = c.constraints.back();
    z = z + (-one) * h1;
    c.constraints.pop_back();
    c.addC(x, y, z);

    c.replace_variable(h0, constants::hashconstant, true);

    return;
}

template <unsigned int height>
void MerkleProof(Circuit &c, string sc, variable root, variable member) {
    for (int i = 0; i < height; ++i) {
        variable arg1 = c.get(sc + "." + to_string(i) + ".arg1");
        variable arg2 = c.get(sc + "." + to_string(i) + ".arg2");
        variable out;
        if (i == height - 1)
            out = root;
        else
            out = c.get(sc + "." + to_string(i) + ".out");

        CRH(c, sc + ".CRH" + to_string(i), out, arg1, arg2);
    }

    variable arg1 = c.get(sc + ".0.arg1");
    variable arg2 = c.get(sc + ".1.arg2");

    auto lc1 = (one * member + (-one) * arg1);
    auto lc2 = (one * member + (-one) * arg2);

    c.addC(lc1, lc2, LC());

    for (int i = 1; i < height; ++i) {
        variable lastout = c.get(sc + "." + to_string(i - 1) + ".out");
        variable arg1 = c.get(sc + "." + to_string(i) + ".arg1");
        variable arg2 = c.get(sc + "." + to_string(i) + ".arg2");

        auto lc1 = (one * lastout + (-one) * arg1);
        auto lc2 = (one * lastout + (-one) * arg2);

        c.addC(lc1, lc2, LC());
    }

    return;
}

void CoinCommitmentAndOwnership(Circuit &c, string sc, variable cm, variable w,
                                variable rho, variable a_sk, variable sn,
                                variable root) {

    variable a_pk = c.get(sc + ".a_pk");
    variable zero = c.get(sc + ".zero");
    MiMC(c, sc + ".PRFaddr", a_pk, zero, a_sk, constants::PRFaddr);
    c.replace_variable(zero, FieldR::zero(), true);

    variable r = c.get(sc + ".r");
    variable s = c.get(sc + ".s");

    variable k = c.get(sc + ".k");
    variable tmp1 = c.get(sc + ".tmp1");

    CRH(c, sc + ".crh1", tmp1, a_pk, rho);
    CRH(c, sc + ".crh2", k, r, tmp1);

    variable tmp2 = c.get(sc + ".tmp2");
    CRH(c, sc + "crh3", tmp2, w, k);
    CRH(c, sc + ".crh4", cm, s, tmp2);

    MiMC(c, sc + ".PRFsn", sn, rho, a_sk, constants::PRFsn);

    MerkleProof<29>(c, sc + ".merkleproof", root, cm);
}

void TSN(Circuit &c, string sc, variable tsn, variable rho, variable role,
         variable a_sk) {
    variable tmp = c.get(sc + ".tmp");
    CRH(c, sc + ".crh", tmp, rho, role);
    MiMC(c, sc + ".PRFtsn", tsn, tmp, a_sk, constants::PRFtsn);
}

void booleanConstrain(Circuit &c, variable x, int numBits) {
    LC lc = LC();

    FieldR twoPower = one;
    for (int i = 0; i < numBits; ++i) {
        variable b_i = c.get(x + ".bits" + to_string(i));
        auto lc1 = LC() + (one * b_i);
        auto lc2 = LC() + one * c.ONE + -one * b_i;
        c.addC(lc1, lc2, LC());

        lc = lc + twoPower * b_i;
        twoPower *= FieldR("2");
    }

    auto lc2 = LC() + one * c.ONE;
    auto lc3 = LC() + one * x;

    c.addC(lc, lc2, lc3);

    return;
}

void lessThanOrEqual(Circuit &c, string sc, variable l, variable r,
                     int numBits) {
    for (int i = 0; i < numBits; ++i) {
        variable t_i = c.get(sc + ".t_" + to_string(i));
        variable pi_i = c.get(sc + ".pi_" + to_string(i));

        variable d_i = c.get(sc + ".d" + to_string(i));

        variable c_i = c.get(r + ".bits" + to_string(i));
        variable a_i = c.get(l + ".bits" + to_string(i));

        c.addC(LC() + one * pi_i, LC() + one * a_i, LC() + one * d_i);

        if (i == numBits - 1) {
            c.addC(LC() + one * pi_i, LC() + one * t_i, LC() + one * c.ONE);
        } else {
            variable pi_iI = c.get(sc + ".pi_" + to_string(i + 1));
            c.addC(LC() + one * pi_i, LC() + one * t_i, LC() + one * pi_iI);
        }

        auto lc1 = LC() + one * c_i;
        auto lc2 = LC() + one * c.ONE + (-one) * a_i;
        auto lc3 = LC() + one * c.ONE + (-one) * t_i;

        c.addC(lc1, lc2, lc3);

        auto lc4 = LC() + one * c.ONE + (-one) * a_i;

        c.addC(lc4, LC() + one * d_i, LC());
    }

    return;
}

void notEqual(Circuit &c, string sc, variable x, variable y) {
    variable tmp = c.get(sc + ".tmp");
    auto lc1 = LC() + one * x + (-one) * y;
    auto lc2 = LC() + one * tmp;

    c.addC(lc1, lc2, LC() + one * c.ONE);
}

void SNnonmembership(Circuit &c, string sc, variable sn, variable root_sn) {
    variable sncontainer = c.get(sc + ".sncontainer");
    MerkleProof<29>(c, sc + ".snmerkleproof", root_sn, sncontainer);

    variable s_less = c.get(sc + ".sless");
    variable s_plus = c.get(sc + ".splus");

    CRH(c, sc + ".CRH", sncontainer, s_less, s_plus);

    booleanConstrain(c, s_less, 255);
    booleanConstrain(c, s_plus, 255);

    lessThanOrEqual(c, sc + "lessThan1", s_less, sn, 255);
    lessThanOrEqual(c, sc + "lessThan2", sn, s_plus, 255);
    notEqual(c, sc + "nequal1", s_less, sn);
    notEqual(c, sc + ".nequal2", sn, s_plus);
}

void proofSansBinomial(Circuit &c, string sc) {
    variable cm = c.get(sc + ".cm");
    variable w = c.get(sc + ".w");
    variable rho = c.get(sc + ".rho");
    variable a_sk = c.get(sc + ".a_sk");
    variable sn = c.get(sc + ".sn");

    variable root_cm = c.get_pub(sc + ".root_cm");
    variable root_sn = c.get_pub(sc + ".root_sn");
    variable tsn = c.get_pub(sc + ".tsn");
    variable role = c.get_pub(sc + ".role");

    CoinCommitmentAndOwnership(c, sc + ".cmSTUFFS", cm, w, rho, a_sk, sn,
                               root_cm);

    SNnonmembership(c, sc + ".snNonmembership", sn, root_sn);
    TSN(c, sc + "tsnproof", tsn, rho, role, a_sk);
}

void selectionFromBinomial(Circuit &c, string sc, variable w, variable role) {
    variable j = c.get(sc + ".j");
    booleanConstrain(c, j, 60);

    
}

#endif /* Circuit.hpp */
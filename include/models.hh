#pragma once

#include "Chili/Tools/FourVector.hh"
#include "LHAPDF/LHAPDF.h"

#include <complex>
#include <vector>

inline double LDot(const chili::FourVector &mom1, const chili::FourVector &mom2) {
    return mom1*mom2;
}

inline constexpr int LeviCivita(const int i, const int j, const int k, const int l) {
    return (i==j||i==k||i==l||j==k||j==l||k==l) ? 0 : -(i-j)*(i-k)*(i-l)*(j-k)*(j-l)*(k-l)/12;
}

// LC example: LTensor[LeviCivitaE, {p1}, {p2}, {p3}, {p4}]
inline double LeviCivitaE(const chili::FourVector &p1, const chili::FourVector &p2,
                          const chili::FourVector &p3, const chili::FourVector &p4) {
    double result = 0;
    for(size_t mu = 0; mu < 4; ++mu) {
        for(size_t nu = 0; nu < 4; ++nu) {
            for(size_t alpha = 0; alpha < 4; ++alpha) {
                for(size_t beta = 0; beta < 4; ++beta) {
                    int i = static_cast<int>(mu);
                    int j = static_cast<int>(nu);
                    int k = static_cast<int>(alpha);
                    int l = static_cast<int>(beta);
                    result += LeviCivita(i, j, k, l)*p1[mu]*p2[nu]*p3[alpha]*p4[beta];
                }
            }
        }
    }
    return result;
}

namespace LHAPDF { class PDF; }

namespace bcnutau {

class ModelBase {
  public:
    ModelBase(const std::string &pdf, double ecm=13000);
    virtual ~ModelBase() = default;

    double Evaluate(const std::vector<chili::FourVector> &mom);
    virtual double MatrixSquare(const std::vector<chili::FourVector> &mom) = 0;
    
  protected:
    std::complex<double> Propagator(const chili::FourVector &mom, double mass, double width) const {
        return mom.Mass2() - mass*mass + std::complex<double>(0, 1)*mass*width;
    }

    std::unique_ptr<LHAPDF::PDF> m_pdf;
    double m_ecm;
};

class StandardModel : public ModelBase {
  public:
    // TODO: Add handling of parameters. Right now they are hard coded
    StandardModel(const std::string &pdf) : ModelBase(pdf) {}

    double MatrixSquare(const std::vector<chili::FourVector> &mom) override;

  private:
    static constexpr double mass_b = 4.3;
    static constexpr double mass_w = 80.385;
    static constexpr double width_w = 2.035;
    static constexpr double vcb2 = 1;
    static constexpr double alpha = 1.0/137.0;
};

}

#pragma once

#include "Chili/Tools/FourVector.hh"
#include "LHAPDF/LHAPDF.h"

#include <complex>
#include <vector>

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

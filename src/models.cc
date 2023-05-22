#include "models.hh"

// Notes: Conversion of momentum from Rodolfo's notes to Chili convention
// p_a = -mom[0]
// p_b = -mom[1]
// p_1 = mom[2]
// p_2 = mom[3]
// p_3 = mom[4]

using bcnutau::ModelBase;
using bcnutau::StandardModel;

ModelBase::ModelBase(const std::string &pdf, double ecm) : m_ecm{ecm} {
    m_pdf = std::unique_ptr<LHAPDF::PDF>(LHAPDF::mkPDF(pdf, 0));
}

double ModelBase::Evaluate(const std::vector<chili::FourVector> &mom) {
    // TODO: Do we want a dynamic factorization scale? If so, what should it be
    const double scale2 = 91.18*91.18;
    double x1 = -mom[0].E()/m_ecm;
    double x2 = -mom[1].E()/m_ecm;
    if(x1 < 0 || x1 > 1) return 0;
    if(x2 < 0 || x2 > 1) return 0;
    double xfxa1 = m_pdf->xfxQ2(21, x1, scale2); // 21 is the pid for gluon
    double xfxa2 = m_pdf->xfxQ2(4, x2, scale2);  // 4 is the pid for charm
    double xfxb1 = m_pdf->xfxQ2(21, x2, scale2); // 21 is the pid for gluon
    double xfxb2 = m_pdf->xfxQ2(4, x1, scale2);  // 4 is the pid for charm
    return (xfxa1/x1*xfxa2/x2 + xfxb1/x1*xfxb2/x2)*MatrixSquare(mom);
}

double StandardModel::MatrixSquare(const std::vector<chili::FourVector> &mom) {
    // TODO: Do we want a dynamic renormalization scale? If so, what should it be
    const double scale2 = 91.18*91.18;
    const double alphas = m_pdf->alphasQ2(scale2);
    constexpr double sinw2 = 0.23;
    double coeff = 4.0*alphas*alpha/sinw2*alpha/sinw2*pow(M_PI, 3)/6.0*vcb2;
    double prop_b = Propagator(mom[2]-mom[0], mass_b, 0).real();
    double prop_w = std::norm(Propagator(mom[3]+mom[4], mass_w, width_w));
    double s = (mom[0]+mom[1]).Mass2();

    double prefactor1 = coeff/prop_b/prop_b/prop_w;
    double prefactor2 = coeff/s/s/prop_w;
    double prefactor3 = coeff/prop_b/s/prop_w;

    double term1 = -prefactor1*64*(mom[4]*mom[1])*
        ((mom[2]*mom[0])*(mom[3]*mom[0])-mass_b*mass_b*mom[3]*(mom[2]+mom[0]));
    double term2 = -prefactor2*32*s*(mom[2]*mom[3])*(mom[4]*mom[0]);
    double term3 = prefactor3*64*((mom[2]*mom[3])*(s/2*mom[4]*(mom[2]-mom[1])-(mom[2]*mom[1])*(mom[4]*mom[0])
                - (mom[2]*mom[0])*(mom[4]*mom[1]))
            - (mom[4]*mom[1])*((mom[2]*mom[1])*(mom[3]*mom[0])+2*(mom[2]*mom[3])*(mom[2]*mom[1])
                - (mom[2]*mom[0])*(mom[3]*mom[1])));

    return term1+term2+term3;
}

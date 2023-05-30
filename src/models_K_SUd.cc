#include "models.hh"
#include <complex>

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


double BSM_S1d::MatrixSquare(const std::vector<chili::FourVector> &mom) {
    // TODO: Do we want a dynamic renormalization scale? If so, what should it be
    auto pa = -mom[0]
    auto pb = -mom[1]
    auto p1 = mom[2]
    auto p2 = mom[3]
    auto p3 = mom[4]
    const double scale2 = 91.18*91.18;
    const double alphas = m_pdf->alphasQ2(scale2);
    constexpr double sinw2 = 0.23;
    double prop_bsm_w = sinw2(-mass_w*mass_w+2*(mom[3]*mom[4])+ width_w*mass_w); 

    auto coeff = 4.0*alpha*gL33*alphas*(gL33*Vcb + gL32*Vcs)*std::conj(Vcb); 
	    
	    //need to define gL33, gR32, gL32, width_S1 and mass_S1

    auto coeff2 = gL33*gL33*4.0*alphas*((gR32*gR32 + (gL33*Vcb + gL32*Vcs)*(gL33*std::conj(Vcb) 
	        + gL32*std::conj(Vcs))).real();


    auto coeff3 = (gL33*gL33*4.0*alphas*(-2*(gL33*Vcb + gL32*Vcs)*(gL33*std::conj(Vcb) 
	          gL32*std::conj(Vcs));

    auto term13 = coeff*LDot(p3, pb)*(2*(2*mass_b*mass_b + LDot(p1, p3)
		    - LDot(p1, pb))*(LDot(p2, p3) - LDot(p2, pb)) 
	            + LDot(p1, p2)*(mass_b*mass_b + 2*LDot(p3, pb)))/
                    (3*sinw2 *Prop(p2 + p3, mass_w, width_w)*Prop(-p3 + pb, mass_S1, 0)*
			      pow(Prop(-p2 - p3 + pb, mass_b, 0),2));


    auto term14 = coeff*(pow(-2*LDot(p1, pb),2)*LDot(p2, p3)
		     - 2*pow(LDot(p1, p3),2)*LDot(p2, pb) 
		     - (LDot(p1, p2)*(mass_b*mass_b - 2*LDot(p2, p3)) 
		     + 2*(mass_b*mass_b + LDot(p1, p2))*LDot(p2, pb))*LDot(p3, pb)
		     - LDot(p1, p3)*((mass_b*mass_b + 2*LDot(p2, p3) 
		     - 2*LDot(p2, pb))*LDot(p2, pb) - 2*LDot(p1, p2)*LDot(p3, pb) 
		     + (2*std::complex<double>(0, 1))*LeviCivitaE(p1,p2,p3,pb) 
		     + LDot(p1, pb)*(LDot(p2, p3)*(mass_b*mass_b + 2*LDot(p1, p3) + 2*LDot(p2, p3)) 
		     + 2*(LDot(p1, p3) - LDot(p2, p3))*LDot(p2, pb) + 2*LDot(p1, p2)*LDot(p3, pb)
		     + (2*std::complex<double>(0, 1))*LeviCivitaE(p1, p2,p3,pb)
		     - std::complex<double>(0, 1)*(mass_b*mass_b + 2*LDot(p2, p3) - 2*LDot(p2, pb))
		     *LeviCivitaE( p1, p2, p3, pb))/
		     (6*sinw2*Prop(p1 + p2, mass_S1, width_S1)* Prop(p2 + p3, mass_w, width_w)
		     *Prop(p1 + p2 + p3, 0, 0)*Prop(-p2 - p3 + pb, mass_b, width_b));

    
    auto term15 = -coeff*(2*mass_b*mass_b*(-LDot(p2, p3) + LDot(p2, pb)) 
	          + LDot(p1, p2)*(mass_b*mass_b - 2*LDot(p2, p3) + 2*LDot(p2, pb) 
	          - 2*LDot(p3, pb)))*LDot(p3, pb))/
		    (6*sinw2*Prop(p1 + p2, mass_S1,width_S1)*  Prop(p2 + p3, mass_w, width_w)
		  * Prop(-p3 + pb, mass_S1, width_S1)* Prop(-p2 - p3 + pb, mass_b, width_b));




    auto term23 = coeff*pow(-2*LDot(p1, pb,2)*LDot(p2, p3) - 2*pow(LDot(p1, p3),2)*LDot(p2, pb)
		  - (LDot(p1, p2)*(mass_b*mass_b - 2*LDot(p2, p3)) + 2*(mass_b*mass_b + LDot(p1, p2))
		  * LDot(p2, pb))*LDot(p3, pb) + LDot(p1, pb)
		  * (LDot(p2, p3) * (mass_b*mass_b + 2*LDot(p1, p3) + 2*LDot(p2, p3))   
	          + 2*(LDot(p1, p3) - LDot(p2, p3))*LDot(p2, pb) +2*LDot(p1, p2)*LDot(p3, pb)
		  - (2*std::complex<double>(0, 1))*LeviCivitaE(p1,p2,p3,pb) 
		  + LDot(p1, p3)*(LDot(p2, pb)*(-mass_b*mass_b - 2*LDot(p2, p3) 
		  + 2*LDot(p2, pb)) + 2*LDot(p1, p2)*LDot(p3, pb) 
	          + (2*std::complex<double>(0, 1))*LeviCivitaE(p1, p2,p3,pb) 
		  + std::complex<double>(0, 1) *(mass_b*mass_b + 2*LDot(p2, p3) - 2*LDot(p2, pb))
		  * LeviCivitaE(p1, p2, p3, pb)) /
		  (6*sinw2*Prop(p2 + p3, mass_w, width_w)*Prop(p1 + p2 + p3, 0, 0)
		  * Prop(-p3 + pb, mass_S1, width_S1)*Prop(-p2 - p3 + pb, mass_b, width_b));

     
   auto term24 = coeff*LDot(p1, p2)*(2*(LDot(p1, p3) + LDot(p2, p3))*(LDot(p1, pb)
	       + LDot(p2, pb)) - (mass_b*mass_b + 2*LDot(p1, p2))*LDot(p3, pb))) /
                 (3*sinw2*Prop(p1 + p2, mass_S1,width_S1)*Prop(p2 + p3, mass_w,width_w)
	       * pow(Prop(p1 + p2 + p3, 0, 0),2));


   auto term25 = - coeff* *LDot(p1, p2)*(mass_b*mass_b + 2*LDot(p1, p2) 
		 - 2*LDot(p1, p3) - 2*LDot(p2, p3))*LDot(p3, pb)) /
		  (6*sinw2*Prop(p1 + p2, mass_S1,width_S1)*Prop(p2 + p3, mass_w,width_w]
		 * Prop(p1 + p2 + p3, 0, 0)*Prop(-p3 + pb, mass_S1, width_S1)); 


   auto term33 = coeff2*LDot(p3, pb)*(2*(2*mass_b*mass_b + LDot(p1, p3) - LDot(p1, pb))
	       * (LDot(p2, p3) - LDot(p2, pb)) + LDot(p1, p2)*(mass_b*mass_b + 2*LDot(p3, pb)))) /
		 (3*pow(Prop(-p3 + pb, mass_S1, width_S1),2)* pow(Prop(-p2 - p3 + pb, mb, width_b),2));


   auto term34 = coeff3*pow(LDot(p1, pb),2)*LDot(p2, p3) 
	       - 2*(gL33*Vcb + gL32*Vcs)*(gL33* + gL32*std::conj(Vcs))*pow(LDot(p1, p3),2)*LDot(p2, pb)
	       - LDot(p1, p3)* (gR32*gR32*mass_b*mass_b + (gL33*Vcb + gL32*Vcs)*(gL33*std::conj(Vcb) 
	         gL32*std::conj(Vcs))*(mass_b*mass_b + 2*LDot(p2, p3)))*LDot(p2, pb) 
	       + 2*(gL33*Vcb + gL32*Vcs)*(gL33*std::conj(Vcb) + gL32*std::conj(Vcs))
	         LDot(p1, p3)*pow(LDot(p2, pb),2) + 2*(gL33*Vcb + gL32*Vcs)
               * (gL33*std::conj(Vcb) + gL32*std::conj(Vcs))*LDot(p1, p2)*LDot(p1, p3)
	       * LDot(p3, pb) + 2*gR32*gL32*mass_b*mass_b*LDot(p2, p3)*LDot(p3, pb) 
	       - 2*gL33*gL33*mass_b*mass_b*Vcb*std::conj(Vcb)*LDot(p2, pb)*LDot(p3, pb) 
	       - 2*gL32*gL33*mass_b*mass_b*Vcs*std::conj(Vcb)*LDot(p2, pb)*LDot(p3, pb) 
	       - 2*gL32*gL33*mass_b*mass_b*Vcb*std::conj(Vcs)*LDot(p2, pb)*LDot(p3, pb) 
	       - 2*gL32*gL32*mass_b*mass_b*Vcs*std::conj(Vcs)*LDot(p2, pb)*LDot(p3, pb) 
	       - LDot(p1, p2)*(gR32*gR32*(mass_b*mass_b - 4*LDot(p2, p3)) + (gL33*Vcb + gL32*Vcs)
	       * (gL33*std::conj(Vcb) + gL32*std::conj(Vcs))*(mass_b*mass_b - 2*LDot(p2, p3) 
	       + 2*LDot(p2, pb)))*LDot(p3, pb) + LDot(p1, pb)
	       * (2*(gL33*Vcb + gL32*Vcs)*(gL33*std::conj(Vcb) + gL32*std::conj(Vcs))
	       * pow(LDot(p2, p3),2) + LDot(p2, p3)*(gR32*gR32*mass_b*mass_b 
	       + (gL33*Vcb + gL32*Vcs)*(gL33*std::conj(Vcb) + gL32*std::conj(Vcs))
	       * (mass_b*mass_b + 2*LDot(p1, p3) - 2*LDot(p2, pb))) + 2*(gL33*Vcb + gL32*Vcs)
	       * (gL33*std::conj(Vcb)+gL32*std::conj(Vcs))*(LDot(p1, p3)*LDot(p2, pb) 
	       + LDot(p1, p2)*LDot(p3, pb) + std::complex<double>(0, 1)*LeviCivitaE(p1,p2,p3,pb)) 
	       + std::complex<double>(0, 1)*(gR32*gR32*mass_b*mass_b - (gL33*Vcb + gL32*Vcs)
	       * (gL33*std::conj(Vcb) + gL32*std::conj(Vcs))*(mass_b*mass_b+ 2*LDot(p1, p3) 
	       + 2*LDot(p2, p3) - 2*LDot(p2, pb)))*LeviCivitaE(p1,p2,p3,pb)) /
	         (6*Prop(p1 + p2, mass_S1, width_S1)* Prop(p1 + p2 + p3, 0, 0)
	       * Prop(-p3 + pb, mass_S1, width_S1)*Prop(-p2 - p3 + pb, mass_b, width_b));
  


    auto term35 = -coeff2* *(2*mass_b*mass_b*(-LDot(p2, p3) + LDot(p2, pb)) 
	        + LDot(p1, p2)*(mass_b*mass_b - 2*LDot(p2, p3) + 2*LDot(p2, pb) - 2*LDot(p3, pb)))
	        * LDot(p3, pb)) /
		  (6*Prop(p1 + p2, mass_S1, wodth_S1]*pow(Prop(-p3 + pb, mass_S1, width_S1),2)
	        * Prop(-p2 - p3 + pb, mass_b, width_b));


    auto term44 = coeff2* LDot(p1, p2)*(2*(LDot(p1, p3)+ LDot(p2, p3))*(LDot(p1, pb) + LDot(p2, pb)) 
	        - (mass_b*mass_b + 2*LDot(p1, p2))*LDot(p3, pb))) /
                  (3*pow(Prop(p1 + p2, mass_S1, width_S1),2)*pow(Prop(p1 + p2 + p3, 0, 0),2));



    auto term45 = -coeff2*LDot(p1, p2)*(mass_b*mass_b + 2*LDot(p1, p2) 
		  -2*LDot(p1, p3) - 2*LDot(p2, p3))*LDot(p3, pb) /
		   (6*pow(Prop(p1 + p2, mass_S1, width_S1],2)*Prop(p1 + p2 + p3, 0, 0)
	          * Prop(-p3 + pb, mass_S1, width_S1)));
	


     auto term55 = -coeff2*LDot(p1, p2)*(mass_b*mass_b + 2*LDot(p1, p2) 
	         - 2*LDot(p1, p3) + 2*LDot(p1, pb) - 2*LDot(p2, p3) + 2*LDot(p2, pb) 
		 - 2*LDot(p3, pb))*LDot(p3, pb) /
		   (6*pow(Prop(p1 + p2, mass_S1,width_S1],2)* pow(Prop(-p3 + pb, mass_S1, width_S1),2));



    return 2*term13.real()+2*term14.real()+2*term15.real()+2*term23.real()+2*term24.real()+2*term25.real()+term33+2*term34.real()+2*term35.real()+term44+2*term45.real()+term55;


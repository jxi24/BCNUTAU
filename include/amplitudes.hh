#pragma once

#include "Chili/Tools/FourVector.hh"

#include <complex>
#include <vector>

namespace bcnutau {

std::complex<double> amp_d0(const std::vector<chili::FourVector> &mom);
std::complex<double> amp_d1(const std::vector<chili::FourVector> &mom);
std::complex<double> amp_d2(const std::vector<chili::FourVector> &mom);
std::complex<double> amp_d3(const std::vector<chili::FourVector> &mom);
std::complex<double> amp_d4(const std::vector<chili::FourVector> &mom);

}

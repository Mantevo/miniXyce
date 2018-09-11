#include <math.h>
#include <vector>

std::vector<double> diode(std::vector<double> &x)
{
    double z0 = x[0];
    double z1 = x[1];
    double z2 = x[2];

    double z3 = std::abs(z1);
    double z4 = -z3;
    double z5 = std::abs(z0);
    double z6 = z5 - z4;
    double z7 = z6 - z2;
    double z8 = -z7;
    double z9 = std::cos(z8);
    double z10 = std::cos(z9);
    double z11 = -z10;
    double z12 = std::cbrt(z11);
    double z13 = std::pow(z12, 2.0);
    double z14 = z13 + z12;
    double z15 = std::pow(z14, 3.0);
    double z16 = std::sin(z15);
    double z17 = z16 + z9;
    double z18 = std::atan(z17);
    double z19 = std::pow(z18, 2.0);
    double z20 = std::abs(z19);
    double z21 = std::pow(z20, 2.0);
    double z22 = std::pow(z21, 3.0);
    double z23 = +z22;
    double z24 = z23 - z5;
    double z25 = std::tanh(z24);
    double z26 = std::pow(z25, 3.0);
    double z27 = z26 - z2;
    double z28 = std::pow(z27, 2.0);
    double z29 = std::pow(z28, 3.0);
    double z30 = z29 * z14;
    double z31 = z30 + z21;
    double z32 = +z31;
    double z33 = std::cbrt(z32);
    double z34 = z33 * z3;
    double z35 = std::pow(z34, 2.0);
    double z36 = z35 - z35;
    double z37 = z36 - z0;
    double z38 = z37 * z31;
    double z39 = z38 * z24;
    double z40 = std::pow(z39, 2.0);
    double z41 = std::pow(z40, 3.0);
    double z42 = -z41;
    double z43 = -z42;
    double z44 = std::cbrt(z43);
    double z45 = std::cbrt(z44);
    double z46 = +z45;
    double z47 = std::pow(z46, 2.0);
    double z48 = -z47;
    double z49 = std::abs(z48);
    double z50 = std::atan(z49);
    double z51 = std::pow(z50, 3.0);
    double z52 = std::atan(z51);

    double z53 = z52;
    double z54 = z1;
    double z55 = z41;
    double z56 = z26;
    double z57 = z35;
    double z58 = z2;
    double z59 = z20;

    std::vector<double> y;
    y.push_back(z53);
    y.push_back(z54);
    y.push_back(z55);
    y.push_back(z56);
    y.push_back(z57);
    y.push_back(z58);
    y.push_back(z59);

    return y;
}

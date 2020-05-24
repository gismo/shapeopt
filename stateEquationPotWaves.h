#ifndef GSSTATEEQUATIONPOTWAVES_H
#define GSSTATEEQUATIONPOTWAVES_H

#include <math.h>
#include <complex>
#include <sstream>
#include <string>
using namespace gismo;

// FIXIT: Clean up. Make more flexible wrt to topology, e.g. use Piecewise function
class gsStateEquationPotWaves{
public:
    gsStateEquationPotWaves(gsMultiPatch<>* mpin, index_t numRefine);

    gsMultiPatch<> solve();

public:
    gsMultiPatch<>* mp;
    index_t degree = 2;
    bool isRefined = false;
    index_t numRef = 1;

    // Assembler variables
    gsMultiBasis<> dbasis;
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    gsSparseSolver<>::LU solver;

    // Wave parameters
    real_t wave_K       = 1;    // Wave number

    // PML parameters
    real_t pml_n        = 2;    // Exponent
    real_t pml_length   = 1;    // width of pml layer
    real_t pml_start    = 1;    // start of pml layer
    real_t pml_C        = 5;    // constant

    gsFunctionExpr<> strech_x_re;
    gsFunctionExpr<> strech_x_im;

    gsFunctionExpr<> strech_y_re;
    gsFunctionExpr<> strech_y_im;

    gsFunctionExpr<> strech_z_re;
    gsFunctionExpr<> strech_z_im;

    gsFunctionExpr<> dstrech_x_re;
    gsFunctionExpr<> dstrech_x_im;

    gsFunctionExpr<> dstrech_y_re;
    gsFunctionExpr<> dstrech_y_im;

    gsFunctionExpr<> dstrech_z_re;
    gsFunctionExpr<> dstrech_z_im;

    gsFunctionExpr<> dstrech_xinv_re;
    gsFunctionExpr<> dstrech_xinv_im;

    gsFunctionExpr<> dstrech_yinv_re;
    gsFunctionExpr<> dstrech_yinv_im;

    gsFunctionExpr<> dstrech_zinv_re;
    gsFunctionExpr<> dstrech_zinv_im;


};

#endif //GSSTATEEQUATIONPOTWAVES_H

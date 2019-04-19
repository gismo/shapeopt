/** @file main.cpp

    @brief Test of calculating derivatives of detJ

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Asger Limkilde
*/

#include <gismo.h>

using namespace gismo;

int main(int argc, char* argv[])
{

    gsMultiPatch<> mp = gsMultiPatch<>(*gsNurbsCreator<>::BSplineSquare(2));
    mp.basis(0).setDegree(2);

    gsMultiBasis<> dbasis(mp);

    gsExprAssembler<> A(1,1);
    gsExprEvaluator<> ev(A);

    A.setIntegrationElements(dbasis);

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    geometryMap G = A.getMap(mp);

    variable g = A.getCoeff(mp);

    gsFunctionExpr<> ff("x**2","y**2",2);
    variable e = A.getCoeff(ff);

    gsFuncCoordinate<real_t> g1_fc(mp,0);

    variable g1 = A.getCoeff(g1_fc);

    gsFunctionExpr<> ff1("1","1",2);
    variable e1 = A.getCoeff(ff1);

    // gsInfo << ev.integral(e.tr()*g) << "\n";
    gsInfo << ev.integral(hess(g1).sqNorm()) << "\n";



    gsInfo << "\n\n==== DONE ====\n";
    return 0;
}

#include <gismo.h>
#include "gsStateEquationAntenna.h"
#include "iostream"
#include "fstream"
using namespace gismo;


gsMatrix<> gsStateEquation::getRhsZeroBC(index_t realOrImag){
	gsSparseMatrix<> mat;
	gsVector<> rhs;
	getTerm(realOrImag,mat,rhs);
	return rhs;
};

/* gsMatrix<> gsStateEquationAntenna::getKu(gsMultiPatch<> sol){
 	// gsInfo << "f = " << *f << "\n" << std::flush;
 	// gsInfo << "ms = " << *ms << "\n" << std::flush;

 	gsFunctionExpr<> gN("0.0",2);
 	gsFunctionExpr<> zero("0.0",2);
   //! [Boundary conditions]

 	gsExprAssembler<> A(1,1);
   geometryMap G = A.getMap(*m_mp);

 	A.setIntegrationElements(dbasis);


   space u = A.getSpace(dbasis);
   u.setInterfaceCont(0);
   u.addBc(bcInfoZero.get("Dirichlet"));

 	variable solVar = A.getCoeff(sol);

 	A.initSystem();

 	A.assemble(igrad(u,G)*(jac(G).inv().tr()*fjac(solVar))*meas(G));

   // delete f;
 	// delete ms;
 	// delete df_dx;
 	// delete df_dy;
 	gsInfo << "Return Ku \n";
 	return A.rhs();

 }
*/

gsMatrix<> gsStateEquation::getDerivativeOfU(){
	// gsInfo << "Get derivative of u\n";
	// gsInfo << "Get sol\n";
	gsMultiPatch<> u_real, u_imag;

	solve(u_real,u_imag);

	// gsInfo << "Get deriv of Ku\n";
	gsMatrix<> dK_realu_real = getDerivativeOfAu(0,u_real);
	gsMatrix<> dK_imagu_real = getDerivativeOfAu(1,u_real);
	gsMatrix<> dK_realu_imag = getDerivativeOfAu(0,u_imag);
	gsMatrix<> dK_imagu_imag = getDerivativeOfAu(1,u_imag);
	// gsInfo << "Get deriv of F\n";
	gsMatrix<> dF_real = getDerivativeOfRhsZeroBC(0);
	gsMatrix<> dF_imag = getDerivativeOfRhsZeroBC(1);

	gsMatrix<> dF(dF_real.rows(),dF_real.cols()*2);
	dF << dF_real, -dF_imag;

	gsMatrix<> dAu(dF_real.rows(),dF_real.cols()*2);
	dAu << dK_realu_real - dK_imagu_imag, -dK_imagu_real - dK_realu_imag;

	// gsInfo << "solve\n";
	gsMatrix<> du = solver.solve(dF.transpose() - dAu.transpose());

	// gsInfo << "return du\n";
	return du.transpose();
}

// Get the rhs of the sytem you need to solve to obtain dudc, without solving
gsMatrix<> gsStateEquation::getDerivativeWithoutSolving(gsMultiPatch<> &u_real, gsMultiPatch<> &u_imag){
	// gsInfo << "Get derivative of u\n";
	// gsInfo << "Get sol\n";

	// gsInfo << "Get deriv of Ku\n";
	gsMatrix<> dK_realu_real = getDerivativeOfAu(0,u_real);
	gsMatrix<> dK_imagu_real = getDerivativeOfAu(1,u_real);
	gsMatrix<> dK_realu_imag = getDerivativeOfAu(0,u_imag);
	gsMatrix<> dK_imagu_imag = getDerivativeOfAu(1,u_imag);
	// gsInfo << "Get deriv of F\n";
	gsMatrix<> dF_real = getDerivativeOfRhsZeroBC(0);
	gsMatrix<> dF_imag = getDerivativeOfRhsZeroBC(1);

	gsMatrix<> dF(dF_real.rows(),dF_real.cols()*2);
	dF << dF_real, -dF_imag;

	gsMatrix<> dAu(dF_real.rows(),dF_real.cols()*2);
	dAu << dK_realu_real - dK_imagu_imag, -dK_imagu_real - dK_realu_imag;

	// Return matrix

	return dF.transpose() - dAu.transpose();

}

// Assumes the system is already factorized, e.g. if you have called getDerivativeWithoutSolving firts;
gsVector<> gsStateEquation::solveAdjoint(gsVector<> &rhs){
	return solver.solve(rhs);
}

void gsStateEquation::printMatSize(gsMatrix<> mat, std::string name){
	gsInfo << "Size of " << name << ":\t (" << mat.rows() << ", " << mat.cols() << ")\n";
}

void gsStateEquation::plotSolution(gsMultiPatch<> &sol, std::string name){
	gsExprAssembler<> A(1,1);
	gsExprEvaluator<> ev(A);

  A.setIntegrationElements(dbasis);

  geometryMap G = A.getMap(*m_mp);

	variable u_sol = A.getCoeff(sol);

	gsInfo<<"Plotting " << name << " in Paraview...\n";
	// ev.options().setSwitch("plot.elements", true);
	ev.writeParaview( u_sol   , G, name);
}

void gsStateEquation::plotSolution(std::string name){
	gsMultiPatch<> ur,ui;
	solve(ur,ui);
	plotSolution(ur,name + "_re");
	plotSolution(ui,name + "_im");
}

void gsStateEquation::plotMagnitude(std::string name){
	gsMultiPatch<> u_real, u_imag;
	solve(u_real,u_imag);

    plotMagnitude(u_real, u_imag, name);
}

void gsStateEquation::plotMagnitude(gsMultiPatch<> &u_real, gsMultiPatch<> &u_imag, std::string name){

	gsExprAssembler<> A(1,1);
	gsExprEvaluator<> ev(A);
    A.setIntegrationElements(dbasis);

    geometryMap G = A.getMap(*m_mp);

	variable u_re = A.getCoeff(u_real);
	variable u_im = A.getCoeff(u_imag);

	gsInfo<<"Plotting " << name << " in Paraview...\n";
	// ev.options().setSwitch("plot.elements", true);
	ev.writeParaview( u_re*u_re + u_im*u_im   , G, name);
}

gsMultiPatch<> gsStateEquation::getPieceWiseFunctionOnPatch(index_t nBoxes, index_t patch, real_t val_on_patch, real_t val_elsewhere){
	gsKnotVector<> u_knots(0,1,0,1);
	gsKnotVector<> v_knots(0,1,0,1);


	gsTensorBSplineBasis<2> T_basis(u_knots,v_knots);

	gsMatrix<> cf_patch(1,1);
	cf_patch << val_on_patch;

	gsMatrix<> cf_elsewhere(1,1);
	cf_elsewhere << val_elsewhere;

	gsMultiPatch<> pw;

	for(index_t i = 0; i < nBoxes; i ++){
		if (i == patch){
  		memory::unique_ptr<gsGeometry<> > b = T_basis.makeGeometry(cf_patch);
			pw.addPatch(*b);
		} else {
  		memory::unique_ptr<gsGeometry<> > b = T_basis.makeGeometry(cf_elsewhere);
			pw.addPatch(*b);
		}
	}

	return pw;
}

void gsStateEquation::getSystem(gsSparseMatrix<> &mat, gsVector<> &rhs){
	gsSparseMatrix<> matReal,matImag;
	gsVector<> rhsReal, rhsImag;

	// Get real parts
	getTerm(0,matReal,rhsReal);
	// gsInfo << "Save Ar \n" << std::flush;
	// std::ofstream file1("Ar.txt");
	// file1 << matReal.toDense();
	// file1.close();
	// gsInfo << "norm of matReal: " << matReal.norm() << "\n" ;

	// Get imaginary parts
	getTerm(1,matImag,rhsImag);
	// gsInfo << "Save Ai \n" << std::flush;
	// std::ofstream file("Ai.txt");
	// file << matImag.toDense();
	// file.close();
	// gsInfo << "norm of matImag: " << matImag.norm() << "\n" ;

	index_t rows = matReal.rows();
	index_t cols = matReal.cols();
	// gsInfo << "DoFs :" << rows*2 << "\n";
	gsSparseMatrix<> A(rows*2,cols*2);
	A.reserve(matReal.nonZeros()*2 + matImag.nonZeros()*2);
	for(index_t c = 0; c < cols; ++c)
	{
	    A.startVec(c); // Important: Must be called once for each column before inserting!

	    for(gsSparseMatrix<>::InnerIterator itR(matReal, c); itR; ++itR){
	         A.insertBack(itR.row(), c) = itR.value();
				 }
	    for(gsSparseMatrix<>::InnerIterator itI(matImag, c); itI; ++itI){
	         A.insertBack(itI.row()+rows, c) = -itI.value();
				 }
	}
	// gsInfo << " next loop...sss\n" << std::flush;
	for(index_t c = 0; c < cols; ++c)
	{
	    A.startVec(c+cols); // Important: Must be called once for each column before inserting!

	    for(gsSparseMatrix<>::InnerIterator itI(matImag, c); itI; ++itI){
	         A.insertBack(itI.row(), c + cols) = -itI.value();
				 }
	    for(gsSparseMatrix<>::InnerIterator itR(matReal, c); itR; ++itR){
	         A.insertBack(itR.row() + rows, c + cols) = -itR.value();
				 }
	}
	// gsInfo << "Finalize A \n" << std::flush;
	A.finalize();



	mat = A;
	rhs.setZero(cols*2);
	rhs << rhsReal, -rhsImag;

    /*std::ofstream file_A;
    gsMatrix<> tmp = A;
    file_A.open("A.txt");
    file_A << tmp;
    file_A.close();

    std::ofstream file_f;
    file_f.open("f.txt");
    file_f << rhs;
    file_f.close();
    */

}

void gsStateEquation::solve(gsMultiPatch<> &u_real, gsMultiPatch<> &u_imag){
    gsVector<> tmp1, tmp2;
    solve(u_real, u_imag, tmp1, tmp2);
}


void gsStateEquation::solve(gsMultiPatch<> &u_real, gsMultiPatch<> &u_imag, gsVector<> &solVec_real, gsVector<> &solVec_imag){
		gsSparseMatrix<> mat;
		gsVector<> rhs;
		getSystem(mat,rhs);

		// gsInfo << "norm of mat: " << mat.norm() << "\n" ;
		// gsInfo << "norm of rhs: " << rhs.norm() << "\n" ;

		// gsInfo << "Save A \n" << std::flush;
		// std::ofstream file("A.txt");
		// file << mat.toDense();
		// file.close();

		gsInfo << "." << std::flush;
		solver.compute(mat);
		solVector = solver.solve(rhs);


		// Plot
		gsExprAssembler<> A(1,1);
		gsExprEvaluator<> ev(A);


		geometryMap G = A.getMap(*m_mp);

		A.setIntegrationElements(dbasis);

		space u = A.getSpace(dbasis);
		u.setInterfaceCont(0);
		// u.addBc(bcInfo.get("Dirichlet"));

		A.initSystem();

		gsMatrix<> solVector_Real = solVector.block(0,0,mat.cols()/2,1);
		gsMatrix<> solVector_Imag = solVector.block(mat.cols()/2,0,mat.cols()/2,1);

		solution u_sol_real = A.getSolution(u,solVector_Real);

		u_sol_real.extract(u_real);

		solution u_sol_imag = A.getSolution(u,solVector_Imag);

		u_sol_imag.extract(u_imag);

        solVec_real = solVector_Real;
        solVec_imag = solVector_Imag;

        // Save solutions. FIXIT: avoid copying data
        m_ur = u_real;
        m_ui = u_imag;
}

gsMatrix<> gsStateEquation::getU(index_t realOrImag){
		gsSparseMatrix<> mat;
		gsVector<> rhs;
		getSystem(mat,rhs);

		gsInfo << "." << std::flush;
		solver.compute(mat);
		solVector = solver.solve(rhs);

		index_t n = mat.cols()/2;

		gsMatrix<> U_Real = solVector.block(0,0,n,1);
		gsMatrix<> U_Imag = solVector.block(n,0,n,1);

		if (realOrImag == 0){
			return U_Real;
		} else {
			return U_Imag;
		}
}

gsMatrix<> gsStateEquation::getU(){
		gsSparseMatrix<> mat;
		gsVector<> rhs;
		getSystem(mat,rhs);

		gsInfo << "." << std::flush;
		solver.compute(mat);
		solVector = solver.solve(rhs);

		return solVector;
}

gsMatrix<> gsStateEquation::getAu(index_t realOrImag, gsMatrix<> U){
		gsSparseMatrix<> mat;
		gsVector<> rhs;
		getSystem(mat,rhs);

		index_t n = mat.cols()/2;

		if (realOrImag == 0){
			return mat.block(0,0,n,n)*U;
		} else {
			return mat.block(n,0,n,n)*U;
		}
}

void gsStateEquation::writeToFile(gsMatrix<> mat, std::string name) const{
    gsInfo << "WRITING to " << name << "\n";
    std::ofstream f(name);
    for(index_t i = 0; i < mat.rows(); i++){
        for(index_t j = 0; j < mat.cols(); j++){
            f << std::setprecision(20) << mat(i,j) << " ";
        }
        f << "\n";
    }
}

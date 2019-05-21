#include <gismo.h>
#include <stdio.h>
#include <math.h>       /* pow */
#include <iostream>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <ctime>

#include "gsDetJacConstraint.h"
#include "gsParamMethod.h"
#include "gsSpringMethod.h"
#include "gsModLiao.h"
#include "gsWinslow.h"
#include "gsAffineOptParamMethod.h"
#include "gsIpOptSparseMatrix.h"
#include "gsShapeOptProblem.h"
#include "gsOptAntenna.h"
#include "gsShapeOptLog.h"

using namespace gismo;

void readFromTxt(std::string name, gsMatrix<> &matrix){
	std::ifstream infile;
	infile.open(name);
	gsInfo << "Loading from " << name << "\n";
	for(int i = 0; i < matrix.rows(); i++){
		for(int j = 0; j < matrix.cols(); j++){
			infile >> matrix(i,j);
		}
	}
	infile.close();
}

void loadCoefs(gsMultiPatch<> &mp, size_t i, std::string folder, std::string name){

	gsMatrix<> cc = mp.patch(i).coefs();
	readFromTxt(BASE_FOLDER "/parametrizations/"+folder+name, cc);

	mp.patch(i).setCoefs(cc);
}

void saveVec(gsVector<> &vec, std::string name){
        gsInfo << "Save to " << name << "\n";
		std::ofstream file (name);
		for(index_t i = 0; i < vec.rows(); i++){
			file << vec[i];
			file << "\n";
		}
		file.close();
}

void saveMat(gsMatrix<> mat, std::string name){
        gsInfo << "Save to " << name << "\n";
		std::ofstream file (name);
		for(index_t i = 0; i < mat.rows(); i++){
            for (index_t j = 0; j < mat.cols(); j++){
			    file << mat(i,j);
                file << " ";
            }
			file << "\n";
		}
		file.close();
}

void convergenceTestOfDetJJacobian(gsOptParamMethod &pM){
	// gsVector<> result = dJC.generateDResultVector();
	// dJC.getDvectors(result);
	gsVector<> result = pM.m_dJC.evalCon();

    saveMat(pM.m_dJC.getJacobian().asDense(),BASE_FOLDER "/../results/J1.txt");

    gsIpOptSparseMatrix J = pM.mapMatrix(pM.m_dJC.space_mapper(),pM.m_dJC.getJacobian());
    saveMat(J.asDense(),BASE_FOLDER "/../results/J2.txt");
	gsMatrix<> Jac = J.asDense();
	gsInfo << "\n Size of D vector : " << result.size() << "\n";
	gsInfo << " Size of Jacobian : ( " << Jac.rows() << ", " << Jac.cols() << ")\n";

	gsVector<> des = pM.getFree();
	gsInfo << "\n Size of design vector : " << des.size() << "\n";

	index_t n = 10;
	gsVector<> Eps(n);
	gsVector<> Error0(n);
	gsVector<> Error1(n);

	gsVector<> ran;
	std::srand((unsigned int) std::time(0));
	ran.setRandom(des.size());


	for(index_t i = 0; i < n; i++){
		// Generate pertubation
		real_t eps = pow(2,-i);
		Eps[i] = eps;

		gsVector<> perturp;
		perturp.setZero(des.size());
		for(index_t i = 0; i < des.size(); i++){
			perturp[i] = ran[i];
		}

		perturp /= perturp.norm();

		perturp *= eps;

		gsVector<> newDes = des + perturp;


		pM.updateFree(newDes);


		// gsVector<> newres = dJC.generateDResultVector();
		// dJC.getDvectors(newres);
		gsVector<> newres = pM.m_dJC.evalCon();

		gsVector<> guess = result + Jac*perturp;

		real_t error0 = (result - newres).norm();
		real_t error1 = (guess - newres).norm();

		Error0[i] = error0;
		Error1[i] = error1;
	}

	gsVector<> rate;
	rate.setZero(n);
	rate.segment(1,n-1) = log10(Error1.segment(1,n-1).array()/Error1.segment(0,n-1).array())/log10(2);
	gsMatrix<> disp(n,4);
	disp << Eps,Error0,Error1,rate;
	gsInfo << "eps \tErr0 \tErr1 \trate\n";
	gsInfo << disp << "\n";

}
//
void convergenceTestOfParaJacobian(gsOptParamMethod &lOP){
	real_t liao = lOP.evalObj();
	gsVector<> grad = lOP.gradObj();

    gsMatrix<> tmp;
	gsMatrix<> hess = lOP.hessObj(tmp);

	gsInfo << "\n" << std::setprecision(10) << liao << "\n";
	// gsInfo << "\n" << grad << "\n";
    gsDebugVar(hess.rows());
    gsDebugVar(hess.cols());

	// std::srand((unsigned int) std::time(0));
	gsVector<> ran;
	ran.setRandom(lOP.n_free);

	gsVector<> des = lOP.getFree();
	gsInfo << "\n Size of design vector : " << des.size() << "\n";

	index_t beg = 0;
	index_t n = 20;
	gsVector<> Eps(n);
	gsVector<> Error0(n);
	gsVector<> Error1(n);
	gsVector<> Error2(n);

	for(index_t i = 0; i < n; i++){
		// Generate pertubation
		real_t eps = pow(2,-beg-i);
		Eps[i] = eps;

		gsVector<> perturp;
		perturp.setZero(lOP.n_free);

		for(index_t i = 0; i < lOP.n_free; i++){
			perturp[i] = ran[i];
		}

		perturp /= perturp.norm();

		perturp *= eps;

		gsVector<> newDes = des + perturp;


		lOP.updateFree(newDes);

		real_t newLiao = lOP.evalObj();
		real_t guess0 = liao;
		real_t guess1 = liao + grad.transpose()*perturp;
		real_t guess2 = liao + grad.transpose()*perturp + 0.5*perturp.transpose()*hess*perturp;

		gsInfo << newLiao <<" " << guess1 << " " << guess2 << "\n";
		// gsInfo << guess0 <<" " << guess1 << " " << newLiao << "\n";

		real_t error0 = std::abs(guess0 - newLiao);
		real_t error1 = std::abs(guess1 - newLiao);
		real_t error2 = std::abs(guess2 - newLiao);

		Error0[i] = error0;
		Error1[i] = error1;
		Error2[i] = error2;
	}

	gsVector<> rate;
	rate.setZero(n);
	rate.segment(1,n-1) = log10(Error1.segment(1,n-1).array()/Error1.segment(0,n-1).array())/log10(2);
	gsVector<> rate2;
	rate2.setZero(n);
	rate2.segment(1,n-1) = log10(Error2.segment(1,n-1).array()/Error2.segment(0,n-1).array())/log10(2);
	gsMatrix<> disp(n,6);
	disp << Eps,Error0,Error1,rate,Error2,rate2;
	gsInfo << "eps \tErr0 \tErr1 \trate \tErr2 \trate\n";
	gsInfo << disp << "\n";

}
//
// void convergenceTestOfParaJacobian(harmonicOptProblem &lOP){
// 	real_t liao = lOP.evalObj();
// 	gsVector<> grad = lOP.gradientObj();
// 	// gsMatrix<> hess = lOP.hessianObj();
//
// 	gsInfo << "\n" << std::setprecision(10) << liao << "\n";
// 	// gsInfo << "\n" << grad << "\n";
//
// 	// for (index_t i = 0; i < lOP.numDesignVars(); i++){
// 	// 	for (index_t j = 0; j < lOP.numDesignVars(); j++){
// 	// 		if (hess(i,j) == 0){ gsInfo << " ";
// 	// 	} else {
// 	// 		gsInfo << 1;
// 	// 	}
// 	// 		;
// 	// 	}
// 	// 	gsInfo <<"\n";
// 	// }
//
// 	// std::srand((unsigned int) std::time(0));
// 	gsVector<> ran;
// 	ran.setRandom(lOP.numDesignVars());
//
// 	gsVector<> des = lOP.dJC.getDesignVariables();
// 	gsInfo << "\n Size of design vector : " << des.size() << "\n";
//
// 	index_t beg = 0;
// 	index_t n = 20;
// 	gsVector<> Eps(n);
// 	gsVector<> Error0(n);
// 	gsVector<> Error1(n);
// 	gsVector<> Error2(n);
//
// 	for(index_t i = 0; i < n; i++){
// 		// Generate pertubation
// 		real_t eps = pow(2,-beg-i);
// 		Eps[i] = eps;
//
// 		gsVector<> perturp;
// 		perturp.setZero(lOP.numDesignVars());
//
// 		gsVector<> lb = lOP.desLowerBounds();
// 		gsVector<> ub = lOP.desUpperBounds();
//
// 		for(index_t i = 0; i < lOP.numDesignVars(); i++){
// 			if( ub[i] != lb[i]){
// 				perturp[i] = ran[i];
// 			}
// 		}
//
// 		perturp /= perturp.norm();
//
// 		perturp *= eps;
//
// 		gsVector<> newDes = des + perturp;
//
//
// 		lOP.dJC.updateDesignVariables(newDes);
//
// 		real_t newLiao = lOP.evalObj();
// 		real_t guess0 = liao;
// 		real_t guess1 = liao + grad.transpose()*perturp;
// 		real_t guess2 = liao + grad.transpose()*perturp ;//+ 0.5*perturp.transpose()*hess*perturp;
//
// 		// gsInfo << newLiao <<" " << guess1 << " " << guess2 << "\n";
// 		// gsInfo << guess0 <<" " << guess1 << " " << newLiao << "\n";
//
// 		real_t error0 = std::abs(guess0 - newLiao);
// 		real_t error1 = std::abs(guess1 - newLiao);
// 		// real_t error2 = std::abs(guess2 - newLiao);
//
// 		Error0[i] = error0;
// 		Error1[i] = error1;
// 		// Error2[i] = error2;
// 	}
//
// 	gsVector<> rate;
// 	rate.setZero(n);
// 	rate.segment(1,n-1) = log10(Error1.segment(1,n-1).array()/Error1.segment(0,n-1).array())/log10(2);
// 	gsVector<> rate2;
// 	rate2.setZero(n);
// 	rate2.segment(1,n-1) = log10(Error2.segment(1,n-1).array()/Error2.segment(0,n-1).array())/log10(2);
// 	gsMatrix<> disp(n,6);
// 	disp << Eps,Error0,Error1,rate,Error2,rate2;
// 	gsInfo << "eps \tErr0 \tErr1 \trate \tErr2 \trate\n";
// 	gsInfo << disp << "\n";
//
// }
//
void convergenceTestOfJacobian(gsShapeOptProblem &sOP){
	// std::srand((unsigned int) std::time(0));
	gsVector<> ran;
	ran.setRandom(sOP.numDesignVars());

	gsVector<> des = sOP.getDesignVariables();
	gsInfo << "\n Size of design vector : " << des.size() << "\n";

	sOP.updateDesignVariables(des);
	real_t obj = sOP.evalObj();
	gsVector<> grad = sOP.gradObj();

	index_t beg = 0;
	index_t n = 20;
	gsVector<> Eps(n);
	gsVector<> Error0(n);
	gsVector<> Error1(n);

	for(index_t i = 0; i < n; i++){
		// Generate pertubation
		real_t eps = pow(2,-beg-i);
		Eps[i] = eps;

		gsVector<> perturp;
		perturp.setZero(sOP.numDesignVars());

		index_t k = sOP.numDesignVars()/8;
		for(index_t i = 0; i < sOP.numDesignVars(); i++){
			perturp[i] = ran[i];
		}

		perturp /= perturp.norm();

		perturp *= eps;

		gsVector<> newDes = des + perturp;

		// sOP.resetParametrizationToReference();
		sOP.updateDesignVariables(newDes);
		// if (i == 1) sOP.writeToFile(sOP.getFlat(),"cpsPerturb.txt");

		real_t newObj = sOP.evalObj();
		real_t guess0 = obj;
		real_t guess1 = obj + grad.transpose()*perturp;

		gsInfo << guess0 <<" " << guess1 << " " << newObj << "\n";

		real_t error0 = std::abs(guess0 - newObj);
		real_t error1 = std::abs(guess1 - newObj);

		Error0[i] = error0;
		Error1[i] = error1;
	}

	gsVector<> rate;
	rate.setZero(n);
	rate.segment(1,n-1) = log10(Error1.segment(1,n-1).array()/Error1.segment(0,n-1).array())/log10(2);
	gsMatrix<> disp(n,4);
	disp << Eps,Error0,Error1,rate;
	gsInfo << "\neps \tErr0 \tErr1 \trate";
	gsInfo << disp << "\n";

}
//
// void convergenceTestOfParaJacobianToFile(paraOptProblem &lOP, std::string name){
//
// 	std::ofstream file;
// 	file.open(name);
//
// 	real_t liao = lOP.evalObj();
// 	gsVector<> grad = lOP.gradientObj();
// 	gsMatrix<> hess = lOP.hessianObj();
//
// 	// for (index_t i = 0; i < lOP.numDesignVars(); i++){
// 	// 	for (index_t j = 0; j < lOP.numDesignVars(); j++){
// 	// 		if (hess(i,j) == 0){ gsInfo << " ";
// 	// 	} else {
// 	// 		gsInfo << 1;
// 	// 	}
// 	// 		;
// 	// 	}
// 	// 	gsInfo <<"\n";
// 	// }
//
// 	std::srand((unsigned int) std::time(0));
// 	gsVector<> ran;
// 	ran.setRandom(lOP.numDesignVars());
//
// 	gsVector<> des = lOP.dJC.getDesignVariables();
// 	gsInfo << "\n Size of design vector : " << des.size() << "\n";
//
// 	index_t beg = 0;
// 	index_t n = 20;
// 	gsVector<> Eps(n);
// 	gsVector<> Error0(n);
// 	gsVector<> Error1(n);
// 	gsVector<> Error2(n);
//
// 	for(index_t i = 0; i < n; i++){
// 		// Generate pertubation
// 		real_t eps = pow(2,-beg-i);
// 		Eps[i] = eps;
//
// 		gsVector<> perturp;
// 		perturp.setZero(lOP.numDesignVars());
//
// 		gsVector<> lb = lOP.desLowerBounds();
// 		gsVector<> ub = lOP.desUpperBounds();
//
// 		for(index_t i = 0; i < lOP.numDesignVars(); i++){
// 			if( ub[i] != lb[i]){
// 				perturp[i] = ran[i];
// 			}
// 		}
//
// 		perturp /= perturp.norm();
//
// 		perturp *= eps;
//
// 		gsVector<> newDes = des + perturp;
//
//
// 		lOP.dJC.updateDesignVariables(newDes);
//
// 		real_t newLiao = lOP.evalObj();
// 		real_t guess0 = liao;
// 		real_t guess1 = liao + grad.transpose()*perturp;
// 		real_t guess2 = liao + grad.transpose()*perturp + 0.5*perturp.transpose()*hess*perturp;
//
// 		// gsInfo << newLiao <<" " << guess1 << " " << guess2 << "\n";
// 		file << guess0 <<" " << guess1 << " " << newLiao << "\n";
//
// 		real_t error0 = std::abs(guess0 - newLiao);
// 		real_t error1 = std::abs(guess1 - newLiao);
// 		real_t error2 = std::abs(guess2 - newLiao);
//
// 		Error0[i] = error0;
// 		Error1[i] = error1;
// 		Error2[i] = error2;
// 	}
//
// 	gsVector<> rate;
// 	rate.setZero(n);
// 	rate.segment(1,n-1) = log10(Error1.segment(1,n-1).array()/Error1.segment(0,n-1).array())/log10(2);
// 	gsVector<> rate2;
// 	rate2.setZero(n);
// 	rate2.segment(1,n-1) = log10(Error2.segment(1,n-1).array()/Error2.segment(0,n-1).array())/log10(2);
// 	gsMatrix<> disp(n,6);
// 	disp << Eps,Error0,Error1,rate,Error2,rate2;
// 	file << "eps \tErr0 \tErr1 \trate \tErr2 \trate\n";
// 	file << disp << "\n";
// 	file.close();
//
// }
//
// // void convergenceTestOfStateEquation(stateEquation &SE, std::string name, bool residualONOFF = false){
// // 	real_t err;
// // 	index_t ndofs, nelems;
// //
// // 	index_t n = 4;
// // 	gsVector<real_t> Errs, Rate;
// // 	gsVector<index_t> NumDofs;
// // 	Errs.setZero(n);
// // 	Rate.setZero(n);
// // 	NumDofs.setZero(n);
// //
// // 	index_t ref = 0;
// // 	// std::ofstream file(name);
// // 	for(index_t i = 0; i < n; i++){
// // 		gsInfo << i << "..\n";
// // 		SE.isRefined = false;
// // 		if (residualONOFF){
// // 			err = SE.getResidual();
// // 		} else {
// // 			SE.solve(err, ndofs, nelems, ref);
// // 		}
// // 		ref = 1;
// // 		Errs[i] = err;
// // 		NumDofs[i] = ndofs;
// // 		if (i > 0){ Rate[i] = log10(err/Errs[i-1])/log10(2.0);}
// // 		// file << sqrt(nelems) << " " << NumDofs[i] << " " << Errs[i] << "\n";
// // 	}
// // 	gsInfo << "Errs: " << Errs.transpose() << "\n";
// // 	gsInfo << "Rate: " << Rate.transpose() << "\n";
// // 	// file.close();
// // }
//
// // void convergenceTestOfSERhsJacobian(paraOptProblem &lOP, stateEquation &SE){
// // 	gsMatrix<> rhs = SE.getRhsZeroBC();
// // 	gsMatrix<> drhs = SE.getDerivativeOfRhsZeroBC();
// //
// // 	// for (index_t i = 0; i < lOP.numDesignVars(); i++){
// // 	// 	for (index_t j = 0; j < lOP.numDesignVars(); j++){
// // 	// 		if (hess(i,j) == 0){ gsInfo << " ";
// // 	// 	} else {
// // 	// 		gsInfo << 1;
// // 	// 	}
// // 	// 		;
// // 	// 	}
// // 	// 	gsInfo <<"\n";
// // 	// }
// //
// // 	std::srand((unsigned int) std::time(0));
// // 	gsVector<> ran;
// // 	ran.setRandom(lOP.numDesignVars());
// //
// // 	gsVector<> des = lOP.dJC.getDesignVariables();
// // 	gsInfo << "\n Size of design vector : " << des.size() << "\n";
// //
// // 	index_t n = 10;
// // 	gsVector<> Eps(n);
// // 	gsVector<> Error0(n);
// // 	gsVector<> Error1(n);
// // 	gsVector<> Error2(n);
// //
// // 	Error2.setZero(n);
// //
// // 	for(index_t i = 0; i < n; i++){
// // 		// Generate pertubation
// // 		real_t eps = pow(2,-i);
// // 		Eps[i] = eps;
// //
// // 		gsVector<> perturp;
// // 		perturp.setZero(lOP.numDesignVars());
// // 		for(index_t i = 0; i < lOP.numDesignVars(); i++){
// // 			perturp[i] = ran[i];
// // 		}
// //
// // 		perturp /= perturp.norm();
// //
// // 		perturp *= eps;
// //
// // 		gsVector<> newDes = des + perturp;
// //
// //
// // 		lOP.updateDesignVariables(newDes);
// //
// // 		gsMatrix<> newRhs = SE.getRhsZeroBC();
// // 		gsMatrix<> guess0 = rhs;
// // 		gsMatrix<> guess1 = rhs + drhs.transpose()*perturp;
// // 		// real_t guess2 = liao + grad.transpose()*perturp + 0.5*perturp.transpose()*hess*perturp;
// //
// // 		if (i == 5 or false){
// // 			gsMatrix<> disp1;
// // 			disp1.setZero(rhs.rows(),3);
// // 			disp1 << guess0, guess1, newRhs;
// // 			gsInfo << disp1;
// // 		}
// // 		// gsInfo << newRhs <<" " << guess1 << " " << guess2 << "\n";
// //
// // 		real_t error0 = (guess0 - newRhs).norm();
// // 		real_t error1 = (guess1 - newRhs).norm();
// // 		// real_t error2 = std::abs(guess2 - newLiao);
// //
// // 		Error0[i] = error0;
// // 		Error1[i] = error1;
// // 		// Error2[i] = error2;
// // 	}
// //
// // 	gsVector<> rate;
// // 	rate.setZero(n);
// // 	rate.segment(1,n-1) = log10(Error1.segment(1,n-1).array()/Error1.segment(0,n-1).array())/log10(2);
// // 	gsVector<> rate2;
// // 	rate2.setZero(n);
// // 	rate2.segment(1,n-1) = log10(Error2.segment(1,n-1).array()/Error2.segment(0,n-1).array())/log10(2);
// // 	gsMatrix<> disp(n,6);
// // 	disp << Eps,Error0,Error1,rate,Error2,rate2;
// // 	gsInfo << "eps \tErr0 \tErr1 \trate \tErr2 \trate\n";
// // 	gsInfo << disp << "\n";
// //
// // }
// //
// // void convergenceTestOfSEKuJacobian(paraOptProblem &lOP, stateEquation &SE){
// // 	gsMultiPatch<> u = SE.solve();
// // 	gsMatrix<> Ku = SE.getKu(u);
// // 	gsMatrix<> dKu = SE.getDerivativeOfKu(u);
// // 	gsInfo << "dKu is obtained\n";
// //
// // 	// std::srand((unsigned int) std::time(0));
// // 	gsVector<> ran;
// // 	ran.setRandom(lOP.numDesignVars());
// //
// // 	gsVector<> des = lOP.dJC.getDesignVariables();
// // 	gsInfo << "\n Size of design vector : " << des.size() << "\n";
// //
// // 	index_t n = 10;
// // 	gsVector<> Eps(n);
// // 	gsVector<> Error0(n);
// // 	gsVector<> Error1(n);
// // 	gsVector<> Error2(n);
// //
// // 	Error2.setZero(n);
// //
// // 	for(index_t i = 0; i < n; i++){
// // 		// Generate pertubation
// // 		real_t eps = pow(2,-i);
// // 		Eps[i] = eps;
// //
// // 		gsVector<> perturp;
// // 		perturp.setZero(lOP.numDesignVars());
// //
// // 		gsVector<> lb = lOP.desLowerBounds();
// // 		gsVector<> ub = lOP.desUpperBounds();
// //
// // 		for(index_t i = 0; i < lOP.numDesignVars(); i++){
// // 				perturp[i] = ran[i];
// // 		}
// //
// // 		perturp /= perturp.norm();
// //
// // 		perturp *= eps;
// //
// // 		gsVector<> newDes = des + perturp;
// //
// //
// // 		lOP.updateDesignVariables(newDes);
// //
// // 		gsMatrix<> newKu = SE.getKu(u);
// // 		gsMatrix<> guess0 = Ku;
// // 		gsMatrix<> guess1 = Ku + dKu.transpose()*perturp;
// // 		// real_t guess2 = liao + grad.transpose()*perturp + 0.5*perturp.transpose()*hess*perturp;
// //
// // 		if (i == 2 or false){
// // 			gsMatrix<> disp1;
// // 			disp1.setZero(Ku.rows(),3);
// // 			disp1 << guess0, guess1, newKu;
// // 			gsInfo << disp1;
// // 		}
// //
// // 		real_t error0 = (guess0 - newKu).norm();
// // 		real_t error1 = (guess1 - newKu).norm();
// // 		// real_t error2 = std::abs(guess2 - newLiao);
// //
// // 		Error0[i] = error0;
// // 		Error1[i] = error1;
// // 		// Error2[i] = error2;
// // 	}
// //
// // 	gsVector<> rate;
// // 	rate.setZero(n);
// // 	rate.segment(1,n-1) = log10(Error1.segment(1,n-1).array()/Error1.segment(0,n-1).array())/log10(2);
// // 	gsVector<> rate2;
// // 	rate2.setZero(n);
// // 	rate2.segment(1,n-1) = log10(Error2.segment(1,n-1).array()/Error2.segment(0,n-1).array())/log10(2);
// // 	gsMatrix<> disp(n,6);
// // 	disp << Eps,Error0,Error1,rate,Error2,rate2;
// // 	gsInfo << "eps \tErr0 \tErr1 \trate \tErr2 \trate\n";
// // 	gsInfo << disp << "\n";
// //
// // }
// //
// // void convergenceTestOfUJacobian(paraOptProblem &lOP, stateEquation &SE){
// // 	gsMatrix<> u = SE.getU();
// // 	gsMatrix<> du = SE.getDerivativeOfU();
// // 	gsInfo << "du is obtained\n";
// //
// // 	// std::srand((unsigned int) std::time(0));
// // 	gsVector<> ran;
// // 	ran.setRandom(lOP.numDesignVars());
// //
// // 	gsVector<> des = lOP.dJC.getDesignVariables();
// // 	gsInfo << "\n Size of design vector : " << des.size() << "\n";
// //
// // 	index_t n = 10;
// // 	gsVector<> Eps(n);
// // 	gsVector<> Error0(n);
// // 	gsVector<> Error1(n);
// // 	gsVector<> Error2(n);
// //
// // 	Error2.setZero(n);
// //
// // 	for(index_t i = 0; i < n; i++){
// // 		// Generate pertubation
// // 		real_t eps = pow(2,-i);
// // 		Eps[i] = eps;
// //
// // 		gsVector<> perturp;
// // 		perturp.setZero(lOP.numDesignVars());
// //
// // 		gsVector<> lb = lOP.desLowerBounds();
// // 		gsVector<> ub = lOP.desUpperBounds();
// //
// // 		for(index_t i = 0; i < lOP.numDesignVars(); i++){
// // 			if( ub[i] != lb[i]){
// // 				perturp[i] = ran[i];
// // 			}
// // 		}
// //
// // 		perturp /= perturp.norm();
// //
// // 		perturp *= eps;
// //
// // 		gsVector<> newDes = des + perturp;
// //
// //
// // 		lOP.updateDesignVariables(newDes);
// //
// // 		gsMatrix<> newU = SE.getU();
// // 		gsMatrix<> guess0 = u;
// // 		gsMatrix<> guess1 = u + du.transpose()*perturp;
// // 		// real_t guess2 = liao + grad.transpose()*perturp + 0.5*perturp.transpose()*hess*perturp;
// //
// // 		if (i == 2 or false){
// // 			gsMatrix<> disp1;
// // 			disp1.setZero(u.rows(),3);
// // 			disp1 << guess0, guess1, newU;
// // 			gsInfo << disp1;
// // 		}
// //
// // 		real_t error0 = (guess0 - newU).norm();
// // 		real_t error1 = (guess1 - newU).norm();
// // 		// real_t error2 = std::abs(guess2 - newLiao);
// //
// // 		Error0[i] = error0;
// // 		Error1[i] = error1;
// // 		// Error2[i] = error2;
// // 	}
// //
// // 	gsVector<> rate;
// // 	rate.setZero(n);
// // 	rate.segment(1,n-1) = log10(Error1.segment(1,n-1).array()/Error1.segment(0,n-1).array())/log10(2);
// // 	gsVector<> rate2;
// // 	rate2.setZero(n);
// // 	rate2.segment(1,n-1) = log10(Error2.segment(1,n-1).array()/Error2.segment(0,n-1).array())/log10(2);
// // 	gsMatrix<> disp(n,6);
// // 	disp << Eps,Error0,Error1,rate,Error2,rate2;
// // 	gsInfo << "eps \tErr0 \tErr1 \trate \tErr2 \trate\n";
// // 	gsInfo << disp << "\n";
// //
// // }
//
// void convergenceTestOfSEAuJacobian(index_t realOrImag,paraOptProblem &lOP, stateEquationAntenna &SE){
//     gsMultiPatch<> u,u_real,u_imag;
//     SE.solve(u_real,u_imag);
//     if (realOrImag == 0){
//         u = u_real;
//     } else {
//         u = u_imag;
//     }
//     gsMatrix<> U = SE.getU(realOrImag);
//
//     gsMatrix<> Ku = SE.getAu(realOrImag,U);
//
//     gsMatrix<> dKu = SE.getDerivativeOfAu(realOrImag,u);
//     gsInfo << "dKu is obtained\n";
//
//     // std::srand((unsigned int) std::time(0));
//     gsVector<> ran;
//     ran.setRandom(lOP.numDesignVars());
//
//     gsVector<> des = lOP.dJC.getDesignVariables();
//     gsInfo << "\n Size of design vector : " << des.size() << "\n";
//
//     index_t n = 15;
//     gsVector<> Eps(n);
//     gsVector<> Error0(n);
//     gsVector<> Error1(n);
//     gsVector<> Error2(n);
//
//     Error2.setZero(n);
//
//     for(index_t i = 0; i < n; i++){
//         // Generate pertubation
//         real_t eps = pow(2,-i);
//         Eps[i] = eps;
//
//         gsVector<> perturp;
//         perturp.setZero(lOP.numDesignVars());
//
//         gsVector<> lb = lOP.desLowerBounds();
//         gsVector<> ub = lOP.desUpperBounds();
//
//         for(index_t i = 0; i < lOP.numDesignVars(); i++){
//             perturp[i] = ran[i];
//         }
//
//         perturp /= perturp.norm();
//
//         perturp *= eps;
//
//         gsVector<> newDes = des + perturp;
//
//         lOP.updateDesignVariables(newDes);
//
//         gsMatrix<> newKu = SE.getAu(realOrImag,U);
//         gsMatrix<> guess0 = Ku;
//         gsMatrix<> guess1 = Ku + dKu.transpose()*perturp;
//         // real_t guess2 = liao + grad.transpose()*perturp + 0.5*perturp.transpose()*hess*perturp;
//
//         if (i == 10 or false){
//             gsMatrix<> disp1;
//             disp1.setZero(Ku.rows(),3);
//             disp1 << guess0, guess1, newKu;
//             // gsInfo << disp1;
//         }
//
//         real_t error0 = (guess0 - newKu).norm();
//         real_t error1 = (guess1 - newKu).norm();
//         // real_t error2 = std::abs(guess2 - newLiao);
//
//         Error0[i] = error0;
//         Error1[i] = error1;
//         // Error2[i] = error2;
//     }
//
//     gsVector<> rate;
//     rate.setZero(n);
//     rate.segment(1,n-1) = log10(Error1.segment(1,n-1).array()/Error1.segment(0,n-1).array())/log10(2);
//     gsVector<> rate2;
//     rate2.setZero(n);
//     rate2.segment(1,n-1) = log10(Error2.segment(1,n-1).array()/Error2.segment(0,n-1).array())/log10(2);
//     gsMatrix<> disp(n,6);
//     disp << Eps,Error0,Error1,rate,Error2,rate2;
//     gsInfo << "eps \tErr0 \tErr1 \trate \tErr2 \trate\n";
//     gsInfo << disp << "\n";
//
// }
//
// void convergenceTestOfSERhsJacobian(index_t realOrImag, paraOptProblem &lOP, stateEquationAntenna &SE){
//     gsMatrix<> rhs = SE.getRhsZeroBC(realOrImag);
//     gsMatrix<> drhs = SE.getDerivativeOfRhsZeroBC(realOrImag);
//
//     // for (index_t i = 0; i < lOP.numDesignVars(); i++){
//     // 	for (index_t j = 0; j < lOP.numDesignVars(); j++){
//     // 		if (hess(i,j) == 0){ gsInfo << " ";
//     // 	} else {
//     // 		gsInfo << 1;
//     // 	}
//     // 		;
//     // 	}
//     // 	gsInfo <<"\n";
//     // }
//
//     std::srand((unsigned int) std::time(0));
//     gsVector<> ran;
//     ran.setRandom(lOP.numDesignVars());
//
//     gsVector<> des = lOP.dJC.getDesignVariables();
//     gsInfo << "\n Size of design vector : " << des.size() << "\n";
//
//     index_t n = 10;
//     gsVector<> Eps(n);
//     gsVector<> Error0(n);
//     gsVector<> Error1(n);
//     gsVector<> Error2(n);
//
//     Error2.setZero(n);
//
//     for(index_t i = 0; i < n; i++){
//         // Generate pertubation
//         real_t eps = pow(2,-i);
//         Eps[i] = eps;
//
//         gsVector<> perturp;
//         perturp.setZero(lOP.numDesignVars());
//         for(index_t i = 0; i < lOP.numDesignVars()/2; i++){
//             perturp[i] = ran[i];
//         }
//
//         perturp /= perturp.norm();
//
//         perturp *= eps;
//
//         gsVector<> newDes = des + perturp;
//
//
//         lOP.updateDesignVariables(newDes);
//
//         gsMatrix<> newRhs = SE.getRhsZeroBC(realOrImag);
//         gsMatrix<> guess0 = rhs;
//         gsMatrix<> guess1 = rhs + drhs.transpose()*perturp;
//         // real_t guess2 = liao + grad.transpose()*perturp + 0.5*perturp.transpose()*hess*perturp;
//
//         if (i == 5 and false){
//             gsMatrix<> disp1;
//             disp1.setZero(rhs.rows(),3);
//             disp1 << guess0, guess1, newRhs;
//             gsInfo << disp1;
//         }
//         // gsInfo << newRhs <<" " << guess1 << " " << guess2 << "\n";
//
//         real_t error0 = (guess0 - newRhs).norm();
//         real_t error1 = (guess1 - newRhs).norm();
//         // real_t error2 = std::abs(guess2 - newLiao);
//
//         Error0[i] = error0;
//         Error1[i] = error1;
//         // Error2[i] = error2;
//     }
//
//     gsVector<> rate;
//     rate.setZero(n);
//     rate.segment(1,n-1) = log10(Error1.segment(1,n-1).array()/Error1.segment(0,n-1).array())/log10(2);
//     gsVector<> rate2;
//     rate2.setZero(n);
//     rate2.segment(1,n-1) = log10(Error2.segment(1,n-1).array()/Error2.segment(0,n-1).array())/log10(2);
//     gsMatrix<> disp(n,6);
//     disp << Eps,Error0,Error1,rate,Error2,rate2;
//     gsInfo << "eps \tErr0 \tErr1 \trate \tErr2 \trate\n";
//     gsInfo << disp << "\n";
//
// }
//
// void convergenceTestOfUJacobian(paraOptProblem &lOP, stateEquationAntenna &SE){
//     gsMatrix<> u = SE.getU();
//     gsMatrix<> du = SE.getDerivativeOfU();
//     gsInfo << "du is obtained\n";
//
//     // std::srand((unsigned int) std::time(0));
//     gsVector<> ran;
//     ran.setRandom(lOP.numDesignVars());
//
//     gsVector<> des = lOP.dJC.getDesignVariables();
//     gsInfo << "\n Size of design vector : " << des.size() << "\n";
//
//     index_t n = 20;
//     gsVector<> Eps(n);
//     gsVector<> Error0(n);
//     gsVector<> Error1(n);
//     gsVector<> Error2(n);
//
//     Error2.setZero(n);
//
//     for(index_t i = 0; i < n; i++){
//         // Generate pertubation
//         real_t eps = pow(2,-i);
//         Eps[i] = eps;
//
//         gsVector<> perturp;
//         perturp.setZero(lOP.numDesignVars());
//
//         gsVector<> lb = lOP.desLowerBounds();
//         gsVector<> ub = lOP.desUpperBounds();
//
//         for(index_t i = 0; i < lOP.numDesignVars(); i++){
//             if( ub[i] != lb[i]){
//                 perturp[i] = ran[i];
//             }
//         }
//
//         perturp /= perturp.norm();
//
//         perturp *= eps;
//
//         gsVector<> newDes = des + perturp;
//
//
//         lOP.updateDesignVariables(newDes);
//
//         gsMatrix<> newU = SE.getU();
//         gsMatrix<> guess0 = u;
//         gsMatrix<> guess1 = u + du.transpose()*perturp;
//         // real_t guess2 = liao + grad.transpose()*perturp + 0.5*perturp.transpose()*hess*perturp;
//
//         if (i == 2 or false){
//             gsMatrix<> disp1;
//             disp1.setZero(u.rows(),3);
//             disp1 << guess0, guess1, newU;
//             gsInfo << disp1;
//         }
//
//         real_t error0 = (guess0 - newU).norm();
//         real_t error1 = (guess1 - newU).norm();
//         // real_t error2 = std::abs(guess2 - newLiao);
//
//         Error0[i] = error0;
//         Error1[i] = error1;
//         // Error2[i] = error2;
//     }
//
//     gsVector<> rate;
//     rate.setZero(n);
//     rate.segment(1,n-1) = log10(Error1.segment(1,n-1).array()/Error1.segment(0,n-1).array())/log10(2);
//     gsVector<> rate2;
//     rate2.setZero(n);
//     rate2.segment(1,n-1) = log10(Error2.segment(1,n-1).array()/Error2.segment(0,n-1).array())/log10(2);
//     gsMatrix<> disp(n,6);
//     disp << Eps,Error0,Error1,rate,Error2,rate2;
//     gsInfo << "eps \tErr0 \tErr1 \trate \tErr2 \trate\n";
//     gsInfo << disp << "\n";
//
// }
//
// void convergenceTestOfDesignUpdateJacobian(shapeOptProblem &sOP){
// 	std::srand((unsigned int) std::time(0));
// 	gsVector<> ran;
// 	ran.setRandom(sOP.numDesignVars());
//
// 	gsVector<> des = sOP.getDesignVariables();
// 	gsInfo << "\n Size of design vector : " << des.size() << "\n";
//
// 	sOP.updateDesignVariables(des);
// 	gsVector<> obj = sOP.getUpdateToCps(des); //OP.pOP.getDesignVariables();
// 	gsMatrix<> grad = sOP.derivativeOfDesignUpdate();
//
// 	index_t beg = -5;
// 	index_t n = 20;
// 	gsVector<> Eps(n);
// 	gsVector<> Error0(n);
// 	gsVector<> Error1(n);
//
// 	for(index_t i = 0; i < n; i++){
// 		// Generate pertubation
// 		gsInfo << "i = " << i << "\n";
// 		real_t eps = pow(2,-beg-i);
// 		Eps[i] = eps;
//
// 		gsVector<> perturp;
// 		perturp.setZero(sOP.numDesignVars());
//
// 		index_t k = sOP.numDesignVars()/8;
// 		for(index_t i = 0; i < sOP.numDesignVars(); i++){
// 			if(i != 0 && i != k && i != 2*k && i != 3*k ){
// 				perturp[i] = ran[i];
// 			}
// 		}
//
// 		perturp /= perturp.norm();
//
// 		perturp *= eps;
//
// 		gsVector<> newDes = des + perturp;
//
// 		// sOP.resetParametrizationToReference();
// 		// sOP.updateDesignVariables(newDes);
// 		// if (i == 1) sOP.writeToFile(sOP.dJC.getDesignVariables(),"cpsPerturb.txt");
//
// 		gsVector<> newObj = sOP.getUpdateToCps(newDes);//OP.pOP.getDesignVariables();
// 		gsVector<> guess0 = obj;
// 		// gsInfo << "size obj (" << obj.rows() << ", " << obj.cols() << ")\n";
// 		// gsInfo << "size grad (" << grad.rows() << ", " << grad.cols() << ")\n";
// 		gsVector<> guess1 = obj + grad*perturp;
//
// 		// gsInfo << guess0 <<" " << guess1 << " " << newObj << "\n";
// 		if (i == 7 or false){
// 			gsMatrix<> disp1;
// 			disp1.setZero(obj.rows(),3);
// 			disp1 << guess0, guess1, newObj;
// 			gsInfo << disp1;
// 		}
//
// 		real_t error0 = (guess0 - newObj).norm();
// 		real_t error1 = (guess1 - newObj).norm();
//
// 		Error0[i] = error0;
// 		Error1[i] = error1;
// 	}
//
// 	gsVector<> rate;
// 	rate.setZero(n);
// 	rate.segment(1,n-1) = log10(Error1.segment(1,n-1).array()/Error1.segment(0,n-1).array())/log10(2);
// 	gsMatrix<> disp(n,4);
// 	disp << Eps,Error0,Error1,rate;
// 	gsInfo << "\neps \tErr0 \tErr1 \trate";
// 	gsInfo << disp << "\n";
//
// }
//
// void checkGradientsWithFD(modLiaoOptProblem &lOP){
// 	gsVector<> grad = lOP.gradientObj();
// 	gsVector<> gradFD = lOP.gradObjFD(lOP.getDesignVariables());
//
// 	gsVector<> diff = grad - gradFD;
//
// 	gsInfo << "Check gradients\n";
// 	for(index_t i = 0; i < grad.rows(); i++){
// 		real_t ad = abs(diff[i]);
// 		if (ad > 1e-12){
// 			gsInfo << "Error of size " << ad << " at index " << i << "\n";
// 		}
// 	}
// }

gsVector<> loadVec(index_t n,std::string name){
		gsVector<> vec;
		vec.setZero(n);
		std::ifstream file (name);
		real_t val;
		for(index_t i = 0; i < n; i++){
			file >> val;
			vec[i] = val;
		}
		return vec;

}


inline bool exists (const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

gsVector<> getCoefficients(gsMultiPatch<> &mp){
  // Count the total number of controlpoints
  index_t n_controlpoints = 0;
  for(index_t i = 0; i < mp.nBoxes(); i++){
    n_controlpoints += mp.patch(i).coefsSize();
  }

  gsVector<> cf(n_controlpoints);
  index_t j = 0;
  for(index_t i = 0; i < mp.nBoxes(); i++){
    for(index_t k = 0; k < mp.patch(i).coefsSize(); k++){
      cf[j] = mp.patch(i).coef(k,0);
      j++;
    }
  }

  return cf;
}

// void testErrorOnDesigns(errorOptProblem &eOP){
//
// 	index_t i = 0;
// 	real_t err;
// 	gsMultiPatch<> sol;
// 	std::string folder = "errorWithGradient1/";
// 	std::string begin = folder + "design_";
// 	std::string end = ".txt";
//
// 	std::ofstream file;
// 	std::string nameErr = folder + "errorsSqMoreQuadPts.txt";
// 	file.open(nameErr);
// 	gsInfo << "Saving err to " << nameErr << "\n";
//
// 	while (true) {
// 		std::string name = begin + std::to_string(i) + end;
// 		gsInfo << "Loading from " << name << "\n";
//
// 		if (not exists(name)){ break; }
//
// 		gsVector<> des = loadVec(eOP.numDesignVars(),name);
//
// 		eOP.updateDesignVariables(des);
//
// 		sol = eOP.SE->solve();
// 		err = eOP.evalObj(sol);
// 		file << err << "\n";
//
// 		gsVector<> solCoeffs = getCoefficients(sol);
//
// 		std::string nameCoeffs = folder + "solCoeffs_" + std::to_string(i) + ".txt";
// 		gsInfo << "Saving to " << nameCoeffs << "\n";
// 		saveVec(solCoeffs,nameCoeffs);
//
// 		i++;
// 		// if(i == 100){break;}
//
// 	}
// 	file.close();
// }
//
// void testConvOfGradOnDesigns(errorOptProblem &eOP, index_t maxiter){
//
// 	index_t i = 0;
// 	real_t err;
// 	gsMultiPatch<> sol;
// 	std::string folder = "errorWithGradient2/";
// 	std::string begin = folder + "design_";
// 	std::string end = ".txt";
//
//
// 	while (true) {
// 		std::string name = begin + std::to_string(i) + end;
// 		gsInfo << "Loading from " << name << "\n";
//
// 		if (not exists(name)){ break; }
//
// 		gsVector<> des = loadVec(eOP.numDesignVars(),name);
//
// 		std::string out = folder + "convTest_" + std::to_string(i) + end;
//
// 		eOP.updateDesignVariables(des);
//
// 		convergenceTestOfParaJacobianToFile(eOP,out);
//
// 		i++;
// 		if(i == maxiter){break;}
//
// 	}
// }

gsMultiPatch<> getGeometry(index_t n, index_t m, index_t degree){

    std::string folder;
    if (n == 4 && m == 4 && degree == 2){
        folder = "/parametrizations/para_x4_y4_p2_q2/";
    } else if (n == 5 && m == 4 && degree == 2){
        folder = "/parametrizations/para_x5_y4_p2_q2/";
    } else if (n == 6 && m == 4 && degree == 2){
        folder = "/parametrizations/para_x6_y4_p2_q2/";
    } else if (n == 7 && m == 4 && degree == 2){
        folder = "/parametrizations/para_x7_y4_p2_q2/";
    } else if (n == 8 && m == 4 && degree == 2){
        folder = "/parametrizations/para_x8_y4_p2_q2/";
    } else if (n == 4 && m == 5 && degree == 2){
        folder = "/parametrizations/para_x4_y5_p2_q2/";
    } else if (n == 5 && m == 5 && degree == 2){
        folder = "/parametrizations/para_x5_y5_p2_q2/";
    } else if (n == 6 && m == 5 && degree == 2){
        folder = "/parametrizations/para_x6_y5_p2_q2/";
    } else if (n == 7 && m == 5 && degree == 2){
        folder = "/parametrizations/para_x7_y5_p2_q2/";
    } else if (n == 8 && m == 5 && degree == 2){
        folder = "/parametrizations/para_x8_y5_p2_q2/";
    }  else if (n == 6 && m == 6 && degree == 2){
        folder = "/parametrizations/para_x6_y6_p2_q2/";
    } else {
        GISMO_ERROR("Parametrization is not generated with these parameters..\n");
    }



	gsInfo << "----------------------\n\n"
	<< "n: " << n << "\n\n"
	<< "m: " << m << "\n\n"
	<< "degree: " << degree << "\n\n"
	<< "----------------------\n\n";

	// 1. construction of a knot vector for each direction
	gsKnotVector<> kv1(0, 1, n - degree - 1, degree + 1);
	gsKnotVector<> kv2(0, 1, m - degree - 1, degree + 1);
	// 2. construction of a basis
	gsTensorBSplineBasis<2, real_t> basis(kv1, kv2);
	// 3. construction of a coefficients
	gsMatrix<> greville = basis.anchors();
	gsMatrix<> coefs (greville.cols(), 2);

	readFromTxt(BASE_FOLDER + folder + "l.txt", coefs);

	// gsInfo << coefs;

	// 4. putting basis and coefficients toghether
	gsTensorBSpline<2, real_t>  left(basis, coefs);

	//! [Geometry data]
	// Define Geometry, must be a gsMultiPatch object

	gsMultiPatch<> patches = gsMultiPatch<>(left);

	readFromTxt(BASE_FOLDER + folder + "b.txt", coefs);

	// 4. putting basis and coefficients toghether
	gsTensorBSpline<2, real_t>  bottom(basis, coefs);
	patches.addPatch(bottom);

	readFromTxt(BASE_FOLDER + folder + "r.txt", coefs);

	// 4. putting basis and coefficients toghether
	gsTensorBSpline<2, real_t>  right(basis, coefs);
	patches.addPatch(right);

	gsTensorBSplineBasis<2, real_t> basisMid(kv1, kv1);
	// 3. construction of a coefficients
	gsMatrix<> grevilleMid = basisMid.anchors();
	gsMatrix<> coefsMid (grevilleMid.cols(), 2);

	readFromTxt(BASE_FOLDER + folder + "m.txt", coefsMid);
	gsTensorBSpline<2, real_t>  middle(basisMid, coefsMid);
	patches.addPatch(middle);

	readFromTxt(BASE_FOLDER + folder + "t.txt", coefs);

	// 4. putting basis and coefficients toghether
	gsTensorBSpline<2, real_t>  top(basis, coefs);
	patches.addPatch(top);

	double tol = 1e-2;
	patches.computeTopology(tol,true);
	patches.closeGaps(tol);

	std::string out = "Geometry";
	gsInfo << "Writing the gsMultiPatch to a paraview file: " << out
	<< "\n\n";
	gsWriteParaview(patches, out);

	return patches;
	// GISMO_ERROR("stop...");

}

gsMultiPatch<> get3DGeometry(){
    gsMultiPatch<> patches;
    gsFileData<> data("geometries3D/cube.xml");

    gsInfo  <<"* There is "<< data.count< gsGeometry<> >() <<" "
            <<data.type< gsGeometry<> >()<<" "<< data.tag< gsGeometry<> >()
            <<" in the file.\n";
    for(index_t i = 1; i <= data.count< gsGeometry<> >(); i++){
        gsGeometry<>::uPtr o = data.getId< gsGeometry<> >(i);
        o->degreeElevate(1);
        o->uniformRefine();
        patches.addPatch(*o);
    }
    patches.computeTopology();

    return patches;
}

int main(int argc, char* argv[]){
gsInfo <<  "Hello G+Smo.\n";

// Parse command line
std::string output(".");
int degree = 2;
int nx = 5;
int ny = 4;

int numRefine = 1;
int maxiter = 10;

int param = 0; // 0: spring, 1: modLiao, 2: winslow

bool plotDesign = false;
bool plotMagnitude = false;
bool plotSolution = false;
bool saveCps = false;
bool useDJC = true;

int quA = 1;
int quB = 1;

int startDes = -1;

gsCmdLine cmd("A test of lumped mass matricies");
cmd.addInt("p", "degree", "Degree of B-Splines.", degree);
cmd.addInt("a", "param", "Parametrization: 0: spring, 1: modLiao, 2: winslow", param);
cmd.addInt("r", "numberRefine", "Number of refinements", numRefine);
cmd.addInt("n", "nx", "Number of splines in first direction", nx);
cmd.addInt("m", "ny", "Number of splines in second direction", ny);
cmd.addInt("i", "maxIter", "Maximal number of reparametrizations", maxiter);

cmd.addInt("A", "quA", "quA", quA);
cmd.addInt("B", "quB", "quB", quB);

cmd.addString("o","output","Name of the output folder (relative to BASE_FOLDER)",output);

cmd.addSwitch("plot", "Create a ParaView visualization file of designs", plotDesign);
cmd.addSwitch("noDJC", "Dont use det jac constraints", useDJC);
cmd.addSwitch("plotmag", "Create a ParaView visualization file of magnitude (obj function)", plotMagnitude);
cmd.addSwitch("plotsol", "Create a ParaView visualization file of solution (real and imag)", plotSolution);
cmd.addSwitch("savecps", "save controlpoinst for design during optimization", saveCps);

cmd.addInt("s", "startDes", "design number to start from", startDes);

cmd.getValues(argc,argv);

char buffer [50];
std::sprintf(buffer,"p = %d \n n = %d\n",degree,numRefine);
gsInfo << buffer;

gsMultiPatch<> patches = getGeometry(nx,ny,degree);

// gsMultiPatch<> patches = get3DGeometry();

// gsMultiPatch<> patches = gsMultiPatch<>(*gsNurbsCreator<>::BSplineSquare(2));
// patches.basis(0).setDegree(degree);
//
// for(int i = 0; i < numRefine; i++){
// 	patches.uniformRefine();
// }
//
// int len = patches[0].coefsSize();
// //
// int pm = 1;
// setupAntennaDomain(patches,degree,sqrt(len),pm);
//
// std::string out = "Geometry";
// gsInfo << "Writing the gsMultiPatch to a paraview file: " << out
// << "\n\n";
// gsWriteParaview(patches, out);

// gsInfo << "Domain is a: \n" << std::flush;
// gsInfo << "The domain is a "<< patches <<"\n";

// gsInfo << "patch 0: " << patches.patch(0) << "\n";

// gsModLiao modLiao(&patches,useDJC);

gsShapeOptLog slog(output);
// gsOptAntenna optA1(&patches,numRefine,&slog,param,quA,quB);
//
// // Start from design at ../results/paramTest/cps_1_0.txt;
// gsVector<> flat = loadVec(optA1.n_flat,BASE_FOLDER "/../results/ParamTests/Winslow/cps_0_80.txt");
// optA1.m_paramMethod->updateFlat(flat);

gsOptAntenna optA(&patches,numRefine,&slog,param,quA,quB);

//
// gsInfo << "obj " << optA.evalObj() << "\n";
// gsInfo << "min d " << optA.m_dJC.evalCon().maxCoeff() << "\n";
//
// gsWriteParaview(patches,BASE_FOLDER "/../results/ParamTests/issue");

// Start from spring
// gsSpringMethod spring(&patches,optA.m_mappers);
// spring.update();

if (param == 0){
    slog << "Using Spring Method for parametrization \n";
    optA.solve();
} else if (param == 1) {
    slog << "Using modLiao for parametrization \n";
    optA.runOptimization(maxiter);
} else if (param == 2) {
    slog << "Using winslow for parametrization \n";
    optA.runOptimization(maxiter);
}


// optA.solve();
// optA.runOptimization(maxiter);
// gsInfo << "obj = " << optA.evalObj() << "\n";

// convergenceTestOfJacobian(optA);

// gsDetJacConstraint dJC(&patches);
// gsInfo << "det J: " << dJC.evalCon().maxCoeff();
//
// gsModLiao modLiao(&patches,useDJC);
// modLiao.print();
// modLiao.update();
// convergenceTestOfDetJJacobian(modLiao);
// convergenceTestOfParaJacobian(modLiao);

// modLiao.solve();
// gsWriteParaview(patches,BASE_FOLDER "/../results/modL");

// gsAffineOptParamMethod affModLiao(&modLiao);
// affModLiao.computeMap();
// affModLiao.update();
// gsVector<> ran;
// ran.setRandom(modLiao.n_tagged);
// ran *= 0.05;
// gsInfo << ran;
// affModLiao.update(modLiao.getTagged() + ran);
// gsWriteParaview(patches,BASE_FOLDER "/../results/affModL");

// gsMatrix<> hOT;
// modLiao.hessObj(hOT);
// gsDebugVar(hOT.rows());
// gsDebugVar(hOT.cols());
//
//
// sM.updateTagged(sM.getTagged() + pert)
// sM.update();
// sM.updateDesignVariables(des);
// gsInfo << "det J: " << dJC.evalCon().maxCoeff();

// gsInfo << "\n ==== DONE ==== \n";
return 0;
}

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
#include "gsNewDetJacConstraint.h"
#include "gsAggregatedConstraint.h"
#include "gsAggFun.h"

#include "gsParamMethod.h"
#include "gsSpringMethod.h"
#include "gsModLiao.h"
#include "gsWinslow.h"
// #include "gsKnupp.h"
#include "gsLiao.h"
#include "gsHarmonic.h"
#include "gsMaxDetJac.h"

#include "gsAffineOptParamMethod.h"
#include "gsIpOptSparseMatrix.h"
#include "gsShapeOptProblem.h"
#include "gsShapeOptWithReg.h"
#include "gsOptAntenna.h"
#include "gsOptParam.h"
#include "gsShapeOptLog.h"

#include "gsStateEquationPotWaves.h"

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
	readFromTxt(BASE_FOLDER "/parametrizations/" + folder + name, cc);

	mp.patch(i).setCoefs(cc);
}

void saveVec(gsVector<> &vec, std::string name){
        gsInfo << "Save to " << name << "\n";
		std::ofstream file (name);
		for(index_t i = 0; i < vec.rows(); i++){
			file << std::setprecision(12) << vec[i];
			file << "\n";
		}
		file.close();
}

void saveVec(gsVector<index_t> &vec, std::string name){
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
			file << std::setprecision(12) << "\n";
		}
		file.close();
}

void saveSparseMat(gsSparseMatrix<> mat, std::string name){
        gsInfo << "Save to " << name << "\n";
        std::ofstream file (name);
        for (int k=0; k<mat.outerSize(); ++k){
            for (gsSparseMatrix<real_t>::InnerIterator it(mat,k); it; ++it)
            {
                file << std::setprecision(12) << it.row();   // row index
                file << " ";
                file << it.col();   // col index (here it is equal to k)
                file << " ";
                file << it.value();
                file << "\n";
            }
        }
    }

void convergenceTestOfDetJJacobian(gsOptParamMethod &pM){
	// gsVector<> result = dJC.generateDResultVector();
	// dJC.getDvectors(result);
	gsVector<> result = pM.m_dJC->evalCon();

    saveMat(pM.m_dJC->getJacobian().asDense(),BASE_FOLDER "/../results/J1.txt");

    gsIpOptSparseMatrix J = pM.mapMatrix(pM.m_dJC->space_mapper(),pM.m_dJC->getJacobian());
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
		gsVector<> newres = pM.m_dJC->evalCon();

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

void convergenceTestOfConstraint(gsOptParamMethod &pM, gsConstraint &con){
	// gsVector<> result = dJC.generateDResultVector();
	// dJC.getDvectors(result);
	gsVector<> result = con.evalCon();

    gsIpOptSparseMatrix J = pM.mapMatrix(con.space_mapper(),con.getJacobian());
	gsMatrix<> Jac = J.asDense();
	gsInfo << "\n Size of D vector : " << result.size() << "\n";
	gsInfo << " Size of Jacobian : ( " << Jac.rows() << ", " << Jac.cols() << ")\n";

	gsVector<> des = pM.getFree();
	gsInfo << "\n Size of design vector : " << des.size() << "\n";

	index_t n = 20;
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
		gsVector<> newres = con.evalCon();

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

void convergenceTestOfParaLagrangianJacobian(gsOptParamMethod &lOP){
	real_t liao = lOP.evalLagrangian();
	gsVector<> grad = lOP.gradLagrangian();

    gsMatrix<> tmp;
	gsMatrix<> hess = lOP.hessLagrangian(tmp);

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

		real_t newLiao = lOP.evalLagrangian();
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

void convergenceTestOfJacobian(gsShapeOptProblem &sOP){
	std::srand((unsigned int) std::time(0));
	gsVector<> ran;
	ran.setRandom(sOP.numDesignVars());

	gsVector<> des = sOP.getDesignVariables();
	gsInfo << "\n Size of design vector : " << des.size() << "\n";

	sOP.updateDesignVariables(des);
	real_t obj = sOP.evalObj();
	gsVector<> grad = sOP.gradObj();

	index_t beg = 3;
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

void convergenceTestOfParamMethodJacobian(gsParamMethod &pM){
	// gsVector<> result = dJC.generateDResultVector();
	// dJC.getDvectors(result);
    gsVector<> flat = pM.getFlat();
    gsVector<> tag = pM.getTagged();

    pM.update();
	gsVector<> free_cps = pM.getFree();

	gsMatrix<> Jac = pM.jacobUpdate(tag);
    pM.updateFlat(flat);

    index_t beg = 5;
	index_t n = 8;
	gsVector<> Eps(n);
	gsVector<> Error0(n);
	gsVector<> Error1(n);

	gsVector<> ran;
	std::srand((unsigned int) std::time(0));
	ran.setRandom(tag.size());

	for(index_t i = 0; i < n; i++){
		// Generate pertubation
		real_t eps = pow(2,-beg-i);
		Eps[i] = eps;

		gsVector<> perturp;
		perturp.setZero(tag.size());
		for(index_t i = 0; i < tag.size(); i++){
			perturp[i] = ran[i];
		}

		// perturp /= perturp.norm();

		gsVector<> newTag = tag + eps*perturp;

		// gsVector<> newres = dJC.generateDResultVector();
		// dJC.getDvectors(newres);
        pM.update(newTag);
		gsVector<> newFree = pM.getFree();
        pM.updateFlat(flat);

        gsVector<> FD = (newFree - free_cps)/eps;

        gsVector<> deriv = Jac*perturp;

        gsInfo << FD.norm() << " \t" << deriv.norm() << "\n";
        if (false && i == 7){
            gsInfo << "FD \t jac \n";
            for (index_t i = 0; i < FD.size(); i++){
                gsInfo << FD[i] << " \t" << deriv[i] << "\n";
            }
        }

		real_t error1 = (FD - deriv).norm();

		Error1[i] = error1;
	}

	gsVector<> rate;
	rate.setZero(n);
	rate.segment(1,n-1) = log10(Error1.segment(1,n-1).array()/Error1.segment(0,n-1).array())/log10(2);
	gsMatrix<> disp(n,3);
	disp << Eps,Error1,rate;
	gsInfo << "eps \tErr1 \trate\n";
	gsInfo << disp << "\n";

}

void convergenceTestOfParaHessian(gsOptParamMethod &lOP){

    gsMatrix<> Hcx;

	gsMatrix<> hess = lOP.hessObj(Hcx);

	std::srand((unsigned int) std::time(0));
	gsVector<> ran_free;
	ran_free.setRandom(lOP.n_free);

	gsVector<> ran_tagged;
	ran_tagged.setRandom(lOP.n_tagged);

	gsVector<> des = lOP.getFree();
	gsVector<> tag = lOP.getTagged();
    gsVector<> flat = lOP.getFlat();

	index_t beg = 7;
	index_t n = 8;
	gsVector<> Eps(n);
	gsVector<> ErrorT(n);
	gsVector<> ErrorF(n);

	for(index_t i = 0; i < n; i++){
		// Generate pertubation
		real_t eps = pow(2,-beg-i);
		Eps[i] = eps;

		gsVector<> p;
		p.setZero(lOP.n_free);
		gsVector<> q;
		q.setZero(lOP.n_tagged);

		for(index_t i = 0; i < lOP.n_free; i++){
			p[i] = ran_free[i];
		}
		for(index_t i = 0; i < lOP.n_tagged; i++){
			q[i] = ran_tagged[i];
		}

		gsVector<> newDes = des + eps*p;
		gsVector<> newTag = tag + eps*q;

        // W(c + eps*p,x + eps*q)
        lOP.updateFlat(flat);
        lOP.updateFree(newDes);
        lOP.updateTagged(newTag);

        real_t Wpq = lOP.evalObj();

        // W(c + eps*p,x)
        lOP.updateFlat(flat);
        lOP.updateFree(newDes);
        lOP.updateTagged(tag);

        real_t Wp = lOP.evalObj();

        // W(c,x + eps*q)
        lOP.updateFlat(flat);
        lOP.updateFree(des);
        lOP.updateTagged(newTag);

        real_t Wq = lOP.evalObj();

        // W(c,x)
        lOP.updateFlat(flat);
        lOP.updateFree(des);
        lOP.updateTagged(tag);

        // W(c + eps*p,x - eps*q)
        lOP.updateFlat(flat);
        lOP.updateFree(newDes);
        lOP.updateTagged(tag - eps*q);

        real_t Wpmq = lOP.evalObj();

        // W(c - eps*p,x + eps*q)
        lOP.updateFlat(flat);
        lOP.updateFree(des - eps*p);
        lOP.updateTagged(newTag);

        real_t Wmpq = lOP.evalObj();

        // W(c - eps*p,x - eps*q)
        lOP.updateFlat(flat);
        lOP.updateFree(des - eps*p);
        lOP.updateTagged(tag - eps*q);

        real_t Wmpmq = lOP.evalObj();

        // W(c,x)
        lOP.updateFlat(flat);
        lOP.updateFree(des);
        lOP.updateTagged(tag);

        real_t W = lOP.evalObj();

        real_t FDtagged2 = (Wpq - Wq - Wp + W)*1/(eps*eps);
        real_t FDtagged = (Wpq - Wpmq - Wmpq + Wmpmq)/(4*eps*eps);

        // W(c - eps*p,x)
        lOP.updateFlat(flat);
        lOP.updateFree(des - eps*p);
        lOP.updateTagged(tag);

        real_t Wmp = lOP.evalObj();

        real_t FDfree = (Wp - 2*W + Wmp)*1/(eps*eps);

		real_t errorT = FDtagged - p.transpose()*Hcx*q;
        gsInfo << FDtagged << " and " << FDtagged2 << ", " << p.transpose()*Hcx*q << "\n";
		real_t errorF = FDfree - p.transpose()*hess*p;
        // gsInfo << FDfree << ", " << p.transpose()*hess*p << "\n";

		ErrorT[i] = errorT;
		ErrorF[i] = errorF;
	}

	gsVector<> rateT;
	rateT.setZero(n);
	rateT.segment(1,n-1) = log10(ErrorT.segment(1,n-1).array()/ErrorT.segment(0,n-1).array())/log10(2);
	gsVector<> rateF;
	rateF.setZero(n);
	rateF.segment(1,n-1) = log10(ErrorF.segment(1,n-1).array()/ErrorF.segment(0,n-1).array())/log10(2);
	gsMatrix<> disp(n,5);
	disp << Eps,ErrorT,rateT,ErrorF,rateF;
	gsInfo << "eps \t\tErrT \trateT \t\tErrF \trateF\n";
	gsInfo << disp << "\n";

}

void convergenceTestOfParaGrad(gsWinslow &lOP){

    gsVector<> gradTagged;
	gsVector<> grad = lOP.gradObj(gradTagged);
    gsVector<> flat = lOP.getFlat();

	gsVector<> ran_free;
	ran_free.setRandom(lOP.n_free);

	gsVector<> ran_tagged;
	ran_tagged.setRandom(lOP.n_tagged);

	gsVector<> des = lOP.getFree();
	gsVector<> tag = lOP.getTagged();

	index_t beg = 7;
	index_t n = 10;
	gsVector<> Eps(n);
	gsVector<> ErrorT(n);
	gsVector<> ErrorF(n);

	for(index_t i = 0; i < n; i++){
		// Generate pertubation
		real_t eps = pow(2,-beg-i);
		Eps[i] = eps;

		gsVector<> p;
		p.setZero(lOP.n_free);
		gsVector<> q;
		q.setZero(lOP.n_tagged);

		for(index_t i = 0; i < lOP.n_free; i++){
			p[i] = ran_free[i];
		}
		for(index_t i = 0; i < lOP.n_tagged; i++){
			q[i] = ran_tagged[i];
		}

		gsVector<> newDes = des + eps*p;
		gsVector<> newTag = tag + eps*q;


        // W(c + eps*p,x)
        lOP.updateFlat(flat);
        lOP.updateFree(newDes);
        // lOP.updateTagged(tag);

        real_t Wp = lOP.evalObj();

        // W(c,x + eps*q)
        lOP.updateFlat(flat);
        lOP.updateFree(des);
        lOP.updateTagged(newTag);

        real_t Wq = lOP.evalObj();

        // W(c,x)
        lOP.updateFlat(flat);
        lOP.updateFree(des);
        // lOP.updateTagged(tag);

        real_t W = lOP.evalObj();

        real_t FDtagged = (Wq - W)/eps;

        real_t FDfree = (Wp -  W)/eps;

		real_t errorT = FDtagged - q.transpose()*gradTagged;
		real_t errorF = FDfree - p.transpose()*grad;
        gsInfo << FDfree << ", " << p.transpose()*grad << ", ";
        gsInfo << FDtagged << ", " << q.transpose()*gradTagged << "\n";

		ErrorT[i] = errorT;
		ErrorF[i] = errorF;

        // Test
	}

	gsVector<> rateT;
	rateT.setZero(n);
	rateT.segment(1,n-1) = log10(ErrorT.segment(1,n-1).array()/ErrorT.segment(0,n-1).array())/log10(2);
	gsVector<> rateF;
	rateF.setZero(n);
	rateF.segment(1,n-1) = log10(ErrorF.segment(1,n-1).array()/ErrorF.segment(0,n-1).array())/log10(2);
	gsMatrix<> disp(n,5);
	disp << Eps,ErrorT,rateT,ErrorF,rateF;
	gsInfo << "eps \tErrT \trateT \tErrF \trateF\n";
	gsInfo << disp << "\n";

}

void convergenceTestOfParaJacobianAll(gsWinslow &lOP){
	real_t liao = lOP.evalObj();
	gsVector<> grad = lOP.gradAll();

	gsMatrix<> hess = lOP.hessAll();

	gsInfo << "\n" << std::setprecision(10) << liao << "\n";
	// gsInfo << "\n" << grad << "\n";
    gsDebugVar(hess.rows());
    gsDebugVar(hess.cols());

	// std::srand((unsigned int) std::time(0));
	gsVector<> ran;
	ran.setRandom(lOP.n_flat);

	gsVector<> des = lOP.getFlat();
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
		perturp.setZero(lOP.n_flat);

		for(index_t i = 0; i < lOP.n_flat; i++){
			perturp[i] = ran[i];
		}

		perturp /= perturp.norm();

		perturp *= eps;

		gsVector<> newDes = des + perturp;


		lOP.updateFlat(newDes);

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

void convergenceTestOfParaHessianAll(gsWinslow &lOP){

    gsVector<> flat = lOP.getFlat();
	gsMatrix<> hess = lOP.hessAll();

	std::srand((unsigned int) std::time(0));
	gsVector<> ran1;
	ran1.setRandom(lOP.n_flat);
	gsVector<> ran2;
	ran2.setRandom(lOP.n_flat);

	index_t beg = 7;
	index_t n = 10;
	gsVector<> Eps(n);
	gsVector<> ErrorA(n);

	for(index_t i = 0; i < n; i++){
        lOP.updateFlat(flat);
		// Generate pertubation
		real_t eps = pow(2,-beg-i);
		Eps[i] = eps;

		gsVector<> p;
		p.setZero(lOP.n_flat);
		gsVector<> q;
		q.setZero(lOP.n_flat);

		for(index_t i = 0; i < lOP.n_flat; i++){
			p[i] = ran1[i];
			q[i] = ran2[i];
		}

        // W(c + eps*p)
        lOP.updateFlat(flat + eps*p);
        real_t Wp = lOP.evalObj();

        // W(c,x)
        lOP.updateFlat(flat);
        real_t W = lOP.evalObj();

        // W(c - eps*p)
        lOP.updateFlat(flat - eps*p);

        real_t Wmp = lOP.evalObj();

        real_t FDall = (Wp - 2*W + Wmp)/(eps*eps);

		real_t errorA = FDall - p.transpose()*hess*p;
        gsInfo << FDall << ", " << p.transpose()*hess*p << "\n";

		ErrorA[i] = errorA;
	}

	gsVector<> rateA;
	rateA.setZero(n);
	rateA.segment(1,n-1) = log10(ErrorA.segment(1,n-1).array()/ErrorA.segment(0,n-1).array())/log10(2);
	gsMatrix<> disp(n,3);
	disp << Eps,ErrorA,rateA;
	gsInfo << "eps \tErrA \trateA \n";
	gsInfo << disp << "\n";

}

void convergenceTestOfParaGradAll(gsWinslow &lOP){

	gsVector<> grad = lOP.gradAll();

	std::srand((unsigned int) std::time(0));
	gsVector<> ran_flat;
	ran_flat.setRandom(lOP.n_flat);

	gsVector<> ran_free;
	ran_free.setRandom(lOP.n_free);

	gsVector<> ran_tagged;
	ran_tagged.setRandom(lOP.n_tagged);

	gsVector<> des = lOP.getFree();
	gsVector<> tag = lOP.getTagged();
	gsVector<> flat = lOP.getFlat();

	index_t beg = 7;
	index_t n = 30;
	gsVector<> Eps(n);
	gsVector<> ErrorA(n);

	for(index_t i = 0; i < n; i++){
		// Generate pertubation
		real_t eps = pow(2,-beg-i);
		Eps[i] = eps;

		gsVector<> p;
		p.setZero(lOP.n_free);
		gsVector<> q;
		q.setZero(lOP.n_tagged);
		gsVector<> pA;
		pA.setZero(lOP.n_flat);


		for(index_t i = 0; i < lOP.n_flat; i++){
			pA[i] = ran_flat[i];
		}
		for(index_t i = 0; i < lOP.n_free; i++){
			p[i] = ran_free[i];
		}
		for(index_t i = 0; i < lOP.n_tagged; i++){
			q[i] = ran_tagged[i];
		}


		gsVector<> newFlat = flat + eps*pA;
		gsVector<> newDes = des + eps*p;
		gsVector<> newTag = tag + eps*q;

        // W(c + eps*p,x + eps*q)
        // lOP.updateFree(newDes);
        // lOP.updateTagged(newTag);
        lOP.updateFlat(newFlat);

        real_t Wpq = lOP.evalObj();

        // W(c,x)
        // lOP.updateFree(des);
        // lOP.updateTagged(tag);
        lOP.updateFlat(flat);

        real_t W = lOP.evalObj();

        real_t FDall = (Wpq - W)/eps;

        // lOP.updateTagged(q);
        // lOP.updateFree(p);
        // gsVector<> flat = lOP.getFlat();

		real_t errorA = FDall - pA.transpose()*grad;
        gsInfo << FDall << ", " << pA.transpose()*grad << "\n";

		ErrorA[i] = errorA;
	}

	gsVector<> rateA;
	rateA.setZero(n);
	rateA.segment(1,n-1) = log10(ErrorA.segment(1,n-1).array()/ErrorA.segment(0,n-1).array())/log10(2);
	gsMatrix<> disp(n,3);
	disp << Eps,ErrorA,rateA;
	gsInfo << "eps \tErrA \trateA \n";
	gsInfo << disp << "\n";

}

gsVector<> loadVec(index_t n,std::string name){
        gsInfo << "load from " << name << "\n";
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

void changeSignOfDetJ(gsGeometry<> & geom){//, index_t n, index_t m){
    geom.scale(-1,0);
    // gsMatrix<> cc = geom.coefs();
    //
    // gsInfo << "cc before: " << cc << "\n\n";
    // for (index_t d = 0; d < geom.targetDim(); d++){
    //     gsMatrix<> cc_d = reshape(cc.col(d),n,m);
    //     cc.col(d) = reshapeBack(cc_d.rowwise().reverse());
    // }
    // geom.setCoefs(cc);
    // gsInfo << "cc after: " << cc << "\n";

}

gsMultiPatch<> getGeometry(index_t n, index_t m, index_t degree){

    std::string folder;
    if (n == 4 && m == 4 && degree == 2){
        folder = "/parametrizations/para_x4_y4_p2_q2/";
    } else if (n == 5 && m == 4 && degree == 2){
        folder = "/parametrizations/para_x5_y4_p2_q2_circ/";
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

	gsInfo << coefs;

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

    // Change sign of determinant.
    for (index_t i = 0; i < patches.nBoxes(); i++){
        changeSignOfDetJ(patches.patch(i));
    }

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

gsMultiPatch<> getJigSaw(index_t i){

    index_t n = 10;
    index_t m = 10;
    index_t degree = 2;

    std::string folder = "/parametrizations/JigSaw/";

	// 1. construction of a knot vector for each direction
	gsKnotVector<> kv1(0, 1, n - degree - 1, degree + 1);
	gsKnotVector<> kv2(0, 1, m - degree - 1, degree + 1);
	// 2. construction of a basis
	gsTensorBSplineBasis<2, real_t> basis(kv1, kv2);
	// 3. construction of a coefficients
	gsMatrix<> greville = basis.anchors();
	gsMatrix<> coefs (greville.cols(), 2);

    gsInfo << BASE_FOLDER + folder + "Jigsaw1.txt \n";
    if (i == 1){
	   readFromTxt(BASE_FOLDER + folder + "Jigsaw1.txt", coefs);
    } else {
	   readFromTxt(BASE_FOLDER + folder + "Jigsaw2.txt", coefs);
    }
    // gsDebugVar(coefs.rows());
    // gsDebugVar(coefs.cols());
    // gsDebugVar(coefs);
    // exit(0);
	// gsInfo << coefs;

	// 4. putting basis and coefficients toghether
	gsTensorBSpline<2, real_t>  left(basis, coefs);
    changeSignOfDetJ(left);

	//! [Geometry data]
	// Define Geometry, must be a gsMultiPatch object

	gsMultiPatch<> patches = gsMultiPatch<>(left);

    // Change sign of determinant.
    // for (index_t i = 0; i < patches.nBoxes(); i++){
        // changeSignOfDetJ(patches.patch(i));
    // }

	double tol = 1e-2;
	patches.computeTopology(tol,true);
	patches.closeGaps(tol);

	return patches;

}

gsMultiPatch<> getJigSaw3D(index_t i){

    index_t n;
    if (i == 1)
        n = 10;
    if (i == 2)
        n = 6;
    if (i == 3)
        n = 5;
    index_t degree = 2;

    std::string folder = "/parametrizations/JigSaw/";

	// 1. construction of a knot vector for each direction
	gsKnotVector<> kv1(0, 1, n - degree - 1, degree + 1);
	gsKnotVector<> kv2(0, 1, n - degree - 1, degree + 1);
	gsKnotVector<> kv3(0, 1, n - degree - 1, degree + 1);
	// 2. construction of a basis
	gsTensorBSplineBasis<3, real_t> basis(kv1, kv2, kv3);
	// 3. construction of a coefficients
	gsMatrix<> greville = basis.anchors();
	gsMatrix<> coefs (greville.cols(), 3);

    gsInfo << BASE_FOLDER + folder + "jigsaw3d.txt \n";
    if (i == 1)
	   readFromTxt(BASE_FOLDER + folder + "jigsaw3d.txt", coefs);
    if (i == 2)
	   readFromTxt(BASE_FOLDER + folder + "jigsaw3d_small.txt", coefs);
    if (i == 3)
	   readFromTxt(BASE_FOLDER + folder + "example.txt", coefs);
    // gsDebugVar(coefs.rows());
    // gsDebugVar(coefs.cols());
    // gsDebugVar(coefs);
    // exit(0);
	// gsInfo << coefs;

	// 4. putting basis and coefficients toghether
	gsTensorBSpline<3, real_t>  left(basis, coefs);
    // changeSignOfDetJ(left);

	//! [Geometry data]
	// Define Geometry, must be a gsMultiPatch object

	gsMultiPatch<> patches = gsMultiPatch<>(left);

    // Change sign of determinant.
    // for (index_t i = 0; i < patches.nBoxes(); i++){
    //     changeSignOfDetJ(patches.patch(i));
    // }

	double tol = 1e-2;
	patches.computeTopology(tol,true);
	patches.closeGaps(tol);

	return patches;

}

gsMultiPatch<> getSeastar(){

    index_t n = 8;
    index_t m = 8;
    index_t degree = 3;

    std::string folder = "/parametrizations/Seastar/";

	// 1. construction of a knot vector for each direction
	gsKnotVector<> kv1(0, 1, n - degree - 1, degree + 1);
	gsKnotVector<> kv2(0, 1, m - degree - 1, degree + 1);
	// 2. construction of a basis
	gsTensorBSplineBasis<2, real_t> basis(kv1, kv2);
	// 3. construction of a coefficients
	gsMatrix<> greville = basis.anchors();
	gsMatrix<> coefs (greville.cols(), 2);

    gsInfo << BASE_FOLDER + folder + "seastar.txt \n";
	readFromTxt(BASE_FOLDER + folder + "seastar.txt", coefs);

    // gsDebugVar(coefs.rows());
    // gsDebugVar(coefs.cols());
    // gsDebugVar(coefs);
    // exit(0);
	// gsInfo << coefs;

	// 4. putting basis and coefficients toghether
	gsTensorBSpline<2, real_t>  left(basis, coefs);
    // changeSignOfDetJ(left);

	//! [Geometry data]
	// Define Geometry, must be a gsMultiPatch object

	gsMultiPatch<> patches = gsMultiPatch<>(left);

    // Change sign of determinant.
    // for (index_t i = 0; i < patches.nBoxes(); i++){
        // changeSignOfDetJ(patches.patch(i));
    // }

	double tol = 1e-2;
	patches.computeTopology(tol,true);
	patches.closeGaps(tol);

	return patches;

}

gsMatrix<> reshape(gsVector<> vec, index_t n, index_t m){
    gsMatrix<> out(n,m);

    index_t c = 0;
    for(index_t i = 0; i < n; i++){
        for(index_t j = 0; j < m; j++){
            gsInfo << "c " << c;
            out(i,j) = vec[c++];
        }
        gsInfo << "\n";
    }
    return out;
}

gsVector<> reshapeBack(gsMatrix<> mat){
    index_t n = mat.rows();
    index_t m = mat.cols();
    gsVector<> out(n*m);

    index_t c = 0;
    for(index_t i = 0; i < n; i++){
        for(index_t j = 0; j < m; j++){
            out[c++] = mat(i,j);
        }
    }
    return out;
}

void testOfParametrizations(std::string output, index_t quA, index_t quB){

    index_t numRefine = 1;
    index_t n_tests = 6;

    // Setup tests
    std::vector< gsMultiPatch<> > P(n_tests);
    std::vector< std::string > names(n_tests);

    // Unit square
    P[0] = gsMultiPatch<>(*gsNurbsCreator<>::BSplineSquare(2));
    names[0] = "/UnitSq/";
    for(int i = 0; i < numRefine; i++){
    	P[0].uniformRefine();
    }

    // Rotated rectangle
    P[1] = gsMultiPatch<>(*gsNurbsCreator<>::BSplineSquare(2));
    P[1].patch(0).scale(2,1);
    P[1].patch(0).rotate(M_PI_4);
    gsVector<> v(2);
    v << 1,2;
    P[1].patch(0).translate(v);

    names[1] = "/Rectangle/";
    for(int i = 0; i < numRefine; i++){
    	P[1].uniformRefine();
    }


    // JigSaws
    P[2] = getJigSaw(1);
    P[2].patch(0).scale(0.1);
    names[2] = "/JigSaw1/";
    P[3] = getJigSaw(2);
    P[3].patch(0).scale(0.1);
    names[3] = "/JigSaw2/";

    P[4] = getSeastar();
    names[4] = "/Seastar/";

    // Initial design
    P[5] = getGeometry(5,4,2);
    names[5] = "/InitDesign/";


    // For all test but the last one
    for (index_t t = 0; t < n_tests - 1; t++){

        gsShapeOptLog slog(output + names[t]);

        gsInfo << P[t] << "\n";
        gsInfo << P[t].patch(0).basis() << "\n";

        gsDetJacConstraint dJC(&P[t]);

        slog << "min d : " << dJC.evalCon().minCoeff() << "\n";
        slog << "max d : " << dJC.evalCon().maxCoeff() << "\n";

        slog << "Spring Method: \n";
        gsInfo << "Spring Method: \n";

        gsSpringMethod sM(&P[t]);
        sM.update();

        std::string name = "Spring";
        slog.plotInParaview(P[t],name);
        slog.saveVec(sM.getFlat(),name);

        slog << "min d : " << dJC.evalCon().minCoeff() << "\n";
        slog << "max d : " << dJC.evalCon().maxCoeff() << "\n";

        slog << "MaxDetJac Method: \n";
        gsInfo << "MaxDetJac Method: \n";

        gsMaxDetJac mDJ(&P[t]);
        mDJ.update();

        name = "MaxDetJac";
        slog.plotInParaview(P[t],name);
        slog.saveVec(sM.getFlat(),name);

        slog << "min d : " << dJC.evalCon().minCoeff() << "\n";
        slog << "max d : " << dJC.evalCon().maxCoeff() << "\n";

        slog << "Winslow: \n";
        gsInfo << "Winslow: \n";

        sM.update();

        gsWinslow winslow(&P[t],true);
        winslow.setQuad(quA,quB);
        winslow.update();

        name = "Winslow";
        slog.plotInParaview(P[t],name);
        slog.saveVec(winslow.getFlat(),name);

        slog << "min d : " << dJC.evalCon().minCoeff() << "\n";
        slog << "max d : " << dJC.evalCon().maxCoeff() << "\n";

        sM.update();
        slog << "Harmonic 11: \n";
        gsInfo << "Harmonic 11: \n";

        gsHarmonic harmonic(&P[t],false);
        harmonic.setQuad(quA,quB);
        harmonic.update();

        slog << "min d : " << dJC.evalCon().minCoeff() << "\n";
        slog << "max d : " << dJC.evalCon().maxCoeff() << "\n";

        name = "Harmonic_11";
        slog.plotInParaview(P[t],name);
        slog.saveVec(harmonic.getFlat(),name);

        sM.update();
        slog << "Harmonic 00: \n";
        gsInfo << "Harmonic 00: \n";
        harmonic.setLambdas(0,0);
        harmonic.update();

        slog << "min d : " << dJC.evalCon().minCoeff() << "\n";
        slog << "max d : " << dJC.evalCon().maxCoeff() << "\n";

        name = "Harmonic_00";
        slog.plotInParaview(P[t],name);
        slog.saveVec(harmonic.getFlat(),name);
    }

    index_t t = n_tests - 1;

    gsShapeOptLog slog(output + names[t]);
    gsOptAntenna optA(&P[t],1,&slog,0,quA,quB);

    gsInfo << P[t] << "\n";
    gsInfo << P[t].patch(0).basis() << "\n";

    gsDetJacConstraint dJC(&P[t]);

    slog << "min d : " << dJC.evalCon().minCoeff() << "\n";
    slog << "max d : " << dJC.evalCon().maxCoeff() << "\n";

    slog << "Spring Method: \n";

    gsSpringMethod sM(&P[t],optA.mappers());
    sM.update();

    std::string name = "Spring";
    slog.plotInParaview(P[t],name);
    slog.saveVec(sM.getFlat(),name);

    slog << "MaxDetJac Method: \n";
    gsInfo << "MaxDetJac Method: \n";

    gsMaxDetJac mDJ(&P[t],optA.mappers());
    mDJ.update();

    name = "MaxDetJac";
    slog.plotInParaview(P[t],name);
    slog.saveVec(sM.getFlat(),name);

    slog << "min d : " << dJC.evalCon().minCoeff() << "\n";
    slog << "max d : " << dJC.evalCon().maxCoeff() << "\n";

    slog << "min d : " << dJC.evalCon().minCoeff() << "\n";
    slog << "max d : " << dJC.evalCon().maxCoeff() << "\n";

    slog << "Winslow: \n";
    gsWinslow winslow(&P[t],optA.mappers(),true);
    winslow.setQuad(quA,quB);
    winslow.update();

    name = "Winslow";
    slog.plotInParaview(P[t],name);
    slog.saveVec(winslow.getFlat(),name);

    slog << "min d : " << dJC.evalCon().minCoeff() << "\n";
    slog << "max d : " << dJC.evalCon().maxCoeff() << "\n";

    sM.update();
    slog << "Harmonic: \n";

    gsHarmonic harmonic(&P[t],optA.mappers(),false);
    harmonic.setQuad(quA,quB);
    harmonic.update();

    slog << "min d : " << dJC.evalCon().minCoeff() << "\n";
    slog << "max d : " << dJC.evalCon().maxCoeff() << "\n";

    name = "Harmonic_11";
    slog.plotInParaview(P[t],name);
    slog.saveVec(harmonic.getFlat(),name);

    sM.update();
    harmonic.setLambdas(0,0);
    harmonic.update();

    slog << "min d : " << dJC.evalCon().minCoeff() << "\n";
    slog << "max d : " << dJC.evalCon().maxCoeff() << "\n";

    name = "Harmonic_00";
    slog.plotInParaview(P[t],name);
    slog.saveVec(harmonic.getFlat(),name);



}

void testOfParametrizations_winslow_inf(std::string output, index_t quA, index_t quB){

    index_t numRefine = 1;
    index_t n_tests = 6;

    // Setup tests
    std::vector< gsMultiPatch<> > P(n_tests);
    std::vector< std::string > names(n_tests);

    // Unit square
    P[0] = gsMultiPatch<>(*gsNurbsCreator<>::BSplineSquare(2));
    names[0] = "/UnitSq/";
    for(int i = 0; i < numRefine; i++){
    	P[0].uniformRefine();
    }

    // Rotated rectangle
    P[1] = gsMultiPatch<>(*gsNurbsCreator<>::BSplineSquare(2));
    P[1].patch(0).scale(2,1);
    P[1].patch(0).rotate(M_PI_4);
    gsVector<> v(2);
    v << 1,2;
    P[1].patch(0).translate(v);

    names[1] = "/Rectangle/";
    for(int i = 0; i < numRefine; i++){
    	P[1].uniformRefine();
    }

    P[2] = getJigSaw(2);
    P[2].patch(0).scale(0.1);
    names[2] = "/Jigsaw2/";

    P[3] = getSeastar();
    names[3] = "/Seastar/";

    // JigSaws
    P[4] = getJigSaw(1);
    P[4].patch(0).scale(0.1);
    names[4] = "/Jigsaw1/";
    P[4].uniformRefine();


    // Initial design
    P[5] = getGeometry(5,4,2);
    names[5] = "/InitDesign/";

    std::string name;

    // For all test but the last one
    for (index_t t = 0; t < n_tests - 1; t++){
        continue; exit(0);
        gsShapeOptLog slog(output + names[t]);

        gsInfo << P[t] << "\n";
        gsInfo << P[t].patch(0).basis() << "\n";

        gsDetJacConstraint dJC(&P[t]);

        slog << "MaxDetJac Method: \n";
        gsInfo << "MaxDetJac Method: \n";

        gsMaxDetJac mDJ(&P[t],true);
        mDJ.update();

        name = "MaxDetJac";
        slog.plotInParaview(P[t],name);
        gsVector<> mdj = mDJ.getFlat();
        slog.saveVec(mdj,name);

        slog << "min d : " << dJC.evalCon().minCoeff() << "\n";
        slog << "max d : " << dJC.evalCon().maxCoeff() << "\n";

        slog << "\n Winslow all constraints: \n";
        gsInfo << "\n --- Winslow all constraints: ---- \n\n";

        mDJ.updateFlat(mdj);

        gsWinslow winslow(&P[t],true,true);
        winslow.setQuad(quA,quB);
        winslow.update();

        name = "Winslow_allConstraints";
        slog.plotInParaview(P[t],name);
        slog.saveVec(winslow.getFlat(),name);

        slog << "min d : " << dJC.evalCon().minCoeff() << "\n";
        slog << "max d : " << dJC.evalCon().maxCoeff() << "\n";

        // slog << "\n Winslow no constraints: \n";
        // gsInfo << "\n --- Winslow no constraints: ---- \n\n";
        //
        // mDJ.updateFlat(mdj);
        //
        // gsWinslow winslow_noCon(&P[t],false);
        // winslow_noCon.setQuad(quA,quB);
        // winslow_noCon.update();
        //
        // name = "Winslow_noConstraints";
        // slog.plotInParaview(P[t],name);
        // slog.saveVec(winslow_noCon.getFlat(),name);
        //
        // slog << "min d : " << dJC.evalCon().minCoeff() << "\n";
        // slog << "max d : " << dJC.evalCon().maxCoeff() << "\n";

        slog << "\n Winslow, check for inf, all constraints: \n";
        gsInfo << "Winslow, check for inf, all constraints: \n";

        mDJ.updateFlat(mdj);

        gsWinslow winslow_check(&P[t],true,true,true,0);
        winslow_check.setQuad(quA,quB);
        winslow_check.update();

        name = "Winslow_check_allConstraints";
        slog.plotInParaview(P[t],name);
        slog.saveVec(winslow_check.getFlat(),name);

        slog << "min d : " << dJC.evalCon().minCoeff() << "\n";
        slog << "max d : " << dJC.evalCon().maxCoeff() << "\n";

        slog << "\n Winslow, check for inf, no constraints: \n";
        gsInfo << "Winslow, check for inf, no constraints: \n";

        mDJ.updateFlat(mdj);

        gsWinslow winslow_check_noCon(&P[t],false,false,true,0);
        winslow_check_noCon.setQuad(quA,quB);
        winslow_check_noCon.update();

        name = "Winslow_check_noConstraints";
        slog.plotInParaview(P[t],name);
        slog.saveVec(winslow_check_noCon.getFlat(),name);

        slog << "min d : " << dJC.evalCon().minCoeff() << "\n";
        slog << "max d : " << dJC.evalCon().maxCoeff() << "\n";
    }

    index_t t = n_tests - 1;


    gsShapeOptLog slog(output + names[t]);
    gsOptAntenna optA(&P[t],1,&slog,0,quA,quB);

    gsInfo << P[t] << "\n";
    gsInfo << P[t].patch(0).basis() << "\n";

    gsDetJacConstraint dJC(&P[t]);

    slog << "MaxDetJac Method: \n";
    gsInfo << "MaxDetJac Method: \n";

    gsMaxDetJac mDJ(&P[t],optA.mappers());
    mDJ.update();

    name = "MaxDetJac";
    slog.plotInParaview(P[t],name);
    gsVector<> mdj = mDJ.getFlat();
    slog.saveVec(mdj,name);

    slog << "min d : " << dJC.evalCon().minCoeff() << "\n";
    slog << "max d : " << dJC.evalCon().maxCoeff() << "\n";

    slog << "\n Winslow all constraints: \n";
    gsInfo << "Winslow all constraints: \n";

    mDJ.updateFlat(mdj);

    gsWinslow winslow(&P[t],optA.mappers(),true,true);
    winslow.setQuad(quA,quB);
    winslow.update();

    name = "Winslow_allConstraints";
    slog.plotInParaview(P[t],name);
    slog.saveVec(winslow.getFlat(),name);

    slog << "min d : " << dJC.evalCon().minCoeff() << "\n";
    slog << "max d : " << dJC.evalCon().maxCoeff() << "\n";

    // slog << "\n Winslow no constraints: \n";
    // gsInfo << "Winslow no constraints: \n";
    //
    // mDJ.updateFlat(mdj);
    //
    // gsWinslow winslow_noCon(&P[t],optA.mappers(),false);
    // winslow_noCon.setQuad(quA,quB);
    // winslow_noCon.update();
    //
    // name = "Winslow_noConstraints";
    // slog.plotInParaview(P[t],name);
    // slog.saveVec(winslow_noCon.getFlat(),name);
    //
    // slog << "min d : " << dJC.evalCon().minCoeff() << "\n";
    // slog << "max d : " << dJC.evalCon().maxCoeff() << "\n";

    slog << "\n Winslow, check for inf, all constraints: \n";
    gsInfo << "Winslow, check for inf, all constraints: \n";

    mDJ.updateFlat(mdj);

    gsWinslow winslow_check(&P[t],optA.mappers(),true,true,true,0);
    winslow_check.setQuad(quA,quB);
    winslow_check.update();

    name = "Winslow_check_allConstraints";
    slog.plotInParaview(P[t],name);
    slog.saveVec(winslow_check.getFlat(),name);

    slog << "min d : " << dJC.evalCon().minCoeff() << "\n";
    slog << "max d : " << dJC.evalCon().maxCoeff() << "\n";

    slog << "\n Winslow, check for inf, no constraints: \n";
    gsInfo << "Winslow, check for inf, no constraints: \n";

    mDJ.updateFlat(mdj);

    gsWinslow winslow_check_noCon(&P[t],optA.mappers(),false,false,true,0);
    winslow_check_noCon.setQuad(quA,quB);
    winslow_check_noCon.update();

    name = "Winslow_check_noConstraints";
    slog.plotInParaview(P[t],name);
    slog.saveVec(winslow_check_noCon.getFlat(),name);

    slog << "min d : " << dJC.evalCon().minCoeff() << "\n";
    slog << "max d : " << dJC.evalCon().maxCoeff() << "\n";

}


void testOfAggregatedConstraints(std::string output, index_t quA, index_t quB){

    index_t numRefine = 1;
    index_t n_tests = 6;

    index_t n_alpha = 100;
    gsVector<> Alpha(n_alpha);
    Alpha.setLinSpaced(n_alpha,-n_alpha+1,0);

    // Setup tests
    std::vector< gsMultiPatch<> > P(n_tests);
    std::vector< std::string > names(n_tests);

    // Unit square
    P[0] = gsMultiPatch<>(*gsNurbsCreator<>::BSplineSquare(2));
    names[0] = "/UnitSq/";
    for(int i = 0; i < numRefine; i++){
    	P[0].uniformRefine();
    }

    // Rotated rectangle
    P[1] = gsMultiPatch<>(*gsNurbsCreator<>::BSplineSquare(2));
    P[1].patch(0).scale(2,1);
    P[1].patch(0).rotate(M_PI_4);
    gsVector<> v(2);
    v << 1,2;
    P[1].patch(0).translate(v);

    names[1] = "/Rectangle/";
    for(int i = 0; i < numRefine; i++){
    	P[1].uniformRefine();
    }

    // JigSaws
    P[2] = getJigSaw(1);
    P[2].patch(0).scale(0.1);
    names[2] = "/Jigsaw1/";
    P[3] = getJigSaw(2);
    P[3].patch(0).scale(0.1);
    names[3] = "/Jigsaw2/";

    P[4] = getSeastar();
    names[4] = "/Seastar/";

    // Initial design
    P[5] = getGeometry(5,4,2);
    names[5] = "/InitDesign/";


    for (index_t t = 0; t < n_tests; t++){

        gsShapeOptLog slog(output + names[t]);

        gsInfo << P[t] << "\n";
        gsInfo << P[t].patch(0).basis() << "\n";

        gsDetJacConstraint dJC(&P[t]);

        std::string name = "d.txt";
        slog.saveVec(dJC.evalCon(),name);

        slog << "alpha min_d agg\n";
        real_t mind = dJC.evalCon().minCoeff();

        for (real_t alpha: Alpha){
            slog << alpha << " ";
            slog << mind << " ";
            gsAggFun fun(dJC.n_constraints,alpha);
            gsAggregatedConstraint aC(&P[t], &dJC, &fun);
            slog << aC.evalCon()[0] << "\n";
        }

    }


}

void testOfAggregatedConstraints2nd(std::string output, index_t quA, index_t quB){

    index_t numRefine = 1;
    index_t n_tests = 5;

    index_t n_alpha = 4;
    gsVector<> Alpha(n_alpha);
    Alpha << -4, -8, -16, -32;

    // Setup tests
    std::vector< gsMultiPatch<> > P(n_tests);
    std::vector< std::string > names(n_tests);

    // Unit square
    P[0] = gsMultiPatch<>(*gsNurbsCreator<>::BSplineSquare(2));
    names[0] = "/UnitSq/";
    for(int i = 0; i < numRefine; i++){
    	P[0].uniformRefine();
    }

    // Rotated rectangle
    P[1] = gsMultiPatch<>(*gsNurbsCreator<>::BSplineSquare(2));
    P[1].patch(0).scale(2,1);
    P[1].patch(0).rotate(M_PI_4);
    gsVector<> v(2);
    v << 1,2;
    P[1].patch(0).translate(v);

    names[1] = "/Rectangle/";
    for(int i = 0; i < numRefine; i++){
    	P[1].uniformRefine();
    }

    // JigSaws
    // P[2] = getJigSaw(1);
    // P[2].patch(0).scale(0.1);
    // names[2] = "/Jigsaw1/";
    P[3] = getJigSaw(2);
    P[3].patch(0).scale(0.1);
    names[3] = "/Jigsaw2/";

    P[2] = getSeastar();
    names[2] = "/Seastar/";

    // Initial design
    P[4] = getGeometry(5,4,2);
    names[4] = "/InitDesign/";


    for (index_t t = 1; t < n_tests; t++){


        gsInfo << P[t] << "\n";
        gsInfo << P[t].patch(0).basis() << "\n";


        gsMultiPatch<> start(P[t]);

        gsInfo << output + names[t] << "\n\n";
        gsShapeOptLog slog(output + names[t]);

        gsVector<> flat;
        if (t == 4){
            gsOptAntenna optA(&start,1,&slog,0,quA,quB);

            gsMaxDetJac mdj(&start,optA.mappers());
            mdj.update();

            flat = mdj.getFlat();
            std::string name = "cps.txt";
            slog.saveVec(flat,name);
            name = "startGuess";
            slog.plotInParaview(start,name);

            gsModLiao winslow(&start,optA.mappers());
            winslow.update();
        } else {
            gsMaxDetJac mdj(&start);
            mdj.update();

            flat = mdj.getFlat();
            std::string name = "cps.txt";
            slog.saveVec(flat,name);
            name = "startGuess";
            slog.plotInParaview(start,name);

            gsWinslow winslow(&start);
            winslow.update();
        }

        gsWinslow win(&start);

        std::string name = "cps";
        slog.saveVec(win.getFlat(),name);
        name = "mp";
        slog.plotInParaview(start,name);

        continue;
        exit(0);
        for (real_t alpha: Alpha){
            gsMultiPatch<> mp(P[t]);

            std::ostringstream strs;
            strs << alpha << "/";
            std::string str = strs.str();


            gsDetJacConstraint dJC(&mp);

            gsInfo << output + names[t] + str << "\n\n";
            gsShapeOptLog slog(output + names[t] + str);

            gsAggFun fun(dJC.n_constraints,alpha);
            gsAggregatedConstraint aC(&mp, &dJC, &fun);

            if (t == 4){
                gsOptAntenna optA(&mp,1,&slog,0,quA,quB);

                gsWinslow winslow(&mp,optA.mappers(),&aC);
                winslow.updateFlat(flat);

                slog << dJC.evalCon().minCoeff() << " " << winslow.m_dJC->evalCon()[0] << "\n";

                winslow.setLog(&slog, &dJC);
                winslow.update();
            } else {
                gsWinslow winslow(&mp,&aC);
                winslow.updateFlat(flat);

                slog << dJC.evalCon().minCoeff() << " " << winslow.m_dJC->evalCon()[0] << "\n";

                winslow.setLog(&slog, &dJC);
                winslow.update();
            }

            gsWinslow winslow(&mp);

            std::string name = "cps.txt";
            slog.saveVec(winslow.getFlat(),name);
            name = "mp";
            slog.plotInParaview(mp,name);
        }

    }


}

void testOfAggregatedConstraints3rd(std::string output, index_t quA, index_t quB, index_t param, index_t numref, index_t maxiter, index_t n, index_t m, index_t degree){

    index_t n_tests = 5;

    index_t n_alpha = 2;
    gsVector<> Alpha(n_alpha);
    Alpha << -32;

    for (real_t alpha: Alpha){
        gsMultiPatch<> mp = getGeometry(n,m,degree);

        std::ostringstream strs;
        strs << alpha << "/";
        std::string str = strs.str();

        gsDetJacConstraint dJC(&mp);

        gsInfo << output + str << "\n\n";
        gsShapeOptLog slog(output + str,true,true,false);

        if (alpha == 1234)
        {
            gsOptAntenna optA(&mp,numref,&slog,param,quA,quB);

            if (param == 0) // Spring method
            {
                continue;
                optA.solve();
            } else {
                optA.runOptimization(maxiter);
            }
        } else {
            gsAggFun fun(dJC.n_constraints,alpha);
            gsAggregatedConstraint aC(&mp, &dJC, &fun);

            gsOptAntenna optA(&mp,numref,&slog,&aC,param,quA,quB);

            if (param == 0) // Spring method
            {
                optA.solve();
            } else {
                optA.runOptimization_aggregatedConstraints(maxiter);
            }
        }
    }


}

gsMultiPatch<> get3DGeometry(){
    gsMultiPatch<> patches;
    gsFileData<> data("geometries3D/halfCube.xml");

    gsInfo  <<"* There is "<< data.count< gsGeometry<> >() <<" "
            <<data.type< gsGeometry<> >()<<" "<< data.tag< gsGeometry<> >()
            <<" in the file.\n";
    // we start at id == 2 since we want a hollow cube. For a patch in the middle
    for(index_t i = 1; i <= data.count< gsGeometry<> >(); i++){
        gsGeometry<>::uPtr o = data.getId< gsGeometry<> >(i);
        o->degreeElevate(1);
        o->uniformRefine();
        patches.addPatch(*o);
    }

    // FIXIT make this part faster, it should be unecessary to construct a gsDetJacConstraint to do this
    gsDetJacConstraint dJC(&patches);

    for (index_t p = 0; p < patches.nBoxes(); p++){
        index_t sign = dJC.getSignOfPatch(p);
        gsInfo << "Sign of patch " << p << " is " << sign << "\n";
        if (sign < 0)
            changeSignOfDetJ(patches.patch(p));
    }

    patches.computeTopology();

    return patches;
}

void generateData02953Project(gsMultiPatch<> & mp, std::string folder){
    // Setup classes
    gsDetJacConstraint dJC(&mp);
    gsWinslow winslow(&mp,false);         // We use the default mappers

    // eval once to setup solver
    gsInfo << " min d : " << dJC.evalCon().minCoeff() << "\n";
    gsInfo << " max d : " << dJC.evalCon().maxCoeff() << "\n";

    // winslow.setQuad(9,10);
    gsInfo << " winslow : " << winslow.evalObj() << "\n";
    // gsInfo << "d \n" << dJC.evalCon() << "\n";
    // exit(0);
    // Write M to file
    gsSparseMatrix<> M = dJC.getMassMatrix(0);
    saveSparseMat(M,BASE_FOLDER + folder + "M.txt");

    // Save the m_mappers in file
    gsVector<index_t> vecx = winslow.mappers()[0].asVector();
    saveVec(vecx, BASE_FOLDER + folder + "px.txt");
    gsVector<index_t> vecy = winslow.mappers()[1].asVector();
    saveVec(vecy, BASE_FOLDER + folder + "py.txt");

    // Save all cps
    gsVector<> cps = winslow.getControlPoints();
    saveVec(cps, BASE_FOLDER + folder + "cps_global.txt");
    cps = winslow.getFlat();
    saveVec(cps, BASE_FOLDER + folder + "cps_flat.txt");
    cps = winslow.getFree();
    saveVec(cps, BASE_FOLDER + folder + "cps_free.txt");
    cps = winslow.getTagged();
    saveVec(cps, BASE_FOLDER + folder + "cps_tagged.txt");


    // Generate Q_i and write to file
    // gsTensorBSplineBasis<2,real_t> bas = dJC.m_detJacBasis.basis(0);
    index_t size = dJC.m_detJacBasis.basis(0).size();
    gsInfo << "size of detJacBasis : " << size << "\n";

    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;

    for (index_t i = 0; i < size; i++){
        gsExprAssembler<> A(1,1);

        gsMultiBasis<> dbasis(mp);

        gsBasisFun<> fun = dJC.m_detJacBasis.basis(0).function(i);

        A.setIntegrationElements(dJC.m_detJacBasis);
        gsExprEvaluator<> ev(A);

        space u = A.getSpace(dbasis);

        variable f = A.getCoeff(fun);

        gsFunctionExpr<> ff("y","-x",2);
        variable var = A.getCoeff(ff);
        auto mat = fjac(var);

        A.initSystem();
        A.assemble(f.val()*grad(u) * mat * grad(u).tr());
        gsInfo << "nnz in A " << A.matrix().nonZeros() << "\n";

        saveSparseMat(A.matrix(),BASE_FOLDER + folder + "Q_" + std::to_string(i) + ".txt");
    }


    // Write initial design
}

int main(int argc, char* argv[]){
gsInfo <<  "Hello G+Smo.\n";

// Parse command line
std::string output("");
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
bool useLag = false;

int quA = 2;
int quB = 2;

real_t lambda_1 = 1;
real_t lambda_2 = 1;

int startDes = -1;

real_t alpha = 0; // Parameter for constraint aggregation
real_t eps = 0; // Parameter for constraint aggregation

gsCmdLine cmd("A test of lumped mass matricies");
cmd.addInt("p", "degree", "Degree of B-Splines.", degree);
cmd.addInt("a", "param", "Parametrization: 0: spring, 1: modLiao, 2: winslow, 3: liao, 4: harmonic, 5: use regularization and perform optimization on all cps", param);
cmd.addInt("r", "numberRefine", "Number of refinements", numRefine);
cmd.addInt("n", "nx", "Number of splines in first direction", nx);
cmd.addInt("m", "ny", "Number of splines in second direction", ny);
cmd.addInt("i", "maxIter", "Maximal number of reparametrizations", maxiter);

cmd.addInt("A", "quA", "quA", quA);
cmd.addInt("B", "quB", "quB", quB);

cmd.addReal("1", "lambda1", "lambda1", lambda_1);
cmd.addReal("2", "lambda2", "lambda2", lambda_2);

cmd.addString("o","output","Name of the output folder (relative to BASE_FOLDER)",output);

cmd.addSwitch("plot", "Create a ParaView visualization file of designs", plotDesign);
cmd.addSwitch("noDJC", "Dont use det jac constraints", useDJC);
cmd.addSwitch("useLag", "Use Lagrangian for linearizations", useLag);
cmd.addSwitch("plotmag", "Create a ParaView visualization file of magnitude (obj function)", plotMagnitude);
cmd.addSwitch("plotsol", "Create a ParaView visualization file of solution (real and imag)", plotSolution);
cmd.addSwitch("savecps", "save controlpoinst for design during optimization", saveCps);

cmd.addInt("s", "startDes", "design number to start from", startDes);

cmd.addReal("l", "alpha", "parameter for smoothed maximum in constraint aggregation", alpha);
cmd.addReal("e", "eps", "parameter for Regularization with winslow", eps);

cmd.getValues(argc,argv);

char buffer [50];
std::sprintf(buffer,"p = %d \n n = %d\n",degree,numRefine);
gsInfo << buffer;

gsShapeOptLog slog(output);

// Test 3D PML and water waves
if (false) {
    gsMultiPatch<> mp = get3DGeometry();

    gsStateEquationPotWaves SE(&mp,3);

    std::string name = "/../results/PML3D/mp";
    slog.plotInParaview(mp,name);
    exit(0);

}

gsMultiPatch<> mp = getGeometry(nx,ny,degree);
// gsMultiPatch<> mp = getJigSaw(1);
// gsMultiPatch<> mp = getSeastar();
// gsMultiPatch<> mp = getJigSaw3D();

// gsMultiPatch<> patches;
// for (index_t pn = 0; pn < mp.nBoxes(); pn ++){
//     gsTensorBSpline<2,real_t> geom = dynamic_cast<gsTensorBSpline<2,real_t>&>(mp.patch(pn));
//     gsTHBSpline<2> Hgeom(geom);
//     patches.addPatch(Hgeom);
// }
// patches.computeTopology();
// mp = patches;

gsMultiPatch<> patches = mp;

// gsMultiPatch<> patches3d = get3DGeometry();
// gsMultiPatch<> patches(patches3d.patch(2));

// gsMultiPatch<> patches = seastar;

// patches.patch(0).scale(0.1);

// gsMultiPatch<> patches = gsMultiPatch<>(*gsNurbsCreator<>::BSplineSquare(degree));
//
// for(int i = 0; i < numRefine; i++){
	// patches.uniformRefine();
// }

gsInfo << "The domain is a "<< patches <<"\n";

// gsStateEquationAntenna SE(&patches, 1);
// SE.printConstants();

// gsGeometry<>::uPtr gg = gsNurbsCreator<>::BSplineSquare();
// gsHBSplineBasis<2> hb(gg->basis());
// gsMatrix<> mm(hb.size(),1); mm.setZero();
// gsHBSpline<2> bb(hb,mm);
// // gsHBSpline<2> bb(gg->basis(),mm);
// exit(0);

// test winslow derivatives
if (false) {
    gsWinslow winslow(&mp,false);
    convergenceTestOfParaJacobian(winslow);
    exit(0);
}

// test harmonic derivatives
if (false) {
    gsHarmonic harmonic(&mp,false);
    convergenceTestOfParaJacobian(harmonic);
    exit(0);
}

if (startDes == 11) {
    gsMultiPatch<> jigsaw = getJigSaw3D(2); //
    mp = jigsaw;

    gsWinslowWithDeriv win(&mp,false,false,true,0);

    std::stringstream stream_in;
    stream_in << BASE_FOLDER << output << "cps" << "_0_" << maxiter << ".txt";
    std::string str_in = stream_in.str();
    gsMatrix<> mat = win.getFlat();
    readFromTxt(str_in,mat);
    // gsInfo << mat;

    win.updateFlat(mat);

    gsDetJacConstraint dJC(&mp,true);
    dJC.evalCon();
    gsShapeOptLog slog1(output,true,false,false);

    // gsOptAntenna optA(&mp,numRefine,&slog1,0,quA,quB);

    std::stringstream stream;
    stream << BASE_FOLDER << output << "detJ";
    std::string str = stream.str();
    dJC.plotDetJ(str);
    exit(0);



}

// calculate detJ constraints for different designs
if (startDes == 10) {
    // gsMultiPatch<> jigsaw = getJigSaw3D(3); //
    // mp = jigsaw;

    gsDetJacConstraint dJC(&mp,true);
    gsShapeOptLog slog1(output,true,false,false);

    // gsOptAntenna optA(&mp,numRefine,&slog1,0,quA,quB);

    // gsWinslowWithDeriv win(&mp,optA.mappers(),false,false,true,0);
    gsWinslowWithDeriv win(&mp,false,false,true,0);
    win.setQuad(quA,quB);

    real_t minD,minDgauss;
    index_t neededRefSteps;

    gsMatrix<> out(maxiter,4);
    out.setZero(maxiter,4);

    for (index_t i = 0; i < maxiter; i += 10)
    {
        std::stringstream stream;
        stream << BASE_FOLDER << output << "cps" << "_0_" << i << ".txt";
        std::string str = stream.str();
        gsMatrix<> mat = win.getFlat();
        readFromTxt(str,mat);
        // gsInfo << mat;

        win.updateFlat(mat);

        minD = dJC.provePositivityOfDetJ_TP(neededRefSteps, 3);

        minDgauss = win.minDetJInGaussPts();

        out(i,0) = minD;
        out(i,1) = neededRefSteps;
        out(i,2) = minDgauss;
        out(i,3) = win.minDetJInGaussPts(15);
    }

    std::stringstream stream;
    stream << BASE_FOLDER << output << "detJ.txt";
    std::string str = stream.str();
    saveMat(out,str);
    exit(0);



}

if (false ) {
    // Note on what to try:
    //  How about if I run without constraints on detJ and afterwards refine m_mp
    //  where the negative coefficients are..?
    gsShapeOptLog slog1(output,true,false,false);
    gsOptAntenna optA(&patches,numRefine,&slog1,param,quA,quB,true);
    gsWinslow winslow(&patches,optA.mappers(),false,false,true,0);

    winslow.setQuad(quA,quB);

    winslow.updateFlat(loadVec(winslow.n_flat,BASE_FOLDER + output));

    gsDetJacConstraint dJC(&patches);
    gsInfo << "mind = " << dJC.evalCon().minCoeff() << "\n";

    dJC.refineUntilPositive(maxiter);

    gsInfo << "mind = " << dJC.evalCon().minCoeff() << "\n";
    exit(0);
}

// Test of refine and plotting det J surfaces
if (false) {
    // Note on what to try:
    //  How about if I run without constraints on detJ and afterwards refine m_mp
    //  where the negative coefficients are..?
    gsShapeOptLog slog1(output,true,false,false);
    gsOptAntenna optA(&patches,numRefine,&slog1,param,quA,quB,true);
    gsWinslow winslow(&patches,optA.mappers(),false,false,true,0);

    winslow.setQuad(quA,quB);

    winslow.updateFlat(loadVec(winslow.n_flat,BASE_FOLDER + output));
    // winslow.update();

    gsDetJacConstraint dJC(&patches);

    gsMultiPatch<> tmp(patches);
    gsInfo << "tmp size " << tmp.patch(0).coefsSize() << "\n";

    // winslow.update();

    std::string name = "winslow";
    gsWriteParaview(patches,name,10000,true);

    gsMultiPatch<> djsurf = dJC.getDetJSurface();

    name = "/../results/testOfDetJ/test2_detJSurface";
    gsWriteParaview(djsurf,BASE_FOLDER + name,10000,false,true);

    gsMultiPatch<> zero = dJC.getDetJSurface(true);

    name = "/../results/testOfDetJ/test2_zero";
    gsWriteParaview(zero,BASE_FOLDER + name,10000,false,false);

    patches.uniformRefine();
    dJC.setup();
    gsInfo << "mind = " << dJC.evalCon().minCoeff() << "\n";

    djsurf = dJC.getDetJSurface();
    dJC.testSplineDetJ();

    name = "/../results/testOfDetJ/test2_detJSurface_uni1";
    gsWriteParaview(djsurf,BASE_FOLDER + name,10000,false,true);

    patches.uniformRefine();
    dJC.setup();
    gsInfo << "mind = " << dJC.evalCon().minCoeff() << "\n";

    djsurf = dJC.getDetJSurface();
    dJC.testSplineDetJ();
    gsInfo << "\n";

    name = "/../results/testOfDetJ/test2_detJSurface_uni2";
    gsWriteParaview(djsurf,BASE_FOLDER + name,10000,false,true);

    patches.uniformRefine();
    dJC.setup();
    gsInfo << "mind = " << dJC.evalCon().minCoeff() << "\n";

    djsurf = dJC.getDetJSurface();
    dJC.testSplineDetJ();
    gsInfo << "\n";

    name = "/../results/testOfDetJ/test2_detJSurface_uni4";
    gsWriteParaview(djsurf,BASE_FOLDER + name,10000,false,true);

    /// ======== WITH ADAPTIVITY ======= //
    // RESET PATCHES
    gsInfo << "\n .......... with adaptivity ............\n";
    patches = tmp;
    gsInfo << "patches size " << patches.patch(0).coefsSize() << "\n";
    dJC.setup();

    djsurf = dJC.getDetJSurface();

    real_t mind = dJC.refineDetJSurfaceUntilPositive(1,djsurf);
    // winslow.refineBasedOnDetJ(0);
    gsInfo << "mind = " << mind << "\n\n";

    name = "/../results/testOfDetJ/test2_detJSurface_adapt1";
    gsWriteParaview(djsurf,BASE_FOLDER + name,10000,true,true);

    // Adapt Twice
    dJC.setup();

    djsurf = dJC.getDetJSurface();

    mind = dJC.refineDetJSurfaceUntilPositive(2,djsurf);
    // winslow.refineBasedOnDetJ(0);
    gsInfo << "mind = " << mind << "\n\n";

    name = "/../results/testOfDetJ/test2_detJSurface_adapt2";
    gsWriteParaview(djsurf,BASE_FOLDER + name,10000,true,true);

    // Adapt Thrice
    dJC.setup();

    djsurf = dJC.getDetJSurface();

    mind = dJC.refineDetJSurfaceUntilPositive(3,djsurf);
    // winslow.refineBasedOnDetJ(0);
    gsInfo << "mind = " << mind << "\n\n";

    name = "/../results/testOfDetJ/test2_detJSurface_adapt3";
    gsWriteParaview(djsurf,BASE_FOLDER + name,10000,true,true);

    exit(0);
}

// Test of runtime
if (true) {
    gsInfo << "test mb\n";

    gsShapeOptLog slog1(output,true,false,false);
    // gsInfo << optWR.evalObj();
    // gsInfo << optWR.gradObj();
    // optWR.solve();

    if (param == 5) // Use regularization
    {
        gsOptAntenna optA(&mp,numRefine,&slog1,0,quA,quB);
        gsShapeOptWithReg optWR(&mp,&optA,numRefine,&slog1,quA,quB,eps);
        optWR.solve();
    } else {
        gsOptAntenna optA(&mp,numRefine,&slog1,param,quA,quB,true);
        // optA.m_paramMethod->updateFlat(loadVec(optA.n_flat,BASE_FOLDER + output + "cps_start.txt"));
        // optA.solve();
        optA.runOptimization(maxiter);
    }

    exit(0);


}

// Test of gsOptParam with reg
if (false) {
    // gsMultiPatch<> jigsaw = getJigSaw(startDes);
    gsMultiPatch<> jigsaw = getJigSaw3D(3); //

    for( index_t r = 0; r < numRefine; r++){
        jigsaw.uniformRefine();
    }

    changeSignOfDetJ(jigsaw.patch(0)); // Use for 3D small

    // Create startguess

    //-------- 3D JIGSAW ------
    index_t size    = jigsaw.patch(0).coefsSize();
    index_t n       = 5;

    gsVector<> vec;
    vec.setLinSpaced(n,-2.5,2.5);

    gsMatrix<> coefs_init(size,3);
    for (index_t i = 0; i < n; i++){
        for(index_t j = 0; j < n; j++){
            for(index_t k = 0; k < n; k++){
                coefs_init(k + j*n + i*n*n,0) = vec[j];
                coefs_init(k + j*n + i*n*n,1) = vec[k];
                coefs_init(k + j*n + i*n*n,2) = vec[i];
            }
        }
    }


    //-------- 2D JIGSAW ------
    //
    // index_t size    = jigsaw.patch(0).coefsSize();
    // index_t n       = sqrt(size);
    //
    // gsVector<> vec;
    // vec.setLinSpaced(n,0,-10);
    //
    // gsVector<> vec2;
    // vec2.setLinSpaced(n,0,10);
    //
    // gsMatrix<> coefs_init(size,2);
    // for (index_t i = 0; i < n; i++){
    //     for(index_t j = 0; j < n; j++){
    //         coefs_init(i + j*n,0) = vec[j];
    //         coefs_init(i + j*n,1) = vec2[i];
    //     }
    // }

    gsMultiPatch<> mp_init(jigsaw);
    mp_init.patch(0).setCoefs(coefs_init);
    changeSignOfDetJ(mp_init.patch(0)); // Use for 3D small

    // gsMultiPatch<>(*gsNurbsCreator<>::BSplineSquare(2));

    gsShapeOptLog slog1(output,true,false,false);

    std::string name = "init";
    slog1.plotInParaview(mp_init,name);

    gsSpringMethod spring(&mp_init);
    gsSpringMethod spring2(&jigsaw);

    // gsInfo << (spring.getFlat() - spring2.getFlat()).norm() << "\n";
    // gsInfo << (mp_init.patch(0).coefs() - jigsaw.patch(0).coefs()).norm() << "\n";
    // mat << spring.getTagged(), spring2.getTagged();

    // gsInfo << spring.n_tagged << "\n";
    gsOptParam optP(&mp_init,&jigsaw,&slog1,param);
    // gsMatrix<> v = optP.currentDesign()*0.99999;
    // optP.m_paramMethod->update();
    // gsInfo << optP.evalObj(v) << "\n";
    // gsInfo << optP.currentDesign() << "\n";
    // gsInfo << optP.numDesignVars() << "\n";
    // optP.solve();
    gsShapeOptWithReg optWR(&mp_init,&optP,numRefine,&slog1,quA,quB,eps);
    // optWR.gradObj();
    optWR.solve();
    exit(0);

}

// Test Regularization with winslow for antenna problem
if (false) {
    gsInfo << output << "\n";
    gsShapeOptLog slog1(output,true,false,false);
    gsOptAntenna optA(&mp,numRefine,&slog1,0,quA,quB);

    // gsSpringMethod spring(&mp,optA.mappers());
    // spring.update();
    // gsWinslowWithDeriv winslow(&mp,optA.mappers(),true,true,false,0);
    // winslow.update();
    //
    // std::stringstream stream;
    //
    // stream << BASE_FOLDER << output << "des";
    // std::string name = stream.str();
    //
    // gsWriteParaview(mp,name,100000,true);
    //
    // gsDetJacConstraint dJC(&mp,true);
    //
    // index_t nrs;
    // gsInfo << "PROVE DETJ>0 : " << dJC.provePositivityOfDetJ_TP(nrs,5) << "\n";
    //
    // gsInfo << "min detJ in pts: " << winslow.minDetJInGaussPts(10) << "\n";

    gsShapeOptWithReg optWR(&mp,&optA,numRefine,&slog1,quA,quB,eps);
    // gsInfo << optWR.evalObj();
    // gsInfo << optWR.gradObj();
    // optWR.solve();
    optWR.runOptimization(maxiter);

    // convergenceTestOfJacobian(optA);
    exit(0);

}

// calculate detJ constraints for different designs
if (false) {
    // gsMultiPatch<> jigsaw = getJigSaw(1);
    gsMultiPatch<> jigsaw = getJigSaw3D(2); //
    gsSpringMethod spring(&jigsaw);
    gsDetJacConstraint dJC(&jigsaw,true);

    gsMatrix<> minDetJ(maxiter,3);

    for (index_t i = 0; i < maxiter; i++){

        std::stringstream stream;
        stream << BASE_FOLDER << output << "cps" << "_0_" << i << ".txt";
        std::string str = stream.str();
        gsMatrix<> mat = spring.getFlat();
        readFromTxt(str,mat);
        // gsInfo << mat;

        spring.updateFlat(mat);

        minDetJ(i,0) = dJC.evalCon().minCoeff();

        gsMultiPatch<> ref(jigsaw);
        ref.uniformRefine();
        gsDetJacConstraint dJCr(&ref,true);

        minDetJ(i,1) = dJCr.evalCon().minCoeff();
        //
        ref.uniformRefine();
        gsDetJacConstraint dJCr2(&ref,true);

        minDetJ(i,2) = dJCr2.evalCon().minCoeff();
        //
        // ref.uniformRefine();
        // gsDetJacConstraint dJCr3(&ref,true);
        //
        // minDetJ(i,3) = dJCr3.evalCon().minCoeff();

    }

    std::stringstream stream;
    stream << BASE_FOLDER << output << "detJ.txt";
    std::string str = stream.str();
    saveMat(minDetJ,str);
    exit(0);



}

// Test of gsOptParam
if (false) {
    // gsMultiPatch<> jigsaw = getJigSaw(1);
    gsMultiPatch<> jigsaw = getJigSaw3D(1); //

    for( index_t r = 0; r < numRefine; r++){
        jigsaw.uniformRefine();
    }

    changeSignOfDetJ(jigsaw.patch(0));

    // Create startguess
    index_t size    = jigsaw.patch(0).coefsSize();
    index_t n       = 10;

    //-------- 3D JIGSAW ------
    gsVector<> vec;
    vec.setLinSpaced(n,-2.5,2.5);

    gsMatrix<> coefs_init(size,3);
    for (index_t i = 0; i < n; i++){
        for(index_t j = 0; j < n; j++){
            for(index_t k = 0; k < n; k++){
                coefs_init(k + j*n + i*n*n,0) = vec[j];
                coefs_init(k + j*n + i*n*n,1) = vec[k];
                coefs_init(k + j*n + i*n*n,2) = vec[i];
            }
        }
    }


    //-------- 2D JIGSAW ------
    //
    // index_t size    = jigsaw.patch(0).coefsSize();
    // index_t n       = 10;
    //
    // gsVector<> vec;
    // vec.setLinSpaced(n,0,-10);
    //
    // gsVector<> vec2;
    // vec2.setLinSpaced(n,0,10);
    //
    // gsMatrix<> coefs_init(size,2);
    // for (index_t i = 0; i < n; i++){
    //     for(index_t j = 0; j < n; j++){
    //         coefs_init(i + j*n,0) = vec[j];
    //         coefs_init(i + j*n,1) = vec2[i];
    //     }
    // }

    gsMultiPatch<> mp_init(jigsaw);
    mp_init.patch(0).setCoefs(coefs_init);
    changeSignOfDetJ(mp_init.patch(0));

    // gsMultiPatch<>(*gsNurbsCreator<>::BSplineSquare(2));

    gsShapeOptLog slog1(output,true,false,false);

    std::string name = "init";
    slog1.plotInParaview(mp_init,name);

    gsSpringMethod spring(&mp_init);
    gsSpringMethod spring2(&jigsaw);

    // gsInfo << (spring.getFlat() - spring2.getFlat()).norm() << "\n";
    // gsInfo << (mp_init.patch(0).coefs() - jigsaw.patch(0).coefs()).norm() << "\n";
    // mat << spring.getTagged(), spring2.getTagged();

    // gsInfo << spring.n_tagged << "\n";
    gsOptParam optP(&mp_init,&jigsaw,&slog1,param);
    // gsMatrix<> v = optP.currentDesign()*0.99999;
    // optP.m_paramMethod->update();
    // gsInfo << optP.evalObj(v) << "\n";
    // gsInfo << optP.currentDesign() << "\n";
    // gsInfo << optP.numDesignVars() << "\n";
    optP.solve();

    exit(0);




}

// Test hessian and jacobUpdate of winslow
if (false) {
    gsShapeOptLog slog1(output,true,false,false);
    gsOptAntenna optA(&patches,numRefine,&slog1,param,quA,quB);

    gsSpringMethod sM(&patches,optA.mappers());
    sM.computeMap();
    sM.update();

    gsWinslow win(&patches,optA.mappers(),false,false,true,0);
    // win.update();

    gsOptAntenna optA2(&patches,numRefine,&slog1,param,quA,quB);

    // convergenceTestOfParaJacobianAll(win);
    // convergenceTestOfParaGradAll(win);
    // convergenceTestOfParaGrad(win);
    // convergenceTestOfParaHessianAll(win);
    // convergenceTestOfParaHessian(win);
    convergenceTestOfParamMethodJacobian(*optA2.m_paramMethod);

    // Save Pc and Px
    exit(0);
}

// Test shape optimization with new winslow
if (false) {
    gsShapeOptLog slog1(output,true,false,false);
    gsOptAntenna optA(&patches,numRefine,&slog1,param,quA,quB);

    gsSpringMethod sM(&patches,optA.mappers());
    sM.computeMap();
    sM.update();

    gsWinslow win(&patches,optA.mappers(),false,false,true,0);
    win.update();

    // optA.m_paramMethod->update();

    // gsAsVector<> result(optA.m_dJC->evalCon().data(),optA.numConstraints());

    // optA.evalCon_into(ctmp,result);

    // optA.updateDesignVariables(sM.getTagged());
    // convergenceTestOfJacobian(optA);
    // convergenceTestOfParamMethodJacobian(*optA.m_paramMethod);

    gsOptAntenna optA2(&patches,numRefine,&slog1,param,quA,quB);
    optA2.solve();

    // optA2.runOptimization(maxiter);


    exit(0);
}

// Test Winslow with and without checking for inf
if (false){
    gsMultiPatch<> js = getGeometry(nx,ny,degree);
    // js.uniformRefine(1);
    //
    // gsMaxDetJac mDJ(&js,true);
    // mDJ.update();
    // mDJ.updateFlat(loadVec(mDJ.n_flat,BASE_FOLDER + output + "Jigsaw1/MaxDetJac.txt"));    //
    gsWinslow win(&js,false,false,true,0);
    win.update();
    // testOfParametrizations_winslow_inf(output,quA,quB);
    exit(0);
}

// Test aggregated constraints with unit square
if (false) {
    gsMultiPatch<> mp = getJigSaw(2);


	std::string out = "winslow";
	gsInfo << "Writing the gsMultiPatch to a paraview file: " << out << "\n\n";
	gsWriteParaview(mp, out);

    gsDetJacConstraint dJC(&mp);
    gsInfo << dJC.evalCon().minCoeff() << "\n";

    gsAggFun fun(dJC.n_constraints,alpha);
    gsAggregatedConstraint aC(&mp, &dJC, &fun);

    gsWinslow winslow(&mp,&aC);
    winslow.update();

	out = "modLiao_aC";
	gsInfo << "Writing the gsMultiPatch to a paraview file: " << out << "\n\n";
	gsWriteParaview(mp, out);
    exit(0);
}

// Test aggregated constraints
if (false) {
    testOfAggregatedConstraints2nd(output,quA,quB);
    //testOfAggregatedConstraints3rd(output,quA,quB,param,numRefine,maxiter,nx,ny,degree);
    exit(0);
}

// Plot magnitude
if (false) {
    gsOptAntenna optA(&patches,numRefine,&slog,param,quA,quB);
    gsWinslow winslow(&patches, optA.mappers(), useDJC);

    gsStateEquationAntenna SE(&patches, 1);

    winslow.updateFlat(loadVec(winslow.n_flat,BASE_FOLDER + output + "Spring/cps_0_178.txt"));

    gsInfo << "OBJ = " << optA.evalObj() << "\n";
    SE.plotMagnitude(BASE_FOLDER  + output + "Spring/mag_0_178");

    winslow.updateFlat(loadVec(winslow.n_flat,BASE_FOLDER + output + "Winslow_Lag/cps_4_117.txt"));
    SE.plotMagnitude(BASE_FOLDER  + output + "Winslow_Lag/mag_4_117");

    winslow.updateFlat(loadVec(winslow.n_flat,BASE_FOLDER + output + "Harmonic_Lag/cps_4_87.txt"));
    SE.plotMagnitude(BASE_FOLDER  + output + "Harmonic_Lag/mag_4_87");

}

// Test constraint aggregation
if (false) {
    gsDetJacConstraint dJC(&patches);
    gsOptAntenna optA(&patches,numRefine,&slog,param,quA,quB);
    gsWinslow winslow(&patches, optA.mappers(), useDJC);

    gsVector<> flat = loadVec(winslow.n_flat,BASE_FOLDER "/../results/TestOfParamAdaptive/spring1161.txt");
    winslow.updateFlat(flat);

    gsMaxDetJac mDJ(&patches, optA.mappers());
    mDJ.update();

    gsAggFun fun(dJC.n_constraints,alpha);

    gsMatrix<> res = fun.eval(dJC.evalCon());

    gsInfo << "\nminimal value of fun: " << dJC.evalCon().minCoeff() << "\n";
    gsInfo << "aggFun: " << res(0,0) << "\n\n";

    gsAggregatedConstraint aC(&patches, &dJC, &fun);

    convergenceTestOfConstraint(winslow,aC);

    exit(0);
}

// Test new det Jac constraints
if (false) {
    gsDetJacConstraint dJ(&mp);

    gsNewDetJacConstraint nDJ(&mp);
    exit(0);
}

// Test Liao derivatives
if (false) {
    gsOptAntenna optA(&patches,numRefine,&slog,param,quA,quB);
    gsHarmonic harmonic(&patches, optA.mappers(), useDJC);
    // gsInfo << "\nLiao \t:\t" << liao.evalObj() << "\n\n";
    harmonic.setQuad(quA,quB);
    convergenceTestOfParaJacobian(harmonic);
    // liao.update();
    exit(0);

}

// Test gsAffineParamMethod with Lagrangian
if (false) {
    gsOptAntenna optA(&patches,numRefine,&slog,param,quA,quB);
    gsWinslow winslow(&patches, optA.mappers(), useDJC);
    winslow.update();

}

// Test gsAffineParamMethod with Lagrangian on full optimization
if (false) {
    // gsOptAntenna optA1(&patches,numRefine,&slog,param,quA,quB,useLag);
    // gsWinslow winslow(&patches, optA1.mappers(), useDJC);

    // gsVector<> flat = loadVec(winslow.n_flat,BASE_FOLDER "/../results/TestOfLagrangian2nd/Spring/cps_0_178.txt");
    // gsInfo << flat << "\n";
    // gsInfo << winslow.n_flat << "\n";
    // winslow.updateFlat(flat);

    gsOptAntenna optA(&patches,numRefine,&slog,param,quA,quB,useLag);

    if (param == 0)
        optA.solve();
    else
        optA.runOptimization(maxiter);

}

// Test lagrange multipliers
if (false) {
    gsOptAntenna optA(&patches,numRefine,&slog,param,quA,quB);
    gsWinslow winslow(&patches, optA.mappers(), useDJC);
    // gsVector<> flat = loadVec(winslow.n_flat,BASE_FOLDER "/../results/TestOfParamAdaptive/spring1161.txt");
    // gsInfo << flat << "\n";
    // gsInfo << winslow.n_flat << "\n";
    // flat.segment(0,winslow.n_flat/2) *= -1;
    //
    winslow.setQuad(quA,quB);
    // winslow.updateFlat(flat);
    // // I dont know how to change gsOptProblem.h file and get it to work here..
    // winslow.update();
    //
    // gsInfo << "lagrangian: " << winslow.evalLagrangian() << "\n";
    // gsInfo << "gradLagrangian: " << winslow.gradLagrangian() << "\n";

    convergenceTestOfParaJacobian(winslow);
    convergenceTestOfParaLagrangianJacobian(winslow);
}

// Test winslow with 3D jigsaw
if (false) {
    // gsSpringMethod spring(&patches);
    // spring.update();

    gsWinslow winslow(&patches,useDJC);
    winslow.print();
    winslow.update();


    std::string name = "jigsaw3d";
    slog.plotInParaview(patches,name);
    exit(0);
}

// Test of active constraint refinement for jigsaw1
if (false){
    // gsOptAntenna optA(&patches,numRefine,&slog,param,quA,quB);
    // gsMaxDetJac mDJ(&patches);
    // mDJ.update();

    gsWinslow winslow(&patches, useDJC);

    // Max 10 iterations
    for (index_t iter = 0; iter < 10; iter++){
        gsInfo << "Iter " << iter << "\n";
        winslow.update();

        std::string iter_str = std::to_string(iter);

        std::string name = "jigsaw_iter" + iter_str;
        slog.plotInParaview(patches,name);

        gsDetJacConstraint dJC(&patches);

        std::vector<bool> elMarked;
        real_t tol = lambda_1;
        dJC.markElements(elMarked,tol);

        for (index_t i = 0; i < elMarked.size(); i++){
            slog << elMarked[i] << " ";
        }
        slog << "\n";

        // gsInfo << "\n" << dJC.evalCon() << "\n";

        dJC.plotDetJ(BASE_FOLDER + output + "/detJ_iter" + iter_str);
        dJC.plotActiveConstraints(elMarked,BASE_FOLDER + output + "/active_iter" + iter_str,tol);

        name = "d_iter" + iter_str;
        slog.saveVec(dJC.evalCon(),name);

        // Stop loop
        if (accumulate(elMarked.begin(),elMarked.end(),0) == 0) {
            break;
        }

        winslow.refineElements(elMarked); // also resets stuff
        winslow.recreateMappers();      // Recreate m_mappers, see gsParamMethod.h for further information
        winslow.m_dJC->setup();          // Reset the gsDetJacConstraint
        dJC.setup();
        winslow.setupOptParameters();   // Reset optimization parameters
        winslow.print();

        dJC.plotDetJ(BASE_FOLDER + output + "/detJ_justRefined_iter" + iter_str);
        dJC.plotActiveConstraints(elMarked,BASE_FOLDER + output + "/active_justRefined_iter" + iter_str,tol);

    }
    exit(0);
    exit(0);
}

// Test of active constraint refinement for 1161 design
if (false){
    gsOptAntenna optA(&patches,numRefine,&slog,param,quA,quB);
    gsWinslow winslow(&patches, optA.mappers(), useDJC);

    gsVector<> flat = loadVec(winslow.n_flat,BASE_FOLDER "/../results/TestOfParamAdaptive/spring1161.txt");
    // gsInfo << flat << "\n";
    // gsInfo << winslow.n_flat << "\n";
    flat.segment(0,winslow.n_flat/2) *= -1;

    winslow.setQuad(quA,quB);
    winslow.updateFlat(flat);

    // Max 10 iterations
    for (index_t iter = 0; iter < 10; iter++){
        gsInfo << "Iter " << iter << "\n";
        winslow.update();

        std::string iter_str = std::to_string(iter);

        std::string name = "1161_iter" + iter_str;
        slog.plotInParaview(patches,name);

        gsDetJacConstraint dJC(&patches);

        std::vector<bool> elMarked;
        real_t tol = 0.001;
        dJC.markElements(elMarked,tol);

        for (index_t i = 0; i < elMarked.size(); i++){
            slog << elMarked[i] << " ";
        }
        slog << "\n";

        // gsInfo << "\n" << dJC.evalCon() << "\n";

        dJC.plotDetJ(BASE_FOLDER + output + "/detJ_iter" + iter_str);
        dJC.plotActiveConstraints(elMarked,BASE_FOLDER + output + "/active_iter" + iter_str,tol);

        name = "d_iter" + iter_str;
        slog.saveVec(dJC.evalCon(),name);

        // Stop loop
        if (accumulate(elMarked.begin(),elMarked.end(),0) == 0) {
            break;
        }

        winslow.refineElements(elMarked); // also resets stuff
        winslow.recreateMappers();      // Recreate m_mappers, see gsParamMethod.h for further information
        winslow.m_dJC->setup();          // Reset the gsDetJacConstraint
        dJC.setup();          // Reset the gsDetJacConstraint
        winslow.setupOptParameters();   // Reset optimization parameters
        winslow.print();

        dJC.plotDetJ(BASE_FOLDER + output + "/detJ_justRefined_iter" + iter_str);
        dJC.plotActiveConstraints(elMarked,BASE_FOLDER + output + "/active_justRefined_iter" + iter_str,tol);

    }
    exit(0);
}

// Test d vector
if (false){
    gsWinslow winslow(&patches, useDJC);
    // winslow.update();

    gsDetJacConstraint dJC(&patches,true);

    const clock_t begin_time = clock();
    for (index_t i = 0; i < 20; i++)
        dJC.evalCon();
    gsInfo << float( clock () - begin_time ) /  CLOCKS_PER_SEC;

    std::string name = "/../results/d";
    slog.saveVec(dJC.evalCon(),name);
    exit(0);
}

// Test a lot of convergence stuff
if (false){
    gsWinslow winslow(&patches, useDJC);

    gsInfo << "\n\nConvergence test of det jac: \n";
    convergenceTestOfDetJJacobian(winslow);

    gsInfo << "\n\nConvergence test of detJac: \n";
    convergenceTestOfParaJacobian(winslow);

    winslow.update();

    std::string name = "/../results/test/winslow";
    slog.plotInParaview(patches,name);

    exit(0);
}

return 0;
}

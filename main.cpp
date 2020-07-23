//#undef __cplusplus
//#define __cplusplus 201402L

#include <gismo.h>
#include <stdio.h>
#include <math.h>       /* pow */
#include <iostream>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <ctime>

#include "gsMyExpressions.h"

#include "gs2NormConstraints.h"
#include "gsDetJacConstraint.h"
#include "gsNewDetJacConstraint.h"
#include "gsAggregatedConstraint.h"
#include "gsAggFun.h"

#include "gsParamMethod.h"
#include "gsSpringMethod.h"
#include "gsSpringMethod2nd.h"
#include "gsModLiao.h"
#include "gsWinslow.h"
#include "gsWinslowPow.h"
// #include "gsKnupp.h"
#include "gsLiao.h"
#include "gsHarmonic.h"
#include "gsMaxDetJac.h"

#include "gsAffineOptParamMethod.h"
#include "gsIpOptSparseMatrix.h"
#include "gsShapeOptProblem.h"
#include "gsShapeOptWithReg.h"
#include "gsOptPotWavesWithReg.h"
#include "gsOptPotWavesWithReg2nd.h"
#include "gsOptAntenna.h"
#include "gsOptPotWaves.h"
#include "gsOptParam.h"
#include "gsOptParamFull.h"
#include "gsOptInit.h"
#include "gsOptInit2nd.h"
#include "gsOptInit3rd.h"
#include "gsOptInit4th.h"
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
		for(index_t j = 0; j < des.size(); j++){
			perturp[j] = ran[j];
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

void convergenceTestOfConstraintFree(gsOptParamMethod &pM, gsConstraint &con){
	// gsVector<> result = dJC.generateDResultVector();
	// dJC.getDvectors(result);
	gsVector<> result = con.evalCon();

    gsIpOptSparseMatrix J = con.getJacobian();
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
		for(index_t j = 0; j < des.size(); j++){
			perturp[j] = ran[j];
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
		for(index_t j = 0; j < des.size(); j++){
			perturp[j] = ran[j];
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

	std::srand((unsigned int) std::time(0));
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

		for(index_t j = 0; j < lOP.n_free; j++){
			perturp[j] = ran[j];
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

	lOP.updateFree(des);

}

void convergenceTestOfFunction(gsFunctionExpr<> fun, gsFunctionExpr<> grad, gsFunctionExpr<> hess ){
    index_t dim = 3;
    gsMatrix<> u(dim,1);
    u << 0,0,0;

	gsMatrix<> fun_eval_mat = fun.eval(u);
    real_t fun_eval = -fun_eval_mat(0,0);
	gsVector<> fun_grad = grad.eval(u);
	gsMatrix<> fun_hess = hess.eval(u);

    gsDebugVar(fun_eval);
    gsDebugVar(fun_grad);
    gsDebugVar(fun_hess);

	std::srand((unsigned int) std::time(0));
	gsVector<> ran;
	ran.setRandom(dim);

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
		perturp.setZero(dim);

		for(index_t i = 0; i < dim; i++){
			perturp[i] = ran[i];
		}

		perturp /= perturp.norm();

		perturp *= eps;

		gsVector<> u_new = u + perturp;

		gsMatrix<> fun_eval_new = fun.eval(u_new);
        real_t fun_new = -fun_eval_new(0,0);

        gsAsConstMatrix<> fun_hessC(fun_hess.data(),dim,dim);

		real_t guess0 = fun_eval;
		real_t guess1 = fun_eval + fun_grad.transpose()*perturp;
		real_t guess2 = fun_eval + fun_grad.transpose()*perturp + 0.5*perturp.transpose()*fun_hessC*perturp;

		gsInfo << fun_new <<" " << guess1 << " " << guess2 << "\n";
		// gsInfo << guess0 <<" " << guess1 << " " << newLiao << "\n";

		real_t error0 = std::abs(guess0 - fun_new);
		real_t error1 = std::abs(guess1 - fun_new);
		real_t error2 = std::abs(guess2 - fun_new);

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

		for(index_t j = 0; j < lOP.n_free; j++){
			perturp[j] = ran[j];
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

    gsDebugVar(sOP.evalObj());
	sOP.updateDesignVariables(des);
	real_t obj = sOP.evalObj();
    gsDebugVar(obj);

    // Time
    clock_t start, stop;
	double totalTime;

	start = clock();
	gsVector<> grad = sOP.gradObj();
	stop = clock();
	totalTime = (stop - start) / (double)CLOCKS_PER_SEC;
    gsInfo << "\n Time to compute grad: " << totalTime << "\n";

	

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
		for(index_t j = 0; j < sOP.numDesignVars(); j++){
			perturp[j] = ran[j];
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

void convergenceTestOfJacobianAll(gsShapeOptProblem &sOP){
	//std::srand((unsigned int) std::time(0));
	gsVector<> ran;
	ran.setRandom(sOP.n_flat);

	gsVector<> des = sOP.getFlat();
	gsInfo << "\n Size of design vector : " << des.size() << "\n";

	sOP.updateFlat(des);
	real_t obj = sOP.evalObj();
	gsVector<> grad = sOP.gradAll();

    //gsDebugVar(grad.norm());

	index_t beg = 0;
	index_t n = 8;
	gsVector<> Eps(n);
	gsVector<> Error0(n);
	gsVector<> Error1(n);

	for(index_t i = 0; i < n; i++){
		// Generate pertubation
		real_t eps = pow(2,-beg-i);
		Eps[i] = eps;

		gsVector<> perturp;
		perturp.setZero(sOP.n_flat);

        index_t i_flat = 0;
        for( index_t d = 0; d < sOP.m_mp->targetDim(); d++)
        {
            for( index_t p = 0; p < sOP.m_mp->nBoxes(); p++)
            {
                for (index_t i = 0; i < sOP.m_mp->patch(p).coefsSize(); i++)
                {
                    if(sOP.isCpsInDomain(i,p,d))
                    {
                        perturp[i_flat] = ran[i_flat];
                    }

                    i_flat++;
                }
            }
        }

		perturp /= perturp.norm();

		perturp *= eps;

		gsVector<> newDes = des + perturp;

		// sOP.resetParametrizationToReference();
		sOP.updateFlat(newDes);
		// if (i == 1) sOP.writeToFile(sOP.getFlat(),"cpsPerturb.txt");
        //

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
	gsInfo << "\neps \tErr0 \tErr1 \trate\n";
	gsInfo << disp << "\n";

}

void convergenceTestOfJacobian(gsOptInit4th &sOP){
	//std::srand((unsigned int) std::time(0));
	gsVector<> ran;
	ran.setRandom(sOP.numDesignVars());

	gsMatrix<> des = sOP.getDesignVars();
	gsInfo << "\n Size of design vector : " << des.size() << "\n";

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

		gsMatrix<> newDes = des + perturp;

		// sOP.resetParametrizationToReference();
		// if (i == 1) sOP.writeToFile(sOP.getFlat(),"cpsPerturb.txt");
        //

        sOP.updateFromDesignVars(newDes);

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
	gsInfo << "\neps \tErr0 \tErr1 \trate\n";
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
		for(index_t j = 0; j < tag.size(); j++){
			perturp[j] = ran[j];
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
            for (index_t j = 0; j < FD.size(); i++){
                gsInfo << FD[j] << " \t" << deriv[j] << "\n";
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

		for(index_t j = 0; j < lOP.n_free; j++){
			p[j] = ran_free[j];
		}
		for(index_t j = 0; j < lOP.n_tagged; j++){
			q[j] = ran_tagged[j];
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

		for(index_t j = 0; j < lOP.n_free; j++){
			p[j] = ran_free[j];
		}
		for(index_t j = 0; j < lOP.n_tagged; j++){
			q[j] = ran_tagged[j];
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

	std::srand((unsigned int) std::time(0));
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

		for(index_t j = 0; j < lOP.n_flat; j++){
			perturp[j] = ran[j];
		}

		perturp /= perturp.norm();

		perturp *= eps;

		gsVector<> newDes = des + perturp;


		lOP.updateFlat(newDes);

		real_t newLiao = lOP.evalObj();
		real_t guess0 = liao;
		real_t guess1 = liao + grad.transpose()*perturp;
		real_t guess2 = liao + grad.transpose()*perturp + 0.5*perturp.transpose()*hess*perturp;

		// gsInfo << newLiao <<" " << guess1 << " " << guess2 << "\n";
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

	lOP.updateFlat(des);

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

		for(index_t j = 0; j < lOP.n_flat; j++){
			p[j] = ran1[j];
			q[j] = ran2[j];
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


		for(index_t j = 0; j < lOP.n_flat; j++){
			pA[j] = ran_flat[j];
		}
		for(index_t j = 0; j < lOP.n_free; j++){
			p[j] = ran_free[j];
		}
		for(index_t j = 0; j < lOP.n_tagged; j++){
			q[j] = ran_tagged[j];
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
        //gsInfo << "load from " << name << "\n";
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
	gsMatrix<> coefs (basis.size(), 2);

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
	gsMatrix<> coefsMid (basisMid.size(), 2);

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
	gsMatrix<> coefs (basis.size(), 2);

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
    if (i == 1 || i == 5 || i == 6 || i == 7)
        n = 10;
    if (i == 2)
        n = 6;
    if (i == 3)
        n = 5;
	if (i == 4 || i == 7)
	{
        std::ostringstream strs;
        if (i == 4)
		    strs << BASE_FOLDER << "/parametrizations/water_passage_3p.xml";

        if (i == 7)
		    strs << BASE_FOLDER << "/parametrizations/water_passage_fixed.xml";
		
		std::string fn(strs.str());

		gsFileData<> fd(fn);

		gsMultiPatch<>::uPtr mp_ptr;
	    mp_ptr = fd.getFirst< gsMultiPatch<> > ();

		return *mp_ptr;
	}
    index_t degree = 2;

    std::string folder = "/parametrizations/JigSaw/";

	// 1. construction of a knot vector for each direction
	gsKnotVector<> kv1(0, 1, n - degree - 1, degree + 1);
	gsKnotVector<> kv2(0, 1, n - degree - 1, degree + 1);
	gsKnotVector<> kv3(0, 1, n - degree - 1, degree + 1);
	// 2. construction of a basis
	gsTensorBSplineBasis<3, real_t> basis(kv1, kv2, kv3);
	// 3. construction of a coefficients
	gsMatrix<> coefs (n*n*n, 3);

    if (i == 1)
	   readFromTxt(BASE_FOLDER + folder + "jigsaw3d_1.txt", coefs);
    if (i == 2)
	   readFromTxt(BASE_FOLDER + folder + "jigsaw3d_small.txt", coefs);
    if (i == 3)
	   readFromTxt(BASE_FOLDER + folder + "example.txt", coefs);
    if (i == 5)
	   readFromTxt(BASE_FOLDER + folder + "jigsaw3d_2.txt", coefs);
    if (i == 6)
	   readFromTxt(BASE_FOLDER + folder + "jigsaw3d_3.txt", coefs);
    if (i == 7)
	   readFromTxt(BASE_FOLDER + folder + "jigsaw3d_4.txt", coefs);
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
	gsMatrix<> coefs (basis.size(), 2);

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

gsMultiPatch<> getTorus(){

    gsMultiPatch<> out(getGeometry(5,4,2).patch(4));

    gsMatrix<> coefs = out.patch(0).coefs();

    patchSide ps = out.boundaries()[0];
    gsVector< unsigned > boundaryDofs_lhs = out.basis(0).boundary(ps);

    ps = out.boundaries()[1];
    gsVector< unsigned > boundaryDofs_rhs = out.basis(0).boundary(ps);

    ps = out.boundaries()[2];
    gsVector< unsigned > boundaryDofs_top = out.basis(0).boundary(ps);

    ps = out.boundaries()[3];
    gsVector< unsigned > boundaryDofs_bot = out.basis(0).boundary(ps);

    for (index_t i = 0; i < boundaryDofs_top.size(); i++)
    {
        coefs.row(boundaryDofs_bot[i]) = coefs.row(boundaryDofs_top[i])*0.33;
    }

    gsVector<> left = coefs.row(boundaryDofs_rhs[0]);
    gsVector<> right = coefs.row(boundaryDofs_bot[4]);

    for (index_t i = 0; i < boundaryDofs_rhs.size(); i++)
    {
        gsDebugVar(coefs.row(boundaryDofs_rhs[i]));
        coefs.row(boundaryDofs_rhs[i]) = (1 - i/3.0)* left + i/3.0 * right;
        gsDebugVar(coefs.row(boundaryDofs_rhs[i]));
    }

    left = coefs.row(boundaryDofs_lhs[0]);
    right = coefs.row(boundaryDofs_bot[0]);

    for (index_t i = 0; i < boundaryDofs_lhs.size(); i++)
    {
        coefs.row(boundaryDofs_lhs[i]) = (1 - i/3.0)* left + i/3.0 * right;

    }


    out.patch(0).setCoefs(coefs);

	return out;

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

void extractCornersAndSaveXML(index_t jigtype, index_t dim)
{

	gsMultiPatch<> jigsaw;
	if (dim == 2) {
		jigsaw = getJigSaw(jigtype);
	} else if (dim == 3) {
		jigsaw = getJigSaw3D(jigtype); 
	} else {
		GISMO_ERROR("dim should be 2 or 3!");
	}
	
	gsMultiPatch<> corner = jigsaw.uniformSplit();

	gsWriteParaview(corner,"corner",10000,true,true);

	gsFileData<> fd;
	fd << corner; 

	index_t num = 1;
	if (dim == 3)
    {
		if (jigtype == 5)
			num = 3;
		if (jigtype == 6)
			num = 2;
		if (jigtype == 7)
			num = 4;
	} 
    else if (dim == 2)
	{
		if (jigtype == 1)
			num = 1;
		if (jigtype == 2)
			num = 2;
	}
	
    std::ostringstream strs;
    if (jigtype == 2 && dim == 3)
    {
        strs << BASE_FOLDER << "/parametrizations/JigsawXML/jig_small" << dim << "d_small";
    }
    else
    {
        strs << BASE_FOLDER << "/parametrizations/JigsawXML/jig" << dim << "d_" << num;
    }
    std::string str = strs.str();

	fd.save(str);
}

gsMultiPatch<> loadCorner(index_t jigtype, index_t dim, std::vector< gsDofMapper > &mappers, index_t numRefine )
{
    std::ostringstream strs;
	strs << BASE_FOLDER << "/parametrizations/JigsawXML/jig" << dim << "d_" << jigtype << ".xml";
	
	std::string fn(strs.str());

	gsFileData<> fd(fn);

	gsMultiPatch<>::uPtr mp_ptr;
	mp_ptr = fd.getFirst< gsMultiPatch<> > ();

	gsMultiPatch<> mp(*mp_ptr);

	for (index_t r = 0; r < numRefine; r++)
		mp.uniformRefine();

	index_t p = 0;
	gsMultiPatch<> out(mp.patch(p));

    if (jigtype == 5)
        changeSignOfDetJ(out.patch(0));

    // Get mappers from multibasis with interfaces glued
    gsMultiBasis<> geoBasis(out);
    for (index_t d = 0; d < out.targetDim(); d++)
    {
        geoBasis.getMapper(iFace::glue,mappers[d],false); // False means that we do not finalize
    }

    // Fix coordinates of boundaries
    for (index_t i = 0; i < mp.nBoundary(); i++){
        // Get boundary local indices
        patchSide ps = mp.boundaries()[i];

		if (ps.patch != p) continue; // Skip other patches

        gsVector<unsigned> boundaryDofs = mp.basis(ps.patch).boundary(ps);

        for (index_t d = 0; d < out.targetDim(); d++)
        {
            mappers[d].markBoundary(ps.patch,boundaryDofs);
        }
    }

    for (index_t i = 0; i < mp.nInterfaces(); i++){
        // Get boundary local indices
		boundaryInterface interface = mp.bInterface(i);
		patchSide ps = interface.first();

		if (ps.patch != p)
		{
			ps = interface.second();
			if (ps.patch != p) 
				continue;
		}

        gsVector<unsigned> boundaryDofs = mp.basis(ps.patch).boundary(ps);

        index_t size = boundaryDofs.size();

		index_t fixDir = 0;
        for (index_t d = 0; d < out.targetDim(); d++)
		{
			real_t diff1 = out.patch(0).coef(boundaryDofs[0],d) - out.patch(0).coef(boundaryDofs[1],d);
			real_t diff2 = out.patch(0).coef(boundaryDofs[0],d) - out.patch(0).coef(boundaryDofs[size/2],d);
			if (diff1 == 0 && diff2 == 0)
			{
				gsInfo << "Fix direction " << d << "\n";
				fixDir = d;
				break;
			}

		}

        mappers[fixDir].markBoundary(ps.patch,boundaryDofs);
    }

    // Finalize mappers
    for (index_t d = 0; d < mp.targetDim(); d++)
    {
        mappers[d].finalize();
    }

    // Tag coordinates of boundaries
    for (index_t i = 0; i < mp.nBoundary(); i++){
        // Get boundary local indices
        patchSide ps = mp.boundaries()[i];
        gsVector<unsigned> boundaryDofs = mp.basis(ps.patch).boundary(ps);

		if (ps.patch != p) continue;

        for (index_t d = 0; d < mp.targetDim(); d++)
        {
            // Tag the controlpoint
            for (index_t j = 0; j < boundaryDofs.size(); j ++){
                mappers[d].markTagged(boundaryDofs[j],ps.patch);
            }
        }
    }

    // Tag interfaces
    for (index_t i = 0; i < mp.nInterfaces(); i++){
        // Get boundary local indices
		boundaryInterface interface = mp.bInterface(i);
		patchSide ps = interface.first();

		if (ps.patch != p)
		{
			ps = interface.second();
			if (ps.patch != p) 
				continue;
		}

        gsVector<unsigned> boundaryDofs = mp.basis(ps.patch).boundary(ps);

        index_t size = boundaryDofs.size();

		index_t fixDir = 0;
        for (index_t d = 0; d < out.targetDim(); d++)
		{
			real_t diff1 = out.patch(0).coef(boundaryDofs[0],d) - out.patch(0).coef(boundaryDofs[1],d);
			real_t diff2 = out.patch(0).coef(boundaryDofs[0],d) - out.patch(0).coef(boundaryDofs[size/2],d);
			if (diff1 == 0 && diff2 == 0)
			{
				gsInfo << "Fix direction " << d << "\n";
				fixDir = d;
				break;
			}

		}

            for (index_t j = 0; j < boundaryDofs.size(); j ++){
                mappers[fixDir].markTagged(boundaryDofs[j],ps.patch);
            }
    }

	return out;

	

}

/*
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
    P[1].patch(0).krotate(M_PI_4);
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


/*
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
*/

gsMultiPatch<> get3DGeometry(){
    gsMultiPatch<> patches;
	gsMultiPatch<>::Ptr patches_ptr = memory::make_shared_not_owned(&patches);
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
    gsDetJacConstraint dJC(patches_ptr);

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
	gsMultiPatch<>::Ptr mp_ptr(&mp);
    // Setup classes
    gsDetJacConstraint dJC(mp_ptr);
    gsWinslow winslow(mp_ptr,false);         // We use the default mappers

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

bool testValidity(gsMultiPatch<> &mp, std::string name, real_t quA, index_t quB, std::string output){
	gsMultiPatch<>::Ptr mp_ptr = memory::make_shared_not_owned(&mp);

    gsDetJacConstraint dJC(mp_ptr,true);

    gsWinslowWithDeriv win(mp_ptr,false,false,true,0);
    win.setQuad(quA,quB);

    real_t minD,minDgauss;
    index_t neededRefSteps;

    gsVector<> out(4);
    out.setZero(4);

    real_t minGaussD = win.minDetJInGaussPts(15);
    out(2) = minGaussD;

    if (minGaussD < 0)
    {
        gsInfo << "detJ < 0 in gauss point!!! \n";
        gsInfo << " value was: " << minGaussD << "\n";
        return 0;
    }

    minD = dJC.provePositivityOfDetJ_TP(neededRefSteps, 5);

    out(0) = minD;
    out(1) = neededRefSteps;

    std::stringstream stream;
    stream << BASE_FOLDER << output << name;
    std::string str = stream.str();
    saveMat(out,str);
	
	return minD > 0;
}

gsMultiPatch<> getInitGuess3d(gsMultiPatch<> &goal)
{
	index_t n = 2;
	index_t degree = 1;

	index_t dim = goal.targetDim();

	gsMultiPatch<> out;
	for( index_t p = 0; p < goal.nBoxes(); p++)
	{	
		gsKnotVector<> kv1(0, 1, n - degree - 1, degree + 1);
		gsKnotVector<> kv2(0, 1, n - degree - 1, degree + 1);
		gsKnotVector<> kv3(0, 1, n - degree - 1, degree + 1);
		// 2. construction of a basis
		gsTensorBSplineBasis<3, real_t> basis(kv1, kv2, kv3);
		// 3. construction of a coefficients
		gsMatrix<> coefs (n*n*n, 3);
		gsTensorBSpline<3, real_t>  tbsout(basis, coefs);

		gsTensorBSpline<3, real_t> tbs = static_cast< gsTensorBSpline<3, real_t>& > (goal.patch(p));

		for( boxCorner bc = boxCorner::getFirst(dim); bc < boxCorner::getEnd(dim); ++bc)
		{
			gsMatrix<>::RowXpr coef = tbs.coefAtCorner(bc);

			// Plug it back into the new one
			unsigned i = tbsout.basis().functionAtCorner(bc);

			for (index_t d = 0; d < dim; d++)
			{
				tbsout.coef(i,d) = coef[d];
			}
			
		}

		for (index_t d = 0; d < dim; d++)	
		{
			tbsout.degreeElevate(tbs.degree(d)-degree,d);

			gsKnotVector<> kv = tbs.knots(d);

			for(gsKnotVector<>::iterator it = kv.begin() + tbs.degree(d) + 1; it != kv.end() - tbs.degree(d) - 1; it++)
			{
				tbsout.insertKnot(*it,d);
			}


		}
		
		out.addPatch(tbsout);
	}
	out.computeTopology();
	//out.closeGaps();
	return out;

}

gsMultiPatch<> getInitGuess3d_unitSq(gsMultiPatch<> &goal)
{
	index_t n = 2;
	index_t degree = 1;

	index_t dim = goal.targetDim();

	gsMultiPatch<> out;
	for( index_t p = 0; p < goal.nBoxes(); p++)
	{	
		gsKnotVector<> kv1(0, 1, n - degree - 1, degree + 1);
		gsKnotVector<> kv2(0, 1, n - degree - 1, degree + 1);
		gsKnotVector<> kv3(0, 1, n - degree - 1, degree + 1);
		// 2. construction of a basis
		gsTensorBSplineBasis<3, real_t> basis(kv1, kv2, kv3);
		// 3. construction of a coefficients
		gsMatrix<> coefs = basis.anchors().transpose();
		gsTensorBSpline<3, real_t>  tbsout(basis, coefs);

		gsTensorBSpline<3, real_t> tbs = static_cast< gsTensorBSpline<3, real_t>& > (goal.patch(p));

		for (index_t d = 0; d < dim; d++)	
		{
			tbsout.degreeElevate(tbs.degree(d)-degree,d);

			gsKnotVector<> kv = tbs.knots(d);

			for(gsKnotVector<>::iterator it = kv.begin() + tbs.degree(d) + 1; it != kv.end() - tbs.degree(d) - 1; it++)
			{
				tbsout.insertKnot(*it,d);
			}


		}
		
		out.addPatch(tbsout);
	}
	//out.computeTopology();
	//out.closeGaps();
	return out;

}

gsMultiPatch<> getInitGuess3d_full(gsMultiPatch<> goal)
{
    gsMultiPatch<> out(goal);
   
    for (index_t p = 0; p < goal.nBoxes(); p++)
    {
        gsMultiPatch<> sp1(goal.patch(p));

        gsMultiPatch<> sp_init = getInitGuess3d_unitSq(sp1);
        sp_init.computeTopology();

        gsMultiPatch<>::Ptr sp_ptr = memory::make_shared_not_owned(&sp1);
        gsMultiPatch<>::Ptr sp_init_ptr = memory::make_shared_not_owned(&sp_init);

        gsInfo << *sp_ptr;
        gsInfo << *sp_init_ptr;
        gsOptInit4th OI(sp_init_ptr, sp_ptr);
        OI.print();

        OI.solve();
        OI.updateMP();

        out.patch(p).setCoefs( sp_init.patch(0).coefs() );
    }

    return out;

}

gsMultiPatch<> getInitGuess2d_unitSq(gsMultiPatch<> &goal)
{
	index_t n = 2;
	index_t degree = 1;

	index_t dim = goal.targetDim();

	gsMultiPatch<> out;
	for( index_t p = 0; p < goal.nBoxes(); p++)
	{	
		gsKnotVector<> kv1(0, 1, n - degree - 1, degree + 1);
		gsKnotVector<> kv2(0, 1, n - degree - 1, degree + 1);
		// 2. construction of a basis
		gsTensorBSplineBasis<2, real_t> basis(kv1, kv2);
		// 3. construction of a coefficients
		gsMatrix<> coefs = basis.anchors().transpose();
		gsTensorBSpline<2, real_t>  tbsout(basis, coefs);

		gsTensorBSpline<2, real_t> tbs = static_cast< gsTensorBSpline<2, real_t>& > (goal.patch(p));

		for (index_t d = 0; d < dim; d++)	
		{
			tbsout.degreeElevate(tbs.degree(d)-degree,d);

			gsKnotVector<> kv = tbs.knots(d);

			for(gsKnotVector<>::iterator it = kv.begin() + tbs.degree(d) + 1; it != kv.end() - tbs.degree(d) - 1; it++)
			{
				tbsout.insertKnot(*it,d);
			}


		}
		
		out.addPatch(tbsout);
	}
	//out.computeTopology();
	//out.closeGaps();
	return out;

}

gsMultiPatch<> getInitGuess2d_full(gsMultiPatch<> goal)
{
    gsMultiPatch<> out(goal);
   
    for (index_t p = 0; p < goal.nBoxes(); p++)
    {
        gsMultiPatch<> sp1(goal.patch(p));

        gsMultiPatch<> sp_init = getInitGuess2d_unitSq(sp1);
        sp_init.computeTopology();

        gsMultiPatch<>::Ptr sp_ptr = memory::make_shared_not_owned(&sp1);
        gsMultiPatch<>::Ptr sp_init_ptr = memory::make_shared_not_owned(&sp_init);

        gsInfo << *sp_ptr;
        gsInfo << *sp_init_ptr;
        gsOptInit4th OI(sp_init_ptr, sp_ptr);
        OI.print();

        OI.solve();
        OI.updateMP();

        out.patch(p).setCoefs( sp_init.patch(0).coefs() );
    }

    return out;

}

gsMultiPatch<> getInitGuess2d(gsMultiPatch<> &goal)
{
	index_t n = 2;
	index_t degree = 1;

	index_t dim = goal.targetDim();

	gsMultiPatch<> out;
	for( index_t p = 0; p < goal.nBoxes(); p++)
	{	
		gsKnotVector<> kv1(0, 1, n - degree - 1, degree + 1);
		gsKnotVector<> kv2(0, 1, n - degree - 1, degree + 1);
		// 2. construction of a basis
		gsTensorBSplineBasis<2, real_t> basis(kv1, kv2);
		// 3. construction of a coefficients
		gsMatrix<> coefs (n*n, 2);
		gsTensorBSpline<2, real_t>  tbsout(basis, coefs);

		gsTensorBSpline<2, real_t> tbs = static_cast< gsTensorBSpline<2, real_t>& > (goal.patch(p));

		for( boxCorner bc = boxCorner::getFirst(dim); bc < boxCorner::getEnd(dim); ++bc)
		{
			gsMatrix<>::RowXpr coef = tbs.coefAtCorner(bc);

			// Plug it back into the new one
			unsigned i = tbsout.basis().functionAtCorner(bc);

			for (index_t d = 0; d < dim; d++)
			{
				tbsout.coef(i,d) = coef[d];
			}
			
		}

		for (index_t d = 0; d < dim; d++)	
		{
			tbsout.degreeElevate(tbs.degree(d)-degree,d);

			gsKnotVector<> kv = tbs.knots(d);
			for(gsKnotVector<>::iterator it = kv.begin() + tbs.degree(d) + 1; it != kv.end() - tbs.degree(d) - 1; it++)
			{
				tbsout.insertKnot(*it,d);
			}


		}
		
		out.addPatch(tbsout);
	}
	out.computeTopology();
	out.closeGaps();
	return out;

}

real_t measureSize(gsGeometry<> &geo)
{
    gsExprAssembler<> A(1,1);

    gsMultiPatch<> mp(geo);
    gsMultiBasis<> dbasis(mp);

    A.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);

    gsExprAssembler<>::geometryMap G = A.getMap(mp);

    return ev.integral(meas(G));

}

void scaleMP(gsMultiPatch<> &mp)
{
    real_t size = measureSize(mp.patch(0));
    gsDebugVar(size);

    for (index_t p = 0; p < mp.nBoxes(); p++)
        mp.patch(p).scale(1.0/(pow(size,1.0/mp.targetDim())));

    size = measureSize(mp.patch(0));
    gsDebugVar(size);

}

bool contains(std::string str, std::string substr)
{   return str.find(substr) != std::string::npos; }

void testOfNvDeriv(gsMultiPatch<> mp)
{
    gsMultiBasis<> dbasis(mp);
    
	gsExprAssembler<> A(1,1);
	gsExprEvaluator<> ev(A);

    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;

    geometryMap G = A.getMap(mp);

    A.setIntegrationElements(dbasis);

	gsMultiBasis<> geom_basis(mp);

    space u = A.getSpace(dbasis);
    space v = A.getTestSpace(u,geom_basis,3);

    gsFunctionExpr<> fun("x",3);
    variable ff = A.getCoeff(fun);
    auto vec = fjac(ff);

    gsFunctionExpr<> zerfun("0.0",3);
    variable zer = ev.getVariable(zerfun);
    auto zero_vec = v*vec;     // Perhaps I need to multiply with u;

    gsBoundaryConditions<> bcInfo;
    for (index_t i = 0; i < mp.nBoundary(); i++)
    {
        patchSide ps = mp.boundaries()[i];
        bcInfo.addCondition(ps.patch,ps.index(),condition_type::neumann,zerfun);

    }

    A.initSystem();
    A.assemble( -2*nvDeriv(v,G)*nv(G)*u.tr(), zero_vec );
	//A.assembleLhsRhsBc(nvDeriv(v,G)*nv(G)*u.tr(), zero_vec , bcInfo.neumannSides()); 

}

int main(int argc, char* argv[]){
gsInfo <<  "Hello G+Smo.\n";

// Parse command line
std::string output("/");
int degree = 2;
int nx = 5;
int ny = 4;

int numRefine = 1;
int maxiter = 10;

int param = 0; // 0: spring, 1: modLiao, 2: winslow
int jig = 1; 
int dim = 2;
int patch = -1;

bool plotDesign = false;
bool plotMagnitude = false;
bool plotSolution = false;
bool saveCps = false;
bool useDJC = true;
bool useLag = false;

int quA = 2;
int quB = 2;

int quA_optParam = 2;
int quB_optParam = 2;

real_t lambda_1 = 0;
real_t lambda_2 = 1;

int startDes = -1;

bool startFromFile = false;
std::string startFile("");

bool optParam = false;
bool optParamXML = false;
bool potWave = false;

real_t alpha = 0; // Parameter for constraint aggregation
real_t eps = 0; // Parameter for constraint aggregation

bool decreaseTau = false;
real_t decrTauFactor = 1.0/4.0;

gsCmdLine cmd("A test of lumped mass matricies");
cmd.addInt("p", "degree", "Degree of B-Splines.", degree);

cmd.addInt("a", "param", "Parametrization: \
        0: spring, \
        1: modLiao, \
        2: winslow, \
        3: liao, \
        4: harmonic, \
        5: use regularization and perform optimization on all cps, \
        6: use winslow with out const5raints, \
        7: maxDetJ, \
        8 : Harmonic without constraints", param);

cmd.addInt("r", "numberRefine", "Number of refinements", numRefine);
cmd.addInt("n", "nx", "Number of splines in first direction", nx);
cmd.addInt("m", "ny", "Number of splines in second direction", ny);
cmd.addInt("i", "maxIter", "Maximal number of reparametrizations", maxiter);

cmd.addInt("j","jig","Which domain to try to parametrize",jig);
cmd.addInt("d","dim","2d or 3d",dim);
cmd.addInt("t","patch","patch",patch);

cmd.addInt("A", "quA", "quA", quA);
cmd.addInt("B", "quB", "quB", quB);

cmd.addInt("3", "quA_optParam", "quA for winslow parametrization", quA_optParam);
cmd.addInt("4", "quB_optParam", "quB for winslow parametrization", quB_optParam);

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

cmd.addSwitch("startFromFile", "Start optimization with design from file as starting guess", startFromFile);
cmd.addString("f", "startFile", "design to start from", startFile);

cmd.addSwitch("optParam", "Run optParam code", optParam);
cmd.addSwitch("optParamXML", "Run optParamXML code", optParamXML);
cmd.addSwitch("potWave", "Run potWave code", potWave);

bool paramTestXML = false;
cmd.addSwitch("paramTestXML", "Run paramTestXML code", paramTestXML);

cmd.addSwitch("decreasingTau", "decrease tau", decreaseTau);
cmd.addReal("c","decrTauFactor", "factor with which to decrease tau", decrTauFactor);

bool useCorner = false;
cmd.addSwitch("useCorner", "use only a corner of multipatch", useCorner);

bool usePow = false;
cmd.addSwitch("usePow", "use winslow: JTJ/detJ^(2/d)", usePow);

bool changeSign = false;
cmd.addSwitch("changeSign", "change sign of detJ", changeSign);

bool glueInterfaces = true;
cmd.addSwitch("glueInterfaces", "glue interfaces? If yes we use multi-linear start guess, otherwise we use the affine startguess", glueInterfaces);

std::string plotDir;
cmd.addString("9", "plotDir", "design to start from", plotDir);

cmd.getValues(argc,argv);

output = "/" + output;
startFile = "/" + startFile;

char buffer [50];
std::sprintf(buffer,"p = %d \n n = %d\n",degree,numRefine);
gsInfo << buffer;

gsShapeOptLog slog(output);

// Test 3D PML and water waves
if (potWave) {
    //gsStateEquationPotWaves SE(numRefine, 3 + lambda_1,3 + lambda_1, 2 + 0.5*lambda_1, 2 + lambda_2, 2 + lambda_2,1 + 0.5*lambda_2);

    //SE.pointSourceTest(BASE_FOLDER + output);
    //SE.convergenceTestOnlyPML( maxiter , BASE_FOLDER + output );
    //SE.convergenceTestNoPML_NoCenter( maxiter , BASE_FOLDER + output );

    gsStateEquationPotWaves SE(numRefine);
    
    //SE.convergenceTestBessel( maxiter, BASE_FOLDER + output );
    //real_t timestep = 0.025*2*M_PI/SE.wave_omega;
    //SE.plotVelocityBessel(timestep , BASE_FOLDER + output + "velocity/veloc");
    //return 0;
    
    gsShapeOptLog slog1(output,true,false,false);
	gsShapeOptLog::Ptr slog1_ptr = memory::make_shared_not_owned(&slog1); 

    gsMultiPatch<>::Ptr mp_ptr = SE.m_mp;
    gsMultiPatch<>::Ptr c_ptr = SE.getCenterReflector();

    gsOptPotWaves optPW(mp_ptr,numRefine,slog1_ptr,param,quA,quB);

    //convergenceTestOfJacobianAll(optPW);
    //

    bool useConstraints = false;

    gsOptPotWaves::Ptr optPW_ptr = memory::make_shared_not_owned( &optPW );
    gsOptPotWavesWithReg2nd optWR(mp_ptr,c_ptr,optPW_ptr,numRefine,slog1_ptr,quA,quB,eps,true,true);

    //gsInfo << " Objective : " << optPW_ptr->evalObj() << "\n";

    //std::string name = "wfun";
    //optPW_ptr->plotGoalFunctions(BASE_FOLDER + output + name);

    //optPW_ptr->gradAll();

    //gsDebugVar(optWR.m_winslow->evalObj());

    if (plotDesign)
    {
    
        std::string nameDetJ = plotDir + "detJ.txt";
		std::ofstream file(nameDetJ);
        for( index_t i = 0; i < 3000; i += 50)
        {

            std::string si = std::to_string(i);
            std::string nameUR = plotDir + "xml/ur_" + si + ".xml";
            std::string nameUI = plotDir + "xml/ui_" + si + ".xml";

            if (!exists(nameUR))
            {
                gsInfo << nameUR << "\n";
                gsInfo << "BREAK\n";
                break;
            }

            file << "i = " << i << "\n";

            // Load solution
	        gsMultiPatch<>::uPtr ur_ptr, ui_ptr;

	        gsFileData<> fdr(nameUR);
	        ur_ptr = fdr.getFirst< gsMultiPatch<> > ();

	        gsFileData<> fdi(nameUI);
	        ui_ptr = fdi.getFirst< gsMultiPatch<> > ();

            std::string nameCps = plotDir + "cps_0_" + si + ".txt";
            gsInfo << nameCps << "\n";

            gsWinslowPow winslow(mp_ptr,optWR.mappers(),false, false, true, 0);
            winslow.setQuad(quA,quB);

            winslow.updateFlat(loadVec(winslow.n_flat,nameCps));

            real_t minDJ0 = winslow.minDetJInGaussPts();
            real_t minDJ15 = winslow.minDetJInGaussPts(15);

            // Compute detJ

            file << "min in qauss pts (" << quA << ", " << quB << ") :" << minDJ0 << "\n";
            file << "min in qauss pts (" << quA << ", " << quB + 15 << ") :" << minDJ15 << "\n";

            if (minDJ15 > 0)
            {
                gsDetJacConstraint dJC(mp_ptr, true);
                index_t neededSteps;
                real_t minJacCoeff = dJC.provePositivityOfDetJ_TP(neededSteps, 3);
                file << "minJacCoeff " << minJacCoeff << ": " << minJacCoeff << "\n";
            }

            // Plot |u|^2
            std::string nameAbsU = plotDir + "paraview/absU_" + si;
            SE.plotMagnitude(*ur_ptr, *ui_ptr, nameAbsU);

            // Plot center as slices
            gsVector<> dir(5);
            gsVector<> par(5);
            dir << 0, 0, 1, 2, 1;
            par << 1, 1, 0, 1, 0;
            for (index_t p = 0; p < 5; p++)
            {
                std::string nameSlice = plotDir + "paraview/side_" + si + "_" + std::to_string(p); 
                gsGeometrySlice<> slice = mp_ptr->patch(p).getIsoParametricSlice(dir[p],par[p]);

                gsWriteParaview(slice, nameSlice, 500);
            }

            file << "================ \n\n";
        }
        file.close();
    
    }
    else
    {

        optWR.runOptimization();

    }


    // ========================================
    // =========   Test constraints   =========
    // ========================================
    
    /*
    gs2NormConstraints con2(mp_ptr,optWR.mappers(), 0.2, 0.5);

    gsWinslow win(mp_ptr,optWR.mappers(),false);

    gsInfo<< win.getTagged().transpose() << "\n\n\n";
    convergenceTestOfConstraintFree(win, con2);
    */
    //
    //
    //mp_ptr->closeGaps();

    //con2.evalCon();
    //con2.getJacobian(win.getFree());

    // ===========================k=============
    // =========  Plot velocity field =========
    // ========================================
    
    /*
    SE.setQuad(quA,quB);
 
    gsMultiPatch<> ur, ui;
    SE.solve(ur,ui);
 
    SE.plotSolution(ur,BASE_FOLDER + output + "ur");
    SE.plotSolution(ui,BASE_FOLDER + output + "ui");
 
    real_t timestep = 0.025*2*M_PI/SE.wave_omega;
 
    SE.plotVelocityField(ur,ui,timestep,BASE_FOLDER + output + "velocity/veloc", true);
    SE.plotVelocityField(ur,ui,timestep,BASE_FOLDER + output + "velocity/scatt", false);
    

    gsInfo << "\n -------------- \n";
    */

    return 0;

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

gsMultiPatch<>::Ptr patches_ptr = memory::make_shared_not_owned(&patches); 
gsMultiPatch<>::Ptr mp_ptr = memory::make_shared_not_owned(&mp);

gsMultiBasis<> bas(patches);

// gsMultiPatch<> patches3d = get3DGeometry();

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
//
//

// get Initial guess of seastar and torus
if (false)
{
    for (index_t i = 0; i < 2; i++)
    {
        gsMultiPatch<> seastar;
        std::string nm;

        if (i == 0)
        {
            seastar = getSeastar();
            nm = "seastar";
        } else {
            seastar = getTorus();
            nm = "torus";
        }

        gsMultiPatch<>::Ptr seastar_ptr = memory::make_shared_not_owned( &seastar );

        gsMultiPatch<> mp_init1 = getInitGuess2d(seastar);
        gsMultiPatch<>::Ptr init1_ptr = memory::make_shared_not_owned( &mp_init1 );

        gsMultiPatch<> mp_init2 = getInitGuess2d_full(seastar);
        gsMultiPatch<>::Ptr init2_ptr = memory::make_shared_not_owned( &mp_init2 );

        gsWinslow winS( seastar_ptr, false );
        std::string name = BASE_FOLDER + output + nm;
        gsVector<> flat = winS.getFlat();
        saveVec(flat, name + ".txt");
        gsWriteParaview(seastar,name,true,true);

        gsWinslow win1( init1_ptr, false );
        name = BASE_FOLDER + output + nm + "_init1";
        flat = win1.getFlat();
        saveVec(flat, name + ".txt");
        gsWriteParaview(mp_init1,name,true,true);

        gsWinslow win2( init2_ptr, false );
        name = BASE_FOLDER + output + nm + "_init2";
        flat = win2.getFlat();
        saveVec(flat, name + ".txt");
        gsWriteParaview(mp_init2,name,true,true);
    }

    exit(0);




}

// Print number of constraints
if (false)
{

    gsInfo << "n  \t p=1 \t p=2 \t p=3 \n";
    for (index_t n = 2; n < 17; n *= 2)
    {

        gsInfo << n;
        for (index_t d = 2; d < 4; d++)
        {

            for (index_t p = 1; p < 4; p++)
            { 
	            gsKnotVector<> kv1(0, 1, n , p + 1);
	            gsKnotVector<> kv2(0, 1, n , p + 1);
	            gsKnotVector<> kv3(0, 1, n , p + 1);

                gsBasis<>::uPtr bas_ptr;
                
                if (d == 2) // 2D
                {
	                bas_ptr = memory::make_unique( new gsTensorBSplineBasis<2, real_t>(kv1, kv2) );
                }
                else // 3D
                {
	                bas_ptr = memory::make_unique( new gsTensorBSplineBasis<3, real_t>(kv1, kv2, kv3) );
                }
	            

                gsMultiBasis<> dbas(*bas_ptr);
                dbas.setDegree(d*p-1);
                dbas.reduceContinuity(1);

                gsInfo << "\t&\t" << bas_ptr->size() << ", " <<  dbas.size();
            }
        }
        gsInfo << "\n";
    }



    exit(0);

}

if (false)
{
    std::stringstream stream;
    stream.precision(16);

    real_t wave_K = 4;

    // --------- J0 --------- //
    index_t N = 50;
    stream << "exp( " << wave_K << " * z ) *(";
    for (index_t k = 0; k < N; k++)
    {
        stream << "(-1)^" << k << "*";
        stream << "(1/4.0 * " << wave_K*wave_K << "*(x^2 + y^2))^" << k << "*";

        // 1/(k!)^2
        stream << "1/(1.0" ;
        for (index_t i = 2; i <= k; i++) // k!
        {
            stream << "*" << i ;
        }
        stream << ")^2 ";

        if (k < N-1)
            stream << " + ";

    }
    stream << ")";

    std::string J0str = stream.str();

    // ----------- Y0 ------------ //
    stream.clear();
    stream.str(std::string());

    stream << "exp( " << wave_K << " * z ) *(";

    // First part
    real_t gamma = 0.57721566490153286060;
    stream << 2.0/M_PI << "* ( log( 1/2 * " << wave_K << "*sqrt(x^2 + y^2) ) + " << gamma << ")*";
    stream << "( " << J0str << " )";

    stream << " + ";
    // Second part
    stream << 2.0/M_PI << " * ( ";

    real_t sum = 0;
    for (index_t n = 1; n < N; n++)
    {
        stream << "(-1)^" << n+1 << "*";
        stream << "(1/4.0 * " << wave_K*wave_K << "*(x^2 + y^2))^" << n << "*";

        // 1/(k!)^2
        stream << "1/(1.0" ;
        for (index_t i = 2; i <= n; i++) // k!
        {
            stream << "*" << i ;
        }
        stream << ")^2 ";

        // 1 + 1/2 + ...
        sum += 1.0/n;

        stream << " * " << sum ;

        if (n < N-1)
            stream << " + ";

    }

    stream << " ) )";

    std::string Y0str = stream.str();

    gsFunctionExpr<> J0(J0str,2);

    index_t neval = 2;
    gsMatrix<> u(2,neval);
    u << 0.5,5,
        0.5,5;

	gsMatrix<> fun_eval_mat = J0.eval(u);

    gsDebugVar(wave_K);
    gsDebugVar(pow(u(0,0),2));
    gsDebugVar(pow(u(1,0),2));
    gsDebugVar(sqrt(pow(u(0,0),2) +  pow(u(1,0),2)));

    for (index_t i = 0; i < neval; i++)
    {
        real_t fun_eval = fun_eval_mat(0,i);

        gsInfo << "J0(" << wave_K*sqrt(pow(u(0,i),2) +  pow(u(1,i),2)) << ") = " << fun_eval << "\n";
    }

    gsFunctionExpr<> Y0(Y0str,2);

	gsMatrix<> y_eval_mat = Y0.eval(u);

    for (index_t i = 0; i < neval; i++)
    {
        real_t fun_eval = y_eval_mat(0,i);

        gsInfo << std::setprecision(12) << "Y0(" << u(0,i) << ", " << u(1,i) << ") = " << fun_eval << "\n";
    }


    return 0;
   

}

// Save as XML
if (false)
{
   gsMultiPatch<> jigsaw3d = getJigSaw3D(jig);

   for (index_t i = 0; i < jigsaw3d.nBoxes(); i++)
        changeSignOfDetJ(jigsaw3d.patch(i));

   gsFileData<> fd;//
   fd << jigsaw3d;

   std::string name = BASE_FOLDER + output;
   gsInfo << "Save in " << name << "\n";
   fd.save( name ); 
   return 0;
}


if(paramTestXML)
{
	gsFileData<> fd(BASE_FOLDER + startFile);

	gsMultiPatch<>::uPtr mp_uptr;

    std::vector< gsDofMapper > mappers(dim);
    gsDebugVar(useCorner);
    if (useCorner)
    {
        gsInfo << "USECORNER\n";
	    gsMultiPatch<> corner = loadCorner(jig,dim,mappers,numRefine);
        mp_uptr = memory::make_unique( new gsMultiPatch<>(corner) );

    }
    else
    {
	    mp_uptr = fd.getFirst< gsMultiPatch<> > ();

	    for(int i = 0; i < numRefine; i++){
	    	mp_uptr->patch(0).uniformRefine();
	    }

    }

    scaleMP(*mp_uptr);

    if (changeSign)
    {
        for (index_t i = 0; i < mp_uptr->nBoxes(); i++)
            changeSignOfDetJ(mp_uptr->patch(i));
    }

    gsMultiPatch<> mp(*mp_uptr);
    mp_ptr = memory::make_shared_not_owned(&mp);

    // START FROM SPRING METHOD
    gsSpringMethod tmp(mp_ptr);
    tmp.update();

    // Put mappers in var
    if (!useCorner)
    {
        mappers = tmp.mappers();
    }

    // Get startGuess
    if ( param == 1 || param == 2 || param == 3 || param == 4 || param == 6 )
    {
        gsMaxDetJac mDJ(mp_ptr,mappers,true); // Use TP solver
        mDJ.update();

        std::string name = "init_detJProof.txt";
        bool valid = testValidity(*mp_ptr,name, quA, quB, output);

        if (!valid)
        {
            gsInfo << "\n ------------ !!!!!! --------------- \n Start guess was not found! \n\n";
        }


    }

    std::string name = "init.txt";
    gsVector<> flat = tmp.getFlat();
    saveVec(flat, BASE_FOLDER + output + name);

    gsParamMethod::Ptr pM_ptr;

    switch (param) {
        case 0 : pM_ptr = memory::make_shared( new gsSpringMethod(mp_ptr,mappers));break; 
        case 1 : pM_ptr = memory::make_shared( new gsModLiao(mp_ptr,mappers, true, true));break;   // detJ
        case 2 : pM_ptr = memory::make_shared( new gsWinslow(mp_ptr,mappers, true, true));break;
        case 3 : pM_ptr = memory::make_shared( new gsLiao(mp_ptr,mappers, true, true));break;
        case 4 : pM_ptr = memory::make_shared( new gsHarmonic(mp_ptr,mappers, true, true));break;
        case 6 : pM_ptr = memory::make_shared( new gsWinslow(mp_ptr,mappers, false, false, true, 0));break;        // No detJ
        case 7 : pM_ptr = memory::make_shared( new gsMaxDetJac(mp_ptr,mappers));break;
        case 8 : pM_ptr = memory::make_shared( new gsHarmonic(mp_ptr,mappers,false));break;        // No detJ
    }

    real_t m_eps = 0.00001;
    if ( param == 1 || param == 2 || param == 3 || param == 4 || param == 6 )
    {
        (std::dynamic_pointer_cast< gsOptParamMethod >(pM_ptr))->m_dJC->setEps(m_eps);
        (std::dynamic_pointer_cast< gsOptParamMethod >(pM_ptr))->setupOptParameters();
    }

    pM_ptr->update();

    name = "detJProof.txt";
    testValidity(*mp_ptr, name, quA, quB, output);

    name = "mp";
    gsWriteParaview(*mp_ptr, BASE_FOLDER + output + name,10000,true,true);

    name = "cps";
    flat = pM_ptr->getFlat();
    saveVec(flat, BASE_FOLDER + output + name);


    return 0;
    
}


// Calculate E_{h/4} for a specific design, and write to file
// --startFromFile
// -startFile "path/to/file"
// -o "folder/to/output/result/"
// -r THE REFINEMENT FROM OPTIMIZATION
//
if (false) 
{
    gsShapeOptLog slog1(output,true,false,false);
	gsShapeOptLog::Ptr slog1_ptr = memory::make_shared_not_owned(&slog1); 

	gsWinslow winslow(mp_ptr);
    gsInfo << "\nLOAD FROM " << startFile << "\n\n";
	winslow.updateFlat( loadVec(winslow.n_flat,startFile));

    gsOptAntenna optA(mp_ptr,numRefine + 2,slog1_ptr,param,quA,quB,true);

    gsInfo << "write to " << output + "Eh4.txt" << "\n";
	std::ofstream file (output + "Eh4.txt",std::ofstream::app);
	file << std::setprecision(12);

    real_t obj = optA.evalObj();
    gsInfo << "Eh4 : " << obj << "\n";
    file << startFile << " " << obj << "\n";
    file.close();
    exit(0);
}

// Check number of uniform constraints
if (false)
{
    gsShapeOptLog slog1(output,true,false,false);
	gsShapeOptLog::Ptr slog1_ptr = memory::make_shared_not_owned(&slog1); 

    gsOptAntenna optA(mp_ptr,numRefine,slog1_ptr,param,quA,quB,true);
    exit(0);
    gsWinslow winslow(mp_ptr,optA.mappers(),false,false,true,0);

    winslow.setQuad(quA,quB);

    winslow.updateFlat(loadVec(winslow.n_flat,output));

    std::string nm = "mp_ptr";
	gsWriteParaview(*mp_ptr,nm,10000,true,true);

    gsDetJacConstraint dJC(mp_ptr);
    gsInfo << "mind = " << dJC.evalCon().minCoeff() << "\n";

    dJC.refineUntilPositiveUniformly(maxiter);

    gsInfo << "\n\nmind = " << dJC.evalCon().minCoeff() << "\n";
    gsInfo << "size = " << dJC.evalCon().size() << "\n";
    return 0;

}


// Test of runtim
if (startDes == 20) 
{
    gsInfo << "test mb\n";

    gsShapeOptLog slog1(output,true,false,false);
	gsShapeOptLog::Ptr slog1_ptr = memory::make_shared_not_owned(&slog1); 
    // gsInfo << optWR.evalObj();
    // gsInfo << optWR.gradObj();
    // optWR.solve();

    if (param == 5) // Use regularization
    {
	
        gsOptAntenna optA(mp_ptr,numRefine,slog1_ptr,6,quA,quB);
		gsOptAntenna::Ptr optA_ptr = memory::make_shared_not_owned(&optA);

        gsShapeOptWithReg optWR(mp_ptr,optA_ptr,numRefine,slog1_ptr,quA,quB,eps);
        //optWR.setWinslowQuad(quA_optParam, quB_optParam);
	
        optWR.solve();
    } else if (param == 6) {
        gsOptAntenna optA(mp_ptr,numRefine,slog1_ptr,param,quA,quB,true);
        // optA.m_paramMethod->updateFlat(loadVec(optA.n_flat,BASE_FOLDER + output));
        // gsInfo << optA.evalObj() << "\n";
        // exit(0);
        // optA.solve();
        optA.runOptimization(maxiter,true); // True means that we use uniform refinement
		// convergenceTestOfJacobian(optA);
    } else if (param == 0) {
        gsOptAntenna optA(mp_ptr,numRefine,slog1_ptr,param,quA,quB,true);
        //optA.m_paramMethod->updateFlat(loadVec(optA.n_flat,BASE_FOLDER + output));
        gsInfo << optA.evalObj() << "\n";
    }

    return 0;

}

if(false)
{
	extractCornersAndSaveXML(jig,dim);
	return 0;
}

if(optParamXML)
{
	gsFileData<> fd(BASE_FOLDER + startFile);

	gsMultiPatch<>::uPtr mp_uptr;

    std::vector< gsDofMapper > mappers(dim);
    gsDebugVar(useCorner);
    if (useCorner)
    {
        gsInfo << "USECORNER\n";
	    gsMultiPatch<> corner = loadCorner(jig,dim,mappers,numRefine);
        mp_uptr = memory::make_unique( new gsMultiPatch<>(corner) );

    }
    else
    {
	    mp_uptr = fd.getFirst< gsMultiPatch<> > ();

	    for(int i = 0; i < numRefine; i++){
	    	mp_uptr->patch(0).uniformRefine();
	    }
    }

    if (changeSign)
    {
        for (index_t i = 0; i < mp_uptr->nBoxes(); i++)
            changeSignOfDetJ(mp_uptr->patch(i));
    }

    scaleMP(*mp_uptr);

    //std::stringstream stream;
    //stream << BASE_FOLDER << output << "jig" << dim << "d_" << jig;
    //std::string nm = stream.str();
    //gsWriteParaview(*mp_uptr,nm,10000);

    gsMultiPatch<> mp(*mp_uptr);
    mp_ptr = memory::make_shared_not_owned(&mp);
    
	gsMultiPatch<> mp_init;
	gsMultiPatch<>::Ptr mp_init_ptr = memory::make_shared_not_owned(&mp_init);

	if (dim == 2) {
    	mp_init = getInitGuess2d(*mp_ptr);
	} else if (dim == 3) {
        if (useCorner)
        {
    	    //mp_init = getInitGuess3d_unitSq(*mp_ptr);
    	    mp_init = getInitGuess3d_full(*mp_ptr);
        }
        else
        {
            if (glueInterfaces)
    	        mp_init = getInitGuess3d(*mp_ptr);
            else
                mp_init = getInitGuess3d_full(mp);
        }
        //mp_init = getInitGuess3d_full(mp);
	}

    for (index_t p = 0; p < mp_init.nBoxes(); p++)
    {
        gsMultiPatch<> sp(mp_init.patch(p));
	    testValidity(sp,"detJ" + std::to_string(p) + ".txt", quA, quB, output);
    }

    gsShapeOptLog slog1(output,true,false,false);
	gsShapeOptLog::Ptr slog1_ptr = memory::make_shared_not_owned(&slog1); 

    std::string name = "init";
    slog1.plotInParaview(mp_init,name);

    name = "goal";
    slog1.plotInParaview(*mp_ptr,name);


    gsOptParam::Ptr optP_ptr;
    if (useCorner)
    {
        gsInfo << "USECORNER\n";
	    gsOptParam optP(mp_init_ptr,mp_ptr,mappers,slog1_ptr,param);
        optP_ptr = memory::make_shared( new gsOptParam(mp_init_ptr,mp_ptr,mappers,slog1_ptr,param) );
    }
    else
    {
        if (glueInterfaces)
            optP_ptr = memory::make_shared( new gsOptParam(mp_init_ptr,mp_ptr,slog1_ptr,param) );
        else
            optP_ptr = memory::make_shared( new gsOptParamFull(mp_init_ptr,mp_ptr,slog1_ptr,param) );
        optP_ptr->setupMappers();
    }

	gsInfo << "EPS = " << eps << "\n";
    gsShapeOptWithReg optWR(mp_init_ptr,optP_ptr,numRefine,slog1_ptr,quA,quB,eps, glueInterfaces, usePow); // glue interfaces, usePow

    //optWR.evalObj();
    if (decreaseTau)
    {
        optWR.runOptimizationUntilPosDetJ(maxiter,decrTauFactor,5);
    }
    else
    {
        optWR.solve();
    }

 	name = "jigsaw";
    slog1.plotInParaview(mp,name);

 	name = "mp_final";
    slog1.plotInParaview(mp_init,name);

	//testValidity(mp_init,"detJ.txt", quA, quB, output);
	
 	// Test snaps
	gsWinslowPow win(mp_init_ptr,optP_ptr->mappers(),false,false,true,0); 
	gsWinslowPow win_jig(mp_ptr,optP_ptr->mappers(),false,false,true,0); 
	gsVector<> tagged_jig = win_jig.getTagged();
	win.updateTagged(tagged_jig);

 	name = "mp_snapped";
    slog1.plotInParaview(mp_init,name);

	std::stringstream stream_out;
	stream_out <<  output << "cps_snapped.txt";
	std::string str_out = stream_out.str();
	gsVector<> flat = win.getFlat();
	saveVec(flat,str_out);

    gsInfo << "Test validity of snap\n";
	testValidity(mp_init,"detJ_snap.txt", quA, quB, output);

 	name = "mp_snappedLin";
    slog1.plotInParaview(mp_init,name);
	
	// Snap bnd control points
	gsWinslow::Ptr win_ptr = memory::make_shared_not_owned(&win);
	gsAffineOptParamMethod lin(win_ptr,false);

	gsVector<> free = lin.getUpdate(tagged_jig);
	win.updateFreeAndTagged(free,tagged_jig);

	std::stringstream stream_out2;
    gsInfo << "Test validity of snap lin\n";
	stream_out2 <<  output << "cps_SnapLin.txt";
	str_out = stream_out2.str();
	flat = win.getFlat();
	saveVec(flat,str_out);

	testValidity(mp_init,"detJ_SnapLin.txt", quA, quB, output);

	return 0;
}

// Old optParamXML code
if(false)
{
	gsFileData<> fd(BASE_FOLDER + startFile);

	gsMultiPatch<>::uPtr mp_uptr;

    if ( contains(startFile,"twistedFlat") )
    {
        gsInfo << "I AM HERE\n";
        gsGeometry<>::uPtr geo_ptr;
	    geo_ptr = fd.getFirst< gsGeometry<> > ();

        mp_uptr = memory::make_unique(new gsMultiPatch<>(*geo_ptr));
    }
    else
    {
	    mp_uptr = fd.getFirst< gsMultiPatch<> > ();
    }

    scaleMP(*mp_uptr);

    gsMultiPatch<> mp(*mp_uptr);
    mp_ptr = memory::make_shared_not_owned(&mp);
    
    gsMultiPatch<> corner;
    std::vector< gsDofMapper > mappers(dim);

	//gsMultiPatch<> corner = loadCorner(jig,dim,mappers,numRefine);
	gsMultiPatch<>::Ptr corner_ptr = memory::make_shared_not_owned(&corner);

	gsMultiPatch<> mp_init;
	gsMultiPatch<>::Ptr mp_init_ptr = memory::make_shared_not_owned(&mp_init);

	if (dim == 2) {
    	mp_init = getInitGuess2d(*mp_ptr);
	} else if (dim == 3) {
    	//mp_init = getInitGuess3d(*mp_ptr);
    	mp_init = getInitGuess3d_unitSq(*mp_ptr);
	}

    mp_init = getInitGuess3d_full(mp);

    for (index_t p = 0; p < mp_init.nBoxes(); p++)
    {
        gsMultiPatch<> sp(mp_init.patch(p));
	    testValidity(sp,"detJ" + std::to_string(p) + ".txt", quA, quB, output);
    }

    gsShapeOptLog slog1(output,true,false,false);
	gsShapeOptLog::Ptr slog1_ptr = memory::make_shared_not_owned(&slog1); 


	for(index_t p = 0; p < mp.nBoxes(); p++)
	{
 		//changeSignOfDetJ(mp.patch(p)); // Use for 3D small
   		//changeSignOfDetJ(mp_init.patch(p)); // Use for 3D small

	}


    std::string name = "init";
    slog1.plotInParaview(mp_init,name);

    name = "goal";
    slog1.plotInParaview(*mp_ptr,name);

	gsOptParamFull optP(mp_init_ptr,mp_ptr,slog1_ptr,param);
    optP.setupMappers();
    gsOptParamFull::Ptr optP_ptr = memory::make_shared_not_owned( &optP);
    // gsMatrix<> v = optP.currentDesign()*0.99999;
    // optP.m_paramMethod->update();
    // gsInfo << optP.evalObj(v) << "\n";
    // gsInfo << optP.currentDesign() << "\n";
    // gsInfo << optP.numDesignVars() << "\n";
    // optP.solve();
	gsInfo << "EPS = " << eps << "\n";
    gsShapeOptWithReg optWR(mp_init_ptr,optP_ptr,numRefine,slog1_ptr,quA,quB,eps, false);

    //optWR.evalObj();
    optWR.solve();

    exit(0);
    // ------ //


    /*
	real_t scale = 0.2;
	for(index_t p = 0; p < corner.nBoxes(); p++)
	{
		corner.patch(p).scale(scale);
		mp_init.patch(p).scale(scale);
 		changeSignOfDetJ(corner.patch(p)); // Use for 3D small
   		changeSignOfDetJ(mp_init.patch(p)); // Use for 3D small
	}




	gsOptParam optP(mp_init_ptr,corner_ptr,mappers,slog1_ptr,param);
    gsOptParam::Ptr optP_ptr = memory::make_shared_not_owned( &optP);
    // gsMatrix<> v = optP.currentDesign()*0.99999;
    // optP.m_paramMethod->update();
    // gsInfo << optP.evalObj(v) << "\n";
    // gsInfo << optP.currentDesign() << "\n";
    // gsInfo << optP.numDesignVars() << "\n";
    // optP.solve();
	gsInfo << "EPS = " << eps << "\n";
    gsShapeOptWithReg optWR(mp_init_ptr,optP_ptr,numRefine,slog1_ptr,quA,quB,eps);
    // optWR.gradObj();
    optWR.solve();

	std::string nm = "jig";
	gsWriteParaview(corner,nm,10000,true,true);

	nm = "init";
	gsWriteParaview(mp_init,nm,10000,true,true);

	testValidity(mp_init,"detJ_init.txt", quA, quB, output);

    // gsMultiPatch<>(*gsNurbsCreator<>::BSplineSquare(2));

 	name = "jigsaw";
    slog1.plotInParaview(corner,name);

 	name = "mp_final";
    slog1.plotInParaview(mp_init,name);

	testValidity(mp_init,"detJ.txt", quA, quB, output);
	
 	// Test snaps
	gsWinslow win(mp_init_ptr,false,false,true,0); 
	gsWinslow win_jig(corner_ptr,false,false,true,0); 
	gsVector<> tagged_jig = win_jig.getTagged();
	win.updateTagged(tagged_jig);

	std::stringstream stream_out;
	stream_out <<  output << "cps_snapped.txt";
	std::string str_out = stream_out.str();
	gsVector<> flat = win.getFlat();
	saveVec(flat,str_out);

	testValidity(mp_init,"detJ_snap.txt", quA, quB, output);
	
	// Snap bnd control points
	gsWinslow::Ptr win_ptr = memory::make_shared_not_owned(&win);
	gsAffineOptParamMethod lin(win_ptr,false);

	gsVector<> free = lin.getUpdate(tagged_jig);
	win.updateFreeAndTagged(free,tagged_jig);

	std::stringstream stream_out2;
	stream_out2 <<  output << "cps_SnapLin.txt";
	str_out = stream_out2.str();
	flat = win.getFlat();
	saveVec(flat,str_out);

	testValidity(mp_init,"detJ_SnapLin.txt", quA, quB, output);

	return 0;
	return 0;
    */
}

// Test of opParam3D
if (false) 
{ 
	gsMultiPatch<> jigsaw = getJigSaw3D(2);
	gsMultiPatch<>::Ptr jigsaw_ptr = memory::make_shared_not_owned(&jigsaw);

	for(int i = 0; i < numRefine; i++){
		jigsaw.uniformRefine();
	}

    changeSignOfDetJ(jigsaw.patch(0)); // Use for 3D small
	
	gsMultiPatch<> mp_des(jigsaw);
	gsMultiPatch<>::Ptr mp_des_ptr = memory::make_shared_not_owned(&mp_des);
	    
	gsSpringMethod spring_jig(jigsaw_ptr);
	gsVector<> tagged_jig = spring_jig.getTagged();
	
	gsWinslow win(mp_des_ptr,false,false,true,0); 

	std::ofstream file (output + "stats.txt");
	file << std::setprecision(12);
	
	for (index_t i = 0; i < 3000; i++) {
	
	    std::stringstream stream_in;
	    stream_in <<  output << "cps" << "_0_" << i << ".txt";
	    std::string str_in = stream_in.str();
	
		if (!exists(str_in)) break;
		
		win.updateFlat(loadVec(win.n_flat,str_in));
	
		gsVector<> tagged = win.getTagged();
		gsVector<> diff = tagged - tagged_jig;
	
		file << i << " ";
		file << diff.norm() << " ";
		file << diff.maxCoeff() << " ";
		file << "\n";
	}
	file.close();

	// Snap bnd control points
//	win.updateTagged(tagged_jig);
//
//	std::stringstream stream_out;
//	stream_out <<  output << "cps_snapped.txt";
//	std::string str_out = stream_out.str();
//	gsVector<> flat = win.getFlat();
//	saveVec(flat,str_out);
//
//	testValidity(mp_des,"detJ_snap.txt", quA, quB, output);

	// Snap bnd control points
	gsWinslow::Ptr win_ptr = memory::make_shared_not_owned(&win);
	gsAffineOptParamMethod lin(win_ptr,false);

	gsVector<> free = lin.getUpdate(tagged_jig);
	win.updateFreeAndTagged(free,tagged_jig);

	std::stringstream stream_out;
	stream_out <<  output << "cps_SnapLin.txt";
	std::string str_out = stream_out.str();
	gsVector<> flat = win.getFlat();
	saveVec(flat,str_out);

	testValidity(mp_des,"detJ_SnapLin.txt", quA, quB, output);

	return 0;
}

if (false) 
{
	gsMultiPatch<> jigsaw = getJigSaw3D(4); 
	for(index_t p = 0; p < jigsaw.nBoxes(); p++)
	{
   		//changeSignOfDetJ(jigsaw.patch(p)); // Use for 3D small
		//jigsaw.patch(p).scale(0.01);
	}
	
 // Evaluate detJ in the corners
 
	gsMatrix<> mat(8,3);
	mat << 1,0,0,
		   0,1,0,
		   0,0,1,
		   1,1,0,
		   1,0,1,
		   0,1,1,
		   0,0,0,
		   1,1,1;
 
	std::vector< gsMatrix<> > outVec;


	for(index_t p = 0; p < jigsaw.nBoxes(); p++)
	{
	    gsExprAssembler<> A(1,1);
        gsMultiBasis<> dbasis(jigsaw);
	    A.setIntegrationElements(dbasis);
	
	    gsExprEvaluator<> ev(A);
	
	    typedef gsExprAssembler<>::geometryMap geometryMap;

	    geometryMap G = A.getMap(jigsaw);
		
	
		for (index_t i = 0; i < 8; i++)
		{
			gsVector<> pt = mat.row(i);
			gsAsConstMatrix<> tmp = ev.eval(jac(G).det(),pt,p);
			//gsInfo << "On patch " << p << " Corner (" << pt[0] << "," << pt[1] << "," << pt[2];
			//gsInfo << ") with value " << tmp(0,0) << "\n";

			if (tmp(0,0) > 0)
			{
				gsInfo << "\n -------------------- \n TMP(0,0) > 0 FOR PATCH " << p << "!!\n --------------- \n";
				for (index_t q = 0; q < jigsaw.nBoxes(); q++)
				{
					gsMatrix<> out(4,3);

					gsAsConstMatrix<> tmp2 = ev.eval(G,pt,p);
					gsDebugVar(tmp2);
					out.row(0) = tmp2.transpose();

					// Find the coefficient

					index_t ind = 1;

					for( index_t i = 0; i < jigsaw.patch(q).coefsSize(); i++)
					{
						gsMatrix<>::RowXpr cf = jigsaw.patch(q).coef(i);
						real_t norm = (cf-tmp2.transpose()).norm();

						if (norm == 0)
						{
							gsInfo << "FOUND COEF " << i << " at patch " << q << " of " << jigsaw.patch(q).coefsSize() << "\n";
							
							gsVector<index_t,3> stride; // Tensor strides
							gsVector<index_t,3> size; // Tensor strides

							// Prepare stride and size
							static_cast<gsTensorBSplineBasis<3>&> (jigsaw.patch(q).basis()).stride_cwise(stride);
							static_cast<gsTensorBSplineBasis<3>&> (jigsaw.patch(q).basis()).size_cwise(size);
							
							for (index_t s = 0; s < 3; s++)
							{
								gsInfo << "FOUND COEF " << i << " , stride = " << stride(s) << "\n";
								gsMatrix<> coefs = jigsaw.patch(q).coefs();

								index_t sign;
								if (q == 0) sign = -1;
								if (q == 1) sign = 1;

								if (s == 1)
									out.row(ind++) = (coefs.row(i - stride(s))).transpose();
								else if (s == 0)
									out.row(ind++) = (coefs.row(i + sign*stride(s))).transpose();
								else
									out.row(ind++) = (coefs.row(i + stride(s))).transpose();
							}
						outVec.push_back(out);
						}
					}
				}
			}

		}


		
	}


	index_t p = 1;
	index_t q = 0;

	std::vector< gsMatrix<> > outVec2;

	for (index_t i = 0; i < jigsaw.nInterfaces(); i++)
	{
        boundaryInterface interface = jigsaw.bInterface(i);

		bool term1 = (interface.first().patch == p && interface.second().patch == q);
		bool term2 =  (interface.first().patch == q && interface.second().patch == p);

		if (!term1 && !term2) continue;

		gsInfo << " ============= WE ARE HERE ============== \n";

		gsVector< index_t > vec(2);
		vec << p, q;

		for(index_t j = 0; j < vec.size(); j++)
		{
			patchSide ps;
			if (j == 0) ps = interface.first();
			if (j == 1) ps = interface.second();

        	gsVector<unsigned> boundaryDofs = jigsaw.basis(ps.patch).boundary(ps);

			gsMatrix<> out(boundaryDofs.size(),3*2);

			for( index_t k = 0; k < boundaryDofs.size(); k++)
			{
				index_t ii = boundaryDofs[k];
				gsMatrix<>::RowXpr cf = jigsaw.patch(ps.patch).coef(ii);
				out.row(k).segment(0,3) = cf;

				index_t dir = ps.direction();

				gsVector<index_t,3> stride; // Tensor strides

				// Prepare stride and size
				static_cast<gsTensorBSplineBasis<3>&> (jigsaw.patch(ps.patch).basis()).stride_cwise(stride);
				
				gsMatrix<> coefs = jigsaw.patch(ps.patch).coefs();
				
				index_t sign;
				if (j == 0) sign = 1;
				if (j == 1) sign = -1;

				out.row(k).segment(3,3) = (coefs.row(ii + sign*stride(dir)));
			}

			outVec.push_back(out);

		}

	}

	for(std::vector< gsMatrix<> >::iterator it = outVec.begin(); it != outVec.end(); ++it) 
	{
		gsDebugVar(*it);
	}

	return 0;

}

if (false)
{

	gsMultiPatch<> jigsaw = getJigSaw3D(4); 

	for(index_t p = 0; p < jigsaw.nBoxes(); p++)
	{
   		changeSignOfDetJ(jigsaw.patch(p)); // Use for 3D small
		jigsaw.patch(p).scale(0.01);
	}
    gsWriteParaview(jigsaw,"before",10000,true,true);

	testValidity(jigsaw,"detJ_before.txt", quA, quB, output);

    index_t p = 1;

    gsVector< index_t > cvec(2);
    cvec << 14, 35;

    gsVector< index_t > qvec(2);
    qvec << 14 + 1, 35 + 1;

    for (index_t j = 0; j < 2; j++)
    {
        index_t c = cvec[j];
        index_t q = qvec[j];

		gsMatrix<>::RowXpr P0 = jigsaw.patch(p).coef(c);
		gsMatrix<>::RowXpr Q1 = jigsaw.patch(p).coef(q);

        jigsaw.patch(p).coef(q) = -Q1 + 2*P0;
    }

    gsWriteParaview(jigsaw,"after",10000,true,true);
	testValidity(jigsaw,"detJ_after.txt", quA, quB, output);

    gsFileData<> fd;
    fd << jigsaw;
    std::ostringstream strs;
    strs << BASE_FOLDER << "/parametrizations/JigsawXML/water_passage_fixed"; 
    std::string str = strs.str();

	fd.save(str);

    return 0;

}


if (optParam) 
{
	index_t jigtype = jig;

	gsMultiPatch<> jigsaw, mp_init;
	if (dim == 2) {
		jigsaw = getJigSaw(jigtype);
    	mp_init = getInitGuess2d(jigsaw);
	} else if (dim == 3) {
		jigsaw = getJigSaw3D(jigtype); 
    	mp_init = getInitGuess3d(jigsaw);
	} else {
		GISMO_ERROR("dim should be 2 or 3!");
	}
	
	if (patch >= 0)
	{
		jigsaw = jigsaw.patch(patch);
		mp_init = mp_init.patch(patch);
	}


	gsMultiPatch<>::Ptr jigsaw_ptr = memory::make_shared_not_owned(&jigsaw);
	gsMultiPatch<>::Ptr mp_init_ptr = memory::make_shared_not_owned(&mp_init);

	if (jigtype == 3 || jigtype == 2 || jigtype == 4 || ((jigtype == 1 || jigtype == 5 || jigtype == 6) && dim == 3))
	{
		if( mp_init.targetDim() == 3)
		{
			for(index_t p = 0; p < jigsaw.nBoxes(); p++)
			{
    			changeSignOfDetJ(jigsaw.patch(p)); // Use for 3D small
	    		changeSignOfDetJ(mp_init.patch(p)); // Use for 3D small
			}
		}
	}

    std::string name;
    switch (jigtype) {
        case 1 : name = "/parametrizations/XML/2D/Jigsaw1"; break;
        case 2 : name = "/parametrizations/XML/2D/Jigsaw2"; break;
    }

    scaleMP(*jigsaw_ptr);

    gsInfo << "Writing to " << BASE_FOLDER + name << " \n \n ";

	gsFileData<> fd;
	fd << jigsaw; 
    fd.save(BASE_FOLDER + name);
    exit(0); 


	for (index_t r = 0; r < numRefine; r++)
	{
		jigsaw.uniformRefine();
		mp_init.uniformRefine();
	}

	//testValidity(jigsaw,"detJ_goal.txt", quA, quB, output);

//	gsDetJacConstraint dJC(jigsaw_ptr, true);
//	dJC.evalCon();
//
//	std::string st = BASE_FOLDER + output + "detJ_wp";
//	dJC.plotDetJ(st);
//
//	return 0;

	std::string nm = "jig";
	gsWriteParaview(jigsaw,nm,10000,true,true);

	nm = "init";
	gsWriteParaview(mp_init,nm,10000,true,true);

	testValidity(mp_init,"detJ_init.txt", quA, quB, output);

    // gsMultiPatch<>(*gsNurbsCreator<>::BSplineSquare(2));

    gsShapeOptLog slog1(output,true,false,false);
	gsShapeOptLog::Ptr slog1_ptr = memory::make_shared_not_owned(&slog1); 

    name = "init";
    slog1.plotInParaview(mp_init,name);


    // gsInfo << (spring.getFlat() - spring2.getFlat()).norm() << "\n";
    // gsInfo << (mp_init.patch(0).coefs() - jigsaw.patch(0).coefs()).norm() << "\n";
    // mat << spring.getTagged(), spring2.getTagged();

    // gsInfo << spring.n_tagged << "\n";
	gsOptParam optP(mp_init_ptr,jigsaw_ptr,slog1_ptr,param);
    gsOptParam::Ptr optP_ptr = memory::make_shared_not_owned( &optP);
    // gsMatrix<> v = optP.currentDesign()*0.99999;
    // optP.m_paramMethod->update();
    // gsInfo << optP.evalObj(v) << "\n";
    // gsInfo << optP.currentDesign() << "\n";
    // gsInfo << optP.numDesignVars() << "\n";
    // optP.solve();
	gsInfo << "EPS = " << eps << "\n";
    gsShapeOptWithReg optWR(mp_init_ptr,optP_ptr,numRefine,slog1_ptr,quA,quB,eps);
    // optWR.gradObj();
    optWR.solve();

 	name = "jigsaw";
    slog1.plotInParaview(jigsaw,name);

 	name = "mp_final";
    slog1.plotInParaview(mp_init,name);

	testValidity(mp_init,"detJ.txt", quA, quB, output);
	
 	// Test snaps
	gsWinslow win(mp_init_ptr,false,false,true,0); 
	gsWinslow win_jig(jigsaw_ptr,false,false,true,0); 
	gsVector<> tagged_jig = win_jig.getTagged();
	win.updateTagged(tagged_jig);

	std::stringstream stream_out;
	stream_out <<  output << "cps_snapped.txt";
	std::string str_out = stream_out.str();
	gsVector<> flat = win.getFlat();
	saveVec(flat,str_out);

	testValidity(mp_init,"detJ_snap.txt", quA, quB, output);
	
	// Snap bnd control points
	gsWinslow::Ptr win_ptr = memory::make_shared_not_owned(&win);
	gsAffineOptParamMethod lin(win_ptr,false);

	gsVector<> free = lin.getUpdate(tagged_jig);
	win.updateFreeAndTagged(free,tagged_jig);

	std::stringstream stream_out2;
	stream_out2 <<  output << "cps_SnapLin.txt";
	str_out = stream_out2.str();
	flat = win.getFlat();
	saveVec(flat,str_out);

	testValidity(mp_init,"detJ_SnapLin.txt", quA, quB, output);

	return 0;
}

if (startFromFile) 
{
    gsShapeOptLog slog1(output,true,false,false);
	gsShapeOptLog::Ptr slog1_ptr = memory::make_shared_not_owned(&slog1); 


	gsWinslow winslow(mp_ptr);
	gsVector<> flat_bef = winslow.getFlat();
	gsInfo << "\n ================ \n Load from" << BASE_FOLDER + startFile << " \n ===================== \n";
	winslow.updateFlat( loadVec(winslow.n_flat,BASE_FOLDER + startFile));
	gsVector<> flat_aft = winslow.getFlat();
	gsDebugVar((flat_bef - flat_aft).norm());
	exit(0);


    if (param == 5) // Use regularization
    {
        gsOptAntenna optA(mp_ptr,numRefine,slog1_ptr,6,quA,quB);
		gsOptAntenna::Ptr optA_ptr = memory::make_shared_not_owned(&optA);

        gsWinslow winslow2(mp_ptr,optA.mappers());
        winslow2.update();

        gsShapeOptWithReg optWR(mp_ptr,optA_ptr,numRefine,slog1_ptr,quA,quB,eps);
    //    optWR.setWinslowQuad(quA_optParam, quB_optParam);

        optWR.solve();
    } else if (param == 6) {
        gsOptAntenna optA(mp_ptr,numRefine,slog1_ptr,param,quA,quB,true);

	gsWinslow winslow2(mp_ptr,optA.mappers());
	winslow2.update();

        optA.setOptParamQuad(quA_optParam, quB_optParam);

        optA.runOptimization(maxiter);
    } else if (param == 0) {
        gsOptAntenna optA(mp_ptr,numRefine,slog1_ptr,param,quA,quB,true);
        optA.m_paramMethod->updateFlat(loadVec(optA.n_flat,BASE_FOLDER + output));
        gsInfo << optA.evalObj() << "\n";
    }

	return 0;


}

// test winslow derivatives
if (true) 
{
    mp = getJigSaw3D(4).patch(0);
    changeSignOfDetJ(mp.patch(0));
    mp_ptr = memory::make_shared_not_owned(&mp);
    gsWinslowPow winslow(mp_ptr,false);

    convergenceTestOfParaJacobianAll(winslow);

    //convergenceTestOfParaJacobian(winslow);

    exit(0);
}

// test harmonic derivatives
if (false) 
{
    gsHarmonic harmonic(mp_ptr,false);
    convergenceTestOfParaJacobian(harmonic);
	return 0;
}

// test of optAntenna derivatives
if (false) 
{
    gsWinslow winslow(mp_ptr,false);
	//gsInfo << winslow.getFlat() << "\n";
	gsShapeOptLog slog1(output,true,false,false);
	gsShapeOptLog::Ptr slog1_ptr = memory::make_shared_not_owned(&slog1); 

	gsOptAntenna optA(mp_ptr,numRefine,slog1_ptr,param,quA,quB,true);
	// gsInfo << "obj : " << optA.evalObj() << "\n";
	convergenceTestOfJacobian(optA);
    return 0;
}

if (startDes == 11) 
{
    gsMultiPatch<> jigsaw = getJigSaw3D(2); //
    mp = jigsaw;

    gsWinslowWithDeriv win(mp_ptr,false,false,true,0);

    std::stringstream stream_in;
    stream_in << BASE_FOLDER << output << "cps" << "_0_" << maxiter << ".txt";
    std::string str_in = stream_in.str();
    gsMatrix<> mat = win.getFlat();
    readFromTxt(str_in,mat);
    // gsInfo << mat;

    win.updateFlat(mat);

    gsDetJacConstraint dJC(mp_ptr,true);
    dJC.evalCon();
    gsShapeOptLog slog1(output,true,false,false);

    // gsOptAntenna optA(&mp,numRefine,&slog1,0,quA,quB);

    std::stringstream stream;
    stream << BASE_FOLDER << output << "detJ";
    std::string str = stream.str();
    dJC.plotDetJ(str);
    return 0;



}

if (false ) 
{
    // Note on what to try:
    //  How about if I run without constraints on detJ and afterwards refine m_mp
    //  where the negative coefficients are..?
    gsShapeOptLog slog1(output,true,false,false);
	gsShapeOptLog::Ptr slog1_ptr = memory::make_shared_not_owned(&slog1); 

    gsOptAntenna optA(patches_ptr,numRefine,slog1_ptr,param,quA,quB,true);
    gsWinslow winslow(patches_ptr,optA.mappers(),false,false,true,0);

    winslow.setQuad(quA,quB);

    winslow.updateFlat(loadVec(winslow.n_flat,BASE_FOLDER + output));

    gsDetJacConstraint dJC(patches_ptr);
    gsInfo << "mind = " << dJC.evalCon().minCoeff() << "\n";

    dJC.refineUntilPositive(maxiter);

    gsInfo << "mind = " << dJC.evalCon().minCoeff() << "\n";
    return 0;
}

// Test of refine and plotting det J surfaces
if (false) 
{
    // Note on what to try:
    //  How about if I run without constraints on detJ and afterwards refine m_mp
    //  where the negative coefficients are..?
    gsShapeOptLog slog1(output,true,false,false);
	gsShapeOptLog::Ptr slog1_ptr = memory::make_shared_not_owned(&slog1); 
    
    gsInfo << *patches_ptr << "\n";

    gsOptAntenna optA(patches_ptr,numRefine,slog1_ptr,param,quA,quB,true);
    gsWinslow winslow(patches_ptr,optA.mappers(),false,false,true,0);

    winslow.setQuad(quA,quB);

    winslow.updateFlat(loadVec(winslow.n_flat,BASE_FOLDER + output));
    // winslow.update();
    //

    gsMultiPatch<> sp(patches_ptr->patch(3));
    gsMultiPatch<>::Ptr sp_ptr = memory::make_shared_not_owned(&sp);

    gsDetJacConstraint dJC(sp_ptr);

    gsMultiPatch<> tmp(patches);
    gsInfo << "tmp size " << tmp.patch(0).coefsSize() << "\n";

    // winslow.update();

    std::string name = "winslow";
    gsWriteParaview(sp,name,10000,true);

    gsMultiPatch<> djsurf = dJC.getDetJSurface();

    name = "/../results/testOfDetJ/test2_detJSurface";
    gsWriteParaview(djsurf,BASE_FOLDER + name,10000,false,true);

    gsMultiPatch<> zero = dJC.getDetJSurface(true);

    name = "/../results/testOfDetJ/test2_zero";
    gsWriteParaview(zero,BASE_FOLDER + name + "_zero",10000,true,false);

    patches.uniformRefine();
    dJC.setup();
    gsInfo << "mind = " << dJC.evalCon().minCoeff() << "\n";

    djsurf = dJC.getDetJSurface();
    dJC.testSplineDetJ();

    name = "/../results/testOfDetJ/test2_detJSurface_uni1";
    gsWriteParaview(djsurf,BASE_FOLDER + name,10000,false,true);

    zero = dJC.getDetJSurface(true);
    gsWriteParaview(zero,BASE_FOLDER + name + "_zero",10000,true,false);

    patches.uniformRefine();
    dJC.setup();
    gsInfo << "mind = " << dJC.evalCon().minCoeff() << "\n";

    djsurf = dJC.getDetJSurface();
    dJC.testSplineDetJ();
    gsInfo << "\n";

    name = "/../results/testOfDetJ/test2_detJSurface_uni2";
    gsWriteParaview(djsurf,BASE_FOLDER + name,10000,false,true);
    zero = dJC.getDetJSurface(true);
    gsWriteParaview(zero,BASE_FOLDER + name + "_zero",10000,true,false);

    patches.uniformRefine();
    dJC.setup();
    gsInfo << "mind = " << dJC.evalCon().minCoeff() << "\n";

    djsurf = dJC.getDetJSurface();
    dJC.testSplineDetJ();
    gsInfo << "\n";

    name = "/../results/testOfDetJ/test2_detJSurface_uni4";
    gsWriteParaview(djsurf,BASE_FOLDER + name,10000,false,true);
    zero = dJC.getDetJSurface(true);
    gsWriteParaview(zero,BASE_FOLDER + name + "_zero",10000,true,false);

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
    gsWriteParaview(djsurf,BASE_FOLDER + name,10000,false,true);
    zero = dJC.getDetJSurface(true);
    mind = dJC.refineDetJSurfaceUntilPositive(1,zero);
    gsWriteParaview(zero,BASE_FOLDER + name + "_zero",10000,true,false);

    // Adapt Twice
    dJC.setup();

    djsurf = dJC.getDetJSurface();

    mind = dJC.refineDetJSurfaceUntilPositive(2,djsurf);
    // winslow.refineBasedOnDetJ(0);
    gsInfo << "mind = " << mind << "\n\n";

    name = "/../results/testOfDetJ/test2_detJSurface_adapt2";
    gsWriteParaview(djsurf,BASE_FOLDER + name,10000,false,true);
    zero = dJC.getDetJSurface(true);
    mind = dJC.refineDetJSurfaceUntilPositive(2,zero);
    gsWriteParaview(zero,BASE_FOLDER + name + "_zero",10000,true,false);

    // Adapt Thrice
    dJC.setup();

    djsurf = dJC.getDetJSurface();

    mind = dJC.refineDetJSurfaceUntilPositive(3,djsurf);
    // winslow.refineBasedOnDetJ(0);
    gsInfo << "mind = " << mind << "\n\n";

    name = "/../results/testOfDetJ/test2_detJSurface_adapt3";
    gsWriteParaview(djsurf,BASE_FOLDER + name,10000,false,true);
    zero = dJC.getDetJSurface(true);
    mind = dJC.refineDetJSurfaceUntilPositive(3,zero);
    gsWriteParaview(zero,BASE_FOLDER + name + "_zero",10000,true,false);

    return 0;
}

// Test of gsOptParam with reg
if (optParam) 
{
    // gsMultiPatch<> jigsaw = getJigSaw(startDes);
	index_t jigtype = 2;
    gsMultiPatch<> jigsaw = getJigSaw3D(jigtype); //
	gsMultiPatch<>::Ptr jigsaw_ptr = memory::make_shared_not_owned(&jigsaw);

    for( index_t r = 0; r < numRefine; r++){
        jigsaw.uniformRefine();
    }

    //changeSignOfDetJ(jigsaw.patch(0)); // Use for 3D small

    // Create startguess

    //-------- 3D JIGSAW ------
    index_t size    = jigsaw.patch(0).coefsSize();
	index_t n;
    gsVector<> vec;

	gsMatrix<> coefs_init;

	if (jigtype == 3)
	{
    	n = 5;
    	vec.setLinSpaced(n,-2.5,2.5);

		coefs_init.setZero(size,3);
    	for (index_t i = 0; i < n; i++){
    	    for(index_t j = 0; j < n; j++){
    	        for(index_t k = 0; k < n; k++){
    	            coefs_init(k + j*n + i*n*n,0) = vec[j];
    	            coefs_init(k + j*n + i*n*n,1) = vec[k];
    	            coefs_init(k + j*n + i*n*n,2) = vec[i];
    	        }
    	    }
    	}
	}
	else if (jigtype == 2)
	{
    	n = 6;
    	vec.setLinSpaced(n,-2.5,2.5);

		coefs_init.setZero(size,3);
    	
    	for (index_t i = 0; i < n; i++){
    	    for(index_t j = 0; j < n; j++){
    	        for(index_t k = 0; k < n; k++){
    	            coefs_init(k + j*n + i*n*n,0) = vec[j];
    	            coefs_init(k + j*n + i*n*n,1) = vec[k];
    	            coefs_init(k + j*n + i*n*n,2) = vec[i];
    	        }
    	    }
    	}
	}
	else
		GISMO_ERROR("not implemented for jigtype != 2 or 3");



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
	gsMultiPatch<>::Ptr mp_init_ptr = memory::make_shared_not_owned(&mp_init);

    mp_init.patch(0).setCoefs(coefs_init);
	if (jigtype == 3 || jigtype == 2)
	{
    	changeSignOfDetJ(jigsaw.patch(0)); // Use for 3D small
    	changeSignOfDetJ(mp_init.patch(0)); // Use for 3D small
	}

    // gsMultiPatch<>(*gsNurbsCreator<>::BSplineSquare(2));

    gsShapeOptLog slog1(output,true,false,false);
	gsShapeOptLog::Ptr slog1_ptr = memory::make_shared_not_owned(&slog1); 

    std::string name = "init";
    slog1.plotInParaview(mp_init,name);

    gsSpringMethod spring(mp_init_ptr);
    gsSpringMethod spring2(jigsaw_ptr);

    // gsInfo << (spring.getFlat() - spring2.getFlat()).norm() << "\n";
    // gsInfo << (mp_init.patch(0).coefs() - jigsaw.patch(0).coefs()).norm() << "\n";
    // mat << spring.getTagged(), spring2.getTagged();

    // gsInfo << spring.n_tagged << "\n";
	gsOptParam optP(mp_init_ptr,jigsaw_ptr,slog1_ptr,param);
    gsOptParam::Ptr optP_ptr = memory::make_shared_not_owned( &optP);
    // gsMatrix<> v = optP.currentDesign()*0.99999;
    // optP.m_paramMethod->update();
    // gsInfo << optP.evalObj(v) << "\n";
    // gsInfo << optP.currentDesign() << "\n";
    // gsInfo << optP.numDesignVars() << "\n";
    // optP.solve();
	gsInfo << "EPS = " << eps << "\n";
    gsShapeOptWithReg optWR(mp_init_ptr,optP_ptr,numRefine,slog1_ptr,quA,quB,eps);
    // optWR.gradObj();
    optWR.solve();

 	name = "jigsaw";
    slog1.plotInParaview(jigsaw,name);

 	name = "mp_final";
    slog1.plotInParaview(mp_init,name);

	testValidity(mp_init,"detJ.txt", quA, quB, output);

    return 0;

}



// Test Regularization with winslow for antenna problem
if (false) 
{
    gsInfo << output << "\n";
    gsShapeOptLog slog1(output,true,false,false);
	gsShapeOptLog::Ptr slog1_ptr = memory::make_shared_not_owned(&slog1); 
    gsOptAntenna optA(mp_ptr,numRefine,slog1_ptr,0,quA,quB);
	gsOptAntenna::Ptr optA_ptr = memory::make_shared_not_owned(&optA);

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
    //// gsInfo << "PROVE DETJ>0 : " << dJC.provePositivityOfDetJ_TP(nrs,5) << "\n";
    //
    // gsInfo << "min detJ in pts: " << winslow.minDetJInGaussPts(10) << "\n";

    gsShapeOptWithReg optWR(mp_ptr,optA_ptr,numRefine,slog1_ptr,quA,quB,eps);
    // gsInfo << optWR.evalObj();
    // gsInfo << optWR.gradObj();
    // optWR.solve();
    optWR.runOptimization(maxiter);

    // convergenceTestOfJacobian(optA);
    return 0;

}

/*
// calculate detJ constraints for different designs
if (false) {
    // gsMultiPatch<> jigsaw = getJigSaw(1);
    gsMultiPatch<> jigsaw = getJigSaw3D(2); //
	gsMultiPatch<>::Ptr jigsaw_ptr(&jigsaw);

    gsSpringMethod spring(jigsaw_ptr);
    gsDetJacConstraint dJC(jigsaw_ptr,true);

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
	gsMultiPatch<>::Ptr jigsaw_ptr;

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
	gsMultiPatch<>::Ptr mp_init_ptr(&mp_init);

    mp_init.patch(0).setCoefs(coefs_init);
    changeSignOfDetJ(mp_init.patch(0));

    // gsMultiPatch<>(*gsNurbsCreator<>::BSplineSquare(2));

    gsShapeOptLog slog1(output,true,false,false);
	gsShapeOptLog::Ptr slog1_ptr(&slog1); 

    std::string name = "init";
    slog1.plotInParaview(mp_init,name);

    gsSpringMethod spring(&mp_init);
    gsSpringMethod spring2(&jigsaw);

    // gsInfo << (spring.getFlat() - spring2.getFlat()).norm() << "\n";
    // gsInfo << (mp_init.patch(0).coefs() - jigsaw.patch(0).coefs()).norm() << "\n";
    // mat << spring.getTagged(), spring2.getTagged();

    // gsInfo << spring.n_tagged << "\n";
    gsOptParam optP(mp_init_ptr,jigsaw_ptr,slog1_ptr,param);
    // gsMatrix<> v = optP.currentDesign()*0.99999;
    // optP.m_paramMethod->update();
    // gsInfo << optP.evalObj(v) << "\n";
    // gsInfo << optP.currentDesign() << "\n";
    // gsInfo << optP.numDesignVars() << "\n";
    optP.solve();

    exit(0);




}

/*
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
*/

return 0;
}

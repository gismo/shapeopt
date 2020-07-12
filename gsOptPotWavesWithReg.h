#ifndef GSPOTWAVESWITHREG_H
#define GSPOTWAVESWITHREG_H
using namespace gismo;

#include "gsShapeOptWithReg.h"
#include "gs2NormConstraints.h"

class gsOptPotWavesWithReg: public gsShapeOptWithReg{
public:

    gsOptPotWavesWithReg(memory::shared_ptr<gsMultiPatch<>> mp, memory::shared_ptr<gsShapeOptProblem> sopt, index_t numRefine, memory::shared_ptr<gsShapeOptLog> slog, real_t quA, index_t quB, real_t eps, bool glueInterfaces = true, bool usePow = false);

    // Overload to get 2NormConstraints
    void evalCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const ;

    void jacobCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const ;

    // Method to setup constraints
    // Overloaded to use gs2NormConstraints
    void setupConstraints();

    void computeJacStructure(); // Overloaded from gsOptProblem

    // Method that is called between each optimization iteration.. Can be used to log
    //      or check stuff. If it returns false the optimization will be interrupted.
    bool intermediateCallback();

    // Method that is called between each optimization iteration.. Can be used to log
    //      or check stuff. If it returns false the optimization will be interrupted.
    // Perhaps we want other behaviour than in gsShapeOptWithReg? 
    //   then we should overload it
    
    //bool intermediateCallback();

    // Method to run optimization again and again while decreasing the regularization parameter
    // Maybe we want to overload it?
    void runOptimization();

public:
    gs2NormConstraints::Ptr m_constraint; // Pointer to constraints

    real_t m_a = 0.25;
    real_t m_b = 0.6;


};



# endif //GSPOTWAVESWITHREG_H

#ifndef GSPOTWAVESWITHREG2_H
#define GSPOTWAVESWITHREG2_H
using namespace gismo;

#include "gsShapeOptWithReg.h"
#include "gsIpOptSparseMatrix.h"

class gsOptPotWavesWithReg2nd: public gsShapeOptWithReg{
public:

    gsOptPotWavesWithReg2nd(gsMultiPatch<>::Ptr mp, gsMultiPatch<>::Ptr center, gsShapeOptProblem::Ptr sopt, index_t numRefine, gsShapeOptLog::Ptr slog, real_t quA, index_t quB, real_t eps, bool glueInterfaces = true, bool usePow = false);

    real_t evalObj() const ;

    gsVector<> gradObj() const ;

    void updateDesignVariables( gsVector<> u ) const;

    void evalCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const ;

    void jacobCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const ;

    void computeJacStructure(); // Overloaded from gsOptProblem

    void updateDesignBounds();

    void setupOptParameters();

    void setupConstraints();

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

    void setWinslowQuad(real_t quA, index_t quB);

public:
    real_t m_a = 0.25;
    real_t m_b = 0.6;

    gsMultiPatch<>::Ptr m_center;
    gsWinslow::Ptr m_center_winslow;

    gsMatrix<> m_constraint_matrix;
    mutable gsIpOptSparseMatrix m_J;
    index_t n_flat_center;

    index_t n_constraints;


};



# endif //GSPOTWAVESWITHREG2_H

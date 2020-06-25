#ifndef GSSHAPEOPTWITHREG_H
#define GSSHAPEOPTWITHREG_H
using namespace gismo;

#include "gsShapeOptProblem.h"
#include "gsShapeOptLog.h"
#include "gsWinslow.h"
#include <gsIpopt/gsOptProblem.h>

class gsShapeOptWithReg: public gsOptProblem<real_t>{
public:

    // Constructs from a multipatch, calling setupMappers, and setupOptParameters
    gsShapeOptWithReg(memory::shared_ptr<gsMultiPatch<>> mp, memory::shared_ptr<gsShapeOptProblem> sopt, index_t numRefine, memory::shared_ptr<gsShapeOptLog> slog, real_t quA, index_t quB, real_t eps, bool glueInterfaces = true, bool usePow = false);

    // Evaluation of the objective, using the design contained in m_mp
    real_t evalObj() const ;

    // Evaluation of the gradient wrt design variables (tagged cps)
    gsVector<> gradObj() const ;

    void evalCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const ;

    void jacobCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const ;

    // Mappers for antenna problem,
    // Interfaces are glued, outer boundary is fixed
    void setupMappers();

    // Method to set the optimization parameters such as design bounds, constraint bounds etc   .
    void setupOptParameters();

    // Update parametrization with u and call evalObj()
    real_t evalObj( const gsAsConstVector<real_t> & u ) const;

    // Update parametrization with u and call gradObj()
    void gradObj_into ( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result ) const;

    // Method that is called between each optimization iteration.. Can be used to log
    //      or check stuff. If it returns false the optimization will be interrupted.
    bool intermediateCallback();

    // Method to run optimization again and again while decreasing the regularization parameter
    void runOptimization(index_t maxiter);

    void setWinslowQuad(real_t quA, index_t quB);

    void runOptimizationUntilPosDetJ(index_t maxiter, real_t k, index_t maxRef);



public:
    real_t m_eps = 1; // Regularization parameter

    memory::shared_ptr<gsMultiPatch<>> m_mp;

    memory::shared_ptr<gsShapeOptLog> m_log;

    memory::shared_ptr<gsShapeOptProblem> m_opt; // Handles shape optimization objective and sensitivities
    memory::shared_ptr<gsWinslow> m_winslow; // Handles the Regularization term, and bookkeeping tasks

    std::vector< gsDofMapper > m_mappers; // Mapper for each coordinate

    gsVector<index_t> m_shift_free;
    gsVector<index_t> m_shift_all;
    gsVector<index_t> m_shift_tagged;

    index_t n_free;
    index_t n_flat;
    index_t n_tagged;
    index_t n_cps;

    index_t counter1 = 0; // Counts the number of iterations?
    index_t counter2 = 0; // Counts the number of iterations?

    bool m_glueInterfaces = true;

};



# endif //GSSHAPEOPTWITHREG_H

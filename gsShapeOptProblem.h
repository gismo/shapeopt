/** @file gsShapeOptProblem.h

@brief  Base class for shape optimization problems.

        It points to an instance of
        a gsParamMethod to handle parametrization, gsShapeOptLog for logging stuff
        during optimization, and gsDetJacConstraint to ensure bijectivity of the
        parametrization.

        To define a specific problem you need to overload the evalObj and gradObj
        method.

        There are two ways to solve the problem,
        1.  if you are using a simple parametrization method, e.g. spring method
            you can call solve();
        2.  If you use a parametrization method that has to reset, you can call
            runOptimization with the maximum number of outer loops as input.

This file is part of the G+Smo library.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Author(s): A. Limkilde, A. Mantzaflaris
*/

#ifndef GSSHAPEOPTPROBLEM_H
#define GSSHAPEOPTPROBLEM_H
using namespace gismo;

#include "gsDetJacConstraint.h"
#include "gsIpOptSparseMatrix.h"
#include "gsParamMethod.h"
#include "gsOptParamMethod.h"
#include "gsShapeOptLog.h"
#include <gsIpopt/gsOptProblem.h>

class gsShapeOptProblem: public gsOptProblem<real_t>{
public:

    // FIXIT: NOT IMPLEMENTED YET
    // Constructs from a pointer to a parametrization method, should glue the interfaces together and
    // eliminate boundaries.
    gsShapeOptProblem(memory::shared_ptr<gsMultiPatch<>> mp, memory::shared_ptr<gsShapeOptLog> slog, bool useDetJCons = false);

    // Constructs from a pointer to a parametrization method, should glue the interfaces together and
    // eliminate boundaries.
    //
    // Used constraints from the \a constraint pointer
    gsShapeOptProblem(memory::shared_ptr<gsMultiPatch<>> mp, memory::shared_ptr<gsShapeOptLog> slog, memory::shared_ptr<gsConstraint> constraint);

    // Constructs from list of mappers, one for each dimension
    //      The design variables for the shape optimization should be tagged in the
    //      mappers. The "inner controlpoints", that needs updation by a parametrization
    //      method should be free. The rest (fixed cps) should be eliminated..
    //
    // Uses gsDetJacConstraints
    gsShapeOptProblem(memory::shared_ptr<gsMultiPatch<>> mp, std::vector< gsDofMapper > mappers, memory::shared_ptr<gsShapeOptLog> slog);

    // Constructs from list of mappers, one for each dimension
    //      The design variables for the shape optimization should be tagged in the
    //      mappers. The "inner controlpoints", that needs updation by a parametrization
    //      method should be free. The rest (fixed cps) should be eliminated..
    //
    // Used constraints from the \a constraint pointer
    gsShapeOptProblem(memory::shared_ptr<gsMultiPatch<>> mp, std::vector< gsDofMapper > mappers, memory::shared_ptr<gsShapeOptLog> slog, memory::shared_ptr<gsConstraint> constraint);

    // Method to set the optimization parameters such as design bounds, constraint bounds etc   .
    // these are automaticly generated from mappers and gsDetJacConstraint m_dJC
    void setupOptParameters();
    void setupConstraints();

    // Method overloaded from gsOptProblem<>
    // Uses gsDetJacConstraint as default constraint
    void computeJacStructure();

    // Evaluation of the objective, using the design variables contained in m_mp
    virtual real_t evalObj() const = 0;

    // Evaluation of the gradient wrt design variables (tagged cps)
    // Use for the optimization
    virtual gsVector<> gradObj() const = 0;

    // Gradient with respect to all control points
    virtual gsVector<> gradAll() const { GISMO_NO_IMPLEMENTATION; };

    virtual gsDofMapper mapper_grad() const { GISMO_NO_IMPLEMENTATION; };

    // Mappers for the specific problem, should be overloaded in children
    virtual void setupMappers() = 0;

    // Evaluation of constraits, uses gsDetJacConstraint as default
    gsVector<> evalCon() const ;

    // Update parametrization with u and call evalObj()
    real_t evalObj( const gsAsConstVector<real_t> & u ) const;

    // Update parametrization with u and call gradObj()
    void gradObj_into ( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result ) const;

    // Update parametrization with u and call evalCon()
    void evalCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const ;

    // Evaluates the jacobian of the constraints. Default behaviour calls gsDetJacConstraint
    void jacobCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const ;

    // Get the tagged Cps. Uses the method from gsParamMethod
    gsVector<> getDesignVariables() const;

    // Update the tagged DoFs by the parametrization methods in gsParamMethod
    // Returns false if the update failed
    bool updateDesignVariables(gsVector<> u) const;

    // Computes the jacobian of the method updateDesignVariables
    //      m_paramMethod gives the jacobian wrt the free Cps,
    //      this method concatenates it with jacobian wrt the tagged cps
    //      but this is simply the identity, since the tagged are set to u.
    gsMatrix<> jacobDesignUpdate() const;

    // Computes the jacobian of the gsDetJacConstraints. The computation is done
    //      by m_dJC. This method maps it to the free and tagged dofs.
    gsIpOptSparseMatrix jacobDetJac() const;

    // Maps a vector from mapper_in indexing, to m_mappers.
    // E.g. used to map gradients with respect to all cps to be respect to free
    // and tagged DoFs...
    //
    // It stacks the dofs such that all coordinates of the free comes first, and
    // then all the coordinates of the tagged. I.e.
    // [x_free, y_free, z_free, x_tagged, y_tagged, z_tagged]'
    //
    // FIXIT: maybe find a more elegant way to handle all of this..
    //      Maybe always assemble wrt. specific mapper, and hold a permutation
    //      in a matrix, that we only need to multiply to map...
    gsMatrix<> mapGradient(gsDofMapper mapper_in, gsMatrix<> mat_in) const;

    // Method to run consecutive solves, resetting the parametrization method each time.
    void runOptimization(index_t maxiter, bool uniRef = false);

    void runOptimization_aggregatedConstraints(index_t maxiter);

    // Method that is called between each optimization iteration.. Can be used to log
    //      or check stuff. If it returns false the optimization will be interrupted.
    bool intermediateCallback();

    // Acces to controlpoints using the implementations in m_paramMethod
    gsVector<> getFree() const { return m_paramMethod->getFree(); };
    gsVector<> getTagged() const { return m_paramMethod->getTagged(); };
    gsVector<> getControlPoints() const { return m_paramMethod->getControlPoints(); };
    gsVector<> getFlat() const { return m_paramMethod->getFlat(); };

    // Accessors
    std::vector< gsDofMapper > mappers() { return m_mappers; };

    // Method to set the bounds on the design variables,
    // should be overloaded in child classes
    virtual void setupDesignBounds() = 0;

    // Method to print optimization parameters
    void print();

    // Method to set number of quadrature points for parametrization method
    void setOptParamQuad(real_t quA, index_t quB);

    virtual real_t getQuA(){ GISMO_NO_IMPLEMENTATION };
    virtual real_t getQuB(){ GISMO_NO_IMPLEMENTATION };

public:
    mutable memory::shared_ptr<gsMultiPatch<>> m_mp;
    mutable memory::shared_ptr<gsConstraint> m_dJC;
    mutable memory::shared_ptr<gsParamMethod> m_paramMethod;

    memory::shared_ptr<gsShapeOptLog> m_log;

    std::vector< gsDofMapper > m_mappers; // Mapper for each coordinate

    index_t n_free;
    index_t n_flat;
    index_t n_tagged;
    index_t n_cps;

    index_t counter1 = 0; // Counts the number of iterations?
    index_t counter2 = 0; // Counts the number of times the reference is updated.

    bool m_useOwnCons = false;
    bool m_useDetJCons;
};



# endif //GSSHAPEOPTPROBLEM_H

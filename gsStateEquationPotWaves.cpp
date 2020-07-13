#include <gismo.h>
#include "gsStateEquationPotWaves.h"
#include "gsMyExpressions.h"
#include <algorithm>
#include <regex>
using namespace gismo;

// Class to represent complex expressions
class complexExpr
{
public: 

// Constructors and accesors
    complexExpr(){};
    complexExpr(std::string re, std::string im): real(re), imag(im){};

    gsFunctionExpr<>::uPtr getFunRe( std::map< std::string, real_t > const_map) const
    {
        gsFunctionExpr<>::uPtr tmp = getFun(real, const_map);
        return tmp;
    }

    gsFunctionExpr<>::uPtr getFunIm( std::map< std::string, real_t > const_map ) const
    {
        gsFunctionExpr<>::uPtr tmp = getFun(imag, const_map);
        return tmp;
    }

    std::string Re() const { return real; };
    std::string Im() const { return imag; };

    std::string ReWithVals(std::map< std::string, real_t> const_map) const
    { 
        return insertValues(real, const_map); 
    };

    std::string ImWithVals(std::map< std::string, real_t> const_map)  const
    {
       return insertValues(imag, const_map); 
    };


// Helpers
    gsFunctionExpr<>::uPtr getFun(std::string str, std::map< std::string, real_t > const_map) const
    {
        std::string str_with_val = insertValues(str,const_map);
    
        gsFunctionExpr<>::uPtr tmp = memory::make_unique(new gsFunctionExpr<>(str_with_val,3));
        return tmp;
    
    }

    std::string insertValues(std::string str, std::map< std::string, real_t > map) const
    {
        std::ostringstream os;       // assemble output string
    
        // Stolen of the interwebs. replace ' ,()' with a sequence of the delimiters that you want to use
        std::regex words_regex("([ ,() - +*/^]|[^ ,() - +*/^]+)"); 
        auto words_begin = std::sregex_iterator(str.begin(), str.end(), words_regex);
        auto words_end = std::sregex_iterator();
    
        for (std::sregex_iterator i = words_begin; i != words_end; ++i)
        {
            std::map< std::string, real_t >::iterator it;
            it = map.find(i->str());
            //gsInfo << i->str() << "||";
    
            if ( it != map.end() ) 
                os << it->second;  
            else 
                os << i->str();
        }
    
        return os.str();
    }

    std::string pad(std::string str)
    {
        return std::string("( " ) + str + std::string( " )");
    }

// Operators

    complexExpr & operator = (const complexExpr & rhs)
    {
        if (this != &rhs)
        { 
            real = rhs.real;
            imag = rhs.imag;

        }
        return *this;
    }

    complexExpr operator * (complexExpr rhs)
    {
        complexExpr res;
        // real*rhs.real - imag*rhs.imag 
        res.real = pad(real) + mult + pad(rhs.real) + minus + pad(imag) + mult + pad(rhs.imag);

        // real*rhs.imag + imag*rhs.real
        res.imag = pad(real) + mult + pad(rhs.imag) + plus + pad(imag) + mult + pad(rhs.real);

        return res;
    }

    complexExpr operator * (std::string factor)
    {
        complexExpr res;
        res.real = pad(real) + mult + factor;

        res.imag = pad(imag) + mult + factor;

        return res;
    }

    complexExpr operator + (complexExpr rhs)
    {
        complexExpr res;
        // real + rhs.real
        res.real = real + plus + rhs.real;

        // imag + rhs.imag
        res.imag = imag + plus + rhs.imag;

        return res;
    }

    complexExpr operator - (complexExpr rhs)
    {
        complexExpr res;
        // real - rhs.real 
        res.real = real + minus + pad(rhs.real);

        // imag - rhs.imag
        res.imag = imag + minus + pad(rhs.imag);

        return res;
    }

    complexExpr operator - ()
    {
        complexExpr res;
        // real - rhs.real 
        res.real = minus + pad(real);

        // imag - rhs.imag
        res.imag = minus + pad(imag);

        return res;
    }

    complexExpr inv()
    {
        complexExpr res;
        
        // compute squared norm
        std::string sqn = sqNorm_asStr();

        // real / (real^2 + imag^2)
        res.real = real + div + pad(sqn);

        // imag / (real^2 + imag^2)
        res.imag = minus + imag + div + pad(sqn);

        return res;
    }

    complexExpr norm()
    {
        complexExpr res;
        
        // compute squared norm
        std::string sqn = sqNorm_asStr();

        // sqrt(real^2 + imag^2)
        res.real = "sqrt" + pad(sqn);

        // 0
        res.imag = "0";

        return res;

    }

    // return real^2 + imag^2
    std::string sqNorm_asStr()
    {
        // (real^2 + imag^2)
        return std::string("( ") + real + std::string(" )^2 ") + plus + std::string("( ") + imag + std::string(" )^2 ");
    }

    // return real^2 + imag^2
    complexExpr sqNorm()
    {
        complexExpr res;
        
        // compute squared norm
        std::string sqn = sqNorm_asStr();

        // real / (real^2 + imag^2)
        res.real = sqn;

        // imag / (real^2 + imag^2)
        res.imag = "0";

        return res;
        // (real^2 + imag^2)
    }

    complexExpr adj()
    {
        complexExpr res;
        
        res.real = real;

        res.imag = minus + pad(imag) ;

        return res;
    }

    complexExpr exp()
    {
        complexExpr res;
        
        res.real = "exp" + pad(real)  + " * cos " + pad(imag) ;

        res.imag = "exp" + pad(real)  + " * sin " + pad(imag) ;

        return res;
    }
    

public: 
    std::string real;
    std::string imag;

    std::string mult    = " * ";
    std::string minus   = " - ";
    std::string div     = " / ";
    std::string plus    = " + ";

};

// ============================ //
//       Helper functions       //
// ============================ //

std::ostream &operator<<(std::ostream &os, complexExpr x)
{ os << x.real << " + i " << x.pad(x.imag) ; return os;}

gsFunctionExpr<>::uPtr getSymmVectorFunRe(std::map< std::string, real_t> & map, 
        const complexExpr & xx, 
        const complexExpr & yx, 
        const complexExpr & zx,
        const complexExpr & yy,
        const complexExpr & zy,
        const complexExpr & zz)
{
    std::string xxs = xx.ReWithVals(map);
    std::string yxs = yx.ReWithVals(map);
    std::string zxs = zx.ReWithVals(map);

    std::string yys = yy.ReWithVals(map);
    std::string zys = zy.ReWithVals(map);

    std::string zzs = zz.ReWithVals(map);

    gsFunctionExpr<>::uPtr ptr = memory::make_unique( new gsFunctionExpr<>(xxs,yxs,zxs,yxs,yys,zys,zxs,zys,zzs,3) ); 

    return ptr;
}

gsFunctionExpr<>::uPtr getSymmVectorFunIm(std::map< std::string, real_t> & map, 
        const complexExpr & xx, 
        const complexExpr & yx, 
        const complexExpr & zx,
        const complexExpr & yy,
        const complexExpr & zy,
        const complexExpr & zz)
{
    std::string xxs = xx.ImWithVals(map);
    std::string yxs = yx.ImWithVals(map);
    std::string zxs = zx.ImWithVals(map);

    std::string yys = yy.ImWithVals(map);
    std::string zys = zy.ImWithVals(map);

    std::string zzs = zz.ImWithVals(map);

    gsFunctionExpr<>::uPtr ptr = memory::make_unique( new gsFunctionExpr<>(xxs,yxs,zxs,yxs,yys,zys,zxs,zys,zzs,3) ); 

    return ptr;
}

gsFunctionExpr<>::uPtr getVectorFunRe(std::map< std::string, real_t> & map, const complexExpr & x, const complexExpr & y, const complexExpr & z)
{
    std::string xs = x.ReWithVals(map);
    std::string ys = y.ReWithVals(map);
    std::string zs = z.ReWithVals(map);
    gsFunctionExpr<>::uPtr ptr = memory::make_unique( new gsFunctionExpr<>(xs,ys,zs,3) ); 

    return ptr;
}

gsFunctionExpr<>::uPtr getVectorFunIm(std::map< std::string, real_t> & map, const complexExpr & x, const complexExpr & y, const complexExpr & z)
{
    std::string xs = x.ImWithVals(map);
    std::string ys = y.ImWithVals(map);
    std::string zs = z.ImWithVals(map);
    gsFunctionExpr<>::uPtr ptr = memory::make_unique( new gsFunctionExpr<>(xs,ys,zs,3) ); 

    return ptr;
}

bool checkForZero(gsVector<unsigned> boundaryDofs, gsMatrix<> coefs, index_t coord = 2)
{
    for (index_t j = 0; j < boundaryDofs.size(); j++)
    {
        index_t k = boundaryDofs[j];
        real_t cps_z = coefs(k,coord); // 2 indicates z coordinate

        if (cps_z != 0)
            return false;
    }

    return true; // Only gets to this point if cps_z = 0 for all bnd points

}


bool checkForPML(gsVector<unsigned> boundaryDofs, gsMatrix<> coefs, real_t lbx, real_t lby, real_t lbz)
{
    for (index_t j = 0; j < boundaryDofs.size(); j++)
    {
        index_t k = boundaryDofs[j];
        real_t cps_x = coefs(k,0); // 0 indicates x coordinate
        real_t cps_y = coefs(k,1); // 1 indicates y coordinate
        real_t cps_z = coefs(k,2); // 2 indicates z coordinate

        bool xbool = (cps_x != lbx) && (cps_x != -lbx);
        bool ybool = (cps_y != lby) && (cps_y != 0);
        bool zbool = (cps_z != 0)   && (cps_z != -lbz);

        if (xbool || ybool || zbool)
            return false;
    }

    return true; 

}

bool checkIfInPML(gsVector<unsigned> boundaryDofs, gsMatrix<> coefs, real_t lx, real_t ly, real_t lz)
{
    for (index_t j = 0; j < boundaryDofs.size(); j++)
    {
        index_t k = boundaryDofs[j];
        real_t cps_x = coefs(k,0); // 0 indicates x coordinate
        real_t cps_y = coefs(k,1); // 1 indicates y coordinate
        real_t cps_z = coefs(k,2); // 2 indicates z coordinate

        bool xbool = (cps_x < -lx) || (cps_x > lx);
        bool ybool = (cps_y < -ly) || (cps_y > ly);
        bool zbool = (cps_z < -lz) || (cps_z > lz);
        
        //gsInfo << "x,y,z = " << xbool << ", "<< ybool << ", "<< zbool << "\n ";

        if (xbool || ybool || zbool)
            return true;
    }

    return false; 

}

gsTensorBSpline<3, real_t> getBox(gsMatrix<> &coefs)
{
    index_t degree  = 1;
    index_t n       = 2;
    index_t dim     = 3;

	gsKnotVector<> kv(0, 1, n - degree - 1, degree + 1);
    
	// 2. construction of a basis
	gsTensorBSplineBasis<3, real_t> basis(kv, kv, kv);

	// 3. construction of a coefficients
	gsTensorBSpline<3, real_t>  tbsout(basis, coefs);

    // To insert more knots use this
    /*
	for (index_t d = 0; d < dim; d++)	
	{
		tbsout.degreeElevate(tbs.degree(d)-degree,d);

		gsKnotVector<> kv = tbs.knots(d);

		for(gsKnotVector<>::iterator it = kv.begin() + tbs.degree(d) + 1; it != kv.end() - tbs.degree(d) - 1; it++)
		{
			tbsout.insertKnot(*it,d);
		}


	}
    */

    return tbsout;

}

gsMatrix<> getBoxHelper(real_t xmin,real_t xmax, real_t ymin,real_t ymax, real_t zmin,real_t zmax)
{
    gsMatrix<> out(8,3);

    out <<  xmin   ,ymin  ,zmin,
            xmax   ,ymin  ,zmin,
            xmin   ,ymax  ,zmin,
            xmax   ,ymax  ,zmin,
            xmin   ,ymin  ,zmax,
            xmax   ,ymin  ,zmax,
            xmin   ,ymax  ,zmax,
            xmax   ,ymax  ,zmax;

    return out;
}

index_t getSignOfDetJ(gsMultiPatch<> mp)
{
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::solution    solution;

    gsExprAssembler<> A(1,1);
    gsExprEvaluator<> ev(A);

    // Elements used for numerical integration
    gsMultiBasis<> dbasis(mp);
    A.setIntegrationElements(dbasis);

    geometryMap G = A.getMap(mp);

    real_t a = ev.integral(jac(G).det()); 

    return (0 < a) - (a < 0);

}

// ============================ //
//       Member functions       //
// ============================ //

gsStateEquationPotWaves::gsStateEquationPotWaves(index_t numRefine, real_t Lx, real_t Ly, real_t Lz, real_t lx, real_t ly, real_t lz): 
    zero("0.0",3), 
    m_numRefine(numRefine)
{
    // Set PML sizes
    pml_Lx = Lx;
    pml_Ly = Ly;
    pml_Lz = Lz;
    pml_lx = lx;
    pml_ly = ly;
    pml_lz = lz;
     
    // Get initial domain
    m_mp = getInitialDomain();

    // Setup dbasis 
    setup();

    // Setup the rest
    constructor();

}

gsStateEquationPotWaves::gsStateEquationPotWaves(gsMultiPatch<>::Ptr mp_ptr, index_t numRefine): 
    zero("0.0",3), 
    m_numRefine(numRefine),
    gsStateEquation(mp_ptr,numRefine)
{
    constructor();
}

gsStateEquationPotWaves::gsStateEquationPotWaves(index_t numRefine): 
    zero("0.0",3), 
    m_numRefine(numRefine)
{
    // Get initial domain
    m_mp = getInitialDomain();

    // Setup dbasis 
    setup();

    // Setup the rest
    constructor();
}

void gsStateEquationPotWaves::constructor()
{

    m_isBndGamma_s_vec.setZero(m_mp->nBoundary());

    complexExpr u_expr("0.0","0.0");
    pde_u_dir_re = u_expr.getFunRe(const_map);
    pde_u_dir_im = u_expr.getFunIm(const_map);

    lapl_uEx = u_expr.getFunRe(const_map);


    // Mark boundaries Gamma_f and Gamma_s
    markBoundaries();

    wave_omega = sqrt(wave_K*wave_g);        

    const_map["C"] = pml_C;
    const_map["n"] = pml_n;
    const_map["Lx"] = pml_Lx;
    const_map["Ly"] = pml_Ly;
    const_map["Lz"] = pml_Lz;
    const_map["lx"] = pml_lx;
    const_map["ly"] = pml_ly;
    const_map["lz"] = pml_lz;

    const_map["K"] = wave_K;
    const_map["K2"] = wave_K*wave_K;
    const_map["g"] = wave_g;
    const_map["omega"] = wave_omega;

    gsInfo << "\n Wave : \n";

    complexExpr wave_u_I("  g/ omega * exp(K*z) * cos(K*x)","- g/ omega * exp(K*z) * sin(K*x)");

    wave_u_I_re = wave_u_I.getFunRe(const_map );
    wave_u_I_im = wave_u_I.getFunIm(const_map );

    gsInfo << "wave_u_I_re: " << *wave_u_I_re << "\n";
    gsInfo << "wave_u_I_im: " << *wave_u_I_im << "\n";

    complexExpr wave_du_I_dx("- g / omega * K * exp(K*z) * sin(K*x)","- g / omega * K * exp(K*z) * cos(K*x)");
    complexExpr wave_du_I_dy("0","0");
    complexExpr wave_du_I_dz("  g / omega * K * exp(K*z) * cos(K*x)", "- g / omega * K * exp(K*z) * sin(K*x)");

    pde_F_re = getVectorFunRe(const_map, -wave_du_I_dx, -wave_du_I_dy, -wave_du_I_dz);
    pde_F_im = getVectorFunIm(const_map, -wave_du_I_dx, -wave_du_I_dy, -wave_du_I_dz);

    // Second derivatives of u_I
    complexExpr xx("- g / omega * K2 * exp(K*z) * cos(K*x)","g / omega * K2 * exp(K*z) * sin(K*x)");
    complexExpr yx("0","0");
    complexExpr zx("- g / omega * K2 * exp(K*z) * sin(K*x)", "- g / omega * K2 * exp(K*z) * cos(K*x)");

    complexExpr yy("0","0");
    complexExpr zy("0","0");

    complexExpr zz("  g / omega * K2 * exp(K*z) * cos(K*x)", "- g / omega * K2 * exp(K*z) * sin(K*x)");

    pde_dF_re = getSymmVectorFunRe(const_map, -xx, -yx, -zx, -yy, -zy, -zz);
    pde_dF_im = getSymmVectorFunIm(const_map, -xx, -yx, -zx, -yy, -zy, -zz);

    // Streching functions for PML
    gsInfo << "\n Pml : \n";

    // the factor (-1)^(n-1) makes sure that the sign is (-1) for even n and 1 for odd n
    //    since we need f'(x) > 0 when xtilde = x + if(x).
    std::string f_x_expr =  
        "switch                                                             "
        "{                                                                  "
        "   case (x < - lx)  : C/(omega * (Lx - lx)^n)*(- x - lx)^ ( n );   "
        "   case (x >   lx)  : - C/(omega * (Lx - lx)^n)*(x - lx)^ ( n );     "
        "   default          : 0;                                           "
        "}                                                                  ";

    std::string f_y_expr =  
        "switch                                                             "
        "{                                                                  "
        "   case (y < - ly)  : C/(omega * (Ly - ly)^n)*(- y - ly)^ ( n );   "
        "   case (y >   ly)  : - C/(omega * (Ly - ly)^n)*(y - ly)^ ( n );     "
        "   default          : 0;                                           "
        "}                                                                  ";

    complexExpr strech_x("x",f_x_expr);
    complexExpr strech_y("y",f_y_expr);

    std::string dsx_expr_im = 
        "switch                                                             "
        "{                                                                  "
        "   case (x < - lx)  : - n*C/(omega * (Lx - lx)^n)*( - x - lx)^(n - 1);   "
        "   case (x >   lx)  : - n*C/(omega * (Lx - lx)^n)*(x - lx)^(n - 1);     "
        "   default          : 0;                                           "
        "}                                                                  ";
        
    std::string dsy_expr_im = 
        "switch                                                             "
        "{                                                                  "
        "   case (y >   ly)  : - n*C/(omega * (Ly - ly)^n)*(y - ly)^(n - 1);     "
        "   default          : 0;                                           "
        "}                                                                  ";

    std::string dsz_expr_im = 
        "switch                                                             "
        "{                                                                  "
        "   case (z < - lz)  :  - n*C/(omega * (Lz - lz)^n)*( - z - lz)^(n - 1);   "
        "   default          : 0;                                           "
        "}                                                                  ";
    
    complexExpr dstrech_x("1",dsx_expr_im);
    complexExpr dstrech_y("1",dsy_expr_im);
    complexExpr dstrech_z("1",dsz_expr_im);

    complexExpr detJ = dstrech_x*dstrech_y*dstrech_z;
    complexExpr detJlen = detJ.norm();

    strech_detJ_re = detJ.getFunRe(const_map);
    strech_detJ_im = detJ.getFunIm(const_map);

    complexExpr dstrech_x_inv = dstrech_x.inv();
    complexExpr dstrech_y_inv = dstrech_y.inv();
    complexExpr dstrech_z_inv = dstrech_z.inv();

    // detJ * Jm1 * Jm1
    //complexExpr xcomp = detJ*dstrech_x_inv.sqNorm();//*dstrech_x_inv.adj();
    //complexExpr ycomp = detJ*dstrech_y_inv.sqNorm();//*dstrech_y_inv.adj();
    //complexExpr zcomp = detJ*dstrech_z_inv.sqNorm();//*dstrech_z_inv.adj();

    complexExpr xcomp = detJ*dstrech_x_inv*dstrech_x_inv;//*dstrech_x_inv.adj();
    complexExpr ycomp = detJ*dstrech_y_inv*dstrech_y_inv;//*dstrech_y_inv.adj();
    complexExpr zcomp = detJ*dstrech_z_inv*dstrech_z_inv;//*dstrech_z_inv.adj();

    strech_detJ_Jm1_Jm1_asVec_re = getVectorFunRe(const_map, xcomp, ycomp, zcomp);
    strech_detJ_Jm1_Jm1_asVec_im = getVectorFunIm(const_map, xcomp, ycomp, zcomp);

    //complexExpr strech_x("x1","y1");
    //complexExpr strech_y("x1","y2");

    std::string factor = "( - 1 / (2 * bell_s2) )";

    complexExpr expon = strech_x*strech_x*factor + strech_y*strech_y*factor;

    // Bell curve
    const_map["bell_s2"] = bell_sigma*bell_sigma;
    const_map["bell_A"] = bell_amplitude;
    //complexExpr bell_cE( "- bell_A * exp(- (x*x + y*y)/(2*bell_s2) )", "0");

    complexExpr bell_cE = expon.exp();
    //gsDebugVar(bell_cE.real);
    //gsDebugVar(bell_cE.imag);

    complexExpr bellTerm = detJ*bell_cE;
   
    bell_term_re = bell_cE.getFunRe(const_map);
    bell_term_im = bell_cE.getFunIm(const_map);


}

void gsStateEquationPotWaves::setup()
{
    for (index_t p = 0; p < m_mp->nBoxes(); p++)
    {
        if (isPatchInDomain(p))
        {
            gsInfo << "Added patch " << p << " to domain\n";
            m_mp_domain.addPatch( m_mp->patch(p) );
        }
    }
    m_mp_domain.computeTopology();
    gsDebugVar(m_mp_domain);
    gsDebugVar(m_mp_domain.domainDim());
    gsDebugVar(m_mp_domain.targetDim());

    gsMultiBasis<> bas(*m_mp);
    dbasis = bas;

    gsMultiBasis<> basD(m_mp_domain);
    dbasis_domain = basD;

    while(dbasis.maxCwiseDegree() < degree)
    {
        dbasis.degreeIncrease();
        dbasis_domain.degreeIncrease();
    }

    // gsInfo << "OBS CHECK DISCRETIZATION OF PDE!!\n";
    for (index_t i = 0; i < m_numRefine; i++){
        dbasis.uniformRefine();
        dbasis_domain.uniformRefine();
    }

    gsDebugVar(dbasis.domainDim());
    gsDebugVar(dbasis.targetDim());
    gsDebugVar(dbasis_domain.domainDim());
    gsDebugVar(dbasis_domain.targetDim());
    

}

void gsStateEquationPotWaves::markBoundaries()
{

    for (index_t i = 0; i < m_mp->nBoundary(); i++)
    {
        patchSide ps = m_mp->boundaries()[i];
        gsVector<unsigned> boundaryDofs = m_mp->basis(ps.patch).boundary(ps);

        // Find Gamma_f
        // Goal is to find the boundaries with z=0
        if (checkForZero(boundaryDofs,m_mp->patch(ps.patch).coefs()))
        {
            gsInfo << "GAMMA_F: Add condition at patch " << ps.patch << " bnd " << ps.index() <<"\n";
            bcInfo_Gamma_f.addCondition(ps.patch, ps.index(), condition_type::neumann, zero);

        }

        gsMatrix<> coefsTranslated = m_mp->patch(ps.patch).coefs();

        for (index_t i = 0; i < coefsTranslated.rows(); i++)
        {
            coefsTranslated(i,2) += pml_lz;
        }

        if (checkForZero(boundaryDofs, coefsTranslated))
        {
            //gsInfo << "GAMMA_B: Add condition at patch " << ps.patch << " bnd " << ps.index() <<"\n";
            bcInfo_Gamma_b.addCondition(ps.patch, ps.index(), condition_type::neumann, zero);
        }

        // Find Gamma_s
        // Goal is to find the boundaries inside the PML domain (ie. not on the PML layers)
        bool cond1 = !checkIfInPML(boundaryDofs,m_mp->patch(ps.patch).coefs(),pml_lx,pml_ly,pml_lz);
        bool cond2 = !checkForZero(boundaryDofs,m_mp->patch(ps.patch).coefs(),2);
        bool cond3 = !checkForZero(boundaryDofs,m_mp->patch(ps.patch).coefs(),1);
        if (cond1 and cond2 and cond3)
        {
            //gsInfo << "GAMMA_S: Add condition at patch " << ps.patch << " bnd " << ps.index() << "\n";
            bcInfo_Gamma_s.addCondition(ps.patch, ps.index(), condition_type::neumann, zero);
            m_isBndGamma_s_vec[i] = 1;
        }
    }

    gsInfo << "bcInfo_Gamma_f : \n" << bcInfo_Gamma_f << "\n";
    gsInfo << "bcInfo_Gamma_b \n" << bcInfo_Gamma_b << "\n";
    gsInfo << "bcInfo_Gamma_s : \n" << bcInfo_Gamma_s << "\n";

}

void gsStateEquationPotWaves::markCenterBessel()
{
    for (index_t i = 0; i < m_mp->nBoundary(); i++)
    {
        patchSide ps = m_mp->boundaries()[i];
        gsVector<unsigned> boundaryDofs = m_mp->basis(ps.patch).boundary(ps);

        // Find Gamma_s
        // Goal is to find the boundaries inside the PML domain (ie. not on the PML layers)
        bool cond1 = !checkIfInPML(boundaryDofs,m_mp->patch(ps.patch).coefs(),pml_lx,pml_ly,pml_lz);
        bool cond2 = !checkForZero(boundaryDofs,m_mp->patch(ps.patch).coefs(),2);
        bool cond3 = !checkForZero(boundaryDofs,m_mp->patch(ps.patch).coefs(),1);
        if (cond1 and cond2 and cond3)
        {
            //gsInfo << "GAMMA_S: Add condition at patch " << ps.patch << " bnd " << ps.index() << "\n";
            bcInfo_Dirichlet_re.addCondition(ps.patch, ps.index(), condition_type::dirichlet, *pde_bessel_J0exp);
            bcInfo_Dirichlet_im.addCondition(ps.patch, ps.index(), condition_type::dirichlet, *pde_bessel_mY0exp);
        }
    }

    gsInfo << "bcInfo_Gamma_f : \n" << bcInfo_Gamma_f << "\n";
    gsInfo << "bcInfo_Gamma_s : \n" << bcInfo_Gamma_s << "\n";
    gsInfo << "bcInfo_Gamma_d : \n" << bcInfo_Dirichlet_re << "\n";
    gsInfo << "bcInfo_Gamma_d : \n" << bcInfo_Dirichlet_im << "\n";

}

void gsStateEquationPotWaves::markBoundariesDirichlet()
{

    for (index_t i = 0; i < m_mp->nBoundary(); i++)
    {
        patchSide ps = m_mp->boundaries()[i];
        gsVector<unsigned> boundaryDofs = m_mp->basis(ps.patch).boundary(ps);

        bool cond1 = checkForZero(boundaryDofs,m_mp->patch(ps.patch).coefs(),1);
        bool cond2 = !(checkIfInPML(boundaryDofs,m_mp->patch(ps.patch).coefs(),pml_lx,pml_ly,pml_lz));

        if (cond1 and cond2)
        {
            
            //gsInfo << "GAMMA_D: Add condition at patch TRE" << ps.patch << " bnd " << ps.index() << "\n";
            bcInfo_Dirichlet_re.addCondition(ps.patch, ps.index(), condition_type::dirichlet, *pde_u_dir_re);
            bcInfo_Dirichlet_im.addCondition(ps.patch, ps.index(), condition_type::dirichlet, *pde_u_dir_im);
        }
    }

    gsInfo << "bcInfo_Gamma_f : \n" << bcInfo_Gamma_f << "\n";
    gsInfo << "bcInfo_Gamma_s : \n" << bcInfo_Gamma_s << "\n";
    gsInfo << "bcInfo_Gamma_d : \n" << bcInfo_Dirichlet_re << "\n";
    gsInfo << "bcInfo_Gamma_d : \n" << bcInfo_Dirichlet_im << "\n";

}

void gsStateEquationPotWaves::markBoundariesDirichletNoPML()
{

    for (index_t i = 0; i < m_mp->nBoundary(); i++)
    {
        patchSide ps = m_mp->boundaries()[i];
        gsVector<unsigned> boundaryDofs = m_mp->basis(ps.patch).boundary(ps);

        if (!checkForZero(boundaryDofs,m_mp->patch(ps.patch).coefs(),2))
        {
            //gsInfo << "GAMMA_D: Add condition at patch " << ps.patch << " bnd " << ps.index() << "\n";
            bcInfo_Dirichlet_re.addCondition(ps.patch, ps.index(), condition_type::dirichlet, *pde_u_dir_re);
            bcInfo_Dirichlet_im.addCondition(ps.patch, ps.index(), condition_type::dirichlet, *pde_u_dir_im);
        }
    }

    gsInfo << "bcInfo_Gamma_f : \n" << bcInfo_Gamma_f << "\n";
    gsInfo << "bcInfo_Gamma_s : \n" << bcInfo_Gamma_s << "\n";
    gsInfo << "bcInfo_Gamma_d : \n" << bcInfo_Dirichlet_re << "\n";
    gsInfo << "bcInfo_Gamma_d : \n" << bcInfo_Dirichlet_im << "\n";

}

void gsStateEquationPotWaves::markBoundariesDirichletNoPML_NoCenter()
{

    for (index_t i = 0; i < m_mp->nBoundary(); i++)
    {
        patchSide ps = m_mp->boundaries()[i];
        gsVector<unsigned> boundaryDofs = m_mp->basis(ps.patch).boundary(ps);

        bool cond1 = !checkForZero(boundaryDofs,m_mp->patch(ps.patch).coefs(),2);
        bool cond2 = !(checkIfInPML(boundaryDofs,m_mp->patch(ps.patch).coefs(),pml_lx,pml_ly,pml_lz));
        if (cond1 and cond2)
        {
            //gsInfo << "GAMMA_D: Add condition at patch " << ps.patch << " bnd " << ps.index() << "\n";
            bcInfo_Dirichlet_re.addCondition(ps.patch, ps.index(), condition_type::dirichlet, *pde_u_dir_re);
            bcInfo_Dirichlet_im.addCondition(ps.patch, ps.index(), condition_type::dirichlet, *pde_u_dir_im);
        }
    }

    gsInfo << "bcInfo_Gamma_f : \n" << bcInfo_Gamma_f << "\n";
    gsInfo << "bcInfo_Gamma_s : \n" << bcInfo_Gamma_s << "\n";
    gsInfo << "bcInfo_Gamma_d : \n" << bcInfo_Dirichlet_re << "\n";
    gsInfo << "bcInfo_Gamma_d : \n" << bcInfo_Dirichlet_im << "\n";

}

bool gsStateEquationPotWaves::isBndGamma_s(index_t b)
{
    return m_isBndGamma_s_vec[b];
}

gsMatrix<> gsStateEquationPotWaves::getCoefs(index_t p)
{
    gsMatrix<> out(8,3);

    real_t lbx = init_lbx; // Length of box in the middle of the domain
    real_t lby = init_lby; // Length of box in the middle of the domain
    real_t lbz = init_lbz; // Length of box in the middle of the domain

    switch (p){

        case 0: out <<  
                        -pml_lx     ,0          ,-pml_lz,
                        -lbx        ,0          ,-lbz,
                        -pml_lx     ,pml_ly     ,-pml_lz,
                        -lbx        ,lby        ,-lbz,
                        -pml_lx     ,0          ,0,
                        -lbx        ,0          ,0,
                        -pml_lx     ,pml_ly     ,0,
                        -lbx        ,lby          ,0;
                break;

        case  1: out <<  
                        pml_lx     ,0          ,0,
                        lbx        ,0          ,0,
                        pml_lx     ,pml_ly     ,0,
                        lbx        ,lby          ,0,
                        pml_lx     ,0          ,-pml_lz,
                        lbx        ,0          ,-lbz,
                        pml_lx     ,pml_ly     ,-pml_lz,
                        lbx        ,lby        ,-lbz;
                break;

        case  2: out << -lbx    ,lby        ,-lbz,
                        lbx     ,lby        ,-lbz,
                        -pml_lx ,pml_ly     ,-pml_lz,
                        pml_lx  ,pml_ly     ,-pml_lz,
                        -lbx    ,lby        ,0,
                        lbx     ,lby        ,0,
                        -pml_lx ,pml_ly     ,0,
                        pml_lx  ,pml_ly     ,0;
                 break;

        case  3: out << -pml_lx ,0          ,-pml_lz,
                        pml_lx  ,0          ,-pml_lz,
                        -pml_lx ,pml_ly     ,-pml_lz,
                        pml_lx  ,pml_ly     ,-pml_lz,
                        -lbx    ,0          ,-lbz,
                        lbx     ,0          ,-lbz,
                        -lbx    ,lby        ,-lbz,
                        lbx     ,lby        ,-lbz;
                 break;

        case  4: out = getBoxHelper(-pml_lx,pml_lx,0,pml_ly,-pml_Lz,-pml_lz); break;

        case  5: out = getBoxHelper(-pml_Lx,-pml_lx,0,pml_ly,-pml_lz,0); break;
        case  6: out = getBoxHelper(pml_lx,pml_Lx,0,pml_ly,-pml_lz,0);   break;

        case  7: out = getBoxHelper(-pml_Lx,-pml_lx,pml_ly,pml_Ly,-pml_lz,0);   break; 
        case  8: out = getBoxHelper(pml_lx,pml_Lx,pml_ly,pml_Ly,-pml_lz,0);     break;

        case  9: out = getBoxHelper(-pml_lx,pml_lx,pml_ly,pml_Ly,-pml_lz,0); break;
        case  10: out = getBoxHelper(-pml_Lx,-pml_lx,0,pml_ly,-pml_Lz,-pml_lz); break;
        case  11: out = getBoxHelper(pml_lx,pml_Lx,0,pml_ly,-pml_Lz,-pml_lz); break;

        case  12: out = getBoxHelper(-pml_lx,pml_lx,pml_ly,pml_Ly,-pml_Lz,-pml_lz); break;
        case  13: out = getBoxHelper(-pml_Lx,-pml_lx,pml_ly,pml_Ly,-pml_Lz,-pml_lz); break;
        case  14: out = getBoxHelper(pml_lx,pml_Lx,pml_ly,pml_Ly,-pml_Lz,-pml_lz); break;

        // The center patch. To use if you want to 'fill it with water'
        case  15: out = getBoxHelper(-lbx,lbx,0,lby,-lbz,0); break;

        /*case  4: out <<  0   ,0  ,0,
                        1   ,0  ,0,
                        0   ,1  ,0,
                        1   ,1  ,0,
                        0   ,0  ,1,
                        1   ,0  ,1,
                        0   ,1  ,1,
                        2   ,1  ,1;
                 break;
        */

    }
    return out;
}

gsMultiPatch<>::Ptr gsStateEquationPotWaves::getInitialDomain(bool includeCenter)
{
    index_t nPatches = 15;

    if (includeCenter)
        nPatches = 16;

    gsInfo << "NPATCHES : " << nPatches << "\n";

	gsMultiPatch<>::Ptr out = memory::make_shared( new gsMultiPatch<> );
    for (index_t p = 0; p < nPatches; p++)
    {
        gsMatrix<> coefs = getCoefs(p);
        gsTensorBSpline<3, real_t> tbs = getBox(coefs);
        out->addPatch( tbs );

        //gsInfo << "Sign of detJ on patch " << p << " is " << getSignOfDetJ(tbs) << "\n";
    }

    out->computeTopology();

    out->uniformRefine();
    out->degreeElevate(m_degree-1);
    out->closeGaps();


    return out;

}

gsMultiPatch<>::Ptr gsStateEquationPotWaves::getInitialDomainTmp(bool includeCenter)
{
    index_t nPatches = 7;

    if (includeCenter)
        nPatches = 8;

    gsVector<> vec(8);
    vec << 0,1,2,3,4,9,12,15;

	gsMultiPatch<>::Ptr out = memory::make_shared( new gsMultiPatch<> );
    for (index_t i = 0; i < nPatches; i++)
    {
        index_t p = vec[i];
        gsMatrix<> coefs = getCoefs(p);
        gsTensorBSpline<3, real_t> tbs = getBox(coefs);
        out->addPatch( tbs );
    }

    out->computeTopology();

    return out;

}

gsMultiPatch<>::Ptr gsStateEquationPotWaves::getInitialDomainNoPML(bool includeCenter)
{
    index_t nPatches = 4;

    if (includeCenter)\
        nPatches = 5;

    gsVector<> vec(5);
    vec << 0,1,2,3,15;

	gsMultiPatch<>::Ptr out = memory::make_shared( new gsMultiPatch<> );
    for (index_t i = 0; i < nPatches; i++)
    {
        index_t p = vec[i];
        gsMatrix<> coefs = getCoefs(p);
        gsTensorBSpline<3, real_t> tbs = getBox(coefs);
        out->addPatch( tbs );
    }

    out->computeTopology();

    return out;

}

gsMultiPatch<>::Ptr gsStateEquationPotWaves::getInitialDomainNoPMLInZDir(bool includeCenter)
{
    index_t nPatches = 5;

    if (includeCenter)\
        nPatches++;

    gsVector<> vec(6);
    vec << 0,1,2,3,9,15;

	gsMultiPatch<>::Ptr out = memory::make_shared( new gsMultiPatch<> );
    for (index_t i = 0; i < nPatches; i++)
    {
        index_t p = vec[i];
        gsMatrix<> coefs = getCoefs(p);
        gsTensorBSpline<3, real_t> tbs = getBox(coefs);
        out->addPatch( tbs );
    }

    out->computeTopology();

    return out;

}

void gsStateEquationPotWaves::getUI(gsMultiPatch<> &uI_re, gsMultiPatch<> &uI_im) 
{
    
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::solution    solution;

    gsExprAssembler<> A(1,1);

    A.options().setInt("quB",m_quB);
    A.options().setReal("quA",m_quA);

    // Elements used for numerical integration
    A.setIntegrationElements(dbasis);

    geometryMap G = A.getMap(*m_mp);

    variable uIRe = A.getCoeff(*wave_u_I_re,G);
    variable uIIm = A.getCoeff(*wave_u_I_im,G);

    space u = A.getSpace(dbasis);

    // Get uI_re
    A.initSystem();

    gsMatrix<> solVector;
    solution u_sol = A.getSolution(u, solVector);

    A.assemble(u*u.tr()*meas(G),u*uIRe*meas(G));

    gsSparseSolver<>::LU solver;
    solver.compute(A.matrix());
    solVector = solver.solve(A.rhs());

    u_sol.extract(uI_re);

    // Get uI_im
    A.initSystem();

    A.assemble(u*uIIm*meas(G));

    solVector = solver.solve(A.rhs());

    u_sol.extract(uI_im);

}

void gsStateEquationPotWaves::plotVelocityField(gsMultiPatch<> &ur, gsMultiPatch<> &ui, real_t timestep, std::string outfile, bool includeIncident)
{
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::solution    solution;

    gsExprAssembler<> A(1,1);

    A.options().setInt("quB",m_quB);
    A.options().setReal("quA",m_quA);

    // Elements used for numerical integration
    gsMultiBasis<> dbas(ur);
    A.setIntegrationElements(dbas);
    gsExprEvaluator<> ev(A);

    geometryMap G = A.getMap(*m_mp);

    variable uSRe = A.getCoeff(ur);
    variable uSIm = A.getCoeff(ui);

    variable uIRe = A.getCoeff(*wave_u_I_re,G);
    variable uIIm = A.getCoeff(*wave_u_I_im,G);

    index_t i = 20;
    real_t t = 0;

    //ev.options().setSwitch("plot.elements", true);
    ev.options().setInt("plot.npts", 8000);

    index_t maxIter = 1000;

    while( t < 2*M_PI/wave_omega)
    {
        // Generate name of file
        std::string name = outfile + "_" + std::to_string( i );

        // Calculate velocity field at time t
        real_t coswt = cos( wave_omega * t);
        real_t sinwt = sin( wave_omega * t);

        auto velocity_field_I = wave_A*uIRe*coswt - wave_A*uIIm*sinwt;
        auto velocity_field_S = wave_A*uSRe*coswt - wave_A*uSIm*sinwt;

        // Plot
        gsInfo<< "Plotting " << name << " in Paraview...\n";
        if (includeIncident)
        {
            ev.writeParaview( velocity_field_I + velocity_field_S    , G, name);
        } else {
            ev.writeParaview( velocity_field_S    , G, name);
        }


        // Update time
        t += timestep;
        i++;

        if (i > maxIter)
           return; 
    }


}

void gsStateEquationPotWaves::plotVelocityBessel(real_t timestep, std::string outfile)
{
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::solution    solution;

    gsExprAssembler<> A(1,1);

    A.options().setInt("quB",m_quB);
    A.options().setReal("quA",m_quA);

    // Elements used for numerical integration
    gsMultiBasis<> dbas(*m_mp);
    A.setIntegrationElements(dbas);
    gsExprEvaluator<> ev(A);

    geometryMap G = A.getMap(*m_mp);

    getBessel();
    variable uRe = A.getCoeff(*pde_bessel_J0exp,G);
    variable uIm = A.getCoeff(*pde_bessel_mY0exp,G);

    // Solve for Bessel
    useNeu = false;
    useDirRI = true;

    pde_u_dir_re = memory::make_unique( new gsFunctionExpr<>(*pde_bessel_J0exp ) );
    pde_u_dir_im = memory::make_unique( new gsFunctionExpr<>(*pde_bessel_mY0exp) );

    gsMultiPatch<> markInner = markInnerDomain();
    variable delta = A.getCoeff(markInner);
    markCenterBessel();

    gsMultiPatch<> ur,ui;
    solve(ur,ui);
    // ...

    useNeu = true;
    useDirRI = false;

    variable uSRe = A.getCoeff(ur);
    variable uSIm = A.getCoeff(ui);

    gsInfo << "Diff real: " << ev.integral(delta*(uSRe-uRe).norm()*meas(G)) << "\n";
    gsInfo << "Diff real: " << ev.integral(delta*(uSIm-uIm).norm()*meas(G)) << "\n";

    std::string nm = outfile + "_ur";
    ev.writeParaview( uSRe    , G, nm);
    nm = outfile + "_ui";
    ev.writeParaview( uSIm    , G, nm);



    index_t i = 20;
    real_t t = 0;

    //ev.options().setSwitch("plot.elements", true);
    ev.options().setInt("plot.npts", 8000);

    index_t maxIter = 1000;

    while( t < 2*M_PI/wave_omega)
    {
        // Generate name of file
        std::string name = outfile + "_" + std::to_string( i );

        // Calculate velocity field at time t
        real_t coswt = cos( wave_omega * t);
        real_t sinwt = sin( wave_omega * t);

        auto velocity_field = wave_A*uRe*coswt - wave_A*uIm*sinwt;

        // Plot
        gsInfo<< "Plotting " << name << " in Paraview...\n";
        ev.writeParaview( velocity_field    , G, name);

        auto velocity_field_s = wave_A*uSRe*coswt - wave_A*uSIm*sinwt;

        // Plot
        gsInfo<< "Plotting " << name << " in Paraview...\n";
        ev.writeParaview( velocity_field    , G, name);
        std::string name2 = outfile + "_s" + std::to_string( i );
        ev.writeParaview( velocity_field_s    , G, name2 );


        // Update time
        t += timestep;
        i++;

        if (i > maxIter)
           return; 
    }


}

void gsStateEquationPotWaves::convergenceTest(index_t max_refine, std::string outfolder)
{

    // convergence test against manufactured solution u = cos(Ky)exp(Kz)
    // run with increased refinement until max_refine

    // Save old rhs
    gsFunctionExpr<> tmp_F_re = *pde_F_re;
    gsFunctionExpr<> tmp_F_im = *pde_F_im;

    // Create new rhs:
    complexExpr u_expr("cos(K*y)*exp(K*z)","0.0");
    pde_u_dir_re = u_expr.getFunRe(const_map);
    pde_u_dir_im = u_expr.getFunIm(const_map);

    complexExpr du_dx_expr("0.0","0.0");
    complexExpr du_dy_expr("- K*sin(K*y)*exp(K*z)","0.0");
    complexExpr du_dz_expr("K*cos(K*y)*exp(K*z)","0.0");

    pde_F_re = getVectorFunRe(const_map, du_dx_expr, du_dy_expr, du_dz_expr);
    pde_F_im = getVectorFunIm(const_map, du_dx_expr, du_dy_expr, du_dz_expr);

    convTestHelper(false, true, max_refine, outfolder);

    // Return pde_F to original
    pde_F_re = memory::make_unique(&tmp_F_re);
    pde_F_im = memory::make_unique(&tmp_F_im);

}

void gsStateEquationPotWaves::convergenceTestNoPML_NoCenter(index_t max_refine, std::string outfolder)
{

    m_mp = getInitialDomainNoPML(false); // true : fills center with water

    // convergence test against manufactured solution u = cos(Ky)exp(Kz)
    // run with increased refinement until max_refine
    complexExpr u_expr("cos(K*y)*exp(K*z)","0.0");
    pde_u_dir_re = u_expr.getFunRe(const_map);
    pde_u_dir_im = u_expr.getFunIm(const_map);

    setup();

    gsBoundaryConditions<> tmp;
    bcInfo_Gamma_f = tmp;

    gsBoundaryConditions<> tmps;
    bcInfo_Gamma_s = tmps;

    gsBoundaryConditions<> tmpb;
    bcInfo_Gamma_b = tmpb;

    markBoundaries();
    markBoundariesDirichletNoPML_NoCenter();

    // Save old rhs
    gsFunctionExpr<> tmp_F_re = *pde_F_re;
    gsFunctionExpr<> tmp_F_im = *pde_F_im;

    // Create new rhs:
    complexExpr du_dx_expr("0.0","0.0");
    complexExpr du_dy_expr("- K*sin(K*y)*exp(K*z)","0.0");
    complexExpr du_dz_expr("K*cos(K*y)*exp(K*z)","0.0");

    pde_F_re = getVectorFunRe(const_map, du_dx_expr, du_dy_expr, du_dz_expr);
    pde_F_im = getVectorFunIm(const_map, du_dx_expr, du_dy_expr, du_dz_expr);

    convTestHelper(true, true, max_refine, outfolder);

    // Return pde_F to original
    pde_F_re = memory::make_unique(&tmp_F_re);
    pde_F_im = memory::make_unique(&tmp_F_im);

}

void gsStateEquationPotWaves::convergenceTestOnlyPML(index_t max_refine, std::string outfolder)
{

    m_mp = getInitialDomainNoPMLInZDir(true); // true : fills center with water

    complexExpr u_expr("cos(K*y)*exp(K*z)","0.0");
    pde_u_dir_re = u_expr.getFunRe(const_map);
    pde_u_dir_im = u_expr.getFunIm(const_map);

    setup();

    gsBoundaryConditions<> tmp;
    bcInfo_Gamma_f = tmp;

    gsBoundaryConditions<> tmpb;
    bcInfo_Gamma_b = tmpb;

    markBoundaries();
    markBoundariesDirichlet();

    convTestHelper(true, false, max_refine, outfolder);

}

void gsStateEquationPotWaves::convergenceTestOnlyPMLAllDir(index_t max_refine, std::string outfolder, real_t angle)
{

    real_t angleRad = angle*M_PI/180;
    gsDebugVar(angle);
    gsDebugVar(angleRad);

    m_mp = getInitialDomain(true); // true : fills center with water
    //m_mp = getInitialDomainTmp(true); // true : fills center with water

    real_t cosA = cos(angleRad);
    real_t sinA = sin(angleRad);
    gsDebugVar(cosA);
    gsDebugVar(sinA);

    const_map["cosA"] = cosA;
    const_map["sinA"] = sinA;

    complexExpr u_expr("cos(K*(sinA*x + cosA*y))*exp(K*z)","0.0");

    pde_u_dir_re = u_expr.getFunRe(const_map);
    pde_u_dir_im = u_expr.getFunIm(const_map);

    setup();

    gsBoundaryConditions<> tmp;
    bcInfo_Gamma_f = tmp;

    gsBoundaryConditions<> tmpb;
    bcInfo_Gamma_b = tmpb;

    markBoundaries();
    markBoundariesDirichlet();

    convTestHelper(true, false, max_refine, outfolder);


}

void gsStateEquationPotWaves::convergenceTestBessel(index_t max_refine, std::string outfolder)
{

    getBessel();
    pde_u_dir_re = memory::make_unique( new gsFunctionExpr<>(*pde_bessel_J0exp));
    pde_u_dir_im = memory::make_unique( new gsFunctionExpr<>(*pde_bessel_mY0exp));

    setup();

    gsBoundaryConditions<> tmp;
    bcInfo_Gamma_f = tmp;

    gsBoundaryConditions<> tmpb;
    bcInfo_Gamma_b = tmpb;

    markBoundaries();
    markCenterBessel();

    useDirRI = true;
    convTestHelper(false, false, max_refine, outfolder);

    gsDebugVar(m_ur); 
    gsDebugVar(m_ui); 
    real_t timestep = 0.025*2*M_PI/wave_omega;
    std::string vfold = "velocity/";
    gsInfo << "Does " << outfolder + vfold << " exist? \n\n";
    plotVelocityField(m_ur,m_ui,timestep,outfolder + vfold,false); // false means to only plot 'scattering field'
    

    useDirRI = false;

}

void gsStateEquationPotWaves::pointSourceTest(std::string outfolder)
{

    useNeuBell = true;
    useDir = false;
    useNeu = false;

    m_mp = getInitialDomain(true); // true : fills center with water
    gsInfo << m_mp << "\n";

    setup();

    complexExpr u_expr("cos(K*sqrt(x*x + y*y))*exp(K*z)","0.0");

    pde_u_dir_re = u_expr.getFunRe(const_map);
    pde_u_dir_im = u_expr.getFunIm(const_map);

    complexExpr laplU_expr("K * exp(K*z) * sin(K*sqrt(x*x + y*y)) / sqrt(x*x + y*y)","0.0");

    lapl_uEx = laplU_expr.getFunRe(const_map);

    gsBoundaryConditions<> tmp;
    bcInfo_Gamma_f = tmp;

    gsBoundaryConditions<> tmps;
    bcInfo_Gamma_s = tmps;

    gsBoundaryConditions<> tmpb;
    bcInfo_Gamma_b = tmpb;

    markBoundaries();

    gsExprAssembler<> A(1,1);
    gsExprEvaluator<> ev(A);

    ev.options().setInt("quB",m_quB); // FIXIT lower no points
    ev.options().setReal("quA", m_quA); // FIXIT lower no points

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    geometryMap G = A.getMap(*m_mp);

    gsMultiPatch<> markInner = markInnerDomain();
    variable delta = A.getCoeff(markInner);

    A.setIntegrationElements(dbasis);

    // Use last solution, only works since we solved in convTestHelper
    variable u_r = A.getCoeff(m_ur); 
    variable u_i = A.getCoeff(m_ui);

	variable B_real = A.getCoeff(*bell_term_re,G);
	variable B_imag = A.getCoeff(*bell_term_im,G);

    ev.options().setInt("quB",m_quB+2); // FIXIT lower no points
    ev.options().setReal("quA", m_quA+2); // FIXIT lower no points

    gsFunctionExpr<> ffun("z",3);
    variable ff = A.getCoeff(ffun);
    auto zdir = fjac(ff);

    if (true) // If plot
    {

        ev.options().setInt("plot.npts", 8000);
    
        std::string name = "ur";
        gsInfo<< "Plotting " << name << " in Paraview...\n";
        ev.writeParaview( u_r , G, outfolder + name);
    
        name = "ui";
        gsInfo<< "Plotting " << name << " in Paraview...\n";
        ev.writeParaview( u_i , G, outfolder + name);

        name = "bcreal";
        gsInfo<< "Plotting " << name << " in Paraview...\n";
        ev.writeParaview((wave_K*u_r - fjac(u_r).tr()*jac(G).inv()*zdir)*delta.val(), G, outfolder + name);

        name = "bcimag";
        gsInfo<< "Plotting " << name << " in Paraview...\n";
        ev.writeParaview((wave_K*u_i - fjac(u_i).tr()*jac(G).inv()*nv(G))*delta.val(), G, outfolder + name);

        name = "Breal";
        gsInfo<< "Plotting " << name << " in Paraview...\n";
        ev.writeParaview(-B_real, G, outfolder + name);

        name = "Bimag";
        gsInfo<< "Plotting " << name << " in Paraview...\n";
        ev.writeParaview(-B_imag, G, outfolder + name);

        name = "bcdiff_re";
        gsInfo<< "Plotting " << name << " in Paraview...\n";
        ev.writeParaview((B_real + wave_K*u_r - fjac(u_r).tr()*jac(G).inv()*zdir)*delta.val(), G, outfolder + name);

        name = "bcdiff_im";
        gsInfo<< "Plotting " << name << " in Paraview...\n";
        ev.writeParaview((B_imag + wave_K*u_i - fjac(u_i).tr()*jac(G).inv()*nv(G))*delta.val(), G, outfolder + name);
    }


    // Return to default settings
    useNeuBell = false;
    useNeu = true;


    real_t timestep = 0.025*2*M_PI/wave_omega;
    std::string vfold = "veloc/";
    gsInfo << "Does " << outfolder + vfold << " exist? \n\n";
    plotVelocityField(m_ur,m_ui,timestep,outfolder + vfold,false); // false means to only plot 'scattering field'

}

void gsStateEquationPotWaves::pointSourceTestForce(std::string outfolder)
{

    useForce = true;
    useDir = true;
    useNeu = false;

    m_mp = getInitialDomain(false); // true : fills center with water

    setup();

    complexExpr u_expr("cos(K*sqrt(x*x + y*y))*exp(K*z)","0.0");

    pde_u_dir_re = u_expr.getFunRe(const_map);
    pde_u_dir_im = u_expr.getFunIm(const_map);

    complexExpr laplU_expr("K * exp(K*z) * sin(K*sqrt(x*x + y*y)) / sqrt(x*x + y*y)","0.0");

    lapl_uEx = laplU_expr.getFunRe(const_map);

    gsBoundaryConditions<> tmp;
    bcInfo_Gamma_f = tmp;

    gsBoundaryConditions<> tmps;
    bcInfo_Gamma_s = tmps;

    gsBoundaryConditions<> tmpb;
    bcInfo_Gamma_b = tmpb;

    markBoundaries();
    markBoundariesDirichletNoPML_NoCenter();

    gsExprAssembler<> A(1,1);
    gsExprEvaluator<> ev(A);

    ev.options().setInt("quB",m_quB); // FIXIT lower no points
    ev.options().setReal("quA", m_quA); // FIXIT lower no points

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    geometryMap G = A.getMap(*m_mp);

    gsMultiPatch<> ur, ui;
    solve(ur,ui);

    A.setIntegrationElements(dbasis);

    variable u_r = A.getCoeff(ur);
    variable u_i = A.getCoeff(ui);

    ev.options().setInt("quB",m_quB+2); // FIXIT lower no points
    ev.options().setReal("quA", m_quA+2); // FIXIT lower no points

    if (true) // If plot
    {

        ev.options().setInt("plot.npts", 8000);
    
        std::string name = "ur";
        gsInfo<< "Plotting " << name << " in Paraview...\n";
        ev.writeParaview( u_r , G, outfolder + name);
    
        name = "ui";
        gsInfo<< "Plotting " << name << " in Paraview...\n";
        ev.writeParaview( u_i , G, outfolder + name);
    }


    // Return to default settings
    useNeuBell = true;

}

void gsStateEquationPotWaves::convergenceTestNoPML(index_t max_refine, std::string outfolder)
{

    
    m_mp = getInitialDomainNoPML(true); // true : fills center with water

    complexExpr u_expr("cos(K*y)*exp(K*z)","0.0");
    pde_u_dir_re = u_expr.getFunRe(const_map);
    pde_u_dir_im = u_expr.getFunIm(const_map);

    setup();

    gsBoundaryConditions<> tmp;
    bcInfo_Gamma_f = tmp;

    markBoundaries();
    markBoundariesDirichletNoPML();

    convTestHelper(true, false, max_refine, outfolder);

}

void gsStateEquationPotWaves::convTestHelper(bool useDirichlet, bool useNeumann, index_t max_refine, std::string outfolder)
{
    std::ofstream f;
    f.open(outfolder + "conv.txt");

    gsExprAssembler<> A(1,1);
    gsExprEvaluator<> ev(A);

    ev.options().setInt("quB",m_quB); // FIXIT lower no points
    ev.options().setReal("quA", m_quA); // FIXIT lower no points

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    // We want to only compare the solution in the domain of interest
    gsMultiPatch<> markInner = markInnerDomain();

    geometryMap G = A.getMap(*m_mp);

    useDir = useDirichlet;
    useNeu = useNeumann;

    for (index_t i = 0; i < max_refine; i++)
    {
        gsMultiPatch<> ur, ui;
        solve(ur,ui);

        A.setIntegrationElements(dbasis);

        variable u_er = A.getCoeff(*pde_u_dir_re,G);
        variable u_ei = A.getCoeff(*pde_u_dir_im,G);
        variable u_r = A.getCoeff(ur);
        variable u_i = A.getCoeff(ui);

        variable delta = A.getCoeff(markInner);

        ev.options().setInt("quB",m_quB+2); // FIXIT lower no points
        ev.options().setReal("quA", m_quA+2); // FIXIT lower no points

        real_t errorSq_re = ev.integral( delta*(u_er - u_r).sqNorm()*meas(G));
        real_t errorSq_im = ev.integral( delta*(u_ei - u_i).sqNorm()*meas(G));
        gsInfo << "size : \t " << dbasis.size() << " , Error real \t : " << sqrt(errorSq_re) << " , Error imag \t : " << sqrt(errorSq_im) <<"\n";

        real_t h = 0;
        for (index_t p = 0; p < m_mp->nBoxes(); p++)
        {
            real_t hp = static_cast<gsTensorBSplineBasis<3>&> (dbasis.basis(p)).getMaxCellLength();
            if (hp > h)
                h = hp;

        }

        f << dbasis.size() << " " << h << " " << sqrt(errorSq_re) << " " << sqrt(errorSq_im) << "\n";
        
        if ( true ) // If plot
        {

            ev.options().setInt("plot.npts", 8000);
        
            std::string name = "delta";
            gsInfo<< "Plotting " << name << " in Paraview...\n";
            ev.writeParaview( delta , G, outfolder + name);
            
            name = "ur";
            gsInfo<< "Plotting " << name << " in Paraview...\n";
            ev.writeParaview( u_r , G, outfolder + name);
        
            name = "ui";
            gsInfo<< "Plotting " << name << " in Paraview...\n";
            ev.writeParaview( u_i , G, outfolder + name);
        
            name = "uer";
            gsInfo<< "Plotting " << name << " in Paraview...\n";
            ev.writeParaview( u_er , G, outfolder + name);

            name = "uei";
            gsInfo<< "Plotting " << name << " in Paraview...\n";
            ev.writeParaview( u_ei , G, outfolder + name);
        
            name = "diff_r";
            gsInfo<< "Plotting " << name << " in Paraview...\n";
            ev.writeParaview( delta*(u_er - u_r)  , G, outfolder + name);

            name = "diff_i";
            gsInfo<< "Plotting " << name << " in Paraview...\n";
            ev.writeParaview( delta*(u_ei - u_i)  , G, outfolder + name);
        }

        dbasis.uniformRefine();

    }
    f.close();


    // Return to default settings
    useDir = false;
    useNeu = true;

}

void gsStateEquationPotWaves::solve(gsMultiPatch<> &u_real, gsMultiPatch<> &u_imag)
{
        gsInfo << ".. ASSEMBLE .. " << std::flush;
		gsSparseMatrix<> mat;
		gsVector<> rhs;
		getSystem(mat,rhs);

		// gsInfo << "norm of mat: " << mat.norm() << "\n" ;
		// gsInfo << "norm of rhs: " << rhs.norm() << "\n" ;

		// gsInfo << "Save A \n" << std::flush;
		// std::ofstream file("A.txt");
		// file << mat.toDense();
		// file.close();

		solver.compute(mat);
        gsInfo << ".. SOLVE .. " << std::flush;
		solVector = solver.solve(rhs);
		


		// Plot
		gsExprAssembler<> A(1,1);
		gsExprEvaluator<> ev(A);


		geometryMap G = A.getMap(*m_mp);

		A.setIntegrationElements(dbasis);

        // Real part
		space u_r = A.getSpace(dbasis);
		u_r.setInterfaceCont(0);
		u_r.addBc(bcInfo_Dirichlet_re.get("Dirichlet"));

		A.initSystem();
		gsMatrix<> solVector_Real = solVector.block(0,0,mat.cols()/2,1);
		solution u_sol_real = A.getSolution(u_r,solVector_Real);
		u_sol_real.extract(u_real);

        // Imag part
		space u_i = A.getSpace(dbasis);
		u_i.setInterfaceCont(0);
		u_i.addBc(bcInfo_Dirichlet_im.get("Dirichlet"));

		A.initSystem();
		gsMatrix<> solVector_Imag = solVector.block(mat.cols()/2,0,mat.cols()/2,1);
		solution u_sol_imag = A.getSolution(u_i,solVector_Imag);
		u_sol_imag.extract(u_imag);

        //gsInfo << "Norm of solVector_imag : " << solVector_Imag.norm() << "\n";
        //gsInfo << "Norm of u_sol : " << ev.integral(u_sol_imag.sqNorm()) << "\n";

        //gsExprAssembler<>::variable utest = ev.getVariable(u_imag);
        //ev.options().setInt("plot.npts", 8000);
        //gsInfo << "Norm of utest : " << ev.integral(utest) << "\n";

        //ev.writeParaview( utest, G, "test");
        //

        m_ur = u_real;
        m_ui = u_imag;



}

gsMultiPatch<> gsStateEquationPotWaves::getPieceWiseFunctionOnPatch(gsVector<> vals)
{
	gsKnotVector<> u_knots(0,1,0,1);
	gsKnotVector<> v_knots(0,1,0,1);
	gsKnotVector<> w_knots(0,1,0,1);

	gsTensorBSplineBasis<3> T_basis(u_knots,v_knots,w_knots);


	gsMultiPatch<> pw;

	for(index_t i = 0; i < vals.size(); i++){
        
	    gsMatrix<> cf_patch(1,1);
	    cf_patch << vals[i];

  		gsGeometry<>::uPtr b = T_basis.makeGeometry(cf_patch);
			pw.addPatch(*b);
	}

	return pw;
}

gsMultiPatch<> gsStateEquationPotWaves::markPatches(gsVector< index_t > patches)
{
    index_t n_patches = m_mp->nBoxes();

    gsVector<> mark;
    mark.setZero(n_patches);

    for( index_t i = 0; i < patches.size(); i++)
    {
        index_t p = patches[i];
        if (p < n_patches)
            mark[ p ] = 1;
    }

    return getPieceWiseFunctionOnPatch(mark);
}

gsMultiPatch<> gsStateEquationPotWaves::markInnerDomain()
{
    gsVector< index_t > patches = getVectorWithDomainPatches();

    return markPatches(patches);

}

void gsStateEquationPotWaves::testSplineSpace(real_t k)
{
    gsExprAssembler<> A(1,1);
    gsExprEvaluator<> ev(A);

    ev.options().setInt("quB",m_quB); // FIXIT lower no points
    ev.options().setReal("quA", m_quA); // FIXIT lower no points

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    // We want to only compare the solution in the domain of interest
    gsMultiPatch<> markInner = markInnerDomain();

    geometryMap G = A.getMap(*m_mp);

    std::string strx = "x^" + std::to_string(k);
    std::string stry = "y^" + std::to_string(k);
    std::string strz = "z^" + std::to_string(k);

    gsFunctionExpr<> poly(strx + stry + strz, 3);

    variable ue = A.getCoeff(poly);

    space u = A.getSpace(dbasis);
    
    gsMatrix<> solVec;
    solution usol = A.getSolution(u,solVec);

    A.setIntegrationElements(dbasis);

    A.initSystem();
    A.assemble(u*u.tr(), u*ue);

    gsSparseSolver<>::LU solverLU;
    solverLU.compute(A.matrix());
    solVec = solverLU.solve(A.rhs());

    gsInfo << "Err: " << ev.integral( (ue-usol).sqNorm() ) << "\n";
}

void gsStateEquationPotWaves::getTerm(index_t realOrImag, gsSparseMatrix<> &mat, gsVector<> &rhs)
{
	// Method to generate system matrices and rhs
	// Input 0 for real part of matrix, and 1 for imaginary part of matrix

	gsExprAssembler<> A(1,1);
	gsExprEvaluator<> ev(A);

    A.options().setReal("quA",m_quA);
    A.options().setInt("quB",m_quB);
    
    geometryMap G = A.getMap(*m_mp);

    A.setIntegrationElements(dbasis);

    space u = A.getSpace(dbasis);
    u.setInterfaceCont(0);

    // Setup terms

    // Laplace Term

    variable detJJm1Jm1_asVec_re = A.getCoeff(*strech_detJ_Jm1_Jm1_asVec_re,G);
    variable detJJm1Jm1_asVec_im = A.getCoeff(*strech_detJ_Jm1_Jm1_asVec_im,G);

    auto dJJJ_re = detJJm1Jm1_asVec_re.asDiag();
    auto dJJJ_im = detJJm1Jm1_asVec_im.asDiag();

    auto laplace_term_real = igrad(u,G)*(dJJJ_re*igrad(u,G).tr())*meas(G);
    auto laplace_term_imag = igrad(u,G)*(dJJJ_im*igrad(u,G).tr())*meas(G);

    variable dJ_re = A.getCoeff(*strech_detJ_re,G);
    variable dJ_im = A.getCoeff(*strech_detJ_im,G);

    // Zero
    gsFunctionExpr<> zerfun("0.0",3);
    variable zer = ev.getVariable(zerfun);

	// Helmholtz Term
	auto helmholtz_term_real = - wave_K*u*u.tr()*dJ_re.val()*nv(G).norm();
	auto helmholtz_term_imag = - wave_K*u*u.tr()*dJ_im.val()*nv(G).norm();

	// Rhs
	variable F_real = A.getCoeff(*pde_F_re,G);
	variable F_imag = A.getCoeff(*pde_F_im,G);

    gsFunctionExpr<> ones_fun("1","1","1",3);
    variable ones = A.getCoeff(ones_fun);

	auto rhs_term_real = (F_real.tr()*nv(G)).val()*u; // We omit *nv(G).norm() since we need unitnv
	auto rhs_term_imag = (F_imag.tr()*nv(G)).val()*u; //*nv(G).norm();

    auto zero_matrix = zer.val()*u*u.tr();  // Perhaps I need to multiply with u*u.tr();
    auto zero_vec = zer.val()*u*nv(G).norm();     // Perhaps I need to multiply with u;

    gsMultiPatch<> markInner = markInnerDomain();
    variable delta = A.getCoeff(markInner);

    // FIXIT : Remove these terms (or comment the out) 
    //         so that we don't spend time precomputing B_real, laplU_real, dJ_re and dJ_im
	// Rhs Bell
	variable B_real = A.getCoeff(*bell_term_re,G);
	variable B_imag = A.getCoeff(*bell_term_im,G);
	auto rhs_bell_real = B_real.val()*u*nv(G).norm();
	auto rhs_bell_imag = B_imag.val()*u*nv(G).norm();

	// Rhs Bell
    //
	variable laplU_real = A.getCoeff(*lapl_uEx,G);

	auto rhs_force_real = laplU_real.val()*dJ_re.val()*u*meas(G);

	auto rhs_force_imag = laplU_real.val()*dJ_im.val()*u*meas(G);

	if (realOrImag == 0){

        if (useDir || useDirRI) 
            u.addBc( bcInfo_Dirichlet_re.get("Dirichlet"));

	    A.initSystem();

		A.assemble(laplace_term_real);
	    A.assembleLhsRhsBc(helmholtz_term_real, zero_vec , bcInfo_Gamma_f.neumannSides()); 
		//A.assembleLhsRhsBc(-helmholtz_term_real, zero_vec , bcInfo_Gamma_b.neumannSides());
        
        if (useForce)
            A.assemble(rhs_force_real);

        if (useNeu)
        {
		    A.assembleLhsRhsBc(zero_matrix, rhs_term_real , bcInfo_Gamma_s.neumannSides());
            //gsDebugVar(A.rhs().norm());
        }

        if (useNeuBell)
        {
		    A.assembleLhsRhsBc(zero_matrix, rhs_bell_real , bcInfo_Gamma_f.neumannSides());
        }
        //gsDebugVar(A.rhs().norm());
        
	} else {
        // It is not a bug that we use bcInfo_Dirichlet_re and not ..._im here!
        // The reason is that when we construct the full real system the rhs gets contributions from both combinations of A_real, A_imag and bcInfo..._re and bcInfo..._im
        // But if we assume that u_imag = 0 then this reduces to
        // | Ar  -Ai | | ur | = | fr - Ar_bnd ur_bnd |
        // | -Ai -Ar | | ui | = | -fi + Ai_bnd ur_bnd |
        //
        // Both rhs we have ur_bnd on the right!
        if (useDir) u.addBc( bcInfo_Dirichlet_re.get("Dirichlet")); 
        if (useDirRI) u.addBc( bcInfo_Dirichlet_im.get("Dirichlet")); 


	    A.initSystem();
        //A.assemble(u*u.tr());
		A.assemble(laplace_term_imag); 
        // FIXIT: If we use abs(detJ) we dont need to assemble this term:
		A.assembleLhsRhsBc(helmholtz_term_imag, zero_vec , bcInfo_Gamma_f.neumannSides()); 

        if (useForce)
            A.assemble(rhs_force_imag);

        if (useNeu)
        {
		    A.assembleLhsRhsBc(zero_matrix, rhs_term_imag , bcInfo_Gamma_s.neumannSides());
        }

        if (useNeuBell)
        {
		    A.assembleLhsRhsBc(zero_matrix, rhs_bell_imag , bcInfo_Gamma_f.neumannSides());
        }

        // If we dont useDirRI then return rhs
    }

	mat = A.matrix();
    rhs = A.rhs();
    //gsDebugVar(rhs.rows());

    if (useDirRI)
    {
        // Assemble real part again
        
	    A.initSystem();
	    A.assemble(laplace_term_real);
	    A.assembleLhsRhsBc(helmholtz_term_real, zero_vec , bcInfo_Gamma_f.neumannSides()); 

        rhs = A.rhs();
        //gsDebugVar(rhs.rows());

    }

}


// Methods for calculating the derivative
gsMatrix<> gsStateEquationPotWaves::getDerivativeOfAu(index_t realOrImag, gsMultiPatch<> sol){
    // gsInfo << "f = " << *f << "\n" << std::flush;
    // gsInfo << "ms = " << *ms << "\n" << std::flush;

    gsFunctionExpr<> zero("0.0",3);
    //! [Boundary conditions]

    gsExprAssembler<> A(1,1);
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    A.options().setReal("quA",m_quA);
    A.options().setInt("quB",m_quB);

    geometryMap G = A.getMap(*m_mp);

    gsExprEvaluator<> ev(A);
    A.setIntegrationElements(dbasis);

    space u = A.getSpace(dbasis);
    u.setInterfaceCont(0);

    gsMultiBasis<> geom_basis(*m_mp);

    space v = A.getTestSpace(u,geom_basis,m_dim);
    // v.setInterfaceCont(0);

    // variable ff = ev.getVariable(f,G);
    // variable dff_dx = A.getCoeff(*df_dx, G);
    // variable dff_dy = A.getCoeff(*df_dy, G);

    variable solVar = A.getCoeff(sol);
    variable zer = A.getCoeff(zero,G);

    variable dJ_re = A.getCoeff(*strech_detJ_re,G);
    variable dJ_im = A.getCoeff(*strech_detJ_im,G);

    variable detJJm1Jm1_asVec_re = A.getCoeff(*strech_detJ_Jm1_Jm1_asVec_re,G);
    variable detJJm1Jm1_asVec_im = A.getCoeff(*strech_detJ_Jm1_Jm1_asVec_im,G);

    auto dJJJ_re = detJJm1Jm1_asVec_re.asDiag();
    auto dJJJ_im = detJJm1Jm1_asVec_im.asDiag();

    auto gradSol = fjac(solVar).tr();
    auto igradSol = fjac(solVar).tr()*jac(G).inv();

    gsMultiPatch<> markInner = markInnerDomain();
    variable delta = A.getCoeff(markInner);

    // Laplace terms
    // OBS: Here it is assumed that dJJJ_re = 1 in the Domain of interest
    // OBS: Here it is assumed that dJJJ_im = 0 in the Domain of interest
    auto d_lapl_real_term_1 = matrix_by_space(jac(G).inv(),jac(v)).trace()*igradSol*(igrad(u,G).tr())*meas(G);

    auto d_lapl_real_term_2 = - collapse(fjac(solVar).tr(),matrix_by_space(jac(G).inv(),jac(v))*jac(G).inv()) * (igrad(u,G).tr())*meas(G);
    
    auto d_lapl_real_term_3 = - collapse(igradSol,matrix_by_space_tr(jac(G).inv(),jac(v)))*(igrad(u,G).tr())*meas(G);
    

	// Helmholtz Term
    auto d_helm_real_term1 = - wave_K*nvDeriv(v,G)*nv(G)/nv(G).norm()*solVar.val()*u.tr()*dJ_re.val();
    // FIXIT : if helmholtx term is complex you need to differentiate the imaginary term also :)

    A.initSystem();
    if (realOrImag == 0){
        // Real sys
        A.assemble(d_lapl_real_term_1);
        A.assemble(d_lapl_real_term_2);
        A.assemble(d_lapl_real_term_3);

	    A.assembleLhsRhsBc(d_helm_real_term1, zeroRhs(v) , bcInfo_Gamma_f.neumannSides()); 
    } else {
        // Imag sys
    }

    return A.matrix(); 

}

gsMatrix<> gsStateEquationPotWaves::getDerivativeOfRhsZeroBC(index_t realOrImag){
	// gsInfo << "getDerivativeOfRhsZeroBC\n" << std::flush;
	// Method to generate system matrices and rhs
	// Input 0 for real part of matrix, and 1 for imaginary part of matrix

	gsExprAssembler<> A(1,1);
	gsExprEvaluator<> ev(A);

    A.options().setReal("quA",m_quA);
    A.options().setInt("quB",m_quB);

	geometryMap G = A.getMap(*m_mp);

	gsFunctionExpr<> zero("0.0",m_dim);
	variable zer = A.getCoeff(zero,G);

	A.setIntegrationElements(dbasis);

    space u = A.getSpace(dbasis);
    u.setInterfaceCont(0);

	gsMultiBasis<> geom_basis(*m_mp);
    space v = A.getTestSpace(u,geom_basis,m_dim);

	// Setup terms
	// Rhs
	variable F_real = A.getCoeff(*pde_F_re,G);
	variable F_imag = A.getCoeff(*pde_F_im,G);

	variable grad_F_real_vec = A.getCoeff(*pde_dF_re,G);
	variable grad_F_imag_vec = A.getCoeff(*pde_dF_im,G);

    auto grad_F_real = reshape(grad_F_real_vec,m_dim,m_dim);
    auto grad_F_imag = reshape(grad_F_imag_vec,m_dim,m_dim);

    gsFunctionExpr<> ones_fun("1","1","1",3);
    variable ones = A.getCoeff(ones_fun);

    gsMultiPatch<> markInner = markInnerDomain();
    variable delta = A.getCoeff(markInner);

    // real terms
	auto drhs_real_term1 = (nvDeriv(v,G)*F_real)*u.tr(); 

	auto drhs_real_term2 = (v*grad_F_real*nv(G))*u.tr(); 

	//auto drhs_real_term3 = (nvDeriv(v,G)*nv(G)/nv(G).norm())* (F_real.tr()*nv(G))*u.tr(); 

    // imag terms
	auto drhs_imag_term1 = (nvDeriv(v,G)*F_imag)*u.tr();//*nv(G).norm(); 

	auto drhs_imag_term2 = (v*grad_F_imag*nv(G))*u.tr();//*nv(G).norm(); 

	//auto drhs_imag_term3 = (nvDeriv(v,G)*nv(G)/nv(G).norm())* (F_imag.tr()*nv(G))*u.tr(); 

	A.initSystem();

	if (realOrImag == 0){
		// FIXIT: why including rhs here? TODO avoid this...
		A.assembleLhsRhsBc(drhs_real_term1 + drhs_real_term2, zeroRhs(v) ,bcInfo_Gamma_s.neumannSides());
	} else {
		// FIXIT: why including rhs here? TODO avoid this...
		A.assembleLhsRhsBc(drhs_imag_term1 + drhs_imag_term2, zeroRhs(v) ,bcInfo_Gamma_s.neumannSides());
	}

	return A.matrix();

}

bool gsStateEquationPotWaves::isPatchInDomain(index_t p)
{
    gsVector< index_t > patches = getVectorWithDomainPatches();

    for (index_t i = 0; i < patches.rows(); i++)
    {
        if (patches[i] == p)
            return true;
    }
    return false;

}    

gsVector< index_t > gsStateEquationPotWaves::getVectorWithDomainPatches()
{
    gsVector< index_t > patches(5);
    patches << 0, 1, 2, 3, 15;

    return patches;
}

gsVector< index_t > gsStateEquationPotWaves::getVectorWithPMLPatches()
{
    gsVector< index_t > domPatches = getVectorWithDomainPatches();

    gsVector< index_t > pmlPatches(m_mp->nBoxes() - domPatches.size());

    index_t i = 0;
    index_t j = 0;
    for (index_t p = 0; p < m_mp->nBoxes(); p++)
    {
        if (p == domPatches[i])
        {
            i++;
            continue;
        } 
        else
        {
            pmlPatches[j] = p;
            j++;
        }
    }


    return pmlPatches;
}

void gsStateEquationPotWaves::getObjFunctions(
        gsFunctionExpr<>::uPtr &expKzpiKx_re,
        gsFunctionExpr<>::uPtr &expKzpiKx_im
        ) 
{

    complexExpr fun(" exp( K * z) * cos( K * x ) ", " exp( K * z) * sin( K * x )");

    expKzpiKx_re = fun.getFunRe(const_map);
    expKzpiKx_im = fun.getFunIm(const_map);

}

void gsStateEquationPotWaves::getObjFunctions(
        gsFunctionExpr<>::uPtr &grad_expKz_expiKx_re,
        gsFunctionExpr<>::uPtr &grad_expKz_expiKx_im,
        gsFunctionExpr<>::uPtr &hess_asVec_expKz_expiKx_re,
        gsFunctionExpr<>::uPtr &hess_asVec_expKz_expiKx_im
) 
{
    // exp(Kz) ( cos(Kx) + i sin(Kx) )
    // d/dx = K exp(Kz) ( - sin(Kx) + i cos(Kx) )
    // d/dy = 0.0
    // d/dz = K exp(Kz) ( cos(Kx) + i sin(Kx) )

    complexExpr d_fun_dx(" - K * exp( K * z) * sin( K * x ) ", " K * exp( K * z) * cos( K * x )");
    complexExpr d_fun_dy(" 0.0 ", "0.0");
    complexExpr d_fun_dz(" K * exp(K * z) * cos(K * x) ", " K * exp( K * z) * sin( K * x ) ");
    grad_expKz_expiKx_re = getVectorFunRe(const_map, d_fun_dx, d_fun_dy, d_fun_dz);
    grad_expKz_expiKx_im = getVectorFunIm(const_map, d_fun_dx, d_fun_dy, d_fun_dz);

    // d2/dx2 = K * K * exp(Kz) ( - cos(Kx) - i sin(Kx) )
    // d2/dy2 = 0.0
    // d2/dz2 = K * K * exp(Kz) ( cos(Kx) + i sin(Kx) )

    complexExpr xx(" - K * K * exp( K * z) * cos( K * x ) ", " - K * K * exp( K * z) * sin( K * x )");
    complexExpr yy(" 0.0 ", "0.0");
    complexExpr zz(" K * K * exp( K * z) * cos( K * x ) ", " K * K * exp( K * z) * sin( K * x ) ");

    complexExpr xy(" 0.0 ", "0.0");
    complexExpr xz(" - K * K * exp( K * z) * sin( K * x ) ", " K * K * exp( K * z) * cos( K * x ) ");

    complexExpr yz(" 0.0 ", " 0.0 ");

    hess_asVec_expKz_expiKx_re = getSymmVectorFunRe(const_map, xx, xy, xz, yy, yz, zz);
    hess_asVec_expKz_expiKx_im = getSymmVectorFunIm(const_map, xx, xy, xz, yy, yz, zz);


}

void gsStateEquationPotWaves::getBessel()
{
    std::stringstream stream;
    stream.precision(12);

    // --------- J0 --------- //
    index_t N = 50;
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

    std::string J0str = stream.str();

    // ----------- Y0 ------------ //
    stream.clear();
    stream.str(std::string());


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

    stream << " )";
    std::string Y0str = stream.str();

    // Multiply with exp(Kz)
    stream.clear();
    stream.str(std::string());
    stream << "exp( " << wave_K << " * z ) *(" << J0str << ")";
    std::string J0exp_str = stream.str();

    stream.clear();
    stream.str(std::string());
    stream << "-exp( " << wave_K << " * z ) *(" << Y0str << ")";
    std::string mY0exp_str = stream.str();

    pde_bessel_J0exp = memory::make_unique( new gsFunctionExpr<>(J0exp_str, 3) );
    pde_bessel_mY0exp = memory::make_unique( new gsFunctionExpr<>(mY0exp_str, 3) );

}

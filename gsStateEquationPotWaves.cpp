#include <gismo.h>
#include "gsStateEquationPotWaves.h"
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

// ============================ //
//       Member functions       //
// ============================ //

gsStateEquationPotWaves::gsStateEquationPotWaves(index_t numRefine, real_t Lx, real_t Ly, real_t Lz, real_t lx, real_t ly, real_t lz): 
    zero("0.0",3), 
    m_numRefine(numRefine)
{
    pml_Lx = Lx;
    pml_Ly = Ly;
    pml_Lz = Lz;
    pml_lx = lx;
    pml_ly = ly;
    pml_lz = lz;
    constructor();
}

gsStateEquationPotWaves::gsStateEquationPotWaves(index_t numRefine): 
    zero("0.0",3), 
    m_numRefine(numRefine)
{
    constructor();
}

void gsStateEquationPotWaves::constructor()
{
    m_mp = getInitialDomain();

    complexExpr u_expr("0.0","0.0");
    pde_u_dir_re = u_expr.getFunRe(const_map);
    pde_u_dir_im = u_expr.getFunIm(const_map);

    setup();

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

    // Streching functions for PML
    gsInfo << "\n Pml : \n";

    std::string dsx_expr_im = 
        "switch                                                             "
        "{                                                                  "
        "   case (x < - lx)  : - n*C/(omega * (Lx - lx)^n)*(- x - lx)^(n - 1);   "
        "   case (x >   lx)  : n*C/(omega * (Lx - lx)^n)*(x - lx)^(n - 1);     "
        "   default          : 0;                                           "
        "}                                                                  ";
        
    std::string dsy_expr_im = 
        "switch                                                             "
        "{                                                                  "
        "   case (y < - ly)  : - n*C/(omega * (Ly - ly)^n)*(- y - ly)^(n - 1);   "
        "   case (y >   ly)  : n*C/(omega * (Ly - ly)^n)*(y - ly)^(n - 1);     "
        "   default          : 0;                                           "
        "}                                                                  ";

    std::string dsz_expr_im = 
        "switch                                                             "
        "{                                                                  "
        "   case (z < - lz)  : - n*C/(omega * (Lz - lz)^n)*(- z - lz)^(n - 1);   "
        "   default          : 0;                                           "
        "}                                                                  ";
    
    complexExpr dstrech_x("1",dsx_expr_im);
    complexExpr dstrech_y("1",dsy_expr_im);
    complexExpr dstrech_z("1",dsz_expr_im);

    complexExpr detJ = dstrech_x*dstrech_y*dstrech_z;
    complexExpr detJlen = detJ.norm();

    strech_detJ_len_re = detJlen.getFunRe(const_map);
    strech_detJ_len_im = detJlen.getFunIm(const_map);

    complexExpr dstrech_x_inv = dstrech_x.inv();
    complexExpr dstrech_y_inv = dstrech_y.inv();
    complexExpr dstrech_z_inv = dstrech_z.inv();

    // detJ * Jm1 * Jm1
    complexExpr xcomp = detJ*dstrech_x_inv.sqNorm();//*dstrech_x_inv.adj();
    complexExpr ycomp = detJ*dstrech_y_inv.sqNorm();//*dstrech_y_inv.adj();
    complexExpr zcomp = detJ*dstrech_z_inv.sqNorm();//*dstrech_z_inv.adj();

    strech_detJ_Jm1_Jm1_asVec_re = getVectorFunRe(const_map, xcomp, ycomp, zcomp);
    strech_detJ_Jm1_Jm1_asVec_im = getVectorFunIm(const_map, xcomp, ycomp, zcomp);



}

void gsStateEquationPotWaves::setup()
{
    gsMultiBasis<> bas(*m_mp);

    dbasis = bas;

    while(dbasis.maxCwiseDegree() < degree)
    {
        dbasis.degreeIncrease();
    }

    // gsInfo << "OBS CHECK DISCRETIZATION OF PDE!!\n";
    for (index_t i = 0; i < m_numRefine; i++){
        dbasis.uniformRefine();
    }

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
            gsInfo << "GAMMA_B: Add condition at patch " << ps.patch << " bnd " << ps.index() <<"\n";
            bcInfo_Gamma_b.addCondition(ps.patch, ps.index(), condition_type::neumann, zero);
        }

        // Find Gamma_s
        // Goal is to find the boundaries inside the PML domain (ie. not on the PML layers)
        if (checkForPML(boundaryDofs,m_mp->patch(ps.patch).coefs(),init_lbx,init_lby,init_lbz))
        {
            gsInfo << "GAMMA_S: Add condition at patch " << ps.patch << " bnd " << ps.index() << "\n";
            bcInfo_Gamma_s.addCondition(ps.patch, ps.index(), condition_type::neumann, zero);
        }
    }

    gsInfo << "bcInfo_Gamma_f : \n" << bcInfo_Gamma_f << "\n";
    gsInfo << "bcInfo_Gamma_b \n" << bcInfo_Gamma_b << "\n";
    gsInfo << "bcInfo_Gamma_s : \n" << bcInfo_Gamma_s << "\n";

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
            
            gsInfo << "GAMMA_D: Add condition at patch TRE" << ps.patch << " bnd " << ps.index() << "\n";
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
            gsInfo << "GAMMA_D: Add condition at patch " << ps.patch << " bnd " << ps.index() << "\n";
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
        bool cond2 = !(checkForPML(boundaryDofs,m_mp->patch(ps.patch).coefs(),init_lbx,init_lby,init_lbz));
        if (cond1 and cond2)
        {
            gsInfo << "GAMMA_D: Add condition at patch " << ps.patch << " bnd " << ps.index() << "\n";
            bcInfo_Dirichlet_re.addCondition(ps.patch, ps.index(), condition_type::dirichlet, *pde_u_dir_re);
            bcInfo_Dirichlet_im.addCondition(ps.patch, ps.index(), condition_type::dirichlet, *pde_u_dir_im);
        }
    }

    gsInfo << "bcInfo_Gamma_f : \n" << bcInfo_Gamma_f << "\n";
    gsInfo << "bcInfo_Gamma_s : \n" << bcInfo_Gamma_s << "\n";
    gsInfo << "bcInfo_Gamma_d : \n" << bcInfo_Dirichlet_re << "\n";
    gsInfo << "bcInfo_Gamma_d : \n" << bcInfo_Dirichlet_im << "\n";

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

        case  1: 
                out = getCoefs(0); 
                out.col(0) = -out.col(0);
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

    if (includeCenter)\
        nPatches = 16;

	gsMultiPatch<>::Ptr out = memory::make_shared( new gsMultiPatch<> );
    for (index_t p = 0; p < nPatches; p++)
    {
        gsMatrix<> coefs = getCoefs(p);
        gsTensorBSpline<3, real_t> tbs = getBox(coefs);
        out->addPatch( tbs );
    }

    out->computeTopology();

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

    variable dJ_le_re = A.getCoeff(*strech_detJ_len_re,G);
    variable dJ_le_im = A.getCoeff(*strech_detJ_len_im,G);

	// Helmholtz Term
	auto helmholtz_term_real = - wave_K*u*u.tr()*dJ_le_re.val()*nv(G).norm();
	auto helmholtz_term_imag = - wave_K*u*u.tr()*dJ_le_im.val()*nv(G).norm();

	// Rhs
	variable F_real = A.getCoeff(*pde_F_re,G);
	variable F_imag = A.getCoeff(*pde_F_im,G);

	auto rhs_term_real = (F_real.tr()*(nv(G)/nv(G).norm())).val()*u*nv(G).norm();
	auto rhs_term_imag = (F_imag.tr()*(nv(G)/nv(G).norm())).val()*u*nv(G).norm();

    // Zero
    gsFunctionExpr<> zerfun("0.0",3);
    variable zer = ev.getVariable(zerfun);

    auto zero_matrix = zer.val()*u*u.tr();  // Perhaps I need to multiply with u*u.tr();
    auto zero_vec = zer.val()*u*nv(G).norm();     // Perhaps I need to multiply with u;

	if (realOrImag == 0){

        if (useDir) u.addBc( bcInfo_Dirichlet_re.get("Dirichlet"));

	    A.initSystem();

		A.assemble(laplace_term_real);
		A.assembleLhsRhsBc(helmholtz_term_real, zero_vec , bcInfo_Gamma_f.neumannSides());
		//A.assembleLhsRhsBc(-helmholtz_term_real, zero_vec , bcInfo_Gamma_b.neumannSides());

        if (useNeu)
        {
		    A.assembleLhsRhsBc(zero_matrix, rhs_term_real , bcInfo_Gamma_s.neumannSides());
        }
	} else {
        // It is not a bug that we use bcInfo_Dirichlet_re and not ..._im here!
        // The reason is that when we construct the full real system the rhs gets contributions from both combinations of A_real, A_imag and bcInfo..._re and bcInfo..._im
        // But if we assume that u_imag = 0 then this reduces to
        // | Ar  -Ai | | ur | = | fr - Ar_bnd ur_bnd |
        // | -Ai -Ar | | ui | = | -fi + Ai_bnd ur_bnd |
        //
        // Both rhs we have ur_bnd on the right!
        if (useDir) u.addBc( bcInfo_Dirichlet_re.get("Dirichlet")); 

	    A.initSystem();
		A.assemble(laplace_term_imag);
		A.assembleLhsRhsBc(helmholtz_term_imag, zero_vec , bcInfo_Gamma_f.neumannSides());

        if (useNeu)
        {
		    A.assembleLhsRhsBc(zero_matrix, rhs_term_imag , bcInfo_Gamma_s.neumannSides());
        }
	}

	mat = A.matrix();
	rhs = A.rhs();

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

void gsStateEquationPotWaves::plotVelocityField(gsMultiPatch<> &ur, gsMultiPatch<> &ui, real_t timestep, std::string outfile)
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

    index_t i = 0;
    real_t t = 0;

    //ev.options().setSwitch("plot.elements", true);
    ev.options().setInt("plot.npts", 8000);

    index_t maxIter = 1000;

    while( t < 2*M_PI)
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
        ev.writeParaview( velocity_field_I + velocity_field_S    , G, name);

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

        variable u_e = A.getCoeff(*pde_u_dir_re,G);
        variable u_r = A.getCoeff(ur);
        variable u_i = A.getCoeff(ui);

        variable delta = A.getCoeff(markInner);

        ev.options().setInt("quB",m_quB+2); // FIXIT lower no points
        ev.options().setReal("quA", m_quA+2); // FIXIT lower no points

        real_t errorSq = ev.integral( delta*(u_e - u_r).sqNorm()*meas(G));
        gsInfo << "size : \t " << dbasis.size() << " , Error \t \t: " << sqrt(errorSq) << "\n";

        real_t h = 0;
        for (index_t p = 0; p < m_mp->nBoxes(); p++)
        {
            real_t hp = static_cast<gsTensorBSplineBasis<3>&> (dbasis.basis(p)).getMaxCellLength();
            if (hp > h)
                h = hp;

        }

        f << dbasis.size() << " " << h << " " << sqrt(errorSq) << "\n";
        
        if ( i == max_refine-1) // If plot
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
        
            name = "uexact";
            gsInfo<< "Plotting " << name << " in Paraview...\n";
            ev.writeParaview( u_e , G, outfolder + name);
        
            name = "diff";
            gsInfo<< "Plotting " << name << " in Paraview...\n";
            ev.writeParaview( delta*(u_e - u_r)  , G, outfolder + name);
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

        gsInfo << "Norm of solVector_imag : " << solVector_Imag.norm() << "\n";
        gsInfo << "Norm of u_sol : " << ev.integral(u_sol_imag.sqNorm()) << "\n";

        gsExprAssembler<>::variable utest = ev.getVariable(u_imag);
        ev.options().setInt("plot.npts", 8000);
        gsInfo << "Norm of utest : " << ev.integral(utest) << "\n";

        ev.writeParaview( utest, G, "test");




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
    gsVector< index_t > patches(5);
    patches << 0, 1, 2, 3, m_mp->nBoxes()-1;

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

#include <gismo.h>
#include <fstream>
#include <string>
#include "gsShapeOptLog.h"
using namespace gismo;

gsShapeOptLog& gsShapeOptLog::operator<< (std::string &str)
{
    std::ofstream f;
    f.open(BASE_FOLDER + m_output + m_logName, std::ofstream::out | std::ofstream::app);
    f << str;
    f.close();
    return *this;
}

gsShapeOptLog& gsShapeOptLog::operator<< (const char str[])
{
    std::ofstream f;
    f.open(BASE_FOLDER + m_output + m_logName, std::ofstream::out | std::ofstream::app);
    f << str;
    f.close();
    return *this;
}

gsShapeOptLog& gsShapeOptLog::operator<< (index_t i)
{
    std::ofstream f;
    f.open(BASE_FOLDER + m_output + m_logName, std::ofstream::out | std::ofstream::app);
    f << i;
    f.close();
    return *this;
}

gsShapeOptLog& gsShapeOptLog::operator<< (real_t r)
{
    std::ofstream f;
    f.open(BASE_FOLDER + m_output + m_logName, std::ofstream::out | std::ofstream::app);
    f << r;
    f.close();
    return *this;
}

void gsShapeOptLog::logObj(real_t obj)
{
    std::ofstream f;
    f.open(BASE_FOLDER + m_output + m_logObjName, std::ofstream::out | std::ofstream::app);
    f << obj << "\n";
    f.close();
}

void gsShapeOptLog::logObj(gsVector<> in)
{
    std::ofstream f;
    f.open(BASE_FOLDER + m_output + m_logObjName, std::ofstream::out | std::ofstream::app);
    for (real_t x:in)
        f << x << " ";

    f << "\n";
    f.close();
}

void gsShapeOptLog::resetLog()
{
    std::ofstream f(BASE_FOLDER + m_output + m_logName);
    f << "=== Log === \n\n";
    f.close();
    std::ofstream f2(BASE_FOLDER + m_output + m_logObjName);
    f2.close();
}

void gsShapeOptLog::saveVec(gsVector<> const &vec, std::string &name)
{
    std::ofstream f(BASE_FOLDER + m_output + name + ".txt");
    for (auto &e : vec) f << std::setprecision(12) << e << "\n";
}

void gsShapeOptLog::saveVec(gsVector<> const &vec, std::string &name, index_t i)
{
    std::stringstream stream;
    stream << name << "_" << i;
    std::string str = stream.str();
    saveVec(vec,str);
}

void gsShapeOptLog::saveVec(gsVector<> const &vec, std::string &name, index_t i, index_t j)
{
    std::stringstream stream;
    stream << name << "_" << i << "_" << j;
    std::string str = stream.str();
    saveVec(vec,str);
}

void gsShapeOptLog::plotInParaview(gsMultiPatch<> &mp, std::string &name)
{
    gsWriteParaview(mp, BASE_FOLDER + m_output + name,10000,true);
}

void gsShapeOptLog::plotInParaview(gsMultiPatch<> &mp, std::string &name, index_t i)
{
    std::stringstream stream;
    stream << name << "_" << i;
    std::string str = stream.str();
    plotInParaview(mp,str);
}

void gsShapeOptLog::plotInParaview(gsMultiPatch<> &mp, std::string &name, index_t i, index_t j)
{
    std::stringstream stream;
    stream << name << "_" << i << "_" << j;
    std::string str = stream.str();
    plotInParaview(mp,str);
}

void gsShapeOptLog::plotMultiPatchOnGeometry(gsMultiPatch<> &mp, gsMultiPatch<> &fun, std::string &name)
{
    // gsInfo << "getDvectors\n" << std::flush;
    // gsInfo<<"Active options:\n"<< A.options() <<"\n";
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    gsExprAssembler<> A(1,1);

    // Elements used for numerical integration
    gsMultiBasis<> bas(fun);
    A.setIntegrationElements(bas);
    gsExprEvaluator<> ev(A);

    geometryMap G = A.getMap(mp);

    variable out = A.getCoeff(fun);

    gsInfo<<"Plotting " << BASE_FOLDER + m_output + name << " in Paraview...\n";
    ev.writeParaview( out   , G, BASE_FOLDER + m_output + name);
    ev.options().setSwitch("plot.elements", true);


}

void gsShapeOptLog::plotMultiPatchOnGeometry(gsMultiPatch<> &mp,  gsMultiPatch<> &fun, std::string &name, index_t i)
{
    std::stringstream stream;
    stream << name << "_" << i;
    std::string str = stream.str();
    plotMultiPatchOnGeometry(mp,fun,str);
}

void gsShapeOptLog::plotSignMultiPatchOnGeometry(gsMultiPatch<> &mp, gsMultiPatch<> fun, std::string &name)
{
    // gsInfo << "getDvectors\n" << std::flush;
    // gsInfo<<"Active options:\n"<< A.options() <<"\n";
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    gsExprAssembler<> A(1,1);

    // Elements used for numerical integration
    gsMultiBasis<> bas(fun);
    A.setIntegrationElements(bas);
    gsExprEvaluator<> ev(A);

    geometryMap G = A.getMap(mp);

    for (index_t p = 0; p < fun.nBoxes(); p++)
    {
        for (index_t i = 0; i < fun.patch(p).coefsSize(); i++)
        {
            for (index_t d = 0; d < fun.targetDim(); d++)
            {
                real_t x = fun.patch(p).coef(i,d);
                fun.patch(p).coef(i,d) = (x > 0) ? 1 : ((x < 0) ? -1 : 0);
            }
        }
    }

    variable out = A.getCoeff(fun);

    gsInfo<<"Plotting " << BASE_FOLDER + m_output + name << " in Paraview...\n";
    ev.writeParaview( out   , G, BASE_FOLDER + m_output + name);
    ev.options().setSwitch("plot.elements", true);
}

void gsShapeOptLog::plotSignMultiPatchOnGeometry(gsMultiPatch<> &mp,  gsMultiPatch<> fun, std::string &name, index_t i)
{
    std::stringstream stream;
    stream << name << "_" << i;
    std::string str = stream.str();
    plotSignMultiPatchOnGeometry(mp,fun,str);
}

void gsShapeOptLog::plotActiveMultiPatchOnGeometry(gsMultiPatch<> &mp, gsMultiPatch<> fun, real_t val, std::string &name)
{
    // gsInfo << "getDvectors\n" << std::flush;
    // gsInfo<<"Active options:\n"<< A.options() <<"\n";
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    gsExprAssembler<> A(1,1);

    // Elements used for numerical integration
    gsMultiBasis<> bas(fun);
    A.setIntegrationElements(bas);
    gsExprEvaluator<> ev(A);

    geometryMap G = A.getMap(mp);

    for (index_t p = 0; p < fun.nBoxes(); p++)
    {
        for (index_t i = 0; i < fun.patch(p).coefsSize(); i++)
        {
            for (index_t d = 0; d < fun.targetDim(); d++)
            {
                real_t x = fun.patch(p).coef(i,d);
                // gsInfo << "x = " << x << " , val = " << val << "\n";
                fun.patch(p).coef(i,d) = (x > -val - 0.0001) ? 1 : 0;
            }
        }
    }

    variable out = A.getCoeff(fun);

    gsInfo<<"Plotting " << BASE_FOLDER + m_output + name << " in Paraview...\n";
    ev.writeParaview( out   , G, BASE_FOLDER + m_output + name);
    ev.options().setSwitch("plot.elements", true);
}

void gsShapeOptLog::plotActiveMultiPatchOnGeometry(gsMultiPatch<> &mp,  gsMultiPatch<> fun, real_t val, std::string &name, index_t i)
{
    std::stringstream stream;
    stream << name << "_" << i;
    std::string str = stream.str();
    plotActiveMultiPatchOnGeometry(mp,fun,val,str);
}

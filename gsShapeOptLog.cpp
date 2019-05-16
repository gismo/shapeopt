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

void gsShapeOptLog::resetLog()
{
    std::ofstream f(BASE_FOLDER + m_output + m_logName);
    f << "=== Log === \n\n";
    f.close();
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
    gsWriteParaview(mp, BASE_FOLDER + m_output + name );
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

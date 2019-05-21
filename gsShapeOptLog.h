/** @file gsShapeOptLog

@brief  Implements logging framework, ie. saving stuff during optimization,
        for doing shape optimization

This file is part of the G+Smo library.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Author(s): A. Limkilde, A. Mantzaflaris
*/
#ifndef GSSHAPEOPTLOG_H
#define GSSHAPEOPTLOG_H
using namespace gismo;

class gsShapeOptLog{
public:

    // Construct with default options
    // output should have the a forward slash in the end, e.g. "/../results/"
    gsShapeOptLog(std::string output): m_output(output) { resetLog(); };

    // Construct with options
    gsShapeOptLog(std::string output, bool saveCps, bool plotDesign, bool plotState):
        m_output(output), m_saveCps(saveCps), m_plotDesign(plotDesign), m_plotState(plotState){ resetLog(); };

    // Operator << is used to log text (strings, char*, index_t and real_t)
    gsShapeOptLog& operator << (std::string &str);
    gsShapeOptLog& operator << (const char str[]);
    gsShapeOptLog& operator << (index_t i);
    gsShapeOptLog& operator << (real_t r);

    void logObj(real_t obj); // Logs the objective value

    // Deletes what is already in the log
    void resetLog();

    // Save vector to file with file name
    void saveVec(gsVector<> const &vec, std::string &name);

    // Save vector to file with file name concatenated with an index
    void saveVec(gsVector<> const &vec, std::string &name, index_t i);

    // Save vector to file with file name concatenated with two indicies
    void saveVec(gsVector<> const &vec, std::string &name, index_t i, index_t j);

    // Write gsMultiPatch to Paraview
    void plotInParaview(gsMultiPatch<> &mp, std::string &name);

    // Write gsMultiPatch to ParaView, with name concatenated with an index
    void plotInParaview(gsMultiPatch<> &mp, std::string &name, index_t i);

    // Write gsMultiPatch to ParaView, with name concatenated with two indicies
    void plotInParaview(gsMultiPatch<> &mp, std::string &name, index_t i, index_t j);

    // Write gsMultiPatch to paraview, concatenated with another multipatch (that defines the geometry)
    // mp denotes geometry, fun denotes the multipatch to plot
    void plotMultiPatchOnGeometry(gsMultiPatch<> &mp, gsMultiPatch<> &fun, std::string &name);

    // Write gsMultiPatch to paraview, concatenated with another multipatch (that defines the geometry)
    // mp denotes geometry, fun denotes the multipatch to plot
    // name concatenated with two indicies
    void plotMultiPatchOnGeometry(gsMultiPatch<> &mp,  gsMultiPatch<> &fun, std::string &name, index_t i);

    // Write the sign of a gsMultiPatch to paraview, concatenated with another multipatch (that defines the geometry)
    // mp denotes geometry, fun denotes the multipatch to plot
    void plotSignMultiPatchOnGeometry(gsMultiPatch<> &mp, gsMultiPatch<> fun, std::string &name);

    // Write the sign of a gsMultiPatch to paraview, concatenated with another multipatch (that defines the geometry)
    // mp denotes geometry, fun denotes the multipatch to plot
    // name concatenated with two indicies
    void plotSignMultiPatchOnGeometry(gsMultiPatch<> &mp, gsMultiPatch<> fun, std::string &name, index_t i);

    // Write the values of a specific value of a gsMultiPatch to paraview, concatenated with another multipatch (that defines the geometry)
    // mp denotes geometry, fun denotes the multipatch to plot
    void plotActiveMultiPatchOnGeometry(gsMultiPatch<> &mp, gsMultiPatch<> fun, real_t val, std::string &name);

    // Write the sign of a gsMultiPatch to paraview, concatenated with another multipatch (that defines the geometry)
    // mp denotes geometry, fun denotes the multipatch to plot
    // name concatenated with two indicies
    void plotActiveMultiPatchOnGeometry(gsMultiPatch<> &mp, gsMultiPatch<> fun, real_t val, std::string &name, index_t i);

    // Accessors
    bool saveCps() { return m_saveCps; };
    bool plotDesign() { return m_plotDesign; };
    bool plotState() { return m_plotState; };
    std::string output() { return m_output; };

public:
    bool m_saveCps = true; // Whether to save the controlPoints each
    bool m_plotDesign = true; // Whether to plot the design
    bool m_plotState = true; // Whether to plot the state

    std::string m_output; // The folder where to output to

    std::string m_logName = "log.txt"; // The name of the log file that the operators << log to
    std::string m_logObjName = "obj.txt"; // The name of the log file for objective function


};

# endif //GSSHAPEOPTLOG_H

#include <gismo.h>
#include "gsStateEquationPotWaves.h"
using namespace gismo;

gsStateEquationPotWaves::gsStateEquationPotWaves(memory::shared_ptr<gsMultiPatch<>> mpin, index_t numRefine):
    mp(mpin), dbasis(*mp)
{
    dbasis.setDegree(degree);

    std::string n = std::to_string(pml_n);
    std::string length = std::to_string(pml_length);
    std::string start = std::to_string(pml_start);
    std::string C = std::to_string(pml_C);

    std::string K = std::to_string(wave_K);

    std::string C_over_K = C + "/" + K;
    std::string C_over_Klength = C + "/(" + K + "*" + length + ")";

    strech_x_re = *(new gsFunctionExpr<>("x",3));
    strech_y_re = *(new gsFunctionExpr<>("y",3));
    strech_z_re = *(new gsFunctionExpr<>("z",3));

    std::string strech_x_im_str = "if(x > " + start + ", " + C_over_K + "*((x - " + start + ")/" + length + ")^" + n + ",0)";
    gsInfo << strech_x_im_str << "\n";
    strech_x_im = *(new gsFunctionExpr<>(strech_x_im_str,3));

    std::string strech_y_im_str = "if(y > " + start + ", " + C_over_K + "*((y - " + start + ")/" + length + ")^" + n + ",0)";
    gsInfo << strech_y_im_str << "\n";
    strech_y_im = *(new gsFunctionExpr<>(strech_y_im_str,3));

    std::string strech_z_im_str = "if(z > " + start + ", " + C_over_K + "*((z - " + start + ")/" + length + ")^" + n + ",0)";
    gsInfo << strech_z_im_str << "\n";
    strech_z_im = *(new gsFunctionExpr<>(strech_z_im_str,3));

    dstrech_x_re = *(new gsFunctionExpr<>("1",3));
    dstrech_y_re = *(new gsFunctionExpr<>("1",3));
    dstrech_z_re = *(new gsFunctionExpr<>("1",3));

    std::string dstrech_x_im_str = "if(x > " + start + ", " + C_over_Klength + "*((x - " + start + ")/" + length + ")^(" + n + "-1),0)";
    gsInfo << dstrech_x_im_str << "\n";
    dstrech_x_im = *(new gsFunctionExpr<>(dstrech_x_im_str,3));

    std::string dstrech_y_im_str = "if(y > " + start + ", " + C_over_Klength + "*((y - " + start + ")/" + length + ")^(" + n + "-1),0)";
    gsInfo << dstrech_y_im_str << "\n";
    dstrech_y_im = *(new gsFunctionExpr<>(dstrech_y_im_str,3));

    std::string dstrech_z_im_str = "if(z > " + start + ", " + C_over_Klength + "*((z - " + start + ")/" + length + ")^(" + n + "-1),0)";
    gsInfo << dstrech_z_im_str << "\n";
    dstrech_z_im = *(new gsFunctionExpr<>(dstrech_z_im_str,3));

    std::string dstrech_xinv_im_str = "if(x > " + start + ", " + C_over_Klength + "*((x - " + start + ")/" + length + ")^(" + n + "-1),0)";
    gsInfo << dstrech_xinv_im_str << "\n";
    dstrech_xinv_im = *(new gsFunctionExpr<>(dstrech_xinv_im_str,3));

    std::string dstrech_yinv_im_str = "if(y > " + start + ", " + C_over_Klength + "*((y - " + start + ")/" + length + ")^(" + n + "-1),0)";
    gsInfo << dstrech_yinv_im_str << "\n";
    dstrech_yinv_im = *(new gsFunctionExpr<>(dstrech_yinv_im_str,3));

    std::string dstrech_zinv_im_str = "if(z > " + start + ", " + C_over_Klength + "*((z - " + start + ")/" + length + ")^(" + n + "-1),0)";
    gsInfo << dstrech_zinv_im_str << "\n";
    dstrech_z_im = *(new gsFunctionExpr<>(dstrech_z_im_str,3));

}

gsMultiPatch<> gsStateEquationPotWaves::solve()
{

}

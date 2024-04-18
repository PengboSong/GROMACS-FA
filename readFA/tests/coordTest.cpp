/*
    coordTest.cpp
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2021/08/26
    Description: Test coordinates related functions.
*/

#include "real.h"
#include "rvec.h"
#include "rmat.h"
#include "rmsd.h"
#include "FAmath.h"

int main()
{
    std::vector<std::string> templates;
    templates.resize(23, "");
    templates[0]  = "ATOM      1  N   ALA     1    %8.3f%8.3f%8.3f  1.00  0.00           N\n";
    templates[1]  = "ATOM      2  H1  ALA     1    %8.3f%8.3f%8.3f  1.00  0.00           H\n";
    templates[2]  = "ATOM      3  H2  ALA     1    %8.3f%8.3f%8.3f  1.00  0.00           H\n";
    templates[3]  = "ATOM      4  H3  ALA     1    %8.3f%8.3f%8.3f  1.00  0.00           H\n";
    templates[4]  = "ATOM      5  CA  ALA     1    %8.3f%8.3f%8.3f  1.00  0.00           C\n";
    templates[5]  = "ATOM      6  HA  ALA     1    %8.3f%8.3f%8.3f  1.00  0.00           H\n";
    templates[6]  = "ATOM      7  CB  ALA     1    %8.3f%8.3f%8.3f  1.00  0.00           C\n";
    templates[7]  = "ATOM      8  HB1 ALA     1    %8.3f%8.3f%8.3f  1.00  0.00           H\n";
    templates[8]  = "ATOM      9  HB2 ALA     1    %8.3f%8.3f%8.3f  1.00  0.00           H\n";
    templates[9]  = "ATOM     10  HB3 ALA     1    %8.3f%8.3f%8.3f  1.00  0.00           H\n";
    templates[10] = "ATOM     11  C   ALA     1    %8.3f%8.3f%8.3f  1.00  0.00           C\n";
    templates[11] = "ATOM     12  O   ALA     1    %8.3f%8.3f%8.3f  1.00  0.00           O\n";
    templates[12] = "ATOM     13  N   ALA     2    %8.3f%8.3f%8.3f  1.00  0.00           N\n";
    templates[13] = "ATOM     14  H   ALA     2    %8.3f%8.3f%8.3f  1.00  0.00           H\n";
    templates[14] = "ATOM     15  CA  ALA     2    %8.3f%8.3f%8.3f  1.00  0.00           C\n";
    templates[15] = "ATOM     16  HA  ALA     2    %8.3f%8.3f%8.3f  1.00  0.00           H\n";
    templates[16] = "ATOM     17  CB  ALA     2    %8.3f%8.3f%8.3f  1.00  0.00           C\n";
    templates[17] = "ATOM     18  HB1 ALA     2    %8.3f%8.3f%8.3f  1.00  0.00           H\n";
    templates[18] = "ATOM     19  HB2 ALA     2    %8.3f%8.3f%8.3f  1.00  0.00           H\n";
    templates[19] = "ATOM     20  HB3 ALA     2    %8.3f%8.3f%8.3f  1.00  0.00           H\n";
    templates[20] = "ATOM     21  C   ALA     2    %8.3f%8.3f%8.3f  1.00  0.00           C\n";
    templates[21] = "ATOM     22  O1  ALA     2    %8.3f%8.3f%8.3f  1.00  0.00           O\n";
    templates[22] = "ATOM     23  O2  ALA     2    %8.3f%8.3f%8.3f  1.00  0.00           O\n";

    DVecVec3 coordmat;
    initmat(coordmat, 3, 23);
    coordmat[0][0]  = {15.570, 14.110, 12.960};
    coordmat[0][1]  = {15.630, 13.120, 12.770};
    coordmat[0][2]  = {15.010, 14.610, 12.300};
    coordmat[0][3]  = {16.550, 14.350, 12.830};
    coordmat[0][4]  = {15.110, 14.320, 14.340};
    coordmat[0][5]  = {15.500, 13.490, 14.940};
    coordmat[0][6]  = {15.370, 15.760, 14.780};
    coordmat[0][7]  = {16.420, 16.050, 14.790};
    coordmat[0][8]  = {14.880, 16.500, 14.140};
    coordmat[0][9]  = {14.980, 15.820, 15.790};
    coordmat[0][10] = {13.630, 13.950, 14.470};
    coordmat[0][11] = {13.210, 13.900, 15.620};
    coordmat[0][12] = {12.850, 13.740, 13.410};
    coordmat[0][13] = {13.000, 14.140, 12.500};
    coordmat[0][14] = {11.620, 12.970, 13.440};
    coordmat[0][15] = {11.040, 13.040, 14.360};
    coordmat[0][16] = {10.810, 13.540, 12.280};
    coordmat[0][17] = {10.280, 14.450, 12.570};
    coordmat[0][18] = {11.470, 13.730, 11.440};
    coordmat[0][19] = { 9.950, 12.920, 12.020};
    coordmat[0][20] = {11.880, 11.480, 13.270};
    coordmat[0][21] = {11.770, 10.850, 14.340};
    coordmat[0][22] = {12.130, 11.040, 12.130};
    
    coordmat[1][0]  = {13.886, 15.442, 11.423};
    coordmat[1][1]  = {13.554, 14.515, 11.197};
    coordmat[1][2]  = {13.165, 16.131, 11.264};
    coordmat[1][3]  = {14.645, 15.526, 10.763};
    coordmat[1][4]  = {14.368, 15.486, 12.812};
    coordmat[1][5]  = {15.103, 14.685, 12.742};
    coordmat[1][6]  = {15.149, 16.744, 13.194};
    coordmat[1][7]  = {15.950, 16.932, 12.479};
    coordmat[1][8]  = {14.526, 17.635, 13.275};
    coordmat[1][9]  = {15.637, 16.458, 14.125};
    coordmat[1][10] = {13.473, 14.911, 13.900};
    coordmat[1][11] = {13.510, 15.456, 15.002};
    coordmat[1][12] = {12.634, 13.935, 13.549};
    coordmat[1][13] = {12.669, 13.459, 12.659};
    coordmat[1][14] = {11.709, 13.297, 14.464};
    coordmat[1][15] = {11.395, 13.947, 15.281};
    coordmat[1][16] = {10.383, 13.041, 13.746};
    coordmat[1][17] = {10.473, 12.817, 12.683};
    coordmat[1][18] = { 9.756, 12.339, 14.297};
    coordmat[1][19] = { 9.786, 13.944, 13.873};
    coordmat[1][20] = {12.240, 11.938, 14.896};
    coordmat[1][21] = {11.898, 11.374, 15.958};
    coordmat[1][22] = {12.988, 11.314, 14.113};
    
    coordmat[2][0]  = {12.432, 16.178, 12.330};
    coordmat[2][1]  = {11.480, 15.868, 12.204};
    coordmat[2][2]  = {12.479, 16.548, 13.268};
    coordmat[2][3]  = {12.735, 16.873, 11.662};
    coordmat[2][4]  = {13.415, 15.084, 12.338};
    coordmat[2][5]  = {13.437, 14.701, 11.318};
    coordmat[2][6]  = {14.795, 15.643, 12.687};
    coordmat[2][7]  = {14.945, 16.649, 12.294};
    coordmat[2][8]  = {15.009, 15.585, 13.754};
    coordmat[2][9]  = {15.605, 15.113, 12.186};
    coordmat[2][10] = {12.823, 14.065, 13.301};
    coordmat[2][11] = {11.860, 14.320, 14.021};
    coordmat[2][12] = {13.493, 12.925, 13.483};
    coordmat[2][13] = {14.240, 12.699, 12.842};
    coordmat[2][14] = {13.102, 11.947, 14.479};
    coordmat[2][15] = {12.018, 11.873, 14.403};
    coordmat[2][16] = {13.910, 10.674, 14.225};
    coordmat[2][17] = {14.980, 10.799, 14.386};
    coordmat[2][18] = {13.575,  9.875, 14.886};
    coordmat[2][19] = {13.797, 10.319, 13.200};
    coordmat[2][20] = {13.411, 12.380, 15.906};
    coordmat[2][21] = {12.764, 11.761, 16.778};
    coordmat[2][22] = {14.193, 13.347, 16.027};

    VecVec3& r0 = coordmat[0];
    real rmsd = 0.;
    rvec3 r0center, rrcenter;
    rvec4 lambda;
    rmat3 U;
    rmat4 q;
    VecVec3 g;
    for (int k = 0; k < 3; ++k)
    {
        VecVec3& rr = coordmat[k];
        if (CalRMSD(rmsd, r0, rr, r0center, rrcenter, lambda, q, true, U, false, g) != -1)
        {
            VecVec3 rrfit;
            rotate(U, rrcenter, r0center, rr, rrfit);
            std::cout << "Frame " << k << std::endl;
            std::cout << "RMSD = " << rmsd << std::endl;
            std::cout << "r0 center = " << r0center << std::endl;
            std::cout << "rr center = " << rrcenter << std::endl;
            std::cout << "lambda = " << lambda << std::endl;
            std::cout << "q = " << q << std::endl;
            std::cout << "U = " << U << std::endl;
            std::cout << "Fit coord = " << std::endl;
            for (int i = 0; i < 23; ++i)
                printf(templates[i].c_str(), rrfit[i][0], rrfit[i][1], rrfit[i][2]);
        }
    }
    
    return 0;
}

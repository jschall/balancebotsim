#include <Eigen/Dense>
#include <iostream>
#include <math.h>
#include "output/balanceboteqns.h"

#define SQ(x) ((x)*(x))

using namespace Eigen;
using namespace std;

static Matrix<double,19,1> states;
static Matrix<double,2,1> inputs;

int main() {
    states << 1, 0, 0, 0, // quat
              0, 0, -1, // pos
              0, // pole_theta
              0, 0, // wheel pos
              0, 0, 0, // ang vel
              0, 0, 0, // vel
              0, // pole_omega
              0, 0; //wheel ang vel
    inputs << 0, 0; // wheel torques

    while (true) {
        for (int i=0; i<1000; i++) {
            double dt = 0.00001;

            Matrix<double,19,19> mm = get_mm(states, inputs);
            Matrix<double,19,1> fo = get_fo(states, inputs);
            Matrix<double,19,1> xdot = mm.llt().solve(fo);
            states += xdot * dt;
            double norm = sqrt(SQ(states(0,0))+SQ(states(1,0))+SQ(states(2,0))+SQ(states(3,0)));
            states(0,0) /= norm;
            states(1,0) /= norm;
            states(2,0) /= norm;
            states(3,0) /= norm;
        }

        double qr = states(0,0);
        double qi = states(1,0);
        double qj = states(2,0);
        double qk = states(3,0);

        double norm = sqrt(SQ(states(0,0))+SQ(states(1,0))+SQ(states(2,0))+SQ(states(3,0)));

        double pitch = 180*asin(2*(qr*qj-qk*qi))/M_PI;

        cout << "       pitch: " << pitch << endl;

        cout << "       quat0: " << states(0,0) << endl
             << "       quat1: " << states(1,0) << endl
             << "       quat2: " << states(2,0) << endl
             << "       quat3: " << states(3,0) << endl
             << "        pos0: " << states(4,0) << endl
             << "        pos1: " << states(5,0) << endl
             << "        pos2: " << states(6,0) << endl
             << "  pole_theta: " << states(7,0) << endl
             << "  lwheel_pos: " << states(8,0) << endl
             << "  rwheel_pos: " << states(9,0) << endl
             << "      omega0: " << states(10,0) << endl
             << "      omega1: " << states(11,0) << endl
             << "      omega2: " << states(12,0) << endl
             << "        vel0: " << states(13,0) << endl
             << "        vel1: " << states(14,0) << endl
             << "        vel2: " << states(15,0) << endl
             << "  pole_omega: " << states(16,0) << endl
             << "lwheel_omega: " << states(17,0) << endl
             << "rwhell_omega: " << states(18,0) << endl << endl;
    }

//         cout << "fo:" << endl << fo << endl << endl << "mm:" << endl << mm << endl;


    return 0;
}

//
//  torus.cpp
//
//
//  Created by Protoss Probe on 2017/04/09.
//  Copyright © 2016-2017年 probe. All rights reserved.
//

#include "poly_run.hpp"
#include "poly_grav.hpp"
#include <boost/array.hpp>
#include <cmath>
#include <iostream>
#include <string>

using namespace std;

class Particle {
  public:
    Particle() = default;
    ~Particle() = default;
    Particle(pos position) : r(position[0]), z(position[1]) {}
    Particle(state_type5 state)
        : r(state[0]), z(state[1]), r_dot(state[2]), z_dot(state[3]),
          lam_dot(state[4]) {}
    Particle(state_type5 state, double time)
        : r(state[0]), z(state[1]), r_dot(state[2]), z_dot(state[3]),
          lam_dot(state[4]), t(time) {}

    double r = 0.0, z = 0.0, r_dot = 0.0, z_dot = 0.0, lam_dot = 0.0, t = 0.0;

    state_type5 convert2state() const { return {r, z, r_dot, z_dot, lam_dot}; }
    pos convert2pos() const { return {r, z}; }
    vel convert2vel() const { return {r_dot, z_dot, lam_dot}; }
};

int main(int argc, char *argv[]) {
    string filename(argv[1]);

    PolyGrav poly(filename);

    poly.init();
    cout << poly.vert_n << '\t' << poly.face_n << '\t' << poly.edge_n << endl;
    poly.principle_axes();
    poly.export_3d_txt("assets/" + filename + "_prin.txt", 'd');

    vec3 pos;
    double x = 1.1;
    // double real = -1 / x;
    // poly.co = -1 / (4. / 3. * M_PI) * 0.5;

    const clock_t start = clock();
    pos = {{x, 0, 0}};
    // for (double i = 1.1; i < 3.0; i = i + 0.1) {
    //     pos[0] = i;
    //     cout << endl;
    //     cout << setprecision(12) << -1 / i << endl;
    //     cout << poly.potential(pos) << endl;
    //     cout << (-1 / i - poly.potential(pos)) / (-1 / i) << endl;
    // }
    double value = poly.potential(pos);
    cout << endl
         << "Cpu Time: "
         << static_cast<double>(clock() - start) / CLOCKS_PER_SEC << endl;
    cout << endl << "------------" << endl;
    // cout << "Real : " << real << endl;
    cout << "Poly Value : " << value << endl;
    // cout << "Ratio : " << real / value << endl;
    cout << "------------" << endl;
}

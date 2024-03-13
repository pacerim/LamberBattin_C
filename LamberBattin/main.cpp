#include <iostream>
#include <cmath>
#include "vector.h"
#include "seebatt.h"
#include "seebattk.h"
#include "LAMBERBATTIN.h"

#define TOL_ 10e-14

using namespace std;

int main() {

    double v[] = {1, 1, 1},
            w1[]  = {1, 2, 3},
            w2[] = {4, 5, 6} ,
            sol[3],
            ro[] = {20.0E6, 20.0E6,0},
            r[] = {-20.0E6, 10.0E6,0},
            tof = 1.0 * 86400,
            vo[3];

    int n;

    // --------------------------- TESTS NORM -------------------------------
    if(fabs(norm(v,3)- sqrt(3)) < TOL_){
        cout << "norm : TRUE" << endl;
    }else{
        cout << "norm : FALSE" << endl;
    }

    // --------------------------- TESTS DOT -------------------------------
    if(fabs(dot(w1, 3,w2 , 3)- 32) < TOL_){
        cout << "dot : TRUE" << endl;
    }else{
        cout << "dot : FALSE" << endl;
    }

    // --------------------------- TESTS CROSS -------------------------------
    cross(sol, n, w1, w2, 3);
    if(sol[0] == -3 && sol[1] == 6 && sol[2] == -3)
        cout << "cross: TRUE" << endl;
    else
        cout << "cross: FALSE" << endl;


    // --------------------------- TESTS SEEBATT -------------------------------
    if(fabs(seebatt(-1) - 2.16754395411523) < TOL_){
        cout << "seebatt : TRUE" << endl;
    }else{
        cout << "seebatt : FALSE" << endl;
    }

    if(fabs(seebatt(0) - 5) < TOL_){
        cout << "seebatt : TRUE" << endl;
    }else{
        cout << "seebatt : FALSE" << endl;
    }

    if(fabs(seebatt(1) - 6.06251330587321) < TOL_){
        cout << "seebatt : TRUE" << endl;
    }else{
        cout << "seebatt : FALSE" << endl;
    }


    // --------------------------- TESTS SEEBATTK -------------------------------
    if(fabs(seebattk(-1) - 0.391304347826087) < TOL_){
        cout << "seebattk : TRUE" << endl;
    }else{
        cout << "seebattk : FALSE" << endl;
    }

    if(fabs(seebattk(0) - 0.333333333333333) < TOL_){
        cout << "seebattk : TRUE" << endl;
    }else{
        cout << "seebattk : FALSE" << endl;
    }

    if(fabs(seebattk(1) - 0.290322580645161) < TOL_){
        cout << "seebattk : TRUE" << endl;
    }else{
        cout << "seebattk : FALSE" << endl;
    }

    /*
    LAMBERBATTIN(ro, r, "retro", tof, vo, v);
    cout << "LAMBERBATTIN : TRUE : vo" <<vo[2] << " v :"<< v[2]<< endl;
     */

    // --------------------------- TESTS LAMBERBATTIN -------------------------------
    if(fabs(LAMBERBATTIN(ro, r, "retro", tof, vo, v) - 0.391304347826087) < TOL_){
        cout << "LAMBERBATTIN : TRUE : vo" <<vo << " v :"<< v<< endl;
    }else{
        cout << "LAMBERBATTIN : FALSE" << endl;
    }

    if(fabs(LAMBERBATTIN(ro, r, "pro", tof, vo, v) - 0.333333333333333) < TOL_){
        cout << "LAMBERBATTIN : TRUE" << endl;
    }else{
        cout << "LAMBERBATTIN : FALSE" << endl;
    }

    return 0;
}
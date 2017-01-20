#ifndef N_BODY_H
#define N_BODY_H

#include <string>
#include "H5Cpp.h"
using namespace std;

//typedef float * (* integrator_func_ptr) (float * f, float * previous_fs);

class Nbody {

public:
    Nbody(int nbodies, int nt, float dt, char * output_file, string integrator_function);

    ~Nbody();

    void evolve();

    float potential_energy();

    float kinetic_energy();

    float moment_of_inertia();

private:

    void step(float * q, float * f);

    bool isinteger(float x);

    void write_to_file();

    void contain_particles();

    void print_state();

    void rk3(void (Nbody::*f)(float *, float *), float * q, float *);

    void dirk3(void (Nbody::*f)(float *, float *), float * q, float *);

    void dormand_prince(void (Nbody::*f)(float *, float *), float * q, float *);

    void adams_bashforth(void (Nbody::*f)(float *, float *), float * q, float * previous_fs);

    void adams_moulton(void (Nbody::*f)(float *, float *), float * q, float * previous_fs);

    float r();

    float G = 5.0;
    int nbodies;
    int nt;
    float dt;
    int t;
    float box_size[2];
    float * q;
    float * masses;
    float U;
    float T;
    float I;
    hid_t outFile, dset[4], mem_space[4], file_space[4];
    void (Nbody::*integrator)(void (Nbody::*f)(float *, float *), float *, float *);

};

#endif

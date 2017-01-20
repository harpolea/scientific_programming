#include "n_body.h"
#include <string>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <cmath>

#ifndef H5_NO_NAMESPACE
#ifndef H5_NO_STD
    using std::cout;
    using std::endl;
#endif  // H5_NO_STD
#endif
#include "H5Cpp.h"
using namespace std;

Nbody::Nbody(int n_bodies, int n_t, float d_t, char * output_file, string integrator_function){
    nbodies = n_bodies;
    nt = n_t;
    dt = d_t;
    t = 0;

    box_size[0] = 10.0;
    box_size[1] = 10.0;

    srand(time(0));

    q = new float[nbodies*2*2];
    masses = new float[nbodies];
    for (int n = 0; n < nbodies; n++) {
        for (int i = 0; i < 2; i++) {
            q[(n * 2)*2 + i] = r() - 0.5;
            q[(n * 2 + 1)*2 + i] = 8.0 * r() - 4.0;
        }
        masses[n] = 15.0 * (0.1 + r());
    }

    U = potential_energy();
    T = kinetic_energy();
    I = moment_of_inertia();

    // set up hdf5 file with datasets
    outFile = H5Fcreate(output_file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    int ndims = 4;
    hsize_t dims[] = {hsize_t(nt+1), hsize_t(nbodies), 2, 2};
    file_space[0] = H5Screate_simple(ndims, dims, NULL);

    hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_layout(plist, H5D_CHUNKED);
    hsize_t chunk_dims[] = {1, hsize_t(nbodies), 2, 2};
    H5Pset_chunk(plist, ndims, chunk_dims);

    // create dataset
    dset[0] = H5Dcreate(outFile, "q", H5T_NATIVE_FLOAT, file_space[0], H5P_DEFAULT, plist, H5P_DEFAULT);

    H5Pclose(plist);

    // make a memory dataspace
    mem_space[0] = H5Screate_simple(ndims, chunk_dims, NULL);

    ndims = 1;
    hsize_t dims2[] = {hsize_t(nt+1)};
    file_space[1] = H5Screate_simple(ndims, dims2, NULL);
    file_space[2] = H5Screate_simple(ndims, dims2, NULL);
    file_space[3] = H5Screate_simple(ndims, dims2, NULL);

    plist = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_layout(plist, H5D_CHUNKED);
    hsize_t chunk_dims2[] = {hsize_t(nt+1)};
    H5Pset_chunk(plist, ndims, chunk_dims2);

    // create dataset
    dset[1] = H5Dcreate(outFile, "U", H5T_NATIVE_FLOAT, file_space[1], H5P_DEFAULT, plist, H5P_DEFAULT);
    dset[2] = H5Dcreate(outFile, "T", H5T_NATIVE_FLOAT, file_space[2], H5P_DEFAULT, plist, H5P_DEFAULT);
    dset[3] = H5Dcreate(outFile, "I", H5T_NATIVE_FLOAT, file_space[3], H5P_DEFAULT, plist, H5P_DEFAULT);

    H5Pclose(plist);

    // make a memory dataspace
    mem_space[1] = H5Screate_simple(ndims, chunk_dims2, NULL);
    mem_space[2] = H5Screate_simple(ndims, chunk_dims2, NULL);
    mem_space[3] = H5Screate_simple(ndims, chunk_dims2, NULL);

    write_to_file();

    if (integrator_function == "rk3") {
        integrator = &Nbody::rk3;
    } else if (integrator_function == "dirk3") {
        integrator = &Nbody::dirk3;
    } else if (integrator_function == "dormand-prince") {
        integrator = &Nbody::dormand_prince;
    } else if (integrator_function == "adams-bashforth") {
        integrator = &Nbody::adams_bashforth;
    } else if (integrator_function == "adams-moulton") {
        integrator = &Nbody::adams_moulton;
    }
}

Nbody::~Nbody(){
    delete[] q;
    delete[] masses;
    H5Fclose(outFile);
}

float Nbody::r() {
    // generate random number between 0 and 1
    return (float) rand() / RAND_MAX;

}

void Nbody::step(float * qq, float * f){
    // step simulation through one timestep

    float * r = new float[nbodies];

    // initialise
    for (int i = 0; i < nbodies*2*2; i++) {
        f[i] = 0.0;
    }

    for (int n = 0; n < nbodies; n++) {
        for (int i = 0; i < 2; i++){
            float dist = 0.0;
            for (int m = 0; m < nbodies; m++) {
                r[m] = qq[(n*2+ 1)*2+i] - qq[(m*2+ 1)*2+i];
                dist += r[m]*r[m];
            }
            dist = sqrt(dist);
            if (dist < 1.0e-10) {
                dist = 1.0e12; //hard-core potential
            }
            for (int m = 0; m < nbodies; m++) {
                f[(m*2)*2+i] += G * masses[n] * r[m] / (dist*dist*dist);
            }
            f[(n*2+1)*2+i] = qq[(n*2)*2+i];
        }
    }

    delete[] r;
}
void Nbody::rk3(void (Nbody::*f)(float *, float *), float * qq, float *){
    // Third-order Runge-Kutta

    float * flux = new float[nbodies*2*2];
    float * qtemp = new float[nbodies*2*2];

    f(qq, flux);
    for (int i = 0; i < nbodies*2*2; i++) {
        qtemp[i] = q[i] + flux[i] * dt;
    }

    f(qtemp, flux);
    for (int i = 0; i < nbodies*2*2; i++) {
        qtemp[i] = 0.25 * (3.0 * q[i] + qtemp[i] + flux[i] * dt);
    }

    f(qtemp, flux);
    for (int i = 0; i < nbodies*2*2; i++) {
        qq[i] = (1./3.) * (q[i] + 2. * qtemp[i] + 2.0 * flux[i] * dt);
    }
    delete[] flux;
    delete[] qtemp;
}

bool Nbody::isinteger(float x){
    return numpy.equal(numpy.mod(x, 1), 0);
}

void Nbody::dirk3(void (Nbody::*f)(float *, float *), float * qq, float *){
    // Implicit diagonal third order Runge-Kutta
    // decomposed flux into linear part (0, v)^T and non-linear part
    // (sum (Gm r / |r|^3), 0)^T
    float mu = 0.5 * (1. - 1./sqrt(3.));
    float nu = 0.5 * (sqrt(3.));
    float gamma = 3. / (2. * (3. + sqrt(3.)));
    float lbda = 3. * (1. + sqrt(3.)) / (2. * (3. + sqrt(3.)));

    float * f1 = new float[nbodies*2*2];
    float * f2 = new float[nbodies*2*2];
    float * qtemp = new float[nbodies*2*2];


    for (int i = 0; i < nbodies * 2 * 2; i++) {
        qtemp[i] = q[i] + dt * mu;
    }
    f(qtemp, f1);

    for (int i = 0; i < nbodies * 2 * 2; i++) {
        qtemp[i] = q[i] + dt * (nu + 2.0 * mu);
    }
    f(qtemp, f2);

    f1 = f(self.q + self.dt * mu)[:,:,:]
    f1[:,1,:] = 0.
    f2 = f(self.q + self.dt * (nu + 2. * mu))[:,:,:]
    f2[:,1,:] = 0.
    q_new = numpy.zeros_like(self.q)

    float A[] = {{0., 0.}, {1., 0.}};
    float M[] = {{1.0, 0.0}, {-dt*mu. 1.0}};

    for (int i = 0; i < nbodies; i++){

        b1 = self.q[i,:,:] + self.dt * mu * f1[i, :, :]
        y1 = solve(M, b1)

        b2 = y1 + self.dt * nu * (A * y1 + f1[i, :, :]) + self.dt * mu * f2[i, :, :]
        y2 = solve(M, b2)

        q_new[i,:,:] = (1 - lbda) * self.q[i,:,:] + lbda * y2 + self.dt * gamma * (A * y2 + f2[i,:,:])
    }

    delete[] f1;
    delete[] f2;
    delete[] qtemp;

}
void Nbody::dormand_prince(void (Nbody::*f)(float *, float *), float * qq, float *){
    // see https://en.wikipedia.org/wiki/Dormand%E2%80%93Prince_method
    float ** k = new float*[7];
    for (int i = 0; i < 7; i++) {
        k[i] = new float[nbodies*2*2];
    }
    float * q_temp = new float[nbodies*2*2];
    float b[] = {35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0.};

    f(q, k[0]);
    for (int i = 0; i < nbodies*2*2; i++) {
        q_temp[i] = q[i] + dt * 1/5 * k[0][i];
    }
    f(q_temp, k[1]);
    for (int i = 0; i < nbodies*2*2; i++) {
        q_temp[i] = q[i] + dt * (3/40 * k[0][i] + 9/40 * k[1][i]);
    }
    f(q_temp, k[2]);
    for (int i = 0; i < nbodies*2*2; i++) {
        q_temp[i] = q[i] + dt * (44/45 * k[0][i] - 56/15 * k[1][i] + 32/9 * k[2][i]);
    }
    f(q_temp, k[3]);
    for (int i = 0; i < nbodies*2*2; i++) {
        q_temp[i] = q[i] + dt * (19372/6561 * k[0][i] - 25360/2187 * k[1][i] + 64448/6561 * k[2][i] - 212/729 * k[3][i]);
    }
    f(q_temp, k[4]);
    for (int i = 0; i < nbodies*2*2; i++) {
        q_temp[i] = q[i] + dt * (9017/3168 * k[0][i] - 355/33 * k[1][i] + 46732/5247 * k[2][i] + 49/176 * k[3][i] - 5103/18656 * k[4][i]);
    }
    f(q_temp, k[5]);
    for (int i = 0; i < nbodies*2*2; i++) {
        q_temp[i] = q[i] + dt * (35/384 * k[0][i] + 500/1113 * k[2][i] + 125/192 * k[3][i] - 2187/6784 * k[4][i] + 11/84 * k[5][i]);
    }
    f(q_temp, k[6]);


    for (int j = 0; j < nbodies*2*2; j++) {
        qq[j] = q[j];
        for (int i = 0; i < 7; i++) {
            qq[j] += dt * k[i][j] * b[i];
        }
    }

    for (int i = 0; i < 7; i++) {
        delete[] k[i];
    }
    delete[] q_temp;
}
void Nbody::adams_bashforth(void (Nbody::*f)(float *, float *), float * qq, float * previous_fs){
    // See http://www.math.iit.edu/~fass/478578_Chapter_2.pdf
    // update previous_fs
    for (int i = 0; i < 4; i++){
        for (int j = 0; j < nbodies*2*2; j++) {
            previous_fs[i*nbodies*2*2+j] = previous_fs[(i+1)*nbodies*2*2+j];
        }

    }

    float * flux = new float[nbodies*2*2];

    f(q, flux);
    for (int i = 0; i < nbodies*2*2; i++) {
        previous_fs[4*nbodies*2*2+i] = flux[i];
    }
    for (int i = 0; i < nbodies*2*2; i++) {
        qq[i] = q[i] + dt / 720.0 * (1901 * previous_fs[4*nbodies*2*2+i]
                                  - 2774 * previous_fs[3*nbodies*2*2+i]
                                  + 2616 * previous_fs[2*nbodies*2*2+i]
                                  - 1274 * previous_fs[nbodies*2*2+i]
                                  + 251 * previous_fs[i]);
    }

    delete[] flux;
}
void Nbody::adams_moulton(void (Nbody::*f)(float *, float *), float * qq, float * previous_fs){
    // See http://www.math.iit.edu/~fass/478578_Chapter_2.pdf
    // predictor
    float * q_temp = new float[nbodies*2*2];
    adams_bashforth(f, q_temp, previous_fs);

    // corrector
    f(q_temp, q_temp);
    for (int i = 0; i < nbodies*2*2; i++) {
        previous_fs[i] = q_temp[i];
    }

    for (int i = 0; i < nbodies*2*2; i++) {
        qq[i] = q[i] + dt / 720.0 * (251 * previous_fs[i]
                                  + 646 * previous_fs[4*nbodies*2*2+i]
                                  - 264 * previous_fs[3*nbodies*2*2+i]
                                  + 106 * previous_fs[2*nbodies*2*2+i]
                                  - 19 * previous_fs[1*nbodies*2*2+i]);
    }

    delete[] q_temp;
}
float Nbody::potential_energy(){
    float Uu = 0;
    float * r = new float[nbodies*2];
    for (int n = 0; n < nbodies; n++){
        for (int i = 0; i < 2; i++) {
            float dist = 0;
            for (int m = 0; m < nbodies; m++) {
                r[m*2+i] = q[(n*2+1)*2+i] - q[(m*2+1)*2+i];
                dist += r[m*2+i] * r[m*2+i];
            }
            dist = sqrt(dist);
            if (abs(dist) < 1.0e-10) dist = 1.0e12;
            for (int m = 0; m < nbodies; m++) {
                Uu += G * masses[n] * masses[m] / dist;
            }
        }
    }
    delete[] r;
    return Uu;
}

float Nbody::kinetic_energy(){
    float p_squared = 0;
    for (int n = 0; n < nbodies; n++){
        for (int i = 0; i < 2; i++) {
            p_squared += 0.5 * masses[n] * q[(n*2+1)*2+i] * q[(n*2+1)*2+i];
        }
    }
    return p_squared;
}
float Nbody::moment_of_inertia(){
    float inertia = 0;
    for (int n = 0; n < nbodies; n++) {
        for (int i = 0; i < 2; i++) {
            inertia += masses[n] * q[(n*2+1)*2+i] * q[(n*2+1)*2+i];
        }
    }
    return inertia;
}
void Nbody::write_to_file(){
        self.output_file["q"][self.time,:,:,:] = self.q
        self.output_file["U"][self.time] = self.U
        self.output_file["T"][self.time] = self.T
        self.output_file["I"][self.time] = self.I
}
void Nbody::contain_particles(){
    // Stop particles moving outside of box
    for (int n = 0; n < nbodies; n++) {
        // x-direction
        if (q[(n*2+1)*2] < -0.5*box_size[0]){
            q[(n*2+1)*2] = -box_size[0] - q[(n*2+1)*2];
            q[(n*2)*2] *= -1.;
        } else if (q[(n*2+1)*2] > 0.5*box_size[0]){
            q[(n*2+1)*2] = box_size[0] - q[(n*2+1)*2];
            q[(n*2)*2] *= -1.;
        }
        // y-direction
        if (q[(n*2+1)*2+1] < -0.5*box_size[1]){
            q[(n*2+1)*2+1] = -box_size[1] - q[(n*2+1)*2+1];
            q[(n*2)*2+1] *= -1.;
        } else if (q[(n*2+1)*2+1] > 0.5*box_size[1]) {
            q[(n*2+1)*2+1] = box_size[1] - q[(n*2+1)*2+1];
            q[(n*2)*2+1] *= -1.;
        }
    }

}
void Nbody::evolve(){
        // run simulation

        float * previous_fs = new float[5*nbodies*2*2];

        for (int i = 0; i < nt; i++){
            t += 1;

            if (integrator == adams_moulton || integrator == adams_bashforth){
                if (t > 5){
                    q = integrator(self.step, previous_fs);
                } else {
                    previous_fs[self.time-1,:,:,:] = step(self.q);
                    q = dormand_prince(self.step);
                }
            } else {
                q = integrator(self.step);
            }

            contain_particles();

            U = potential_energy();
            T = kinetic_energy();
            I = moment_of_inertia();

            write_to_file();
        }

        delete[] previous_fs;
}
void Nbody::print_state(){
        // print state to screen
        printf("time = %d,  U = %f,  T = %f,  I = %f", t, U, T, I);
}

int main(){
    int nbodies = 100;
    int nt = 1000;
    float dt = 0.0001;
    char output_file[] = "nbody_data.h5";
    string integrator_function = "rk3";
    Nbody nbody = Nbody(nbodies, nt, dt, output_file, integrator_function);
    nbody.evolve();
}

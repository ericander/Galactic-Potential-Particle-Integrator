#ifndef __PARTICLE_HPP_INCLUDED__
#define __PARTICLE_HPP_INCLUDED__

#include "Galaxy.hpp"
#include "Satellite.hpp"
#include <vector>
class Particle {
    /* Class for particls that that can be integrated in galactic 
     * potentials. The class uses a fourth order Runge-Kutta to 
     * integrate particle positions in the potentials.
     *
     * Instance Variables:
     *      npar
     *          Number of particles
     *      t
     *          Current time for every particle
     *      x,y,z
     *          Vector with positional coordinate
     *      vx,vy,vz
     *          Vector with time-derivative of position
     *
     *****************************************************************/
    // Last modified: Eric Andersson, 11 Jan 2018.
    
    public:
        double npar, t;
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> z;
        std::vector<double> vx;
        std::vector<double> vy;
        std::vector<double> vz;

        Particle(std::vector<double>, std::vector<double>, 
                 std::vector<double>, std::vector<double>,
                 std::vector<double>, std::vector<double>);
        void update(Galaxy, Satellite, double&, bool);
        void rk4(Galaxy, Satellite, double [3], double [3],
                int, double);
        void rkck(Galaxy, Satellite, double [3], double [3],
                int, double&, double&, double);
};
#endif

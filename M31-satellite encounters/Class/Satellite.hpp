#ifndef __SATELLITE_HPP_INCLUDED__
#define __SATELLITE_HPP_INCLUDED__

#include "Galaxy.hpp"
class Satellite {
    /* Class for the satellite galaxy.
     *
     * Class Variables:
     *      M_s
     *          Mass of the satellite.
     *      r_s
     *          Characteristic radius of the satellite.
     *      x, y, z
     *          Position of the satellite galaxy.
     *      vx, vy, vz
     *          Velocity of the satellite galaxy.
     *      t
     *          Current time.
     *****************************************************************/
    // Last modified: Eric Andersson, 11 Jan 2018.
    
    public:
        double M_s, r_s, x, y, z, vx, vy, vz, t;

    Satellite(double, double, double, double, double,
              double, double, double);
    void rk4(Galaxy, Satellite, double [3], double [3], double);
    void rkck(Galaxy, Satellite, double [3], double [3], double&,
                double&);
    void update(Galaxy, Satellite, double&, double&, bool);
    void acceleration(Satellite, double, double, double, double *);
    double potential(Satellite, double, double, double);
};

#endif

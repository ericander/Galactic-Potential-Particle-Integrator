#include "Galaxy.hpp"
#include <cmath>
#include <iostream>
Galaxy::Galaxy ( double init_M_b, double init_r_b,
                      double init_M_d, double init_A, double init_B,
                      double init_r_c, double init_V_h, 
                      double init_vr_bulge, double init_vl_disc,
                      double init_vt_disc) {
    /* Constructor for Galaxy class. Initiates the major galaxy as a 
     * Galaxy object.
     *
     *      Arguments:
     *          init_M_b
     *              Characteristic bulge mass.
     *          init_r_b
     *              Characteristic bulge radius.
     *          init_M_d
     *              Characteristic disc mass.
     *          init_A
     *              Characteristic disc radius.
     *          init_B
     *              Characteristic disc thickness.
     *          init_r_c
     *              Characteristic halo radius.
     *          init_V_h
     *              Rotational velocity in halo at r_c.
     *          init_vr_disc
     *              Visible radius of bulge.
     *          init_vl_disc
     *              Visible length of disc.
     *          init_vt_disc
     *              Visible thickness of disc.
     *
     *****************************************************************/
    // Last modified: Eric Andersson, 11 Jan 2018.

    M_b = init_M_b;
    r_b = init_r_b;
    M_d = init_M_d;
    A = init_A;
    B = init_B;
    r_c = init_r_c;
    V_h = init_V_h;
    vr_bulge = init_vr_bulge;
    vl_disc = init_vl_disc;
    vt_disc = init_vt_disc;
}

void Galaxy::acceleration (Galaxy g,
        double x, double y, double z, double *a ) {
    /* Function that computes the acceleration of a particle in the 
     * galactic potential.
     *
     *      Arguments:
     *          g
     *              Major galaxy object.
     *          x,y,z
     *              Coordinates of particle.
     *          *a
     *              Pointer to first element in a vector [3] to which 
     *              the acceleration is set as (ax, ay, az).
     *
     *****************************************************************/
    // Last modified: Eric Andersson, 11 Jan 2018.

    // Allocate memory for acceleration components
    double a_bulge [3];
    double a_disc [3];
    double a_halo [3];

    // Bulge acceleration.
    g.bulge_acceleration(g, x, y, z, a_bulge); 
    
    // Disc acceleration.
    g.disc_acceleration(g, x, y, z, a_disc);

    // Halo acceleration.
    g.halo_acceleration(g, x, y, z, a_halo);

    // Sum up and assign acceleration.
    *a = a_bulge[0] + a_disc[0] + a_halo[0];
    *(a+1) = a_bulge[1] + a_disc[1] + a_halo[1];
    *(a+2) = a_bulge[2] + a_disc[2] + a_halo[2];
}

void Galaxy::bulge_acceleration(Galaxy g, double x, double y, 
        double z, double *a) {
    /* Function for calculating the bulge component of the
     * acceleration. The bulge potential is approximated with a
     * Hernquist potential given by
     *       
     *   \Phi_{\rm b}(r) = - \frac{GM_b}{r_b + r}.
     *
     * See Hernquist (1990) for details.
     *
     *      Arguments:
     *          g
     *              Major galaxy object.
     *          x,y,z
     *              Coordinates of particle.
     *          *a
     *              Pointer to first element in a vector [3] to which 
     *              the acceleration is set as (ax, ay, az).
     *
     *****************************************************************/
    // Last modified: Eric Andersson, 11 Jan 2018.

    // Constants
    const double PI = 3.14159265359;             // Value of pi
    const double G = 4.558e-13 * pow(PI, 2);     // kpc³ M_sun⁻¹ Myr⁻²

    // Common denominator
    double denominator = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) * 
       pow(g.r_b + sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)), 2);

    // x-direction
    *a = - G * g.M_b * x / denominator;

    // y-direction
    *(a+1) = - G * g.M_b * y / denominator;

    // z-direction
    *(a+2) = - G * g.M_b * z / denominator;

}

void Galaxy::disc_acceleration(Galaxy g, double x, double y, double z,
        double *a) {
    /* Function for calculating the disc component of the
     * acceleration. The potential is approximated by a Miyamoto-Nagai
     * potential given by 
     *
     *   \Phi_{\rm disc}(R,z) = 
     *   - \frac{GM_d}{\sqrt{R^2 + (a + \sqrt{z^2 + b^2})^2}}.
     *
     * See Miamoyo & Nagai (1975) for details.
     *
     *      Arguments:
     *          g
     *              Major galaxy object.
     *          x,y,z
     *              Coordinates of particle.
     *          *a
     *              Pointer to first element in a vector [3] to which 
     *              the acceleration is set as (ax, ay, az).
     *
     *****************************************************************/
    // Last modified: Eric Andersson, 11 Jan 2018.
    
    // Constants
    const double PI = 3.14159265359;             // Value of pi
    const double G = 4.558e-13 * pow(PI, 2);     // kpc³ M_sun⁻¹ Myr⁻²
    
    // Denominator x and y
    double xydenominator = pow(
            ( pow(g.A + sqrt( pow(g.B, 2) + pow(z, 2) ), 2) 
            + pow(x, 2) + pow(y, 2) ), 3./2.);

    // x-direction
    *a = - G * g.M_d * x / xydenominator;

    // y-direction
    *(a+1) = - G * g.M_d * y / xydenominator;

    // Denominator for z
    double zdenominator = xydenominator * sqrt(pow(g.B, 2) + pow(z, 2));

    // z-direction
    *(a+2) = - G * g.M_d * z * (g.A + sqrt( pow(g.B, 2) + pow(z, 2) ) )
        / zdenominator;

}

void Galaxy::halo_acceleration(Galaxy g, double x, double y, double z,
        double *a) {
    /* Function for calculating the halo component of the
     * acceleration. The halo potential is approximated by a 
     * logarithmic potential, given by
     *
     *   \Phi_{\rm halo}(r) = \frac{1}{2}V_h^2\ln(r_h^2 + r^2).
     *
     * See, e.g., Font et al. (2006) for details. Outside the halo the
     * potential is given by the potential of a point mass.
     *
     *      Arguments:
     *          g
     *              Major galaxy object.
     *          x,y,z
     *              Coordinates of particle.
     *          *a
     *              Pointer to first element in a vector [3] to which 
     *              the acceleration is set as (ax, ay, az).
     *
     *****************************************************************/
    // Last modified: Eric Andersson, 11 Jan 2018.
    
    // Constants
    const double PI = 3.14159265359;             // Value of pi
    const double G = 4.558e-13 * pow(PI, 2);     // kpc³ M_sun⁻¹ Myr⁻²
    const double M_M31 = 1.24e12;                // M_sun
    const double K = 0.001022;                   // Conversion factor
    double r_cut = 155.08;                       // kpc
    
    if ( sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) < r_cut) {
        // Acceleration if particle lie within the halo.

        // Common denominator.
        double denominator = pow(g.r_c, 2) + pow(x, 2) + pow(y, 2) 
            + pow(z, 2);
        
        // x-direction
        *a = - pow(g.V_h*K,2) * x / denominator;

        // y-direction
        *(a+1) = - pow(g.V_h*K,2) * y / denominator;

        // z-direction
        *(a+2) = - pow(g.V_h*K,2) * z / denominator;
    
    }else {
        // Acceleration if the particle lie outside the halo.
        
        // Common denominator
        double denominator = pow( 
                pow(x, 2) + pow(y, 2) + pow(z, 2), 3./2.);

        // Halo mass
        double M_halo = M_M31 - (g.M_b + g.M_d);

        // x-direction
        *a = - G * M_halo * x / denominator;

        // y-direction
        *(a+1) = - G * M_halo * y / denominator;

        // z-direction
        *(a+2) = - G * M_halo * z / denominator;

    }

}

double Galaxy::collision(Galaxy g, double x, double y, double z) {
    /* Function that determines whether the satellite has is passing 
     * through visible parts of the major galaxy or not.
     *
     *      Arguments:
     *          g
     *              Galaxy object.
     *          x, y, z
     *              Position of object that can collide with Galaxy.
     *
     ******************************************************************/
    // Last modified: Eric Andersson, 11 Jan 2018.
    
    double collision = 0;
    
    // Check if object passes trough bulge.
    double r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    if(r < g.vr_bulge) {
        collision = 1;
    }

    // Check if object passes trough disc.
    double R = sqrt(pow(x, 2) + pow(y, 2));
    if(R < g.vl_disc && std::abs(z) < g.vt_disc) {
        collision = 1; 
    }
    return collision;
}

#include "Satellite.hpp"
#include "Galaxy.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

Satellite::Satellite(double M_s_init, double r_s_init, 
                    double x0, double y0, double z0, 
                    double vx0, double vy0, double vz0) {
    /* Constructor for Satellite class. Initiates the galaxy as a 
     * Satellite object.
     *
     *      Arguments:
     *          M_s_init
     *              Characteristic dwarf galaxy mass.
     *          r_s_init
     *              Characteristic dwarf galaxy radius.
     *          x0
     *              Initial x position of dwarf galaxy.
     *          y0
     *              Initial y position of dwarf galaxy.
     *          z0
     *              Initial z position of dwarf galaxy.
     *          vx0
     *              Initial x velocity of dwarf galaxy.
     *          vy0
     *              Initial y velocity of dwarf galaxy.
     *          vz0
     *              Initial z velocity of dwarf galaxy.
     *
     *****************************************************************/
    // Last modified: Eric Andersson, 11 Jan 2018.
    
    // Initialize class variables.
    M_s = M_s_init;
    r_s = r_s_init;
    x = x0;
    y = y0;
    z = z0;
    vx = vx0;
    vy = vy0;
    vz = vz0;
    t = 0;
}

void Satellite::rk4(Galaxy g, Satellite s, double r [3], 
        double v [3], double dt) {
    /* Function that integrates the dwarf galaxy trajectory by taking 
     * a fourth-order Runge-Kutta step with stepsize dt.
     *
     *      Arguments:
     *          g
     *              Galaxy object for calling acceleration function.
     *          r
     *              Vector with the positional coordinates (x, y, z)
     *          v
     *              Vector with the velocity (vx, vy, vz)
     *          dt
     *              Step size of the integration in time.
     *
     *****************************************************************/
    // Last modified: Eric Andersson, 11 Jan 2018.
    
    // Initialize increments.
    double dr1 [3];
    double dr2 [3];
    double dr3 [3];
    double dr4 [3];
    double dv1 [3];
    double dv2 [3];
    double dv3 [3];
    double dv4 [3];
    
    // Initialize variables.
    double a [3];       // Total acceleration
    double refSys [3];  // Acceleration of refrence system

    g.acceleration(g, r[0], r[1], r[2], a);
      
    // Account for acceleration of the reference system.
    s.acceleration(s, 0, 0, 0, refSys);    
    for(int i=0; i<3; i++) {
        a[i] -= refSys[i];
    }

    // Compute the first increments.
    for(int i=0; i<3; i++){
        dr1[i] = dt * v[i];
        dv1[i] = dt * a[i];
    }
    
    // Compute the second increments
    g.acceleration(g,
        r[0]+dr1[0]/2.,r[1]+dr1[1]/2.,r[2]+dr1[2]/2.,a);
    
    // Account for acceleration of referece fram.
    s.acceleration(s, -(dr1[0]/2.), -(dr1[1]/2.),
        -(dr1[2]/2.), refSys);    

    for(int i=0; i<3; i++) {
        a[i] -= refSys[i];
    }

    for(int i=0; i<3; i++){
        dr2[i] = dt * (v[i] + dv1[i]/2.);
        dv2[i] = dt * a[i];
    }
    
    // Compute the third increments
    g.acceleration(g,
        r[0]+dr2[0]/2.,r[1]+dr2[1]/2.,r[2]+dr2[2]/2., a);

    // Account for acceleration of referece fram.
    s.acceleration(s, -(dr2[0]/2.), -(dr2[1]/2.), -(dr2[2]/2.),
        refSys);    
    
    for(int i=0; i<3; i++) {
        a[i] -= refSys[i];
    }
        
    for(int i=0; i<3; i++){
        dr3[i] = dt * (v[i] + dv2[i]/2.);
        dv3[i] = dt * a[i];
    }
    
    // Compute the fourth increments
    g.acceleration(g, r[0]+dr3[0], r[1]+dr3[1], r[2]+dr3[2], a);
        
    // Account for acceleration of referece fram.
    s.acceleration(s, -dr3[0], -dr3[1], -dr3[2], refSys);    
        
    for(int i=0; i<3; i++) {
        a[i] -= refSys[i];
    }
    
    for(int i=0; i<3; i++){
        dr4[i] = dt * (v[i] + dv3[i]);
        dv4[i] = dt * a[i];
    }
        
    // Sum up all increments. 
    x += dr1[0]/6. + dr2[0]/3. + dr3[0]/3. + dr4[0]/6.; 
    y += dr1[1]/6. + dr2[1]/3. + dr3[1]/3. + dr4[1]/6.;
    z += dr1[2]/6. + dr2[2]/3. + dr3[2]/3. + dr4[2]/6.;
    vx += dv1[0]/6. + dv2[0]/3. + dv3[0]/3. + dv4[0]/6.; 
    vy += dv1[1]/6. + dv2[1]/3. + dv3[1]/3. + dv4[1]/6.; 
    vz += dv1[2]/6. + dv2[2]/3. + dv3[2]/3. + dv4[2]/6.; 
}

void Satellite::rkck(Galaxy g, Satellite s, double r [3], 
        double v [3], double &dt, double &time) {
    /* Function that integrates the dwarf galaxy trajectory by taking 
     * a Cash-Karp Runge-Kutta step which adapts the time step size 
     * given a maximum allowed error in the integration.
     *
     *      Arguments:
     *          g
     *              Galaxy object for calling acceleration function.
     *          r
     *              Vector with the positional coordinates (x, y, z).
     *          v
     *              Vector with the velocity (vx, vy, vz).
     *          dt
     *              Pointer to timestep value location.
     *          time
     *              Pointer to time value location.
     *
     *****************************************************************/
    // Last modified: Eric Andersson, 11 Jan 2018.
    
    // Maximum error that we accept.
    static const double eps = 1e-6;
    const double SAFETY = 0.9;
    const double ERRCON = 1.89e-4;

    // Cash-Karp Parameters for embedded Runge-Kutta method.
    static const double a2=0.2, a3=0.3, a4=0.6, a5=1.0, a6=0.875,
        b21=0.2, b31=3.0/40.0, b32=9.0/40.0, b41=0.3, b42=-0.9,
        b43=1.2, b51=-11.0/54.0, b52=2.5, b53=-70.0/27.0,
        b54=35.0/27.0, b61=1631.0/55296.0, b62=175.0/512.0,
        b63=575.0/13824.0, b64=44275.0/110592.0, b65=253.0/4096.0,
        c1=37.0/378.0, c2=0, c3=250.0/621.0, c4=125.0/594.0, c5=0,
        c6=512.0/1771.0, dc1=c1-2825.0/27648.0, dc3=c3-18575.0/48384.0,
        dc4=c4-13525.0/55296.0, dc5=-277.0/14336.0, dc6=c6-0.25;
        
    // Initialzing increments.
    double dr1 [3], dr2 [3], dr3 [3], dr4 [3], dr5 [3], dr6 [3];
    double dv1 [3], dv2 [3], dv3 [3], dv4 [3], dv5 [3], dv6 [3];

    // Initialize variables.
    double xtemp, ytemp, ztemp, vxtemp, vytemp, vztemp;
    double a [3];
    double a_s [3];
    double dttemp;
    double errmax;
    double rerr [3], verr [3];
    double refSys [3];

    // Scale the error to each value.
    double rscal [3], vscal [3];
    for (int i=0; i<3; i++) {
        // Avoid infinitly small errors.
        rscal[i] = (r[i] == 0 ? eps : eps * r[i]);
        vscal[i] = (v[i] == 0 ? eps : eps * v[i]);
    }
    
    // Attempt integration.
    dttemp = dt;
    for (;;) {
        // Compute first increment.
        g.acceleration(g, r[0], r[1], r[2], a);
        // Account for acceleration of the reference system.
        s.acceleration(s, 0, 0, 0, refSys);    
        for(int i=0; i<3; i++) {
            a[i] -= refSys[i];
        }

        for(int i=0; i<3; i++) {
            dr1[i] = dttemp * v[i];
            dv1[i] = dttemp * a[i];
        }

        // Compute seccond increment.
        g.acceleration(g, r[0] + b21 * dr1[0],
                          r[1] + b21 * dr1[1], 
                          r[2] + b21 * dr1[2], a);        
 
        // Account for acceleration of referece fram.
        s.acceleration(s, -(b21 * dr1[0]), -(b21 * dr1[1]),
            -(b21 * dr1[2]), refSys);    

        for(int i=0; i<3; i++) {
            a[i] -= refSys[i];
        }

        for(int i=0; i<3; i++) {
            dr2[i] = dttemp * (v[i] + b21 * dv1[i]);
            dv2[i] = dttemp * a[i];
        }
    
        // Compute third increment.
        g.acceleration(g, r[0] + b31 * dr1[0] + b32 * dr2[0],
                          r[1] + b31 * dr1[1] + b32 * dr2[1], 
                          r[2] + b31 * dr1[2] + b32 * dr2[2], a);

        // Account for acceleration of referece fram.
        s.acceleration(s, -(b31 * dr1[0] + b32 * dr2[0]),
            -(b31 * dr1[1] + b32 * dr2[1]),
            -(b31 * dr1[2] + b32 * dr2[2]), refSys);    
        for(int i=0; i<3; i++) {
            a[i] -= refSys[i];
        }

        for(int i=0; i<3; i++) {
            dr3[i] = dttemp * (v[i] + b31 * dv1[i] + b32 * dv2[i]);
            dv3[i] = dttemp * a[i];
        }

        // Compute fourth increment.
        g.acceleration(g, 
                r[0] + b41 * dr1[0] + b42 * dr2[0] + b43 * dr3[0],
                r[1] + b41 * dr1[1] + b42 * dr2[1] + b43 * dr3[1],
                r[2] + b41 * dr1[2] + b42 * dr2[2] + b43 * dr3[2],
                a);
        
        // Account for acceleration of referece fram.
        s.acceleration(s, 
            -(b41 * dr1[0] + b42 * dr2[0] + b43 * dr3[0]),
            -(b41 * dr1[1] + b42 * dr2[1] + b43 * dr3[1]),
            -(b41 * dr1[2] + b42 * dr2[2] + b43 * dr3[2]),
            refSys);    
        
        for(int i=0; i<3; i++) {
            a[i] -= refSys[i];
        }
        
        for(int i=0; i<3; i++) {
            dr4[i] = dttemp * 
                (v[i] + b41*dv1[i] + b42*dv2[i] + b43*dv3[i]);
            dv4[i] = dttemp * a[i];
        }

        // Compute fifth increment.
        g.acceleration(g, 
                r[0] + b51 * dr1[0] + b52 * dr2[0] + b53 * dr3[0]
                + b54 * dr4[0],
                r[1] + b51 * dr1[1] + b52 * dr2[1] + b53 * dr3[1]
                + b54 * dr4[1],
                r[2] + b51 * dr1[2] + b52 * dr2[2] + b53 * dr3[2]
                + b54 * dr4[2],
                a);
 
        // Account for acceleration of referece fram.
        s.acceleration(s, 
            -(b51*dr1[0] + b52*dr2[0] + b53*dr3[0] + b54*dr4[0]),
            -(b51*dr1[1] + b52*dr2[1] + b53*dr3[1] + b54*dr4[1]),
            -(b51*dr1[2] + b52*dr2[2] + b53*dr3[2] + b54*dr4[2]),
            refSys);    

        for(int i=0; i<3; i++) {
            a[i] -= refSys[i];
        }

        for(int i=0; i<3; i++) {
            dr5[i] = dttemp * 
                (v[i] + b51*dv1[i] + b52*dv2[i] + b53*dv3[i] + 
                b54 * dv4[i]);
            dv5[i] = dttemp * a[i];
        }

        // Compute sixth increment.
        g.acceleration(g, 
                r[0] + b61 * dr1[0] + b62 * dr2[0] + b63 * dr3[0]
                + b64 * dr4[0] + b65 * dr5[0],
                r[1] + b61 * dr1[1] + b62 * dr2[1] + b63 * dr3[1]
                + b64 * dr4[1] + b65 * dr5[1],
                r[2] + b61 * dr1[2] + b62 * dr2[2] + b63 * dr3[2]
                + b64 * dr4[2] + b65 * dr5[2],
                a);
        
        // Account for acceleration of referece fram.
        s.acceleration(s, 
            -(b61 * dr1[0] + b62 * dr2[0] + b63 * dr3[0] 
                + b64 * dr4[0] + b65 * dr5[0]),
            -(b61 * dr1[1] + b62 * dr2[1] + b63 * dr3[1]
                + b64 * dr4[1] + b65 * dr5[1]),
            -(b61 * dr1[2] + b62 * dr2[2] + b63 * dr3[2]
                + b64 * dr4[2] + b65 * dr5[2]),
            refSys);    

        for(int i=0; i<3; i++) {
            a[i] -= refSys[i];
        }

        for(int i=0; i<3; i++) {
            dr6[i] = dttemp * 
                (v[i] + b61*dv1[i] + b62*dv2[i] + b63*dv3[i] + 
                b64 * dv4[i] + b65 * dv5[i]);
            dv6[i] = dttemp * a[i];
        }
        
        // Accumalate the results.
        xtemp = x + c1 * dr1[0] + c2 * dr2[0] + c3 * dr3[0] +
                 c4 * dr4[0] + c5 * dr5[0] + c6 * dr6[0];
        ytemp = y + c1 * dr1[1] + c2 * dr2[1] + c3 * dr3[1] +
                 c4 * dr4[1] + c5 * dr5[1] + c6 * dr6[1];
        ztemp = z + c1 * dr1[2] + c2 * dr2[2] + c3 * dr3[2] +
                 c4 * dr4[2] + c5 * dr5[2] + c6 * dr6[2];
        vxtemp = vx + c1 * dv1[0] + c2 * dv2[0] + c3 * dv3[0] +
                  c4 * dv4[0] + c5 * dv5[0] + c6 * dv6[0];
        vytemp = vy + c1 * dv1[1] + c2 * dv2[1] + c3 * dv3[1] +
                  c4 * dv4[1] + c5 * dv5[1] + c6 * dv6[1];
        vztemp = vz + c1 * dv1[2] + c2 * dv2[2] + c3 * dv3[2] +
                  c4 * dv4[2] + c5 * dv5[2] + c6 * dv6[2];

        // Estimate the error of integration attempt.
        for (int i=0; i<3; i++) {
            rerr[i] = dc1*dr1[i] + dc3*dr3[i] + dc4*dr4[i] + dc5*dr5[i]
                + dc6*dr6[i];
            verr[i] = dc1*dv1[i] + dc3*dv3[i] + dc4*dv4[i] + dc5*dv5[i]
                + dc6*dv6[i];
        }
       
        // Find maximum error.
        errmax = 0;
        for (int i=0; i<3; i++) {
            errmax = std::max(errmax, std::fabs(rerr[i]/rscal[i]));
            errmax = std::max(errmax, std::fabs(verr[i]/vscal[i]));
        }

        errmax /= eps;  // Scale relative to required tolerance.
        if (errmax <= 1.0) {
            // Accurate integration successfull.
            x = xtemp;
            y = ytemp;
            z = ztemp;
            vx = vxtemp;
            vy = vytemp;
            vz = vztemp;
            t += dttemp;
            // Break if satellite has caught up to particles.
            if (t >= time) {
                if (t <= 0.999 * time){
                    // Keep track of time divergnce.
                    std::cout << "WARNING! Time-difference is " << 
                                t - time << std::endl;
                }
                break;
            }
        }else {
            // Decrease timestep by some integer value with is no
            // smaller than a tenth of original timestep.
            dttemp = std::max(dt/10., dttemp / 2.0); 
        }
        
        // Keep track of infinitly small time-step exception.
        if (dttemp == 0) {
        std::cout << "Stepsize underflow in Runge-Kutta." << std::endl;
            exit (EXIT_FAILURE);
        } 
    }
}

void Satellite::update(Galaxy g, Satellite s, double &dt,
                        double &time, bool adapt_timestep) {
    /* Function that updates the satellite galaxy position.
     *
     *      Arguments:
     *          g
     *              Galaxy object for calling acceleration function.
     *          s
     *              Satellite object.
     *          dt 
     *              Step size in time.
     *          time
     *              Current time in simulation.
     *          adapt_timestep
     *              True/False for adapive time step integrator.
     *
     *****************************************************************/
    // Last modified: Eric Andersson, 11 Jan 2018.
    
    // Initiate vectors holding phase space coordinates.
    double r [3] = {x, y, z};
    double v [3] = {vx, vy, vz};

    // Check whether to use adaptive timestep or not and integrate
    // accordingly.
    if (adapt_timestep) {
        this->rkck(g, s, r, v, dt, time);
    }else {
        this->rk4(g, s, r, v, dt);
        t += dt;
    }
}

void Satellite::acceleration(Satellite s, 
    double xp, double yp, double zp, double *a){
    /* Function that calculates the acceleration of a particle in the
     * satellite galaxy potential. The satellite potential is given by
     * a Plummer model, given by 
     *
     *   \Phi_p = \frac{GM_{\rm p}}{\sqrt{r_{\rm p}^2 + r^2}}.
     *
     * See Plummer (1911) for details.
     *
     *      Arguments:
     *          s
     *              Satellite object.
     *          xp 
     *              Position of particle in x.
     *          yp 
     *              Position of particle in y.
     *          zp 
     *              Position of particle in z.
     *          a
     *              Pointer to location where acceleration value will
     *              be updated.
     *
     *****************************************************************/
    // Last modified: Eric Andersson, 11 Jan 2018.
    
    // Constants
    const double PI = 3.14159265359;        // the value of pi
    const double G = 4.558e-13 * pow(PI, 2);// kpc³ M_sun⁻¹ Myr⁻²
    
    // Calculate denominator.
    double denominator = pow(
        pow(xp - s.x, 2) + pow(yp - s.y, 2) + 
        pow(zp - s.z, 2) + pow(s.r_s, 2), 3./2.);

    // Calculate the acceleration in each direction.
    *a = - G * s.M_s * (xp - s.x) / denominator;
    *(a+1) = - G * s.M_s * (yp - s.y) / denominator;
    *(a+2) = - G * s.M_s * (zp - s.z) / denominator;
}

double Satellite::potential(Satellite s, double x, double y, double z){
    /* Function that computes the strenth of the potential field at a 
     * given position. The satellite potential is given by
     * a Plummer model, given by 
     *
     *   \Phi_p = -\frac{GM_{\rm p}}{\sqrt{r_{\rm p}^2 + r^2}}.
     *
     * See Plummer (1911) for details.
     *
     *      Arguments:
     *          s
     *              Satellite object.
     *          x 
     *              Position of particle in x.
     *          y 
     *              Position of particle in y.
     *          z 
     *              Position of particle in z.
     *
     *****************************************************************/
    // Last modified: Eric Andersson, 20 Feb 2018.
    
    // Constants
    const double PI = 3.14159265359;        // the value of pi
    const double G = 4.558e-13 * pow(PI, 2);// kpc³ M_sun⁻¹ Myr⁻²
    
    double Phi;

    double r = sqrt(pow(x,2) + pow(y,2) + pow(z,2));
    double rs = sqrt(pow(s.x,2) + pow(s.y,2) + pow(s.z,2));
    
    // Position relative the satellite.
    r -= rs;

    Phi = - s.M_s * G / (sqrt(pow(s.r_s,2) + pow(r,2)));

    return Phi;
}

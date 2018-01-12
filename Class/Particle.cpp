#include "Particle.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <stdlib.h>
Particle::Particle(std::vector<double> x0, std::vector<double> y0,
        std::vector<double> z0, std::vector<double> vx0, 
        std::vector<double> vy0, std::vector<double> vz0) {
    /* Constructor for the Particle object. Initiates a set of 
     * particles. The parameter param for particle n is given by
     * particle.param[n] where particle is the Particle object.
     *
     *      Arguments:
     *          x0,y0,z0
     *              Vector with initial Cartiesian coordinates (x,y,z)
     *              for each particle.
     *          vx0,vy0,vz0
     *              Vector with initial velocities (vx, vy, vz).
     *
     *****************************************************************/
    // Last modified: Eric Andersson, 11 Jan 2018.
    
    // Constants
    double K = 0.001022;            // Conversion factor

    // Initialize all class variables.
    npar = x0.size();
    t = 0;
    x = x0;
    y = y0;
    z = z0;
    vx = vx0;
    vy = vy0;
    vz = vz0;

    // Transform to correct unit
    for (int i = 0; i<npar; i++) {
        vx[i] *= K;
        vy[i] *= K;
        vz[i] *= K;
    }
}

void Particle::rk4(Galaxy g, Satellite s, double r [3], 
        double v [3], int index, double dt) {
    /* Function that integrates the dwarf galaxy trajectory by taking 
     * a fourth-order Runge-Kutta step with stepsize dt.
     *
     *      Arguments:
     *          g
     *              Galaxy object for calling acceleration function.
     *          s
     *              Satellite object for calling acceleration function.
     *          r
     *              Vector with the positional coordinates (x, y, z).
     *          v
     *              Vector with the velocity (vx, vy, vz).
     *          index
     *              Index of integrated particle.
     *          dt
     *              Step size of the integration in time.
     *
     *****************************************************************/
    // Last modified: Eric Andersson, 11 Jan 2018.
    
    // Initiate increments.
    double dr1 [3];
    double dr2 [3];
    double dr3 [3];
    double dr4 [3];
    double dv1 [3];
    double dv2 [3];
    double dv3 [3];
    double dv4 [3];
    
    // Initiate variables.
    double a [3];       // Total acceleration
    double a_s [3];     // Acceleration from satellite
    double refSys [3];  // Acceleration of reference system

    // Acceleration of the reference system.
    s.acceleration(s, 0, 0, 0, refSys);
    
    // Compute the first increments
    g.acceleration(g, r[0], r[1], r[2], a);
    s.acceleration(s, r[0], r[1], r[2], a_s);
    
    for(int i=0; i<3; i++){
        a[i] += a_s[i] - refSys[i];
    }

    for(int i=0; i<3; i++){
        dr1[i] = dt * v[i];
        dv1[i] = dt * a[i];
    }
        
    // Compute the second increments
    g.acceleration(g, r[0]+dr1[0]/2.,r[1]+dr1[1]/2.,r[2]+dr1[2]/2.,a);
    s.acceleration(s, r[0]+dr1[0]/2.,r[1]+dr1[1]/2.,r[2]+dr1[2]/2.,a_s);
    
    for(int i=0; i<3; i++){
        a[i] += a_s[i] - refSys[i];
    }

    for(int i=0; i<3; i++){
        dr2[i] = dt * (v[i] + dv1[i]/2.);
        dv2[i] = dt * a[i];
    }
        
    // Compute the third increments
    g.acceleration(g, r[0]+dr2[0]/2.,r[1]+dr2[1]/2.,r[2]+dr2[2]/2., a);
    s.acceleration(s, r[0]+dr2[0]/2.,r[1]+dr2[1]/2.,r[2]+dr2[2]/2.,
            a_s);
        
    for(int i=0; i<3; i++){
        a[i] += a_s[i] - refSys[i];
    }

    for(int i=0; i<3; i++){
        dr3[i] = dt * (v[i] + dv2[i]/2.);
        dv3[i] = dt * a[i];
    }
    // Compute the fourth increments
    g.acceleration(g, r[0]+dr3[0], r[1]+dr3[1], r[2]+dr3[2], a);
    s.acceleration(s, r[0]+dr3[0], r[1]+dr3[1], r[2]+dr3[2], a_s);
         
    for(int i=0; i<3; i++){
        a[i] += a_s[i] - refSys[i];
    }

    for(int i=0; i<3; i++){
        dr4[i] = dt * (v[i] + dv3[i]);
        dv4[i] = dt * a[i];
    }
    
    // Sum up all increments. 
    x[index] += dr1[0]/6. + dr2[0]/3. + dr3[0]/3. + dr4[0]/6.; 
    y[index] += dr1[1]/6. + dr2[1]/3. + dr3[1]/3. + dr4[1]/6.;
    z[index] += dr1[2]/6. + dr2[2]/3. + dr3[2]/3. + dr4[2]/6.;
    vx[index] += dv1[0]/6. + dv2[0]/3. + dv3[0]/3. + dv4[0]/6.; 
    vy[index] += dv1[1]/6. + dv2[1]/3. + dv3[1]/3. + dv4[1]/6.; 
    vz[index] += dv1[2]/6. + dv2[2]/3. + dv3[2]/3. + dv4[2]/6.; 
}

void Particle::rkck(Galaxy g, Satellite s, double r [3], double v [3],
        int index, double &dt, double &err, double eps) {
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
     *          index
     *              Index of integrated particle.
     *          dt
     *              Adress to timestep value.
     *          err
     *              Adress to maximum error value.
     *          eps
     *              Maximim tolerance on error
     *****************************************************************/
    // Last modified: Eric Andersson, 11 Jan 2018.
    
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
    double a [3];
    double a_s [3];
    double rerr [3], verr [3]; 
    double refSys [3];

    // Scale the desired error to each value.
    double rscal [3], vscal [3];
    for (int i=0; i<3; i++) {
        // Avoid infinitly small errors.
        rscal[i] = (r[i] == 0 ? eps*0.00001 : eps * r[i]);
        vscal[i] = (v[i] == 0 ? eps*0.00001 : eps * v[i]);
    }

    // Acceleration of the reference system.
    s.acceleration(s, 0, 0, 0, refSys);

    // Perform integration.
    // Compute first increment.
    g.acceleration(g, r[0], r[1], r[2], a);
    s.acceleration(s, r[0], r[1], r[2], a_s);
     
    for(int i=0; i<3; i++) {
        a[i] += a_s[i] - refSys[i];
    }
    
    for(int i=0; i<3; i++) {
        dr1[i] = dt * v[i];
        dv1[i] = dt * a[i];
    }

    // Compute seccond increment.
    g.acceleration(g, r[0] + b21 * dr1[0],
                      r[1] + b21 * dr1[1], 
                      r[2] + b21 * dr1[2], a);
    s.acceleration(s, r[0] + b21 * dr1[0], 
                      r[1] + b21 * dr1[1], 
                      r[2] + b21 * dr1[2], a_s);
    
    for(int i=0; i<3; i++) {
        a[i] += a_s[i] - refSys[i];
    }

    for(int i=0; i<3; i++) {
        dr2[i] = dt * (v[i] + b21 * dv1[i]);
        dv2[i] = dt * a[i];
    }
    
    // Compute third increment.
    g.acceleration(g, r[0] + b31 * dr1[0] + b32 * dr2[0],
                      r[1] + b31 * dr1[1] + b32 * dr2[1], 
                      r[2] + b31 * dr1[2] + b32 * dr2[2], a);
    s.acceleration(s, r[0] + b31 * dr1[0] + b32 * dr2[0], 
                      r[1] + b31 * dr1[1] + b32 * dr2[1], 
                      r[2] + b31 * dr1[2] + b32 * dr2[2], a_s);

    for(int i=0; i<3; i++) {
        a[i] += a_s[i] - refSys[i];
    }

    for(int i=0; i<3; i++) {
        dr3[i] = dt * (v[i] + b31 * dv1[i] + b32 * dv2[i]);
        dv3[i] = dt * a[i];
    }

    // Compute fourth increment.
    g.acceleration(g, 
                r[0] + b41 * dr1[0] + b42 * dr2[0] + b43 * dr3[0],
                r[1] + b41 * dr1[1] + b42 * dr2[1] + b43 * dr3[1],
                r[2] + b41 * dr1[2] + b42 * dr2[2] + b43 * dr3[2],
                a);
    s.acceleration(s, 
                r[0] + b41 * dr1[0] + b42 * dr2[0] + b43 * dr3[0],
                r[1] + b41 * dr1[1] + b42 * dr2[1] + b43 * dr3[1],
                r[2] + b41 * dr1[2] + b42 * dr2[2] + b43 * dr3[2],
                a_s);

    for(int i=0; i<3; i++) {
        a[i] += a_s[i] - refSys[i];
    }

    for(int i=0; i<3; i++) {
        dr4[i] = dt * 
                (v[i] + b41*dv1[i] + b42*dv2[i] + b43*dv3[i]);
        dv4[i] = dt * a[i];
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
    s.acceleration(s, 
                r[0] + b51 * dr1[0] + b52 * dr2[0] + b53 * dr3[0]
                + b54 * dr4[0],
                r[1] + b51 * dr1[1] + b52 * dr2[1] + b53 * dr3[1]
                + b54 * dr4[1],
                r[2] + b51 * dr1[2] + b52 * dr2[2] + b53 * dr3[2]
                + b54 * dr4[2],
                      a_s);

    for(int i=0; i<3; i++) {
        a[i] += a_s[i] - refSys[i];
    }

    for(int i=0; i<3; i++) {
        dr5[i] = dt * 
                (v[i] + b51*dv1[i] + b52*dv2[i] + b53*dv3[i] + 
                b54 * dv4[i]);
        dv5[i] = dt * a[i];
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
    s.acceleration(s, 
                r[0] + b61 * dr1[0] + b62 * dr2[0] + b63 * dr3[0]
                + b64 * dr4[0] + b65 * dr5[0],
                r[1] + b61 * dr1[1] + b62 * dr2[1] + b63 * dr3[1]
                + b64 * dr4[1] + b65 * dr5[1],
                r[2] + b61 * dr1[2] + b62 * dr2[2] + b63 * dr3[2]
                + b64 * dr4[2] + b65 * dr5[2],
                a_s);

    for(int i=0; i<3; i++) {
        a[i] += a_s[i] - refSys[i];
    }

    for(int i=0; i<3; i++) {
        dr6[i] = dt * 
                (v[i] + b61*dv1[i] + b62*dv2[i] + b63*dv3[i] + 
                b64 * dv4[i] + b65 * dv5[i]);
        dv6[i] = dt * a[i];
    }
   
    // Accumalate the results.
    x[index] += c1 * dr1[0] + c2 * dr2[0] + c3 * dr3[0] +
                c4 * dr4[0] + c5 * dr5[0] + c6 * dr6[0];
    y[index] += c1 * dr1[1] + c2 * dr2[1] + c3 * dr3[1] +
                c4 * dr4[1] + c5 * dr5[1] + c6 * dr6[1];
    z[index] += c1 * dr1[2] + c2 * dr2[2] + c3 * dr3[2] +
                c4 * dr4[2] + c5 * dr5[2] + c6 * dr6[2];
    vx[index] += c1 * dv1[0] + c2 * dv2[0] + c3 * dv3[0] +
                 c4 * dv4[0] + c5 * dv5[0] + c6 * dv6[0];
    vy[index] += c1 * dv1[1] + c2 * dv2[1] + c3 * dv3[1] +
                 c4 * dv4[1] + c5 * dv5[1] + c6 * dv6[1];
    vz[index] += c1 * dv1[2] + c2 * dv2[2] + c3 * dv3[2] +
                 c4 * dv4[2] + c5 * dv5[2] + c6 * dv6[2];

    // Estimate the errors.
    for (int i=0; i<3; i++) {
        rerr[i] = dc1*dr1[i] + dc3*dr3[i] + dc4*dr4[i] + dc5*dr5[i]
                + dc6*dr6[i];
        verr[i] = dc1*dv1[i] + dc3*dv3[i] + dc4*dv4[i] + dc5*dv5[i]
                + dc6*dv6[i];
    }
       
    // Estimate largest error in integration.
    err = 0;
    for (int i=0; i<3; i++) {
        err = std::max(err, std::fabs(rerr[i]/rscal[i]));
        err = std::max(err, std::fabs(verr[i]/vscal[i]));
    }
}

void Particle::update(Galaxy g, Satellite s, double &dt, 
        bool adapt_timestep) {
    /* Function that updates particle positions.
     *
     *      Arguments:
     *          g
     *              Galaxy object for calling acceleration function.
     *          s
     *              Satellite object.
     *          dt 
     *              Step size in time.
     *          adapt_timestep
     *              True/False for adapive time step integrator.
     *
     *****************************************************************/
    // Last modified: Eric Andersson, 11 Jan 2018.
    
    double r [3], v [3];
    
    if (adapt_timestep){
        // Use adaptive timestep integrator.
        // Parameter for time-step adaptation.
        double err = 0;    // Maximum error for particlular particle.
        double errmax = 0; // Maximum error among all particles.
        static const double eps = 1e-6; // Relative error tolerance.
        const double SAFETY = 0.9;
        const double ERRCON = 1.89e-4;
        const double dtmax = 0.99;
    
        for(int i=0; i<npar; i++){
            r[0] = x[i];
            r[1] = y[i];
            r[2] = z[i];
            v[0] = vx[i];
            v[1] = vy[i];
            v[2] = vz[i];
            rkck(g, s, r, v, i, dt, err, eps);
            errmax = std::max(errmax, err);
        }
        t += dt;
        // Adapt the timestep.
        double dttemp;
        errmax /= eps;  // Scale relative to required tolerance.
        if (errmax <= 1.0) {
            // Increase time-step size if accurate enough.
            if (errmax > ERRCON) dt *= SAFETY*std::pow(errmax, -0.2);
            else std::max(dt *= 2.0, dtmax); // No more than a factor 2 
                                             // increase or dt > dtmax.
        }else {
            // Decrease the time-step to satisfy error requirement.
            dttemp = SAFETY*dt*std::pow(errmax, -0.2);
            
            // Decrease with max factor 10.
            dt = (dt >= 0.0 ? 
                std::max(dttemp, 0.1*dt) : std::min(dttemp, 0.1*dt));
        
            // Keep track of infinitly small time-step exception.
            if (dt == 0) {
                std::cout << "Stepsize underflow in Runge-Kutta." 
                    << std::endl;
                exit (EXIT_FAILURE);
            }
        }
    }else {
        // Use regular fourth-order Runge-Kutta.
        for(int i=0; i<npar; i++){
            r[0] = x[i];
            r[1] = y[i];
            r[2] = z[i];
            v[0] = vx[i];
            v[1] = vy[i];
            v[2] = vz[i];
            rk4(g, s, r, v, i, dt);
        }
        t += dt;
    }
}

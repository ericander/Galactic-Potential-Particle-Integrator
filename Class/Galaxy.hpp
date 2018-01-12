#ifndef __GALAXY_HPP_INCLUDED__
#define __GALAXY_HPP_INCLUDED__

class Galaxy {
    /* Class for the major galaxy. The potential of this galaxy consist
     * of a bulge, disc and halo potential which is fixed at origo of
     * the simulated space.
     *
     * Class variables:
     *      M_b
     *          Bulge mass [M_sun]
     *      r_b
     *          Bulge radius [kpc]
     *      M_d
     *          Disc mass [M_sun]
     *      A 
     *          Disc radial scale [kpc]
     *      B
     *          Disc thickness scale [kpc]
     *      r_c
     *          Characteristic radius of halo [kpc]
     *      V_h
     *          Rotational velocity at r_c [kpc Myr⁻¹]
     *      vr_disc
     *          Visible radius of bulge.
     *      vl_disc
     *          Visible length of disc.
     *      vt_disc
     *          Visible thickness of disc.
     *
     ******************************************************************/
    // Eric Andersson, 11 Jan 2018.
    public:
        double M_b, r_b, M_d, A, B, r_c, V_h, 
               vr_bulge, vl_disc, vt_disc; 
        Galaxy ( double, double, double, double, double, double,
                        double, double, double, double );
        void acceleration ( Galaxy, double, double, double, 
                double* );
        void bulge_acceleration ( Galaxy, double, double, double, 
                double* );
        void disc_acceleration ( Galaxy, double, double, double, 
                double* );
        void halo_acceleration ( Galaxy, double, double, double,
                double* );
        double collision(Galaxy, double, double, double);
};

#endif

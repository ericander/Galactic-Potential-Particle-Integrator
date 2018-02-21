/* ~/main.cpp
 *
 * Discription:
 * Program for simulating motions of globular clusters (GC) in galactic
 * potentials. The program can load multiple test-particles (GC's), 
 * however they are considered to be mass-less and will not interact
 * with each other.
 *
 * Author: Eric Andersson (eric)
 *********************************************************************/
// Eric Andersson, 30 Aug 2017
#include <iostream>
#include <cmath>
#include <vector>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
using namespace std;

#include "Class/Galaxy.hpp"
#include "Class/Particle.hpp"
#include "Class/Satellite.hpp"

int main (int argc, char* argv[]) {
    /* Main function of the program. This is the first function that
     * the program will call when executed. This function will set 
     * all parameters, initiate all variables and create all objects
     * that is used in the simulation. The function will then simulate 
     * the encounter by updating position and velocity the satellite 
     * galaxy and its particle population. This function will also print
     * usefull information during the run in the terminal and store 
     * data in given output files. 
     *
     *****************************************************************/
    // Last modified: Eric Andersson, 11 Jan 2017.

    // Declare and set system parameters.
    bool adaptive_timestep = false; // If true -> use adaptive timestep 
    double dt = 0.1;                // Time-step in Myr.
    double dt_rec = 5;              // Data recorded every dt_rec Myr.
    double time_tot = 14000;        // Time after encounter in Myr.
    double nruns = 1;               // Number of runs unless specified.
    int startRUN;
    int finalRUN;
    double K = 0.001022;            // Conversion factor for output.
    double r0 = 500;                // Start distance of satellite.

    // Set parameters for the potential fields.
    // Major galaxy
    double M_b = 3.3e10;       // (9.2e10)    Bulge mass [M_sun]. 
    double r_b = 0.61;         // (0.7)       Bulge radius [kpc].
    double M_d = 1.034365e11;  // (1.3e11)    Disc mass [M_sun].
    double A = 6.42833;        // (6.5)       Disc scale [kpc].
    double B = 0.26476;        // (0.26)      Disc scale [kpc].
    double r_h = 8.18;         // (12)        Halo radius [kpc].
    double V_h = 186;          //             Halo velocity [km/s].
    double delta = 27e4;       //             Halo denisty paramter.
    double rho = 136;          //             Critical density.
    double r_b_lum = 3;        // Luminous bulge radius [kpc].
    double h_d_lum = 1;        // Luminous disc thickness [kpc].
    double r_d_lum = 15;       // Luninous disc length [kpc].

    // Satellite galaxy
    double M_s = 1e9;          // Satellite mass [M_sun].
    double r_s = 0.8;          // Satellite radius [kpc].

    // Print initial statement.
    cout << "#############################################" << endl;
    cout << "### Running Galactic Potential Integrator ###" << endl;
    cout << "######### Created by Eric Andersson #########" << endl;
    cout << "#############################################" << endl;
    
    cout << " " << endl;
    cout << "Program will now set-up the suit." << endl;
    cout << " " << endl;

    // Set up what encounters to run.
    if (argc == 3){
        startRUN = atoi(argv[1]);
        finalRUN = atoi(argv[2]);
        cout << "Will run encounters RUN" << startRUN << " to RUN" << 
            finalRUN-1 << endl;
    } else{
        cout << "Encounters not specified. Will run all " << nruns << 
            " encounter, starting with RUN000." << endl;
        startRUN = 0;
        finalRUN = nruns;
    }
    
    cout << "Simulation set-up was successful. Will now start runs." 
        << endl;
    cout << " " << endl;
    
    // Start the simulations.
    time_t start = time(0);
    char* starttime = ctime(&start);
    cout << "Started at " << starttime << endl;
    
    for (int run = startRUN; run < finalRUN; run++) {
        cout << "Running simulation number " << run << endl;
        cout << " " << endl;
        
        int terminate = 0;      // Used later. Needs reset each run.
        
        // Set up directories.
        cout << "Creating directories..." << endl;
        stringstream directory;
        directory.str("");
        directory << "RUN";
        directory << setfill('0') << setw(3) << run;
        string dirstr = directory.str();
        stringstream strs;
        strs << "mkdir ./";
        strs << dirstr; 
        strs << "/data";
        string str = strs.str();
        const char *command = str.c_str();
        system(command);
        cout << "...Done." << endl;

        // Initial conditions.
        cout << "Setting up initial conditions for satellite galaxy..."
            << endl;
        
        // Allocate memory for initial conditions.
        vector<double> xp0;
        vector<double> yp0;
        vector<double> zp0;
        vector<double> vxp0;
        vector<double> vyp0;
        vector<double> vzp0;
        double xs0, ys0, zs0, vxs0, vys0, vzs0;
    
        // Read in dwarf initial conditions.
        cout << "Reading data from ./" << dirstr << 
            "/dwarf_IC.txt" << endl;
      
        // Read initial conditions file for satellite.
        int count = 0;
        double x_data, y_data, z_data, vx_data, vy_data, vz_data;
        strs.str("./");
        strs << dirstr;
        strs << "/dwarf_IC.txt";
        string initstr = strs.str();
        const char *filename = initstr.c_str();
        ifstream infile(filename);
        infile >> x_data >> y_data >> z_data >> 
                    vx_data >> vy_data >> vz_data;
        xs0 = x_data;
        ys0 = y_data;
        zs0 = z_data;
        vxs0 = vx_data;
        vys0 = vy_data;
        vzs0 = vz_data;
    
        // Initialize galaxies.
        cout << "Initializing galaxies..." << endl;
        Galaxy M31(M_b, r_b, M_d, A, B, r_h, V_h, delta, rho,  
                r_b_lum, h_d_lum, r_d_lum);
        Satellite S(M_s, r_s, xs0, ys0, zs0, vxs0, vys0, vzs0);

        // Print information about model.
        cout << "################ MODEL #################" << endl;
        cout << "Major galaxy model:" << endl; 
        
        // The Andromeda galaxy
        if (M31.model == 'B'){
            cout << "Using M31 potential model from Bekki (2001)" 
                << endl;
            cout << "Hernqvist potential (bulge)" << endl;
            cout << "Miamoto-Nagai potential (disc)" << endl;
            cout << "Logarithmic potential (halo) with 155 kpc cutoff" 
            << endl;
            cout << " " << endl; 
            cout << "Parameters used:" << endl;
            cout << "M_b = " << M31.M_b << " Msun" << endl;
            cout << "r_b = " << M31.r_b << " kpc" << endl;
            cout << "M_d = " << M31.M_d << " Msun" << endl;
            cout << "A = " << M31.A << " kpc" << endl;
            cout << "B = " << M31.B << " kpc" << endl;
            cout << "r_c = " << M31.r_c << " kpc" << endl;
            cout << "V_h = " << M31.V_h << " km/s" << endl;
            cout << "########################################" << endl;
        }else if (M31.model == 'G'){
            cout << "Using M31 potential model from Geehan (2006)" 
                << endl;
            cout << "Hernqvist potential (bulge)" << endl;
            cout << "Miamoto-Nagai potential (disc)" << endl;
            cout << "NFW potential (halo)" 
            << endl;
            cout << " " << endl; 
            cout << "Parameters used:" << endl;
            cout << "M_b = " << M31.M_b << " Msun" << endl;
            cout << "r_b = " << M31.r_b << " kpc" << endl;
            cout << "M_d = " << M31.M_d << " Msun" << endl;
            cout << "A = " << M31.A << " kpc" << endl;
            cout << "B = " << M31.B << " kpc" << endl;
            cout << "r_c = " << M31.r_c << " kpc" << endl;
            cout << "delta = " << M31.delta << endl;
            cout << "rho = " << M31.rho << " Msun/kpc^3" << endl;
            cout << "########################################" << endl;
        }
        // The dwarf galaxy
        cout << " " << endl;
        cout << "Satellite galaxy model:" << endl;
        cout << "Plummer potential (dSph)" << endl;
        cout << " " << endl;
        cout << "Parameters used:" << endl;
        cout << "M_s = " << S.M_s << " Msun" << endl;
        cout << "r_s = " << S.r_s << " kpc" << endl;
        cout << "########################################" << endl;
        cout << "...Done." << endl;
        cout << " " << endl;

        // Print information about data storing.
        cout << "Setting up files for storing data..." << endl;
        cout << "Parameters used in all encounters stored in info.txt"
            << endl;

        // Save data about parametrization of the encounter set.
        FILE * myInfoFile;
        myInfoFile = fopen("info.txt", "a");
        fprintf(myInfoFile, "%s%i\n", "RUN ", run);
        fprintf(myInfoFile, "%s\n", "Simulation information");
        if (adaptive_timestep) {
            fprintf(myInfoFile, "%s%s\n", 
                    "adaptive_timestep:", "true");
        }else {
            fprintf(myInfoFile, "%s%s\n", 
                    "adaptive_timestep = ", "false");
        }
        fprintf(myInfoFile, "%s%.3f\n", "dt (Myr) = ", dt);
        fprintf(myInfoFile, "%s%.0f\n", "dt_rec (Myr) = ", dt_rec);
        fprintf(myInfoFile, "%s%.0f\n", 
                "time_tot (Myr) = ", time_tot);
        fprintf(myInfoFile, "%s%.0f\n\n", "nruns = ", nruns);
        fprintf(myInfoFile, "%s\n%s\n", "Potential parameters",
                "Major galaxy");
        fprintf(myInfoFile, "%s%.0f\n", "M_b = ", M_b);
        fprintf(myInfoFile, "%s%.0f\n", "r_b = ", r_b);
        fprintf(myInfoFile, "%s%.0f\n", "M_d = ", M_d);
        fprintf(myInfoFile, "%s%.3f\n", "A = ", A);
        fprintf(myInfoFile, "%s%.3f\n", "B = ", B);
        fprintf(myInfoFile, "%s%.0f\n", "r_h = ", r_h);
        fprintf(myInfoFile, "%s%.0f\n", "V_h = ", V_h);
        fprintf(myInfoFile, "%s%.0f\n", "r_b_lum = ", r_b_lum);
        fprintf(myInfoFile, "%s%.0f\n", "h_d_lum = ", h_d_lum);
        fprintf(myInfoFile, "%s%.0f\n", "r_d_lum = ", r_d_lum);
        fprintf(myInfoFile, "%s%.0f\n", "V_h = ", V_h);
        fprintf(myInfoFile, "\n%s\n", "Satellite galaxy");
        fprintf(myInfoFile, "%s%.0f\n", "M_s = ", M_s);
        fprintf(myInfoFile, "%s%.0f\n", "r_s = ", r_s);
        fprintf(myInfoFile, "%s%.2f\n", "r_p = ", 
                sqrt(pow(S.x,2) + pow(S.y,2) + pow(S.z,2)));
        fprintf(myInfoFile, "%s%.5f\n", "v_max = ", 
                sqrt(pow(S.vx,2) + pow(S.vy,2) + pow(S.vz,2)));
        fprintf(myInfoFile, "%s\n\n", "##############################");
        fclose(myInfoFile);
        cout << "Done..." << endl;
        cout << " " << endl;

        // Integrate dwarf to initial position.
        cout << "Will now attempt to integrate dwarf galaxy to an " <<
            "intial distance of " << r0 << " kpc." << endl;
        double r = sqrt(pow(S.x,2) + pow(S.y,2) + pow(S.z,2));
        cout << "Integrating..." << endl;
        while (r < r0){
            S.update(M31, S, dt, S.t, adaptive_timestep);
            r = sqrt(pow(S.x,2) + pow(S.y,2) + pow(S.z,2));
            if (S.t > 14000){
                cout << "TimeError: Time of impact would exceed age "
                    "of Universe. Run was terminated." << endl;
                stringstream terminatestr;
                terminatestr << dirstr;
                terminatestr << "/TERMINATED.out";
                string terminated = terminatestr.str();
                const char *terminatefile = terminated.c_str();
                ofstream outfile (terminatefile);
                outfile.close();
                terminate = 1;
                break;
            }
        }
        if (terminate == 1) continue;
        cout << "...Done." << endl;

        // Reverse trajectory of dwarf.
        double t0 = S.t;
        S.t *= -1;
        S.vx *= -1;
        S.vy *= -1;
        S.vz *= -1;
            
        // Read in and initiate cluster population.
        cout << "Read in and initialize cluster population..." << endl;
        strs.str("./GC_sample.txt");
        string pinitstr = strs.str();
        const char *pfilename = pinitstr.c_str();
        ifstream pinfile(pfilename);
        while (pinfile >> x_data >> y_data >> z_data >> 
                    vx_data >> vy_data >> vz_data) {
            xp0.push_back(x_data + S.x);
            yp0.push_back(y_data + S.y);
            zp0.push_back(z_data + S.z);
            vxp0.push_back(vx_data + S.vx);
            vyp0.push_back(vy_data + S.vy);
            vzp0.push_back(vz_data + S.vz);
        }
        Particle p(xp0, yp0, zp0, vxp0, vyp0, vzp0, S.t);
        cout << "Done..." << endl;
        cout << "Initiated " << p.npar << " clusters." << endl;
        cout << " " << endl;

        // Run simulation of encounter.
        // Initialize variables that will be used later.
        double progress = 0;
        double current_progress = 0;
        double record = S.t;
        double collision = 0;
        double Ek, Ep;
        
        cout << "The program will now integrate the dwarf and " <<
            "cluster trajectories for " << time_tot + t0 << " Myr."
            << endl;
        cout << "The separation will be at a minimum " << t0 << 
            " Myr into the simulation." << endl;
        cout << "Integrating..." << endl;
        while (time_tot > p.t) {
            if (record <= p.t) {
                
                // Save data for particles.
                for(int i=0; i<p.npar; i++){
                    // Compute energy 
                    Ek = 0.5 * (pow(p.vx[i],2) + pow(p.vy[i],2) + 
                        pow(p.vz[i],2));
                    Ep = M31.potential(M31, p.x[i], p.y[i], p.z[i]) + 
                        S.potential(S, p.x[i], p.y[i], p.z[i]);
                    
                    // Set up directories.
                    strs.str("./");
                    strs << dirstr;
                    strs << "/data/particle_";
                    strs << i;
                    strs << ".txt";
                    string pstorestr = strs.str();
                    FILE * myFile;
                    myFile = fopen(pstorestr.c_str(), "a");
                    fprintf( myFile, 
                    "%.2f\t%.5f\t%.5f\t%.5f\t%.8f\t%.8f\t%.8f\t%.5f\t%.5f\n", 
                    p.t, p.x[i], p.y[i], p.z[i], p.vx[i]/K, p.vy[i]/K, 
                    p.vz[i]/K, Ek, Ep);
                    fclose(myFile);
                }

                // Compute satellite energy 
                Ek = 0.5 * (pow(S.vx,2) + pow(S.vy,2) + 
                        pow(S.vz,2));
                Ep = M31.potential(M31, S.x, S.y, S.z);
                
                // Save data for satellite.
                strs.str("./");
                strs << dirstr;
                strs << "/data/satellite.txt";
                string sstorestr = strs.str();
                FILE * myFile;
                myFile = fopen(sstorestr.c_str(), "a");
                fprintf( myFile, 
                "%.2f\t%.5f\t%.5f\t%.5f\t%.8f\t%.8f\t%.8f\t%.0f\t%.5f\t%.5f\n", 
                S.t, S.x, S.y, S.z, S.vx/K, S.vy/K, S.vz/K, collision, Ek, Ep);
                fclose(myFile);
                
                // Save timestep data.
                if (adaptive_timestep) {
                    strs.str("./");
                    strs << dirstr;
                    strs << "/data/timestep.txt";
                    FILE * timeStepData;
                    string dtstr = strs.str();
                    timeStepData = fopen(dtstr.c_str(), "a");
                    fprintf(timeStepData, "%.1f\t%.5f\n", p.t, dt);
                    fclose(timeStepData);
                }

                record += dt_rec;
            }
        
            // Update system.
            p.update(M31, S, dt, adaptive_timestep);            
            S.update(M31, S, dt, p.t, adaptive_timestep);
            
            if (collision == 0) {
                collision = M31.collision(M31, S.x, S.y, S.z);
            }

            // Print progress.
            progress = ((p.t + t0) / (time_tot + t0))*100.0;
            
            if (current_progress < progress){
                cout << "] Progress: " << int(progress) << " %\r";
                cout.flush();
                current_progress += 1;
            }
        }
        cout << endl;
        cout << "... Done." << endl;
        cout << "Simulation ended at t = " << p.t << " Myr" << endl;
    }
    time_t end = time(0);
    char* endtime = ctime(&end);
    cout << "Ended at " << endtime << endl;
}

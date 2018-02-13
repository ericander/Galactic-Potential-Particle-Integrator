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
    double time_tot = 14000;        // Total simulated time in Myr.
    double nruns = 1;               // Number of runs.
    int startRUN;
    int finalRUN;

    // Set parameters for the potential fields.
    // Major galaxy
    double M_b = 9.2e10;            // Bulge mass [M_sun]. 
    double r_b = 0.7;               // Bulge radius [kpc].
    double M_d = 1.3e11;            // Disc mass [M_sun].
    double A = 6.5;                 // Disc scale [kpc].
    double B = 0.26;                // Disc scale [kpc].
    double r_h = 12;                // Halo critical radius [kpc].
    double V_h = 186;               // Halo rotational velocity [km/s].
    double r_b_lum = 3;             // Luminous bulge radius [kpc].
    double h_d_lum = 1;             // Luminous disc thickness [kpc].
    double r_d_lum = 15;            // Luninous disc length [kpc].

    // Satellite galaxy
    double M_s = 1e9;               // Satellite mass [M_sun].
    double r_s = 0.8;               // Satellite radius [kpc].

    // Print initial statement.
    cout << "#############################################" << endl;
    cout << "### Running Galactic Potential Integrator ###" << endl;
    cout << "######### Created by Eric Andersson #########" << endl;
    cout << "#############################################" << endl;
    cout << "Simulations will run for " << time_tot << " Myr." << endl;

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

    double nSteps = time_tot / dt;
    
    for (int run = startRUN; run < finalRUN; run++) {
        // Set up directories.
        stringstream directory;
        directory.str("");
        directory << "RUN";
        directory << setfill('0') << setw(3) << run;
        string dirstr = directory.str();
        cout << "Running simulation number " << run << endl;
        stringstream strs;
        strs << "mkdir ./";
        strs << dirstr; 
        strs << "/data";
        string str = strs.str();
        const char *command = str.c_str();
        system(command);
         
        // Particle and satellite phase-space coordinates.
        vector<double> xp0;
        vector<double> yp0;
        vector<double> zp0;
        vector<double> vxp0;
        vector<double> vyp0;
        vector<double> vzp0;
        double xs0, ys0, zs0, vxs0, vys0, vzs0;
    
        // Read in data.
        cout << "Reading data from ./" << dirstr << 
            "/initial_conditions.txt ..." << endl;
      
        int count = 0;
        double x_data, y_data, z_data, vx_data, vy_data, vz_data;
        strs.str("./");
        strs << dirstr;
        strs << "/initial_conditions.txt";
        string initstr = strs.str();
        const char *filename = initstr.c_str();
        ifstream infile(filename);
        while (infile >> x_data >> y_data >> z_data >> 
                    vx_data >> vy_data >> vz_data) {
            if (count == 0) {
                xs0 = x_data;
                ys0 = y_data;
                zs0 = z_data;
                vxs0 = vx_data;
                vys0 = vy_data;
                vzs0 = vz_data;
            }
            else {
                xp0.push_back(x_data);
                yp0.push_back(y_data);
                zp0.push_back(z_data);
                vxp0.push_back(vx_data);
                vyp0.push_back(vy_data);
                vzp0.push_back(vz_data);
            }
            count++;
        }
        cout << "... Done." << endl; 
    
        // Initialize galaxies and globular clusters.
        cout << "Initializing galaxies and particles..." << endl;
        Galaxy M31(M_b, r_b, M_d, A, B, r_h, V_h, 
                r_b_lum, h_d_lum, r_d_lum);
        Satellite S(M_s, r_s, xs0, ys0, zs0, vxs0, vys0, vzs0);
        Particle p(xp0, yp0, zp0, vxp0, vyp0, vzp0);
        cout << "... Done." << endl;

        // Print information about model.
        cout << "################ MODEL #################" << endl;
        cout << "Initiated with major galaxy, satellite galaxy and "
                << p.npar << " particle(s)." << endl;
        cout << "Major galaxy model:" << endl; 
        cout << "Hernqvist potential (bulge)" << endl;
        cout << "Miamoto-Nagai potential (disc)" << endl;
        cout << "Logarithmic potential (halo) with 300 kpc cutoff" 
            << endl;
        cout << "Major galaxy (Andromeda) parameters:" << endl;
        cout << "M_b = " << M31.M_b << endl;
        cout << "r_b = " << M31.r_b << endl;
        cout << "M_d = " << M31.M_d << endl;
        cout << "A = " << M31.A << endl;
        cout << "B = " << M31.B << endl;
        cout << "r_c = " << M31.r_c << endl;
        cout << "V_h = " << M31.V_h << endl;
        cout << "########################################" << endl;
        cout << "Satellite galaxy model:" << endl;
        cout << "Plummer potential (dSph)" << endl;
        cout << "Satellite galaxy (dSph) parameters:" << endl;
        cout << "M_s = " << S.M_s << endl;
        cout << "r_s = " << S.r_s << endl;
        cout << "########################################" << endl;
        
        // Print information about data storing.
        cout << "Setting up files for storing data." << endl;
        cout << "Parameters used in all encounters stored in info.txt"
            << endl;

        // Save data about parametrization of the encounter set.
        
        FILE * myInfoFile;
        myInfoFile = fopen("info.txt", "a");
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
        fclose(myInfoFile);

        cout << "Data will be stored in ./" << dirstr << 
            "/data/ every " << dt_rec << " Myr." << endl;
        
        // Warm-up if using adaptive timestep
        if (adaptive_timestep) {
            cout << "Warming up simulation to adapt time-step." << endl;
            
            // Warm up to set initial timestep.
            for (int i=0; i<1000; i++){
                cout.setstate(ios_base::failbit);
                p.update(M31, S, dt, adaptive_timestep);            
                S.update(M31, S, dt, p.t, adaptive_timestep);
                cout.clear();
            }
            // Reset time.
            p.t = 0;
            S.t = 0;
        }

        // Run simulations.
        double progress = 0;
        double record = 0;
        double collision = 0;

        cout << "Simulating model ..." << endl;
        while (time_tot > p.t) {
            if (record <= p.t) {
                // Save data for particles.
                for(int i=0; i<p.npar; i++){
                    strs.str("./");
                    strs << dirstr;
                    strs << "/data/particle_";
                    strs << i;
                    strs << ".txt";
                    string pstorestr = strs.str();
                    FILE * myFile;
                    myFile = fopen(pstorestr.c_str(), "a");
                    fprintf( myFile, 
                    "%.2f\t%.5f\t%.5f\t%.5f\t%.8f\t%.8f\t%.8f\n", 
                    p.t, p.x[i], p.y[i], p.z[i], p.vx[i], p.vy[i], 
                    p.vz[i]);
                    fclose(myFile);
                }
                // Save data for satellite.
                strs.str("./");
                strs << dirstr;
                strs << "/data/satellite.txt";
                string sstorestr = strs.str();
                FILE * myFile;
                myFile = fopen(sstorestr.c_str(), "a");
                fprintf( myFile, 
                "%.2f\t%.5f\t%.5f\t%.5f\t%.8f\t%.8f\t%.8f\t%.0f\n", 
                S.t, S.x, S.y, S.z, S.vx, S.vy, S.vz, collision);
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
            progress = p.t / time_tot;    
            cout << "] Progress: " << int(progress * 100.0) << " %\r";
            cout.flush();
        }
        cout << endl;
        cout << "... Done." << endl;
        cout << "Simulation ended at t = " << p.t << " Myr" << endl;
    }
}

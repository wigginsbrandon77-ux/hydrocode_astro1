#ifndef HYDRO_SOLVER_3D_H
#define HYDRO_SOLVER_3D_H

#include <stdio.h>
#include <stdbool.h>
#include "input_deck.h"
#include "gravity_solver_3d.h"

//#define RHO_FLOOR 1.0e-3
#define RHO_FLOOR 1.0e-5
#define P_FLOOR 1.0e-6
#define IDX3D(i, j, k, ny, nz) (((i) * (ny) + (j)) * (nz) + (k))

// Structure for a single cell's conservative state variables (3D)
typedef struct {
    double rho;   // Density
    double rhou;  // x-momentum
    double rhov;  // y-momentum
    double rhow;  // z-momentum
    double E;     // Total Energy
} State3D;

// Structure for a single cell's primitive state variables (3D)
typedef struct {
    double rho; // Density
    double u;   // x-velocity
    double v;   // y-velocity
    double w;   // z-velocity
    double p;   // Pressure
} Primitive3D;

// Main 3D solver structure
typedef struct GodunovHLLC3D {
    int nx, ny, nz;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    double dx, dy, dz;
    double dt;
    double cfl;

    BCType bc_xmin, bc_xmax;
    BCType bc_ymin, bc_ymax;
    BCType bc_zmin, bc_zmax;
    
    ProblemType problem; 
    double gamma;        
    
    State3D *U;       // Current state
    State3D *U_half;  // Half-step state
    State3D *F_x;     // x-direction fluxes
    State3D *F_y;     // y-direction fluxes
    State3D *F_z;     // z-direction fluxes
    
    // Gravity field arrays
    double *gx;
    double *gy;
    double *gz;
    
    // Gravity solver
    GravitySolver3D *gravity_solver;
    bool enable_gravity;  // Flag to enable/disable gravity
    double G_const;

} GodunovHLLC3D;


// --- Public Function Prototypes ---

// Solver creation/destruction
GodunovHLLC3D* solver3d_create(int nx, int ny, int nz, 
                                double xmin, double xmax, 
                                double ymin, double ymax,
                                double zmin, double zmax,
                                double G);
void solver3d_destroy(GodunovHLLC3D *solver);

// Initialization
void set_initial_conditions_3d(GodunovHLLC3D *solver, ProblemType problem);
void set_initial_conditions_3d_from_deck(GodunovHLLC3D *solver, const InputDeck *deck);

// Main simulation loop
void run_simulation_3d(GodunovHLLC3D *solver, double t_end, double output_dt);

// Main simulation loop with restart capability
void run_simulation_3d_restart(GodunovHLLC3D *solver, double t_start, double t_end, 
                                double output_dt, int start_output_count);

// I/O functions
void write_output_to_file_3d(const GodunovHLLC3D *solver, int step);
void write_output_vtk(const GodunovHLLC3D *solver, int step);

// Restart function
int load_from_csv_restart(GodunovHLLC3D *solver, const char *filename);

// Convenience function to create from input deck
GodunovHLLC3D* solver3d_create_from_deck(const InputDeck *deck);

#endif // HYDRO_SOLVER_3D_H

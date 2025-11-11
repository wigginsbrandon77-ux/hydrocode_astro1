#include "hydro_solver_3d.h"
#include "gravity_solver_3d.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <sys/stat.h>
#include <sys/types.h>

// Helper macros
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define IDX3D(i, j, k, ny, nz) (((i) * (ny) + (j)) * (nz) + (k))

// Forward declarations of internal functions
static Primitive3D cons_to_prim_3d(State3D U, const GodunovHLLC3D *solver);
static State3D prim_to_cons_3d(Primitive3D P, const GodunovHLLC3D *solver);
static double sound_speed_3d(double rho, double p, const GodunovHLLC3D *solver);
static State3D compute_flux_x_3d(State3D U, Primitive3D P);
static State3D compute_flux_y_3d(State3D U, Primitive3D P);
static State3D compute_flux_z_3d(State3D U, Primitive3D P);
static State3D hllc_riemann_solver_x_3d(State3D UL, State3D UR, const GodunovHLLC3D *solver);
static State3D hllc_riemann_solver_y_3d(State3D UL, State3D UR, const GodunovHLLC3D *solver);
static State3D hllc_riemann_solver_z_3d(State3D UL, State3D UR, const GodunovHLLC3D *solver);
static void compute_dt_3d(GodunovHLLC3D *solver);
static void apply_boundary_conditions_3d(GodunovHLLC3D *solver);
static void godunov_step_3d(GodunovHLLC3D *solver);

// =============================================================================
// SOLVER CREATION AND DESTRUCTION
// =============================================================================

GodunovHLLC3D* solver3d_create(int nx, int ny, int nz, 
                                double xmin, double xmax, 
                                double ymin, double ymax,
                                double zmin, double zmax,
                                double G) {
    GodunovHLLC3D *solver = (GodunovHLLC3D*) malloc(sizeof(GodunovHLLC3D));
    if (!solver) {
        fprintf(stderr, "Failed to allocate solver structure.\n");
        return NULL;
    }
    
    solver->nx = nx;
    solver->ny = ny;
    solver->nz = nz;
    
    solver->xmin = xmin;
    solver->xmax = xmax;
    solver->ymin = ymin;
    solver->ymax = ymax;
    solver->zmin = zmin;
    solver->zmax = zmax;
    
    solver->dx = (xmax - xmin) / nx;
    solver->dy = (ymax - ymin) / ny;
    solver->dz = (zmax - zmin) / nz;
    
    solver->cfl = 0.4;
    solver->dt = 0.0;
    solver->gamma = 5.0 / 3.0;
    solver->G_const = G;
    solver->enable_gravity = (G > 0.0);  // Enable gravity if G > 0
    
    // Initialize gravity solver if needed
    solver->gravity_solver = NULL;
    if (solver->enable_gravity) {
        solver->gravity_solver = gravity_solver_create(nx, ny, nz,
                                                        xmin, xmax,
                                                        ymin, ymax,
                                                        zmin, zmax,
                                                        G);
        if (!solver->gravity_solver) {
            fprintf(stderr, "Warning: Failed to create gravity solver. Disabling gravity.\n");
            solver->enable_gravity = false;
        } else {
            // Set gravity solver to use isolated boundary conditions
            gravity_set_boundary_type(solver->gravity_solver, GRAV_BC_ISOLATED);
            printf("Gravity solver initialized with G = %.6e\n", G);
        }
    }
    
    // Allocate arrays
    int ncells = nx * ny * nz;
    
    solver->U = (State3D*) malloc(ncells * sizeof(State3D));
    solver->U_half = (State3D*) malloc(ncells * sizeof(State3D));
    solver->F_x = (State3D*) malloc(ncells * sizeof(State3D));
    solver->F_y = (State3D*) malloc(ncells * sizeof(State3D));
    solver->F_z = (State3D*) malloc(ncells * sizeof(State3D));
    
    solver->gx = (double*) calloc(ncells, sizeof(double));
    solver->gy = (double*) calloc(ncells, sizeof(double));
    solver->gz = (double*) calloc(ncells, sizeof(double));
    
    if (!solver->U || !solver->U_half || !solver->F_x || !solver->F_y || !solver->F_z ||
        !solver->gx || !solver->gy || !solver->gz) {
        fprintf(stderr, "Failed to allocate solver arrays.\n");
        solver3d_destroy(solver);
        return NULL;
    }
    
    // Initialize arrays to zero
    memset(solver->U, 0, ncells * sizeof(State3D));
    memset(solver->U_half, 0, ncells * sizeof(State3D));
    
    return solver;
}

void solver3d_destroy(GodunovHLLC3D *solver) {
    if (!solver) return;
    
    free(solver->U);
    free(solver->U_half);
    free(solver->F_x);
    free(solver->F_y);
    free(solver->F_z);
    free(solver->gx);
    free(solver->gy);
    free(solver->gz);
    
    // Destroy gravity solver if it exists
    if (solver->gravity_solver) {
        gravity_solver_destroy(solver->gravity_solver);
    }
    
    free(solver);
}

GodunovHLLC3D* solver3d_create_from_deck(const InputDeck *deck) {
    GodunovHLLC3D *solver = solver3d_create(
        deck->nx, deck->ny, deck->nz,
        deck->xmin, deck->xmax,
        deck->ymin, deck->ymax,
        deck->zmin, deck->zmax,
        deck->G_const
    );
    
    if (!solver) return NULL;
    
    solver->cfl = deck->cfl;
    solver->gamma = deck->gamma;
    solver->problem = deck->problem;
    
    solver->bc_xmin = deck->bc_xmin;
    solver->bc_xmax = deck->bc_xmax;
    solver->bc_ymin = deck->bc_ymin;
    solver->bc_ymax = deck->bc_ymax;
    solver->bc_zmin = deck->bc_zmin;
    solver->bc_zmax = deck->bc_zmax;
    
    // Gravity is automatically enabled if G > 0
    // Can be explicitly controlled via input deck if desired
    
    return solver;
}

// =============================================================================
// STATE CONVERSIONS
// =============================================================================

static Primitive3D cons_to_prim_3d(State3D U, const GodunovHLLC3D *solver) {
    Primitive3D P;
    
    P.rho = MAX(U.rho, RHO_FLOOR);
    P.u = U.rhou / P.rho;
    P.v = U.rhov / P.rho;
    P.w = U.rhow / P.rho;
    
    double ke = 0.5 * P.rho * (P.u * P.u + P.v * P.v + P.w * P.w);
    P.p = (solver->gamma - 1.0) * (U.E - ke);
    P.p = MAX(P.p, P_FLOOR);
    
    return P;
}

static State3D prim_to_cons_3d(Primitive3D P, const GodunovHLLC3D *solver) {
    State3D U;
    
    U.rho = P.rho;
    U.rhou = P.rho * P.u;
    U.rhov = P.rho * P.v;
    U.rhow = P.rho * P.w;
    
    double ke = 0.5 * P.rho * (P.u * P.u + P.v * P.v + P.w * P.w);
    U.E = P.p / (solver->gamma - 1.0) + ke;
    
    return U;
}

static double sound_speed_3d(double rho, double p, const GodunovHLLC3D *solver) {
    return sqrt(solver->gamma * p / rho);
}

// =============================================================================
// FLUX COMPUTATIONS
// =============================================================================

static State3D compute_flux_x_3d(State3D U, Primitive3D P) {
    State3D F;
    
    F.rho = U.rhou;
    F.rhou = U.rhou * P.u + P.p;
    F.rhov = U.rhov * P.u;
    F.rhow = U.rhow * P.u;
    F.E = (U.E + P.p) * P.u;
    
    return F;
}

static State3D compute_flux_y_3d(State3D U, Primitive3D P) {
    State3D F;
    
    F.rho = U.rhov;
    F.rhou = U.rhou * P.v;
    F.rhov = U.rhov * P.v + P.p;
    F.rhow = U.rhow * P.v;
    F.E = (U.E + P.p) * P.v;
    
    return F;
}

static State3D compute_flux_z_3d(State3D U, Primitive3D P) {
    State3D F;
    
    F.rho = U.rhow;
    F.rhou = U.rhou * P.w;
    F.rhov = U.rhov * P.w;
    F.rhow = U.rhow * P.w + P.p;
    F.E = (U.E + P.p) * P.w;
    
    return F;
}

// =============================================================================
// HLLC RIEMANN SOLVERS
// =============================================================================

static State3D hllc_riemann_solver_x_3d(State3D UL, State3D UR, const GodunovHLLC3D *solver) {
    Primitive3D PL = cons_to_prim_3d(UL, solver);
    Primitive3D PR = cons_to_prim_3d(UR, solver);
    
    double cL = sound_speed_3d(PL.rho, PL.p, solver);
    double cR = sound_speed_3d(PR.rho, PR.p, solver);
    
    // Wave speed estimates
    double SL = MIN(PL.u - cL, PR.u - cR);
    double SR = MAX(PL.u + cL, PR.u + cR);
    
    // If the wave structure is entirely left or right, return the appropriate flux
    if (SL >= 0.0) {
        return compute_flux_x_3d(UL, PL);
    }
    if (SR <= 0.0) {
        return compute_flux_x_3d(UR, PR);
    }
    
    // Compute contact wave speed
    double num = PR.p - PL.p + PL.rho * PL.u * (SL - PL.u) - PR.rho * PR.u * (SR - PR.u);
    double den = PL.rho * (SL - PL.u) - PR.rho * (SR - PR.u);
    double S_star = num / den;
    
    State3D F_hllc;
    
    if (S_star >= 0.0) {
        // Use left star state
        double factor = (SL - PL.u) / (SL - S_star);
        State3D UL_star;
        UL_star.rho = PL.rho * factor;
        UL_star.rhou = UL_star.rho * S_star;
        UL_star.rhov = UL_star.rho * PL.v;
        UL_star.rhow = UL_star.rho * PL.w;
        UL_star.E = UL_star.rho * (UL.E / PL.rho + (S_star - PL.u) * (S_star + PL.p / (PL.rho * (SL - PL.u))));
        
        State3D FL = compute_flux_x_3d(UL, PL);
        F_hllc.rho = FL.rho + SL * (UL_star.rho - UL.rho);
        F_hllc.rhou = FL.rhou + SL * (UL_star.rhou - UL.rhou);
        F_hllc.rhov = FL.rhov + SL * (UL_star.rhov - UL.rhov);
        F_hllc.rhow = FL.rhow + SL * (UL_star.rhow - UL.rhow);
        F_hllc.E = FL.E + SL * (UL_star.E - UL.E);
    } else {
        // Use right star state
        double factor = (SR - PR.u) / (SR - S_star);
        State3D UR_star;
        UR_star.rho = PR.rho * factor;
        UR_star.rhou = UR_star.rho * S_star;
        UR_star.rhov = UR_star.rho * PR.v;
        UR_star.rhow = UR_star.rho * PR.w;
        UR_star.E = UR_star.rho * (UR.E / PR.rho + (S_star - PR.u) * (S_star + PR.p / (PR.rho * (SR - PR.u))));
        
        State3D FR = compute_flux_x_3d(UR, PR);
        F_hllc.rho = FR.rho + SR * (UR_star.rho - UR.rho);
        F_hllc.rhou = FR.rhou + SR * (UR_star.rhou - UR.rhou);
        F_hllc.rhov = FR.rhov + SR * (UR_star.rhov - UR.rhov);
        F_hllc.rhow = FR.rhow + SR * (UR_star.rhow - UR.rhow);
        F_hllc.E = FR.E + SR * (UR_star.E - UR.E);
    }
    
    return F_hllc;
}

static State3D hllc_riemann_solver_y_3d(State3D UL, State3D UR, const GodunovHLLC3D *solver) {
    Primitive3D PL = cons_to_prim_3d(UL, solver);
    Primitive3D PR = cons_to_prim_3d(UR, solver);
    
    double cL = sound_speed_3d(PL.rho, PL.p, solver);
    double cR = sound_speed_3d(PR.rho, PR.p, solver);
    
    double SL = MIN(PL.v - cL, PR.v - cR);
    double SR = MAX(PL.v + cL, PR.v + cR);
    
    if (SL >= 0.0) {
        return compute_flux_y_3d(UL, PL);
    }
    if (SR <= 0.0) {
        return compute_flux_y_3d(UR, PR);
    }
    
    double num = PR.p - PL.p + PL.rho * PL.v * (SL - PL.v) - PR.rho * PR.v * (SR - PR.v);
    double den = PL.rho * (SL - PL.v) - PR.rho * (SR - PR.v);
    double S_star = num / den;
    
    State3D F_hllc;
    
    if (S_star >= 0.0) {
        double factor = (SL - PL.v) / (SL - S_star);
        State3D UL_star;
        UL_star.rho = PL.rho * factor;
        UL_star.rhou = UL_star.rho * PL.u;
        UL_star.rhov = UL_star.rho * S_star;
        UL_star.rhow = UL_star.rho * PL.w;
        UL_star.E = UL_star.rho * (UL.E / PL.rho + (S_star - PL.v) * (S_star + PL.p / (PL.rho * (SL - PL.v))));
        
        State3D FL = compute_flux_y_3d(UL, PL);
        F_hllc.rho = FL.rho + SL * (UL_star.rho - UL.rho);
        F_hllc.rhou = FL.rhou + SL * (UL_star.rhou - UL.rhou);
        F_hllc.rhov = FL.rhov + SL * (UL_star.rhov - UL.rhov);
        F_hllc.rhow = FL.rhow + SL * (UL_star.rhow - UL.rhow);
        F_hllc.E = FL.E + SL * (UL_star.E - UL.E);
    } else {
        double factor = (SR - PR.v) / (SR - S_star);
        State3D UR_star;
        UR_star.rho = PR.rho * factor;
        UR_star.rhou = UR_star.rho * PR.u;
        UR_star.rhov = UR_star.rho * S_star;
        UR_star.rhow = UR_star.rho * PR.w;
        UR_star.E = UR_star.rho * (UR.E / PR.rho + (S_star - PR.v) * (S_star + PR.p / (PR.rho * (SR - PR.v))));
        
        State3D FR = compute_flux_y_3d(UR, PR);
        F_hllc.rho = FR.rho + SR * (UR_star.rho - UR.rho);
        F_hllc.rhou = FR.rhou + SR * (UR_star.rhou - UR.rhou);
        F_hllc.rhov = FR.rhov + SR * (UR_star.rhov - UR.rhov);
        F_hllc.rhow = FR.rhow + SR * (UR_star.rhow - UR.rhow);
        F_hllc.E = FR.E + SR * (UR_star.E - UR.E);
    }
    
    return F_hllc;
}

static State3D hllc_riemann_solver_z_3d(State3D UL, State3D UR, const GodunovHLLC3D *solver) {
    Primitive3D PL = cons_to_prim_3d(UL, solver);
    Primitive3D PR = cons_to_prim_3d(UR, solver);
    
    double cL = sound_speed_3d(PL.rho, PL.p, solver);
    double cR = sound_speed_3d(PR.rho, PR.p, solver);
    
    double SL = MIN(PL.w - cL, PR.w - cR);
    double SR = MAX(PL.w + cL, PR.w + cR);
    
    if (SL >= 0.0) {
        return compute_flux_z_3d(UL, PL);
    }
    if (SR <= 0.0) {
        return compute_flux_z_3d(UR, PR);
    }
    
    double num = PR.p - PL.p + PL.rho * PL.w * (SL - PL.w) - PR.rho * PR.w * (SR - PR.w);
    double den = PL.rho * (SL - PL.w) - PR.rho * (SR - PR.w);
    double S_star = num / den;
    
    State3D F_hllc;
    
    if (S_star >= 0.0) {
        double factor = (SL - PL.w) / (SL - S_star);
        State3D UL_star;
        UL_star.rho = PL.rho * factor;
        UL_star.rhou = UL_star.rho * PL.u;
        UL_star.rhov = UL_star.rho * PL.v;
        UL_star.rhow = UL_star.rho * S_star;
        UL_star.E = UL_star.rho * (UL.E / PL.rho + (S_star - PL.w) * (S_star + PL.p / (PL.rho * (SL - PL.w))));
        
        State3D FL = compute_flux_z_3d(UL, PL);
        F_hllc.rho = FL.rho + SL * (UL_star.rho - UL.rho);
        F_hllc.rhou = FL.rhou + SL * (UL_star.rhou - UL.rhou);
        F_hllc.rhov = FL.rhov + SL * (UL_star.rhov - UL.rhov);
        F_hllc.rhow = FL.rhow + SL * (UL_star.rhow - UL.rhow);
        F_hllc.E = FL.E + SL * (UL_star.E - UL.E);
    } else {
        double factor = (SR - PR.w) / (SR - S_star);
        State3D UR_star;
        UR_star.rho = PR.rho * factor;
        UR_star.rhou = UR_star.rho * PR.u;
        UR_star.rhov = UR_star.rho * PR.v;
        UR_star.rhow = UR_star.rho * S_star;
        UR_star.E = UR_star.rho * (UR.E / PR.rho + (S_star - PR.w) * (S_star + PR.p / (PR.rho * (SR - PR.w))));
        
        State3D FR = compute_flux_z_3d(UR, PR);
        F_hllc.rho = FR.rho + SR * (UR_star.rho - UR.rho);
        F_hllc.rhou = FR.rhou + SR * (UR_star.rhou - UR.rhou);
        F_hllc.rhov = FR.rhov + SR * (UR_star.rhov - UR.rhov);
        F_hllc.rhow = FR.rhow + SR * (UR_star.rhow - UR.rhow);
        F_hllc.E = FR.E + SR * (UR_star.E - UR.E);
    }
    
    return F_hllc;
}

// =============================================================================
// TIME STEP COMPUTATION
// =============================================================================

static void compute_dt_3d(GodunovHLLC3D *solver) {
    double max_speed = 1.0e-10;
    
    for (int i = 0; i < solver->nx; i++) {
        for (int j = 0; j < solver->ny; j++) {
            for (int k = 0; k < solver->nz; k++) {
                int idx = IDX3D(i, j, k, solver->ny, solver->nz);
                Primitive3D P = cons_to_prim_3d(solver->U[idx], solver);
                double cs = sound_speed_3d(P.rho, P.p, solver);
                
                double speed = fabs(P.u) + cs;
                max_speed = MAX(max_speed, speed);
                
                speed = fabs(P.v) + cs;
                max_speed = MAX(max_speed, speed);
                
                speed = fabs(P.w) + cs;
                max_speed = MAX(max_speed, speed);
            }
        }
    }
    
    double dt_x = solver->dx / max_speed;
    double dt_y = solver->dy / max_speed;
    double dt_z = solver->dz / max_speed;
    double dt_min = MIN(dt_x, MIN(dt_y, dt_z));
    
    solver->dt = solver->cfl * dt_min;
}

// =============================================================================
// BOUNDARY CONDITIONS
// =============================================================================

static void apply_boundary_conditions_3d(GodunovHLLC3D *solver) {
    // X-direction boundaries
    for (int j = 0; j < solver->ny; j++) {
        for (int k = 0; k < solver->nz; k++) {
            int idx_inner = IDX3D(0, j, k, solver->ny, solver->nz);
            int idx_outer = IDX3D(solver->nx - 1, j, k, solver->ny, solver->nz);
            
            // X-min boundary
            if (solver->bc_xmin == BC_PERIODIC) {
                // Already handled by copying appropriate cells
            } else if (solver->bc_xmin == BC_REFLECTIVE) {
                Primitive3D P = cons_to_prim_3d(solver->U[idx_inner], solver);
                P.u = -P.u; // Reflect x-velocity
                solver->U[idx_inner] = prim_to_cons_3d(P, solver);
            } else if (solver->bc_xmin == BC_OUTFLOW) {
                // Zero-gradient extrapolation: copy from interior
                int idx_boundary = IDX3D(0, j, k, solver->ny, solver->nz);
                int idx_interior = IDX3D(1, j, k, solver->ny, solver->nz);
                solver->U[idx_boundary] = solver->U[idx_interior];
            }
            
            // X-max boundary
            if (solver->bc_xmax == BC_PERIODIC) {
                // Already handled
            } else if (solver->bc_xmax == BC_REFLECTIVE) {
                Primitive3D P = cons_to_prim_3d(solver->U[idx_outer], solver);
                P.u = -P.u;
                solver->U[idx_outer] = prim_to_cons_3d(P, solver);
            } else if (solver->bc_xmax == BC_OUTFLOW) {
                // Zero-gradient extrapolation: copy from interior
                int idx_boundary = IDX3D(solver->nx - 1, j, k, solver->ny, solver->nz);
                int idx_interior = IDX3D(solver->nx - 2, j, k, solver->ny, solver->nz);
                solver->U[idx_boundary] = solver->U[idx_interior];
            }
        }
    }
    
    // Y-direction boundaries
    for (int i = 0; i < solver->nx; i++) {
        for (int k = 0; k < solver->nz; k++) {
            int idx_inner = IDX3D(i, 0, k, solver->ny, solver->nz);
            int idx_outer = IDX3D(i, solver->ny - 1, k, solver->ny, solver->nz);
            
            if (solver->bc_ymin == BC_REFLECTIVE) {
                Primitive3D P = cons_to_prim_3d(solver->U[idx_inner], solver);
                P.v = -P.v;
                solver->U[idx_inner] = prim_to_cons_3d(P, solver);
            } else if (solver->bc_ymin == BC_OUTFLOW) {
                // Zero-gradient extrapolation: copy from interior
                int idx_boundary = IDX3D(i, 0, k, solver->ny, solver->nz);
                int idx_interior = IDX3D(i, 1, k, solver->ny, solver->nz);
                solver->U[idx_boundary] = solver->U[idx_interior];
            }
            
            if (solver->bc_ymax == BC_REFLECTIVE) {
                Primitive3D P = cons_to_prim_3d(solver->U[idx_outer], solver);
                P.v = -P.v;
                solver->U[idx_outer] = prim_to_cons_3d(P, solver);
            } else if (solver->bc_ymax == BC_OUTFLOW) {
                // Zero-gradient extrapolation: copy from interior
                int idx_boundary = IDX3D(i, solver->ny - 1, k, solver->ny, solver->nz);
                int idx_interior = IDX3D(i, solver->ny - 2, k, solver->ny, solver->nz);
                solver->U[idx_boundary] = solver->U[idx_interior];
            }
        }
    }
    
    // Z-direction boundaries
    for (int i = 0; i < solver->nx; i++) {
        for (int j = 0; j < solver->ny; j++) {
            int idx_inner = IDX3D(i, j, 0, solver->ny, solver->nz);
            int idx_outer = IDX3D(i, j, solver->nz - 1, solver->ny, solver->nz);
            
            if (solver->bc_zmin == BC_REFLECTIVE) {
                Primitive3D P = cons_to_prim_3d(solver->U[idx_inner], solver);
                P.w = -P.w;
                solver->U[idx_inner] = prim_to_cons_3d(P, solver);
            } else if (solver->bc_zmin == BC_OUTFLOW) {
                // Zero-gradient extrapolation: copy from interior
                int idx_boundary = IDX3D(i, j, 0, solver->ny, solver->nz);
                int idx_interior = IDX3D(i, j, 1, solver->ny, solver->nz);
                solver->U[idx_boundary] = solver->U[idx_interior];
            }
            
            if (solver->bc_zmax == BC_REFLECTIVE) {
                Primitive3D P = cons_to_prim_3d(solver->U[idx_outer], solver);
                P.w = -P.w;
                solver->U[idx_outer] = prim_to_cons_3d(P, solver);
            } else if (solver->bc_zmax == BC_OUTFLOW) {
                // Zero-gradient extrapolation: copy from interior
                int idx_boundary = IDX3D(i, j, solver->nz - 1, solver->ny, solver->nz);
                int idx_interior = IDX3D(i, j, solver->nz - 2, solver->ny, solver->nz);
                solver->U[idx_boundary] = solver->U[idx_interior];
            }
        }
    }
}

// =============================================================================
// GRAVITY SOLVER
// =============================================================================

static void solve_gravity_field(GodunovHLLC3D *solver) {
    if (!solver->enable_gravity || !solver->gravity_solver) {
        return;
    }
    
    // Extract density field from conservative variables
    int ncells = solver->nx * solver->ny * solver->nz;
    double *rho_field = (double*) malloc(ncells * sizeof(double));
    
    for (int i = 0; i < solver->nx; i++) {
        for (int j = 0; j < solver->ny; j++) {
            for (int k = 0; k < solver->nz; k++) {
                int idx = IDX3D(i, j, k, solver->ny, solver->nz);
                rho_field[idx] = solver->U[idx].rho;
            }
        }
    }
    
    // Solve Poisson equation for gravitational potential
    int status = gravity_solve(solver->gravity_solver, rho_field);
    
    if (status != 0) {
        fprintf(stderr, "Warning: Gravity solver did not converge fully.\n");
    }
    
    // Copy gravitational acceleration to hydro solver arrays
    for (int idx = 0; idx < ncells; idx++) {
        solver->gx[idx] = solver->gravity_solver->gx[idx];
        solver->gy[idx] = solver->gravity_solver->gy[idx];
        solver->gz[idx] = solver->gravity_solver->gz[idx];
    }
    
    free(rho_field);
}

// =============================================================================
// GRAVITY SOURCE TERMS
// =============================================================================

static void apply_gravity_source_terms(GodunovHLLC3D *solver, double dt) {
    if (!solver->enable_gravity) {
        return;
    }
    
    // Apply gravity source terms to momentum and energy
    // d(rho u)/dt += rho * g
    // dE/dt += rho * (u·g)
    
    for (int i = 0; i < solver->nx; i++) {
        for (int j = 0; j < solver->ny; j++) {
            for (int k = 0; k < solver->nz; k++) {
                int idx = IDX3D(i, j, k, solver->ny, solver->nz);
                
                double rho = solver->U[idx].rho;
                double u = solver->U[idx].rhou / rho;
                double v = solver->U[idx].rhov / rho;
                double w = solver->U[idx].rhow / rho;
                
                double gx = solver->gx[idx];
                double gy = solver->gy[idx];
                double gz = solver->gz[idx];
                
                // Update momentum
                solver->U[idx].rhou += dt * rho * gx;
                solver->U[idx].rhov += dt * rho * gy;
                solver->U[idx].rhow += dt * rho * gz;
                
                // Update energy: dE/dt = rho * (v · g)
                double work = rho * (u * gx + v * gy + w * gz);
                solver->U[idx].E += dt * work;
            }
        }
    }
}

// =============================================================================
// GODUNOV TIME STEP
// =============================================================================

static void godunov_step_3d(GodunovHLLC3D *solver) {
    // Solve for gravitational field if gravity is enabled
    if (solver->enable_gravity) {
        solve_gravity_field(solver);
    }
    
    // Compute timestep
    compute_dt_3d(solver);
    
    const double dt = solver->dt;
    const double dx = solver->dx;
    const double dy = solver->dy;
    const double dz = solver->dz;
    
    // =========================================================================
    // STEP 1: Compute x-direction fluxes
    // =========================================================================
    for (int i = 0; i < solver->nx - 1; i++) {
        for (int j = 0; j < solver->ny; j++) {
            for (int k = 0; k < solver->nz; k++) {
                int idx_L = IDX3D(i, j, k, solver->ny, solver->nz);
                int idx_R = IDX3D(i + 1, j, k, solver->ny, solver->nz);
                
                solver->F_x[idx_L] = hllc_riemann_solver_x_3d(solver->U[idx_L], solver->U[idx_R], solver);
            }
        }
    }
    
    // Handle periodic BC for x-direction fluxes
    if (solver->bc_xmin == BC_PERIODIC && solver->bc_xmax == BC_PERIODIC) {
        for (int j = 0; j < solver->ny; j++) {
            for (int k = 0; k < solver->nz; k++) {
                int idx_last = IDX3D(solver->nx - 1, j, k, solver->ny, solver->nz);
                int idx_first = IDX3D(0, j, k, solver->ny, solver->nz);
                solver->F_x[idx_last] = hllc_riemann_solver_x_3d(solver->U[idx_last], solver->U[idx_first], solver);
            }
        }
    }
    
    // =========================================================================
    // STEP 2: Compute y-direction fluxes
    // =========================================================================
    for (int i = 0; i < solver->nx; i++) {
        for (int j = 0; j < solver->ny - 1; j++) {
            for (int k = 0; k < solver->nz; k++) {
                int idx_L = IDX3D(i, j, k, solver->ny, solver->nz);
                int idx_R = IDX3D(i, j + 1, k, solver->ny, solver->nz);
                
                solver->F_y[idx_L] = hllc_riemann_solver_y_3d(solver->U[idx_L], solver->U[idx_R], solver);
            }
        }
    }
    
    // Handle periodic BC for y-direction fluxes
    if (solver->bc_ymin == BC_PERIODIC && solver->bc_ymax == BC_PERIODIC) {
        for (int i = 0; i < solver->nx; i++) {
            for (int k = 0; k < solver->nz; k++) {
                int idx_last = IDX3D(i, solver->ny - 1, k, solver->ny, solver->nz);
                int idx_first = IDX3D(i, 0, k, solver->ny, solver->nz);
                solver->F_y[idx_last] = hllc_riemann_solver_y_3d(solver->U[idx_last], solver->U[idx_first], solver);
            }
        }
    }
    
    // =========================================================================
    // STEP 3: Compute z-direction fluxes
    // =========================================================================
    for (int i = 0; i < solver->nx; i++) {
        for (int j = 0; j < solver->ny; j++) {
            for (int k = 0; k < solver->nz - 1; k++) {
                int idx_L = IDX3D(i, j, k, solver->ny, solver->nz);
                int idx_R = IDX3D(i, j, k + 1, solver->ny, solver->nz);
                
                solver->F_z[idx_L] = hllc_riemann_solver_z_3d(solver->U[idx_L], solver->U[idx_R], solver);
            }
        }
    }
    
    // Handle periodic BC for z-direction fluxes
    if (solver->bc_zmin == BC_PERIODIC && solver->bc_zmax == BC_PERIODIC) {
        for (int i = 0; i < solver->nx; i++) {
            for (int j = 0; j < solver->ny; j++) {
                int idx_last = IDX3D(i, j, solver->nz - 1, solver->ny, solver->nz);
                int idx_first = IDX3D(i, j, 0, solver->ny, solver->nz);
                solver->F_z[idx_last] = hllc_riemann_solver_z_3d(solver->U[idx_last], solver->U[idx_first], solver);
            }
        }
    }
    
    // =========================================================================
    // STEP 4: Update conserved variables
    // =========================================================================
    for (int i = 1; i < solver->nx - 1; i++) {
        for (int j = 1; j < solver->ny - 1; j++) {
            for (int k = 1; k < solver->nz - 1; k++) {
                int idx = IDX3D(i, j, k, solver->ny, solver->nz);
                int idx_im1 = IDX3D(i - 1, j, k, solver->ny, solver->nz);
                int idx_jm1 = IDX3D(i, j - 1, k, solver->ny, solver->nz);
                int idx_km1 = IDX3D(i, j, k - 1, solver->ny, solver->nz);
                
                State3D dU;
                dU.rho = -(dt / dx) * (solver->F_x[idx].rho - solver->F_x[idx_im1].rho)
                         -(dt / dy) * (solver->F_y[idx].rho - solver->F_y[idx_jm1].rho)
                         -(dt / dz) * (solver->F_z[idx].rho - solver->F_z[idx_km1].rho);
                
                dU.rhou = -(dt / dx) * (solver->F_x[idx].rhou - solver->F_x[idx_im1].rhou)
                          -(dt / dy) * (solver->F_y[idx].rhou - solver->F_y[idx_jm1].rhou)
                          -(dt / dz) * (solver->F_z[idx].rhou - solver->F_z[idx_km1].rhou);
                
                dU.rhov = -(dt / dx) * (solver->F_x[idx].rhov - solver->F_x[idx_im1].rhov)
                          -(dt / dy) * (solver->F_y[idx].rhov - solver->F_y[idx_jm1].rhov)
                          -(dt / dz) * (solver->F_z[idx].rhov - solver->F_z[idx_km1].rhov);
                
                dU.rhow = -(dt / dx) * (solver->F_x[idx].rhow - solver->F_x[idx_im1].rhow)
                          -(dt / dy) * (solver->F_y[idx].rhow - solver->F_y[idx_jm1].rhow)
                          -(dt / dz) * (solver->F_z[idx].rhow - solver->F_z[idx_km1].rhow);
                
                dU.E = -(dt / dx) * (solver->F_x[idx].E - solver->F_x[idx_im1].E)
                       -(dt / dy) * (solver->F_y[idx].E - solver->F_y[idx_jm1].E)
                       -(dt / dz) * (solver->F_z[idx].E - solver->F_z[idx_km1].E);
                
                solver->U[idx].rho += dU.rho;
                solver->U[idx].rhou += dU.rhou;
                solver->U[idx].rhov += dU.rhov;
                solver->U[idx].rhow += dU.rhow;
                solver->U[idx].E += dU.E;
                
                // Apply floors
                solver->U[idx].rho = MAX(solver->U[idx].rho, RHO_FLOOR);
                Primitive3D P = cons_to_prim_3d(solver->U[idx], solver);
                P.p = MAX(P.p, P_FLOOR);
                solver->U[idx] = prim_to_cons_3d(P, solver);
            }
        }
    }
    
    // Apply gravity source terms
    if (solver->enable_gravity) {
        apply_gravity_source_terms(solver, dt);
    }
    
    // Handle periodic boundaries
    if (solver->bc_xmin == BC_PERIODIC && solver->bc_xmax == BC_PERIODIC) {
        for (int j = 0; j < solver->ny; j++) {
            for (int k = 0; k < solver->nz; k++) {
                int idx_0 = IDX3D(0, j, k, solver->ny, solver->nz);
                int idx_nx = IDX3D(solver->nx - 1, j, k, solver->ny, solver->nz);
                solver->U[idx_0] = solver->U[idx_nx];
            }
        }
    }
    
    if (solver->bc_ymin == BC_PERIODIC && solver->bc_ymax == BC_PERIODIC) {
        for (int i = 0; i < solver->nx; i++) {
            for (int k = 0; k < solver->nz; k++) {
                int idx_0 = IDX3D(i, 0, k, solver->ny, solver->nz);
                int idx_ny = IDX3D(i, solver->ny - 1, k, solver->ny, solver->nz);
                solver->U[idx_0] = solver->U[idx_ny];
            }
        }
    }
    
    if (solver->bc_zmin == BC_PERIODIC && solver->bc_zmax == BC_PERIODIC) {
        for (int i = 0; i < solver->nx; i++) {
            for (int j = 0; j < solver->ny; j++) {
                int idx_0 = IDX3D(i, j, 0, solver->ny, solver->nz);
                int idx_nz = IDX3D(i, j, solver->nz - 1, solver->ny, solver->nz);
                solver->U[idx_0] = solver->U[idx_nz];
            }
        }
    }
    
    // Apply other boundary conditions
    apply_boundary_conditions_3d(solver);
}

// =============================================================================
// STELLAR PROFILE READER
// =============================================================================

typedef struct {
    int n_points;
    double *r;      // Radius
    double *rho;    // Density
    double *u;      // Specific internal energy
    double R_star;  // Stellar radius
} StellarProfile1D;

static StellarProfile1D* read_stellar_profile(const char *filename) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Error: Could not open stellar profile file '%s'\n", filename);
        return NULL;
    }
    
    // Count lines (skip header lines starting with #)
    int n_points = 0;
    char line[512];
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] != '#' && strlen(line) > 5) {
            n_points++;
        }
    }
    rewind(fp);
    
    printf("Reading stellar profile from %s...\n", filename);
    printf("  Found %d data points\n", n_points);
    
    // Allocate profile structure
    StellarProfile1D *profile = (StellarProfile1D*) malloc(sizeof(StellarProfile1D));
    profile->n_points = n_points;
    profile->r = (double*) malloc(n_points * sizeof(double));
    profile->rho = (double*) malloc(n_points * sizeof(double));
    profile->u = (double*) malloc(n_points * sizeof(double));
    
    // Read data
    int i = 0;
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] != '#' && strlen(line) > 5) {
            sscanf(line, "%lf %lf %lf", &profile->r[i], &profile->rho[i], &profile->u[i]);
            i++;
        }
    }
    fclose(fp);
    
    // Find stellar radius (where density drops significantly)
    double rho_max = profile->rho[0];
    double rho_threshold = 1e-4 * rho_max;
    profile->R_star = profile->r[n_points-1];  // Default to max radius
    
    for (int j = 0; j < n_points; j++) {
        if (profile->rho[j] < rho_threshold) {
            profile->R_star = profile->r[j];
            break;
        }
    }
    
    printf("  Stellar radius: R = %.4f\n", profile->R_star);
    printf("  Central density: ρ_c = %.4e\n", profile->rho[0]);
    printf("  Central internal energy: u_c = %.4e\n", profile->u[0]);
    
    return profile;
}

static void destroy_stellar_profile(StellarProfile1D *profile) {
    if (!profile) return;
    free(profile->r);
    free(profile->rho);
    free(profile->u);
    free(profile);
}

static void interpolate_profile(const StellarProfile1D *profile, double r, 
                                 double *rho_out, double *u_out) {
    // Linear interpolation of 1D profile
    
    if (r <= profile->r[0]) {
        // Inside innermost point - use central value
        *rho_out = profile->rho[0];
        *u_out = profile->u[0];
        return;
    }
    
    if (r >= profile->r[profile->n_points - 1]) {
        // Outside profile - use atmosphere
        *rho_out = 1e-6;
        *u_out = 1e-6;
        return;
    }
    
    // Find bracketing points
    int i_low = 0;
    for (int i = 0; i < profile->n_points - 1; i++) {
        if (r >= profile->r[i] && r < profile->r[i+1]) {
            i_low = i;
            break;
        }
    }
    
    int i_high = i_low + 1;
    
    // Linear interpolation
    double r_low = profile->r[i_low];
    double r_high = profile->r[i_high];
    double f = (r - r_low) / (r_high - r_low);
    
    *rho_out = profile->rho[i_low] + f * (profile->rho[i_high] - profile->rho[i_low]);
    *u_out = profile->u[i_low] + f * (profile->u[i_high] - profile->u[i_low]);
}

static void set_initial_conditions_stellar_profile(GodunovHLLC3D *solver,
                                                    const char *profile_filename) {
    printf("Setting up stellar profile from file...\n");
    
    // Read 1D profile
    StellarProfile1D *profile = read_stellar_profile(profile_filename);
    if (!profile) {
        fprintf(stderr, "Failed to read profile. Using uniform density instead.\n");
        // Fallback: uniform sphere
        for (int i = 0; i < solver->nx; i++) {
            for (int j = 0; j < solver->ny; j++) {
                for (int k = 0; k < solver->nz; k++) {
                    Primitive3D P;
                    P.rho = 1.0;
                    P.u = 0.0;
                    P.v = 0.0;
                    P.w = 0.0;
                    P.p = 1.0;
                    int idx = IDX3D(i, j, k, solver->ny, solver->nz);
                    solver->U[idx] = prim_to_cons_3d(P, solver);
                }
            }
        }
        return;
    }
    
    // Map 1D profile onto 3D grid
    // Assume star is centered at domain center
    double center_x = 0.5 * (solver->xmin + solver->xmax);
    double center_y = 0.5 * (solver->ymin + solver->ymax);
    double center_z = 0.5 * (solver->zmin + solver->zmax);
    
    printf("Mapping profile to 3D grid...\n");
    printf("  Star center: (%.3f, %.3f, %.3f)\n", center_x, center_y, center_z);
    
    for (int i = 0; i < solver->nx; i++) {
        for (int j = 0; j < solver->ny; j++) {
            for (int k = 0; k < solver->nz; k++) {
                double x = solver->xmin + (i + 0.5) * solver->dx;
                double y = solver->ymin + (j + 0.5) * solver->dy;
                double z = solver->zmin + (k + 0.5) * solver->dz;
                
                // Distance from star center
                double r = sqrt(pow(x - center_x, 2) + 
                              pow(y - center_y, 2) + 
                              pow(z - center_z, 2));
                
                // Interpolate density and internal energy from profile
                double rho, u;
                interpolate_profile(profile, r, &rho, &u);
                
                // Set primitive variables
                Primitive3D P;
                P.rho = rho;
                P.u = 0.0;  // No initial velocity
                P.v = 0.0;
                P.w = 0.0;
                
                // Pressure from ideal gas: P = (γ-1) ρ u
                P.p = (solver->gamma - 1.0) * rho * u;
                
                // Ensure floors
                P.rho = MAX(P.rho, RHO_FLOOR);
                P.p = MAX(P.p, P_FLOOR);
                
                int idx = IDX3D(i, j, k, solver->ny, solver->nz);
                solver->U[idx] = prim_to_cons_3d(P, solver);
            }
        }
    }
    
    printf("Stellar profile mapped to grid successfully.\n");
    
    // Clean up
    destroy_stellar_profile(profile);
}

// =============================================================================
// BINARY STAR SYSTEM SETUP WITH TWO DIFFERENT PROFILES
// =============================================================================

static void set_initial_conditions_binary_stars_dual(GodunovHLLC3D *solver,
                                                      const char *profile_filename1,
                                                      const char *profile_filename2,
                                                      double x1, double y1, double z1,
                                                      double x2, double y2, double z2,
                                                      double vx1, double vy1, double vz1,
                                                      double vx2, double vy2, double vz2,
                                                      bool auto_orbit,
                                                      double separation) {
    printf("Setting up binary star system with two different profiles...\n");
    
    // Read both profiles
    StellarProfile1D *profile1 = read_stellar_profile(profile_filename1);
    StellarProfile1D *profile2 = read_stellar_profile(profile_filename2);
    
    if (!profile1 || !profile2) {
        fprintf(stderr, "Failed to read one or both profiles. Cannot create binary system.\n");
        if (profile1) destroy_stellar_profile(profile1);
        if (profile2) destroy_stellar_profile(profile2);
        return;
    }
    
    printf("Profile 1 (%s):\n", profile_filename1);
    printf("  Stellar radius: R1 = %.4f\n", profile1->R_star);
    printf("  Central density: ρ_c,1 = %.4e\n", profile1->rho[0]);
    printf("  Central internal energy: u_c,1 = %.4e\n", profile1->u[0]);
    
    printf("Profile 2 (%s):\n", profile_filename2);
    printf("  Stellar radius: R2 = %.4f\n", profile2->R_star);
    printf("  Central density: ρ_c,2 = %.4e\n", profile2->rho[0]);
    printf("  Central internal energy: u_c,2 = %.4e\n", profile2->u[0]);
    
    // Calculate stellar masses by integrating density profile
    double M_star1 = 0.0;
    for (int i = 0; i < profile1->n_points - 1; i++) {
        double r1 = profile1->r[i];
        double r2 = profile1->r[i+1];
        double rho1 = profile1->rho[i];
        double rho2 = profile1->rho[i+1];
        double dr = r2 - r1;
        double r_mid = 0.5 * (r1 + r2);
        double rho_mid = 0.5 * (rho1 + rho2);
        double dM = 4.0 * M_PI * rho_mid * r_mid * r_mid * dr;
        M_star1 += dM;
    }
    
    double M_star2 = 0.0;
    for (int i = 0; i < profile2->n_points - 1; i++) {
        double r1 = profile2->r[i];
        double r2 = profile2->r[i+1];
        double rho1 = profile2->rho[i];
        double rho2 = profile2->rho[i+1];
        double dr = r2 - r1;
        double r_mid = 0.5 * (r1 + r2);
        double rho_mid = 0.5 * (rho1 + rho2);
        double dM = 4.0 * M_PI * rho_mid * r_mid * r_mid * dr;
        M_star2 += dM;
    }
    
    printf("\n");
    printf("Stellar masses:\n");
    printf("  M1 ≈ %.4e\n", M_star1);
    printf("  M2 ≈ %.4e\n", M_star2);
    printf("  M_total = %.4e\n", M_star1 + M_star2);
    
    // If auto_orbit is enabled, calculate positions and velocities for circular orbit
    if (auto_orbit) {
        printf("\nAuto-calculating orbital parameters for circular orbit...\n");
        
        double M_total = M_star1 + M_star2;
        
        // Center of mass frame: place stars symmetrically about origin
        // Distance from center of mass for each star
        double r1_com = separation * M_star2 / M_total;
        double r2_com = separation * M_star1 / M_total;
        
        // Place stars along x-axis
        x1 = -r1_com;
        y1 = 0.0;
        z1 = 0.0;
        
        x2 = r2_com;
        y2 = 0.0;
        z2 = 0.0;
        
        // Circular orbital velocity: v = sqrt(G * M_total / separation)
        double v_orbital = sqrt(solver->G_const * M_total / separation);
        
        // Each star's velocity in center-of-mass frame
        double v1 = v_orbital * M_star2 / M_total;
        double v2 = v_orbital * M_star1 / M_total;
        
        // Velocities perpendicular to separation (in y-direction)
        vx1 = 0.0;
        vy1 = +v1;
        vz1 = 0.0;
        
        vx2 = 0.0;
        vy2 = -v2;
        vz2 = 0.0;
        
        // Orbital period
        double period = 2.0 * M_PI * separation / v_orbital;
        
        printf("  Calculated separation: a = %.4f\n", separation);
        printf("  Star 1 distance from COM: %.4f\n", r1_com);
        printf("  Star 2 distance from COM: %.4f\n", r2_com);
        printf("  Orbital velocity: v_orb = %.4e\n", v_orbital);
        printf("  Star 1 velocity: v1 = %.4e\n", v1);
        printf("  Star 2 velocity: v2 = %.4e\n", v2);
        printf("  Orbital period: P = %.4e\n", period);
        printf("  Suggested t_end: %.2f - %.2f (1-2 orbits)\n", period, 2.0 * period);
    }
    
    printf("\n");
    printf("Final binary configuration:\n");
    printf("  Star 1 position: (%.3f, %.3f, %.3f)\n", x1, y1, z1);
    printf("  Star 1 velocity: (%.3e, %.3e, %.3e)\n", vx1, vy1, vz1);
    printf("  Star 2 position: (%.3f, %.3f, %.3f)\n", x2, y2, z2);
    printf("  Star 2 velocity: (%.3e, %.3e, %.3e)\n", vx2, vy2, vz2);
    printf("\n");
    
    // Map both stars onto 3D grid
    printf("Mapping binary system to 3D grid...\n");
   
    for (int i = 0; i < solver->nx; i++) {
        for (int j = 0; j < solver->ny; j++) {
            for (int k = 0; k < solver->nz; k++) {
                double x = solver->xmin + (i + 0.5) * solver->dx;
                double y = solver->ymin + (j + 0.5) * solver->dy;
                double z = solver->zmin + (k + 0.5) * solver->dz;
                
                // Distance from each star center
                double r1 = sqrt(pow(x - x1, 2) + pow(y - y1, 2) + pow(z - z1, 2));
                double r2 = sqrt(pow(x - x2, 2) + pow(y - y2, 2) + pow(z - z2, 2));
                
                // Interpolate density and internal energy from profiles
                double rho1, u1, rho2, u2;
                interpolate_profile(profile1, r1, &rho1, &u1);
                interpolate_profile(profile2, r2, &rho2, &u2);
                
                // Assign velocity based on which star we're inside (by radius)
                Primitive3D P;
                
                // Check if we're inside either star's radius
                bool inside_star1 = (r1 < profile1->R_star);
                bool inside_star2 = (r2 < profile2->R_star);
                
                if (inside_star1 && !inside_star2) {
                    // Inside star 1 only
                    P.rho = rho1;
                    P.u = vx1;
                    P.v = vy1;
                    P.w = vz1;
                    P.p = (solver->gamma - 1.0) * rho1 * u1;
                    
                } else if (inside_star2 && !inside_star1) {
                    // Inside star 2 only
                    P.rho = rho2;
                    P.u = vx2;
                    P.v = vy2;
                    P.w = vz2;
                    P.p = (solver->gamma - 1.0) * rho2 * u2;
                    
                } else if (inside_star1 && inside_star2) {
                    // Overlap region - blend both stars
                    // Use density-weighted average
                     double rho_total = rho1 + rho2;
                    double weight1 = rho1 / rho_total;
                    double weight2 = rho2 / rho_total;
                    
                    P.rho = rho_total;
                    P.u = weight1 * vx1 + weight2 * vx2;
                    P.v = weight1 * vy1 + weight2 * vy2;
                    P.w = weight1 * vz1 + weight2 * vz2;
                    
                    double p1 = (solver->gamma - 1.0) * rho1 * u1;
                    double p2 = (solver->gamma - 1.0) * rho2 * u2;
                    P.p = p1 + p2;
                    
                } else {
                    // Outside both stars - ambient medium
                    P.rho = RHO_FLOOR;
                    P.u = 0.0;  // ← Zero velocity in ambient
                    P.v = 0.0;  // ← Zero velocity in ambient
                    P.w = 0.0;  // ← Zero velocity in ambient
                    P.p = P_FLOOR;
                }
                
                // Apply floors
                P.rho = MAX(P.rho, RHO_FLOOR);
                P.p = MAX(P.p, P_FLOOR);
                
                int idx = IDX3D(i, j, k, solver->ny, solver->nz);
                solver->U[idx] = prim_to_cons_3d(P, solver);
            }
        }
    }
    
    printf("Binary star system mapped to grid successfully.\n");
    
    // Clean up
    destroy_stellar_profile(profile1);
    destroy_stellar_profile(profile2);
}

// =============================================================================
// LEGACY BINARY STAR SYSTEM SETUP (same profile for both stars)
// =============================================================================

static void set_initial_conditions_binary_stars(GodunovHLLC3D *solver,
                                                 const char *profile_filename) {
    printf("Setting up binary star system from file...\n");
    
    // Read 1D profile
    StellarProfile1D *profile = read_stellar_profile(profile_filename);
    if (!profile) {
        fprintf(stderr, "Failed to read profile. Cannot create binary system.\n");
        return;
    }
    
    // Calculate stellar mass by integrating density profile
    // M ≈ 4π ∫ ρ(r) r² dr
    double M_star = 0.0;
    for (int i = 0; i < profile->n_points - 1; i++) {
        double r1 = profile->r[i];
        double r2 = profile->r[i+1];
        double rho1 = profile->rho[i];
        double rho2 = profile->rho[i+1];
        
        // Trapezoidal rule for integration
        double dr = r2 - r1;
        double r_mid = 0.5 * (r1 + r2);
        double rho_mid = 0.5 * (rho1 + rho2);
        double dM = 4.0 * M_PI * rho_mid * r_mid * r_mid * dr;
        M_star += dM;
    }
    
    printf("  Stellar mass: M ≈ %.4e\n", M_star);
    
    // Binary orbital parameters
    // For simplicity: circular orbit, equal masses, orbit in x-y plane
    double separation = 3.0 * profile->R_star;  // Initial separation (3 stellar radii)
    double M_total = 2.0 * M_star;
    
    // For circular orbit: v = sqrt(G*M_total / separation)
    // But in center-of-mass frame with equal masses, each star has v/2
    double v_orbital = sqrt(solver->G_const * M_total / separation);
    
    // For equal masses, each moves with v_orbital/2 perpendicular to separation
    double v1 = v_orbital * 0.5;
    double v2 = v_orbital * 0.5;
    
    // Star positions: along x-axis, symmetric about origin
    double x1 = -separation / 2.0;
    double y1 = 0.0;
    double z1 = 0.0;
    
    double x2 = separation / 2.0;
    double y2 = 0.0;
    double z2 = 0.0;
    
    // Velocities: perpendicular to separation (in y-direction), opposite directions
    double vx1 = 0.0;
    double vy1 = +v1;  // Star 1 moves in +y
    double vz1 = 0.0;
    
    double vx2 = 0.0;
    double vy2 = -v2;  // Star 2 moves in -y
    double vz2 = 0.0;
    
    // Orbital period
    double period = 2.0 * M_PI * separation / v_orbital;
    
    printf("\n");
    printf("Binary system parameters:\n");
    printf("  Stellar radius: R = %.4f\n", profile->R_star);
    printf("  Stellar mass: M ≈ %.4e\n", M_star);
    printf("  Separation: a = %.4f (%.2f R)\n", separation, separation / profile->R_star);
    printf("  Orbital velocity: v = %.4e\n", v_orbital);
    printf("  Orbital period: P = %.4e\n", period);
    printf("  Star 1 position: (%.3f, %.3f, %.3f)\n", x1, y1, z1);
    printf("  Star 1 velocity: (%.3f, %.3f, %.3f)\n", vx1, vy1, vz1);
    printf("  Star 2 position: (%.3f, %.3f, %.3f)\n", x2, y2, z2);
    printf("  Star 2 velocity: (%.3f, %.3f, %.3f)\n", vx2, vy2, vz2);
    printf("\n");
    
    // Map both stars onto 3D grid
printf("Mapping binary system to 3D grid...\n");
    
    for (int i = 0; i < solver->nx; i++) {
        for (int j = 0; j < solver->ny; j++) {
            for (int k = 0; k < solver->nz; k++) {
                double x = solver->xmin + (i + 0.5) * solver->dx;
                double y = solver->ymin + (j + 0.5) * solver->dy;
                double z = solver->zmin + (k + 0.5) * solver->dz;
                
                // Distance from each star center
                double r1 = sqrt(pow(x - x1, 2) + pow(y - y1, 2) + pow(z - z1, 2));
                double r2 = sqrt(pow(x - x2, 2) + pow(y - y2, 2) + pow(z - z2, 2));
                
                // Interpolate density and internal energy from profile
                double rho1, u1, rho2, u2;
                interpolate_profile(profile, r1, &rho1, &u1);
                interpolate_profile(profile, r2, &rho2, &u2);
                
                // Assign velocity based on which star we're inside (by radius)
                Primitive3D P;
                
                // Check if we're inside either star's radius
                bool inside_star1 = (r1 < profile->R_star);
                bool inside_star2 = (r2 < profile->R_star);
                
                if (inside_star1 && !inside_star2) {
                    // Inside star 1 only
                    P.rho = rho1;
                    P.u = vx1;
                    P.v = vy1;
                    P.w = vz1;
                    P.p = (solver->gamma - 1.0) * rho1 * u1;
                    
                } else if (inside_star2 && !inside_star1) {
                    // Inside star 2 only
                    P.rho = rho2;
                    P.u = vx2;
                    P.v = vy2;
                    P.w = vz2;
                    P.p = (solver->gamma - 1.0) * rho2 * u2;
                    
                } else if (inside_star1 && inside_star2) {
                    // Overlap region - blend both stars
                    // Use density-weighted average
                    double rho_total = rho1 + rho2;
                    double weight1 = rho1 / rho_total;
                    double weight2 = rho2 / rho_total;
                    
                    P.rho = rho_total;
                    P.u = weight1 * vx1 + weight2 * vx2;
                    P.v = weight1 * vy1 + weight2 * vy2;
                    P.w = weight1 * vz1 + weight2 * vz2;
                    
                    double p1 = (solver->gamma - 1.0) * rho1 * u1;
                    double p2 = (solver->gamma - 1.0) * rho2 * u2;
                    P.p = p1 + p2;
                    
                } else {
                    // Outside both stars - ambient medium
                    P.rho = RHO_FLOOR;
                    P.u = 0.0;  // ← Zero velocity in ambient
                    P.v = 0.0;  // ← Zero velocity in ambient
                    P.w = 0.0;  // ← Zero velocity in ambient
                    P.p = P_FLOOR;
                }
                
                // Apply floors
                P.rho = MAX(P.rho, RHO_FLOOR);
                P.p = MAX(P.p, P_FLOOR);
                
                int idx = IDX3D(i, j, k, solver->ny, solver->nz);
                solver->U[idx] = prim_to_cons_3d(P, solver);
            }
        }
    }
    
    printf("Binary star system mapped to grid successfully.\n");
    printf("Expected orbital period: P ≈ %.3f\n", period);
    printf("Suggested simulation time: t_end = %.1f (1-2 orbits)\n", 2.0 * period);
    
    // Clean up
    destroy_stellar_profile(profile);
}
// =============================================================================
// INITIAL CONDITIONS
// =============================================================================

void set_initial_conditions_3d(GodunovHLLC3D *solver, ProblemType problem) {
    if (problem == KELVIN_HELMHOLTZ_3D) {
        printf("Setting up 3D Kelvin-Helmholtz problem...\n");
        
        // Classic KH setup with two layers of different density/velocity
        const double rho1 = 2.0;  // Upper layer density
        const double rho2 = 1.0;  // Lower layer density
        const double v1 = -0.5;   // Upper layer velocity (in x-direction)
        const double v2 = 0.5;    // Lower layer velocity (in x-direction)
        const double p0 = 2.5;    // Pressure
        const double w0 = 0.1;    // Shear layer width
        const double A = 0.1;     // Perturbation amplitude (increased from 0.01)
        
        for (int i = 0; i < solver->nx; i++) {
            for (int j = 0; j < solver->ny; j++) {
                for (int k = 0; k < solver->nz; k++) {
                    double x = solver->xmin + (i + 0.5) * solver->dx;
                    double y = solver->ymin + (j + 0.5) * solver->dy;
                    double z = solver->zmin + (k + 0.5) * solver->dz;
                    
                    Primitive3D P;
                    
                    // Use y-coordinate for the shear layer at y=0
                    double y_interface = 0.0;
                    
                    // Smooth transition using tanh
                    double f = 0.5 * (1.0 + tanh((y - y_interface) / w0));
                    
                    // Set base state
                    P.rho = rho1 * f + rho2 * (1.0 - f);
                    P.u = v1 * f + v2 * (1.0 - f);  // Shear in x-direction
                    
                    // Add velocity perturbations localized to the shear layer
                    // Perturbations should be perpendicular to the flow (in y and z)
                    // and concentrated near the interface
                    
                    // Localization factor: concentrate perturbations near y=0
                    double envelope = exp(-y * y / (4.0 * w0 * w0));
                    
                    // Single-mode perturbation for clear instability
                    // Use wavelength comparable to domain size
                    double lambda_x = 2.0;  // Wavelength in x
                    double kx = 2.0 * M_PI / lambda_x;
                    
                    // Add random phase in z for 3D character
                    double lambda_z = 2.0;
                    double kz = 2.0 * M_PI / lambda_z;
                    
                    // Velocity perturbation in y-direction (perpendicular to shear)
                    double pert_v = A * sin(kx * x) * cos(kz * z) * envelope;
                    
                    // Small z-perturbation for 3D structure
                    double pert_w = 0.5 * A * cos(kx * x) * sin(kz * z) * envelope;
                    
                    P.v = pert_v;
                    P.w = pert_w;
                    P.p = p0;
                    
                    int idx = IDX3D(i, j, k, solver->ny, solver->nz);
                    solver->U[idx] = prim_to_cons_3d(P, solver);
                }
            }
        }
        
        printf("3D Kelvin-Helmholtz initial conditions set.\n");
        
    } else if (problem == BLAST_WAVE_3D) {
        printf("Setting up 3D Blast Wave problem...\n");
        
        // Sedov-Taylor blast wave
        const double rho_ambient = 1.0;
        const double p_ambient = 1.0e-5;
        const double E_blast = 1.0;  // Total energy in blast
        const double r_blast = 0.1;  // Blast radius
        
        double center_x = 0.5 * (solver->xmin + solver->xmax);
        double center_y = 0.5 * (solver->ymin + solver->ymax);
        double center_z = 0.5 * (solver->zmin + solver->zmax);
        
        // Calculate volume of blast region
        double V_blast = (4.0 / 3.0) * M_PI * r_blast * r_blast * r_blast;
        double p_blast = (solver->gamma - 1.0) * E_blast / V_blast;
        
        for (int i = 0; i < solver->nx; i++) {
            for (int j = 0; j < solver->ny; j++) {
                for (int k = 0; k < solver->nz; k++) {
                    double x = solver->xmin + (i + 0.5) * solver->dx;
                    double y = solver->ymin + (j + 0.5) * solver->dy;
                    double z = solver->zmin + (k + 0.5) * solver->dz;
                    
                    double r = sqrt((x - center_x) * (x - center_x) +
                                    (y - center_y) * (y - center_y) +
                                    (z - center_z) * (z - center_z));
                    
                    Primitive3D P;
                    P.rho = rho_ambient;
                    P.u = 0.0;
                    P.v = 0.0;
                    P.w = 0.0;
                    
                    if (r < r_blast) {
                        P.p = p_blast;
                    } else {
                        P.p = p_ambient;
                    }
                    
                    int idx = IDX3D(i, j, k, solver->ny, solver->nz);
                    solver->U[idx] = prim_to_cons_3d(P, solver);
                }
            }
        }
        
        printf("3D Blast Wave initial conditions set.\n");
    } else if (problem == STELLAR_PROFILE_3D) {
        printf("Setting up 3D stellar profile from file...\n");
        
        // Read profile from file
        // File should be named "stellar_profile.txt" in run directory
        set_initial_conditions_stellar_profile(solver, "stellar_profile.txt");
        
        printf("Stellar profile initial conditions set.\n");
    } else if (problem == BINARY_STAR_3D) {
        printf("Setting up 3D binary star system from file...\n");
        printf("Note: Using legacy mode (same profile for both stars)\n");
        printf("For different profiles, use set_initial_conditions_3d_from_deck()\n");
        
        // Read profile from file and create binary system
        // File should be named "stellar_profile.txt" in run directory
        set_initial_conditions_binary_stars(solver, "stellar_profile.txt");
        
        printf("Binary star initial conditions set.\n");
    }
}

// =============================================================================
// INITIAL CONDITIONS WITH INPUT DECK (NEW)
// =============================================================================

// Helper function: Simple pseudo-random number generator (Linear Congruential Generator)
static unsigned long rng_state = 1;

static void seed_rng(int seed) {
    rng_state = (unsigned long)seed;
}

static double rand_uniform(void) {
    // LCG parameters from Numerical Recipes
    rng_state = (rng_state * 1664525UL + 1013904223UL) & 0xFFFFFFFFUL;
    return (double)rng_state / 4294967296.0;  // Return value in [0, 1)
}

static double rand_normal(void) {
    // Box-Muller transform to generate normal distribution
    double u1 = rand_uniform();
    double u2 = rand_uniform();
    return sqrt(-2.0 * log(u1 + 1e-10)) * cos(2.0 * M_PI * u2);
}

// Helper function: Associated Legendre polynomial P_l^m(x)
// Simplified implementation for low l values (l <= 4)
static double legendre_Plm(int l, int m, double x) {
    if (m < 0 || m > l) return 0.0;
    
    // Ensure x is in valid range
    if (x > 1.0) x = 1.0;
    if (x < -1.0) x = -1.0;
    
    double pmm = 1.0;
    if (m > 0) {
        double somx2 = sqrt((1.0 - x) * (1.0 + x));
        double fact = 1.0;
        for (int i = 1; i <= m; i++) {
            pmm *= -fact * somx2;
            fact += 2.0;
        }
    }
    
    if (l == m) return pmm;
    
    double pmmp1 = x * (2 * m + 1) * pmm;
    if (l == m + 1) return pmmp1;
    
    double pll = 0.0;
    for (int ll = m + 2; ll <= l; ll++) {
        pll = ((2 * ll - 1) * x * pmmp1 - (ll + m - 1) * pmm) / (ll - m);
        pmm = pmmp1;
        pmmp1 = pll;
    }
    
    return pll;
}

// Helper function: Real spherical harmonic Y_l^m(theta, phi)
static double spherical_harmonic_real(int l, int m, double theta, double phi) {
    double cos_theta = cos(theta);
    double Plm = legendre_Plm(l, abs(m), cos_theta);
    
    // Normalization factor
    double norm = sqrt((2.0 * l + 1.0) / (4.0 * M_PI));
    
    // Include factorial ratio in normalization
    // For simplicity, we'll use approximate normalization
    
    if (m > 0) {
        return norm * Plm * cos(m * phi);
    } else if (m < 0) {
        return norm * Plm * sin(abs(m) * phi);
    } else {
        return norm * Plm;
    }
}

// Apply random Fourier mode perturbations to velocity field
static void apply_velocity_perturbations(GodunovHLLC3D *solver, 
                                         double amplitude,
                                         double kmin, double kmax,
                                         int n_modes) {
    printf("  Adding random velocity perturbations:\n");
    printf("    Amplitude: %.3e (fraction of sound speed)\n", amplitude);
    printf("    Wavenumber range: [%.2f, %.2f]\n", kmin, kmax);
    printf("    Number of modes: %d\n", n_modes);
    
    // Generate random Fourier modes
    typedef struct {
        double kx, ky, kz;  // Wave vector
        double ax, ay, az;  // Amplitude for each velocity component
        double phase;       // Random phase
    } FourierMode;
    
    FourierMode *modes = (FourierMode*) malloc(n_modes * sizeof(FourierMode));
    
    for (int n = 0; n < n_modes; n++) {
        // Random wavenumber magnitude
        double k_mag = kmin + (kmax - kmin) * rand_uniform();
        
        // Random direction (uniform on sphere)
        double phi = 2.0 * M_PI * rand_uniform();
        double cos_theta = 2.0 * rand_uniform() - 1.0;
        double sin_theta = sqrt(1.0 - cos_theta * cos_theta);
        
        modes[n].kx = k_mag * sin_theta * cos(phi);
        modes[n].ky = k_mag * sin_theta * sin(phi);
        modes[n].kz = k_mag * cos_theta;
        
        // Random amplitudes (each component independent)
        modes[n].ax = rand_normal();
        modes[n].ay = rand_normal();
        modes[n].az = rand_normal();
        
        // Random phase
        modes[n].phase = 2.0 * M_PI * rand_uniform();
    }
    
    // Apply perturbations to grid
    for (int i = 0; i < solver->nx; i++) {
        for (int j = 0; j < solver->ny; j++) {
            for (int k = 0; k < solver->nz; k++) {
                int idx = IDX3D(i, j, k, solver->ny, solver->nz);
                
                double x = solver->xmin + (i + 0.5) * solver->dx;
                double y = solver->ymin + (j + 0.5) * solver->dy;
                double z = solver->zmin + (k + 0.5) * solver->dz;
                
                // Skip cells outside octant
                if (x < 0.0 || y < 0.0 || z < 0.0) continue;
                
                // Get current state
                Primitive3D P = cons_to_prim_3d(solver->U[idx], solver);
                
                // Calculate local sound speed for scaling
                double c_s = sound_speed_3d(P.rho, P.p, solver);
                
                // Sum contributions from all modes
                double du = 0.0, dv = 0.0, dw = 0.0;
                for (int n = 0; n < n_modes; n++) {
                    double arg = modes[n].kx * x + modes[n].ky * y + modes[n].kz * z + modes[n].phase;
                    double wave = cos(arg);
                    
                    du += modes[n].ax * wave;
                    dv += modes[n].ay * wave;
                    dw += modes[n].az * wave;
                }
                
                // Normalize and scale by amplitude and sound speed
                double norm = sqrt((double)n_modes);
                P.u += (du / norm) * amplitude * c_s;
                P.v += (dv / norm) * amplitude * c_s;
                P.w += (dw / norm) * amplitude * c_s;
                
                // Update conservative variables
                solver->U[idx] = prim_to_cons_3d(P, solver);
            }
        }
    }
    
    free(modes);
    printf("  Velocity perturbations applied.\n");
}

// Apply random Fourier mode perturbations to density field
static void apply_density_perturbations(GodunovHLLC3D *solver,
                                        double amplitude,
                                        double kmin, double kmax,
                                        int n_modes) {
    printf("  Adding random density perturbations:\n");
    printf("    Fractional amplitude: %.3e\n", amplitude);
    printf("    Wavenumber range: [%.2f, %.2f]\n", kmin, kmax);
    printf("    Number of modes: %d\n", n_modes);
    
    // Generate random Fourier modes
    typedef struct {
        double kx, ky, kz;
        double amp;
        double phase;
    } DensityMode;
    
    DensityMode *modes = (DensityMode*) malloc(n_modes * sizeof(DensityMode));
    
    for (int n = 0; n < n_modes; n++) {
        double k_mag = kmin + (kmax - kmin) * rand_uniform();
        
        double phi = 2.0 * M_PI * rand_uniform();
        double cos_theta = 2.0 * rand_uniform() - 1.0;
        double sin_theta = sqrt(1.0 - cos_theta * cos_theta);
        
        modes[n].kx = k_mag * sin_theta * cos(phi);
        modes[n].ky = k_mag * sin_theta * sin(phi);
        modes[n].kz = k_mag * cos_theta;
        modes[n].amp = rand_normal();
        modes[n].phase = 2.0 * M_PI * rand_uniform();
    }
    
    // Apply perturbations
    for (int i = 0; i < solver->nx; i++) {
        for (int j = 0; j < solver->ny; j++) {
            for (int k = 0; k < solver->nz; k++) {
                int idx = IDX3D(i, j, k, solver->ny, solver->nz);
                
                double x = solver->xmin + (i + 0.5) * solver->dx;
                double y = solver->ymin + (j + 0.5) * solver->dy;
                double z = solver->zmin + (k + 0.5) * solver->dz;
                
                if (x < 0.0 || y < 0.0 || z < 0.0) continue;
                
                Primitive3D P = cons_to_prim_3d(solver->U[idx], solver);
                
                // Sum mode contributions
                double delta_rho = 0.0;
                for (int n = 0; n < n_modes; n++) {
                    double arg = modes[n].kx * x + modes[n].ky * y + modes[n].kz * z + modes[n].phase;
                    delta_rho += modes[n].amp * cos(arg);
                }
                
                // Normalize and apply
                double norm = sqrt((double)n_modes);
                double rho_pert = 1.0 + (delta_rho / norm) * amplitude;
                
                // Ensure density doesn't go negative
                if (rho_pert < 0.1) rho_pert = 0.1;
                
                P.rho *= rho_pert;
                
                // Update conservative variables
                solver->U[idx] = prim_to_cons_3d(P, solver);
            }
        }
    }
    
    free(modes);
    printf("  Density perturbations applied.\n");
}

// Apply spherical harmonic mode perturbations
static void apply_spherical_harmonic_perturbation(GodunovHLLC3D *solver,
                                                  int l, int m,
                                                  double amp_rho,
                                                  double amp_vel) {
    printf("  Adding spherical harmonic mode perturbation:\n");
    printf("    Mode: Y_%d^%d\n", l, m);
    printf("    Density amplitude: %.3e (fractional)\n", amp_rho);
    printf("    Velocity amplitude: %.3e (fraction of c_s)\n", amp_vel);
    
    // Physical interpretation of modes:
    if (l == 1 && m == 0) printf("    (Dipole mode - global shift/rotation axis asymmetry)\n");
    else if (l == 1) printf("    (Dipole mode - lateral asymmetry)\n");
    else if (l == 2 && m == 0) printf("    (Quadrupole mode - prolate/oblate deformation)\n");
    else if (l == 2) printf("    (Quadrupole mode - elliptical asymmetry)\n");
    else printf("    (Higher-order mode - complex asymmetry pattern)\n");
    
    for (int i = 0; i < solver->nx; i++) {
        for (int j = 0; j < solver->ny; j++) {
            for (int k = 0; k < solver->nz; k++) {
                int idx = IDX3D(i, j, k, solver->ny, solver->nz);
                
                double x = solver->xmin + (i + 0.5) * solver->dx;
                double y = solver->ymin + (j + 0.5) * solver->dy;
                double z = solver->zmin + (k + 0.5) * solver->dz;
                
                if (x < 0.0 || y < 0.0 || z < 0.0) continue;
                
                double r = sqrt(x*x + y*y + z*z);
                if (r < 1e-10) continue;
                
                // Convert to spherical coordinates
                double theta = acos(z / r);
                double phi = atan2(y, x);
                
                // Compute spherical harmonic
                double Y_lm = spherical_harmonic_real(l, m, theta, phi);
                
                // Get current state
                Primitive3D P = cons_to_prim_3d(solver->U[idx], solver);
                double c_s = sound_speed_3d(P.rho, P.p, solver);
                
                // Apply density perturbation
                P.rho *= (1.0 + amp_rho * Y_lm);
                
                // Apply velocity perturbation (radial and tangential components)
                // Radial component
                double v_r = amp_vel * c_s * Y_lm;
                P.u += v_r * (x / r);
                P.v += v_r * (y / r);
                P.w += v_r * (z / r);
                
                // Update conservative variables
                solver->U[idx] = prim_to_cons_3d(P, solver);
            }
        }
    }
    
    printf("  Spherical harmonic perturbations applied.\n");
}

// Main perturbation application function
static void apply_perturbations_to_octant(GodunovHLLC3D *solver, const InputDeck *deck) {
    printf("Applying perturbations with seed = %d\n", deck->perturbation_seed);
    seed_rng(deck->perturbation_seed);
    
    // Apply random velocity perturbations
    if (deck->perturb_velocity_amplitude > 0.0 && deck->perturb_velocity_modes > 0) {
        printf("\n--- Random Velocity Perturbations ---\n");
        printf("Physical motivation: Seeds convective instabilities and turbulence.\n");
        printf("Astrophysical context: Pre-SN convection, neutrino-driven convection.\n\n");
        
        apply_velocity_perturbations(solver,
                                     deck->perturb_velocity_amplitude,
                                     deck->perturb_velocity_kmin,
                                     deck->perturb_velocity_kmax,
                                     deck->perturb_velocity_modes);
    }
    
    // Apply random density perturbations
    if (deck->perturb_density_amplitude > 0.0 && deck->perturb_density_modes > 0) {
        printf("\n--- Random Density Perturbations ---\n");
        printf("Physical motivation: Composition variations and shell structure.\n");
        printf("Astrophysical context: Burning shell interfaces, composition gradients.\n\n");
        
        apply_density_perturbations(solver,
                                    deck->perturb_density_amplitude,
                                    deck->perturb_density_kmin,
                                    deck->perturb_density_kmax,
                                    deck->perturb_density_modes);
    }
    
    // Apply spherical harmonic perturbation
    if (deck->enable_sph_harmonic_mode) {
        printf("\n--- Spherical Harmonic Mode Perturbation ---\n");
        printf("Physical motivation: Large-scale asymmetries from rotation, magnetic fields.\n");
        printf("Astrophysical context: SASI modes, rotation-induced deformation, binary effects.\n\n");
        
        apply_spherical_harmonic_perturbation(solver,
                                              deck->sph_harmonic_l,
                                              deck->sph_harmonic_m,
                                              deck->sph_harmonic_amp_rho,
                                              deck->sph_harmonic_amp_vel);
    }
    
    printf("\nAll perturbations applied successfully.\n");
}

// Stellar octant explosion setup
static void set_initial_conditions_stellar_octant_explosion(GodunovHLLC3D *solver,
                                                             const char *profile_filename,
                                                             double explosion_energy,
                                                             double explosion_radius,
                                                             const char *explosion_type) {
    printf("Setting up stellar octant with explosion...\n");
    printf("  Profile: %s\n", profile_filename);
    printf("  Explosion energy: %.3e (code units)\n", explosion_energy);
    printf("  Explosion radius: %.3f\n", explosion_radius);
    printf("  Explosion type: %s\n", explosion_type);
    
    // Read stellar profile
    StellarProfile1D *profile = read_stellar_profile(profile_filename);
    if (!profile) {
        fprintf(stderr, "Failed to read stellar profile. Cannot create stellar octant.\n");
        return;
    }
    
    printf("  Stellar radius: R = %.4f\n", profile->R_star);
    
    // Calculate the volume within explosion_radius for energy distribution
    double explosion_volume = 0.0;
    int explosion_cell_count = 0;
    
    // First pass: count cells in explosion region
    for (int i = 0; i < solver->nx; i++) {
        for (int j = 0; j < solver->ny; j++) {
            for (int k = 0; k < solver->nz; k++) {
                double x = solver->xmin + (i + 0.5) * solver->dx;
                double y = solver->ymin + (j + 0.5) * solver->dy;
                double z = solver->zmin + (k + 0.5) * solver->dz;
                
                // Only fill octant where x >= 0, y >= 0, z >= 0
                if (x < 0.0 || y < 0.0 || z < 0.0) continue;
                
                double r = sqrt(x*x + y*y + z*z);
                
                if (r < explosion_radius) {
                    explosion_cell_count++;
                    explosion_volume += solver->dx * solver->dy * solver->dz;
                }
            }
        }
    }
    
    printf("  Explosion region: %d cells, volume = %.4e\n", explosion_cell_count, explosion_volume);
    
    // Calculate energy per cell (for thermal explosion) or velocity magnitude (for kinetic)
    double energy_per_cell = 0.0;
    double velocity_magnitude = 0.0;
    bool is_thermal = (strcmp(explosion_type, "thermal") == 0);
    
    if (is_thermal) {
        energy_per_cell = explosion_energy / (double)explosion_cell_count;
        printf("  Energy per cell: %.3e\n", energy_per_cell);
    } else {
        // For kinetic explosion, estimate velocity from kinetic energy
        // KE = 0.5 * M * v^2, so v = sqrt(2 * KE / M)
        // We'll estimate mass in explosion region
        double total_mass_in_explosion = 0.0;
        
        for (int i = 0; i < solver->nx; i++) {
            for (int j = 0; j < solver->ny; j++) {
                for (int k = 0; k < solver->nz; k++) {
                    double x = solver->xmin + (i + 0.5) * solver->dx;
                    double y = solver->ymin + (j + 0.5) * solver->dy;
                    double z = solver->zmin + (k + 0.5) * solver->dz;
                    
                    if (x < 0.0 || y < 0.0 || z < 0.0) continue;
                    
                    double r = sqrt(x*x + y*y + z*z);
                    
                    if (r < explosion_radius) {
                        double rho, u;
                        interpolate_profile(profile, r, &rho, &u);
                        double cell_mass = rho * solver->dx * solver->dy * solver->dz;
                        total_mass_in_explosion += cell_mass;
                    }
                }
            }
        }
        
        velocity_magnitude = sqrt(2.0 * explosion_energy / total_mass_in_explosion);
        printf("  Total mass in explosion region: %.3e\n", total_mass_in_explosion);
        printf("  Radial velocity magnitude: %.3e\n", velocity_magnitude);
    }
    
    // Set up the grid
    int cells_set = 0;
    for (int i = 0; i < solver->nx; i++) {
        for (int j = 0; j < solver->ny; j++) {
            for (int k = 0; k < solver->nz; k++) {
                int idx = IDX3D(i, j, k, solver->ny, solver->nz);
                
                double x = solver->xmin + (i + 0.5) * solver->dx;
                double y = solver->ymin + (j + 0.5) * solver->dy;
                double z = solver->zmin + (k + 0.5) * solver->dz;
                
                Primitive3D P;
                
                // Only fill octant where x >= 0, y >= 0, z >= 0
                if (x >= 0.0 && y >= 0.0 && z >= 0.0) {
                    double r = sqrt(x*x + y*y + z*z);
                    
                    // Interpolate stellar profile
                    double rho, u;
                    interpolate_profile(profile, r, &rho, &u);
                    
                    P.rho = rho;
                    P.u = 0.0;
                    P.v = 0.0;
                    P.w = 0.0;
                    
                    // Calculate pressure from internal energy
                    // p = rho * u * (gamma - 1)
                    P.p = rho * u * (solver->gamma - 1.0);
                    
                    // Add explosion
                    if (r < explosion_radius) {
                        if (is_thermal) {
                            // Thermal explosion: add energy
                            // E_new = E_old + energy_per_cell
                            // For ideal gas: E = rho * u + 0.5 * rho * v^2
                            // We add to internal energy, so u_new = u + energy_per_cell / rho
                            double u_internal = u;
                            u_internal += energy_per_cell / (rho * solver->dx * solver->dy * solver->dz);
                            P.p = rho * u_internal * (solver->gamma - 1.0);
                        } else {
                            // Kinetic explosion: add radial velocity
                            if (r > 1e-10) {  // Avoid division by zero at center
                                P.u = velocity_magnitude * (x / r);
                                P.v = velocity_magnitude * (y / r);
                                P.w = velocity_magnitude * (z / r);
                            }
                        }
                    }
                    
                    cells_set++;
                } else {
                    // Outside octant: set to low-density atmosphere
                    P.rho = RHO_FLOOR;
                    P.u = 0.0;
                    P.v = 0.0;
                    P.w = 0.0;
                    P.p = P_FLOOR;
                }
                
                solver->U[idx] = prim_to_cons_3d(P, solver);
            }
        }
    }
    
    printf("  Set up %d cells in octant\n", cells_set);
    printf("Stellar octant explosion initial conditions complete.\n");
    
    destroy_stellar_profile(profile);
}

// =============================================================================
// ADD TORUS TO OCTANT SETUP
// =============================================================================

static void add_torus_to_octant(GodunovHLLC3D *solver, const InputDeck *deck) {
    printf("\nAdding torus to octant setup...\n");
    printf("  Major radius R: %.3f\n", deck->torus_major_radius);
    printf("  Minor radius r: %.3f\n", deck->torus_minor_radius);
    printf("  Density: %.3e\n", deck->torus_density);
    printf("  Pressure: %.3e\n", deck->torus_pressure);
    printf("  Rotation velocity: %.3f\n", deck->torus_rotation_velocity);
    
    int torus_cells = 0;
    double R = deck->torus_major_radius;  // Major radius (distance from origin to tube center)
    double r = deck->torus_minor_radius;  // Minor radius (tube thickness)
    
    // Loop over all cells
    for (int i = 0; i < solver->nx; i++) {
        for (int j = 0; j < solver->ny; j++) {
            for (int k = 0; k < solver->nz; k++) {
                int idx = IDX3D(i, j, k, solver->ny, solver->nz);
                
                double x = solver->xmin + (i + 0.5) * solver->dx;
                double y = solver->ymin + (j + 0.5) * solver->dy;
                double z = solver->zmin + (k + 0.5) * solver->dz;
                
                // Only fill octant where x >= 0, y >= 0, z >= 0
                if (x < 0.0 || y < 0.0 || z < 0.0) continue;
                
                // Calculate distance from z-axis in x-y plane
                double rho_cyl = sqrt(x*x + y*y);
                
                // Torus condition: (rho_cyl - R)^2 + z^2 <= r^2
                // This defines a torus in the x-y plane, symmetric about z-axis
                double torus_distance = sqrt((rho_cyl - R) * (rho_cyl - R) + z * z);
                
                if (torus_distance <= r) {
                    // This cell is inside the torus
                    Primitive3D P;
                    P.rho = deck->torus_density;
                    P.p = deck->torus_pressure;
                    
                    // Set rotation velocity (azimuthal direction in x-y plane)
                    if (rho_cyl > 1e-10) {  // Avoid division by zero at axis
                        // Azimuthal velocity: v_phi = v_rot * (tangent vector)
                        // Tangent vector in x-y plane: (-y, x, 0) / rho_cyl
                        P.u = -deck->torus_rotation_velocity * y / rho_cyl;
                        P.v = deck->torus_rotation_velocity * x / rho_cyl;
                        P.w = 0.0;
                    } else {
                        P.u = 0.0;
                        P.v = 0.0;
                        P.w = 0.0;
                    }
                    
                    // Convert to conservative variables and overwrite cell
                    solver->U[idx] = prim_to_cons_3d(P, solver);
                    torus_cells++;
                }
            }
        }
    }
    
    printf("  Torus added: %d cells filled\n", torus_cells);
}

// =============================================================================
// INITIAL CONDITIONS FROM INPUT DECK
// =============================================================================

void set_initial_conditions_3d_from_deck(GodunovHLLC3D *solver, const InputDeck *deck) {
    if (deck->problem == STELLAR_OCTANT_EXPLOSION_3D) {
        printf("Setting up stellar octant with explosion...\n");
        
        set_initial_conditions_stellar_octant_explosion(solver,
                                                        deck->stellar_profile,
                                                        deck->explosion_energy,
                                                        deck->explosion_radius,
                                                        deck->explosion_type);
        
        // Apply perturbations if enabled
        if (deck->enable_perturbations) {
            printf("\nApplying perturbations to seed turbulence...\n");
            apply_perturbations_to_octant(solver, deck);
        }
        
        // Add torus if enabled
        if (deck->enable_torus) {
            add_torus_to_octant(solver, deck);
        }
        
        printf("Stellar octant explosion initial conditions set.\n");
    } else if (deck->problem == BINARY_STAR_3D) {
        printf("Setting up 3D binary star system with dual profiles...\n");
        
        // Use the new dual-profile binary setup
        set_initial_conditions_binary_stars_dual(solver,
                                                 deck->star1_profile,
                                                 deck->star2_profile,
                                                 deck->star1_x, deck->star1_y, deck->star1_z,
                                                 deck->star2_x, deck->star2_y, deck->star2_z,
                                                 deck->star1_vx, deck->star1_vy, deck->star1_vz,
                                                 deck->star2_vx, deck->star2_vy, deck->star2_vz,
                                                 deck->binary_auto_orbit,
                                                 deck->binary_separation);
        
        printf("Binary star initial conditions set.\n");
    } else {
        // For all other problems, use the standard initialization
        set_initial_conditions_3d(solver, deck->problem);
    }
}

// =============================================================================
// RESTART FROM CSV FILE
// =============================================================================

int load_from_csv_restart(GodunovHLLC3D *solver, const char *filename) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Error: Could not open restart file '%s'\n", filename);
        return -1;
    }
    
    printf("Loading restart data from: %s\n", filename);
    
    // Read header line
    char line[1024];
    if (!fgets(line, sizeof(line), fp)) {
        fprintf(stderr, "Error: Could not read header from restart file\n");
        fclose(fp);
        return -1;
    }
    
    // Verify header format (should be: x,y,z,rho,u,v,w,p,gx,gy,gz)
    if (strstr(line, "rho") == NULL || strstr(line, "u") == NULL) {
        fprintf(stderr, "Error: Invalid header format in restart file\n");
        fprintf(stderr, "Expected: x,y,z,rho,u,v,w,p,gx,gy,gz\n");
        fclose(fp);
        return -1;
    }
    
    // Initialize all cells to zero first
    int ncells = solver->nx * solver->ny * solver->nz;
    memset(solver->U, 0, ncells * sizeof(State3D));
    
    // Read data lines
    int points_read = 0;
    double x, y, z, rho, u, v, w, p, gx, gy, gz;
    
    while (fgets(line, sizeof(line), fp)) {
        // Parse CSV line - try full format with gravity first
        int parsed = sscanf(line, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                           &x, &y, &z, &rho, &u, &v, &w, &p, &gx, &gy, &gz);
        
        // If that fails, try format without gravity fields
        if (parsed < 8) {
            parsed = sscanf(line, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                           &x, &y, &z, &rho, &u, &v, &w, &p);
            gx = gy = gz = 0.0;
        }
        
        if (parsed < 8) {
            continue; // Skip malformed lines
        }
        
        // Find the corresponding cell index
        int i = (int)((x - solver->xmin) / solver->dx);
        int j = (int)((y - solver->ymin) / solver->dy);
        int k = (int)((z - solver->zmin) / solver->dz);
        
        // Check if indices are in valid range
        if (i >= 0 && i < solver->nx &&
            j >= 0 && j < solver->ny &&
            k >= 0 && k < solver->nz) {
            
            int idx = IDX3D(i, j, k, solver->ny, solver->nz);
            
            // Convert primitive to conservative variables
            Primitive3D P;
            P.rho = rho;
            P.u = u;
            P.v = v;
            P.w = w;
            P.p = p;
            
            solver->U[idx] = prim_to_cons_3d(P, solver);
            
            // Also restore gravity field if available
            if (solver->enable_gravity && parsed >= 11) {
                solver->gx[idx] = gx;
                solver->gy[idx] = gy;
                solver->gz[idx] = gz;
            }
            
            points_read++;
        }
    }
    
    fclose(fp);
    
    printf("Successfully loaded %d data points from restart file\n", points_read);
    
    // Check if we read enough points
    if (points_read < ncells / 10) {
        fprintf(stderr, "Warning: Only loaded %d points out of %d cells. This may indicate:\n",
                points_read, ncells);
        fprintf(stderr, "  - The restart file was subsampled\n");
        fprintf(stderr, "  - Different grid resolution\n");
        fprintf(stderr, "  - Incomplete or corrupted file\n");
    }
    
    return 0;
}

// =============================================================================
// SIMULATION MAIN LOOP
// =============================================================================

void run_simulation_3d(GodunovHLLC3D *solver, double t_end, double output_dt) {
    double t = 0.0;
    int step_counter = 0;
    int output_count = 0;
    double next_output_time = 0.0;
    
    printf("Starting 3D simulation from t=0.0 to t=%.5f\n", t_end);
    printf("Grid: %d x %d x %d = %d cells\n", solver->nx, solver->ny, solver->nz, 
           solver->nx * solver->ny * solver->nz);
    
    // Write initial conditions
    write_output_to_file_3d(solver, output_count);
    write_output_vtk(solver, output_count);
    output_count++;
    next_output_time = output_dt;
    
    while (t < t_end) {
        godunov_step_3d(solver);
        t += solver->dt;
        step_counter++;
        
        // Check if it's time for output
        if (t >= next_output_time) {
            write_output_to_file_3d(solver, output_count);
            write_output_vtk(solver, output_count);
            printf("Step %d, t = %.5f, dt = %.5e, Output %d\n", 
                   step_counter, t, solver->dt, output_count);
            output_count++;
            next_output_time += output_dt;
        }
        
        if (step_counter % 10 == 0) {
            printf("Step %d, t = %.5f, dt = %.5e\n", step_counter, t, solver->dt);
        }
    }
    
    // Final output
    if (t - (next_output_time - output_dt) > 1e-6) {
        write_output_to_file_3d(solver, output_count);
        write_output_vtk(solver, output_count);
    }
    
    printf("Simulation completed: %d steps, final time: %.5f\n", step_counter, t);
}

void run_simulation_3d_restart(GodunovHLLC3D *solver, double t_start, double t_end, 
                                double output_dt, int start_output_count) {
    double t = t_start;
    int step_counter = 0;
    int output_count = start_output_count;
    double next_output_time = t_start + output_dt;
    
    printf("Restarting 3D simulation from t=%.5f to t=%.5f\n", t_start, t_end);
    printf("Starting output count: %d\n", start_output_count);
    printf("Grid: %d x %d x %d = %d cells\n", solver->nx, solver->ny, solver->nz, 
           solver->nx * solver->ny * solver->nz);
    
    // Write initial output at restart time
    write_output_to_file_3d(solver, output_count);
    write_output_vtk(solver, output_count);
    output_count++;
    
    while (t < t_end) {
        godunov_step_3d(solver);
        t += solver->dt;
        step_counter++;
        
        // Check if it's time for output
        if (t >= next_output_time) {
            write_output_to_file_3d(solver, output_count);
            write_output_vtk(solver, output_count);
            printf("Step %d, t = %.5f, dt = %.5e, Output %d\n", 
                   step_counter, t, solver->dt, output_count);
            output_count++;
            next_output_time += output_dt;
        }
        
        if (step_counter % 10 == 0) {
            printf("Step %d, t = %.5f, dt = %.5e\n", step_counter, t, solver->dt);
        }
    }
    
    // Final output
    if (t - (next_output_time - output_dt) > 1e-6) {
        write_output_to_file_3d(solver, output_count);
        write_output_vtk(solver, output_count);
    }
    
    printf("Simulation completed: %d steps, final time: %.5f\n", step_counter, t);
}

// =============================================================================
// OUTPUT FUNCTIONS
// =============================================================================

void write_output_to_file_3d(const GodunovHLLC3D *solver, int step) {
    char filename[256];
    sprintf(filename, "output/output_%06d.csv", step);
    
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Error: Could not open file %s for writing.\n", filename);
        return;
    }
    
    fprintf(fp, "x,y,z,rho,u,v,w,p,gx,gy,gz\n");  // Add gx,gy,gz
    
    // Write every Nth point to avoid huge files (subsample for large grids)
    int skip = 1;
    if (solver->nx > 128) skip = 2;
    if (solver->nx > 256) skip = 4;
    
    for (int i = 0; i < solver->nx; i += skip) {
        for (int j = 0; j < solver->ny; j += skip) {
            for (int k = 0; k < solver->nz; k += skip) {
                int idx = IDX3D(i, j, k, solver->ny, solver->nz);
                double x = solver->xmin + (i + 0.5) * solver->dx;
                double y = solver->ymin + (j + 0.5) * solver->dy;
                double z = solver->zmin + (k + 0.5) * solver->dz;
                
                Primitive3D P = cons_to_prim_3d(solver->U[idx], solver);
                double gx = solver->enable_gravity ? solver->gx[idx] : 0.0;
                double gy = solver->enable_gravity ? solver->gy[idx] : 0.0;
                double gz = solver->enable_gravity ? solver->gz[idx] : 0.0;
                fprintf(fp, "%.6f,%.6f,%.6f,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n", 
                    x, y, z, P.rho, P.u, P.v, P.w, P.p, gx, gy, gz);
                //fprintf(fp, "%.6f,%.6f,%.6f,%.6e,%.6e,%.6e,%.6e,%.6e\n", 
                //        x, y, z, P.rho, P.u, P.v, P.w, P.p);
            }
        }
    }
    
    fclose(fp);
}

void write_output_vtk(const GodunovHLLC3D *solver, int step) {
    char filename[256];
    sprintf(filename, "output/output_%06d.vtk", step);
    
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Error: Could not open file %s for writing.\n", filename);
        return;
    }
    
    // Write VTK header
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "3D Hydro Simulation Output\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET STRUCTURED_POINTS\n");
    fprintf(fp, "DIMENSIONS %d %d %d\n", solver->nx, solver->ny, solver->nz);
    fprintf(fp, "ORIGIN %f %f %f\n", solver->xmin, solver->ymin, solver->zmin);
    fprintf(fp, "SPACING %f %f %f\n", solver->dx, solver->dy, solver->dz);
    fprintf(fp, "POINT_DATA %d\n", solver->nx * solver->ny * solver->nz);
    
    // Write density
    fprintf(fp, "SCALARS density float 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (int i = 0; i < solver->nx; i++) {
        for (int j = 0; j < solver->ny; j++) {
            for (int k = 0; k < solver->nz; k++) {
                int idx = IDX3D(i, j, k, solver->ny, solver->nz);
                Primitive3D P = cons_to_prim_3d(solver->U[idx], solver);
                fprintf(fp, "%e\n", P.rho);
            }
        }
    }
    
    // Write pressure
    fprintf(fp, "SCALARS pressure float 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (int i = 0; i < solver->nx; i++) {
        for (int j = 0; j < solver->ny; j++) {
            for (int k = 0; k < solver->nz; k++) {
                int idx = IDX3D(i, j, k, solver->ny, solver->nz);
                Primitive3D P = cons_to_prim_3d(solver->U[idx], solver);
                fprintf(fp, "%e\n", P.p);
            }
        }
    }
    
    // Write velocity as vector
    fprintf(fp, "VECTORS velocity float\n");
    for (int i = 0; i < solver->nx; i++) {
        for (int j = 0; j < solver->ny; j++) {
            for (int k = 0; k < solver->nz; k++) {
                int idx = IDX3D(i, j, k, solver->ny, solver->nz);
                Primitive3D P = cons_to_prim_3d(solver->U[idx], solver);
                fprintf(fp, "%e %e %e\n", P.u, P.v, P.w);
            }
        }
    }
    
    fclose(fp);
}
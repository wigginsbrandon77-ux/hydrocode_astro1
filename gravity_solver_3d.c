#include "gravity_solver_3d.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

// =============================================================================
// 2D MODE DETECTION (NEW)
// =============================================================================

/**
 * @brief Check if solver should operate in 2D mode
 * 
 * Returns true if:
 * - nz == 1 (single z-zone), OR
 * - Domain has extreme aspect ratio (pancake geometry)
 * 
 * In 2D mode, the solver:
 * - Solves only 2D Poisson equation (no z-derivatives)
 * - Sets gz = 0 everywhere
 * - Avoids numerical issues from extreme grid anisotropy
 */
static bool is_2d_mode(const GravitySolver3D *solver) {
    // Definitely 2D if only one z-zone
    if (solver->nz == 1) {
        return true;
    }
    
    // Check aspect ratio - if z-dimension is much thinner than x,y
    double Lx = solver->xmax - solver->xmin;
    double Lz = solver->zmax - solver->zmin;
    double aspect_ratio = Lx / Lz;
    
    // Warn if aspect ratio is extreme but not quite 2D
    if (aspect_ratio > 10.0 && solver->nz > 1) {
        return true;  // Quasi-2D mode
    }
    
    return false;
}

// =============================================================================
// SOLVER CREATION AND DESTRUCTION
// =============================================================================

GravitySolver3D* gravity_solver_create(int nx, int ny, int nz,
                                        double xmin, double xmax,
                                        double ymin, double ymax,
                                        double zmin, double zmax,
                                        double G) {
    GravitySolver3D *solver = (GravitySolver3D*) malloc(sizeof(GravitySolver3D));
    if (!solver) {
        fprintf(stderr, "Error: Failed to allocate gravity solver structure.\n");
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
    
    solver->G = G;
    
    // Default SOR parameters
    solver->omega = 1.5;           // Relaxation parameter
    solver->max_iterations = 10000;
    solver->tolerance = 1.0e-6;
    
    solver->bc_type = GRAV_BC_ISOLATED;  // Default for astrophysical problems
    
    // Allocate arrays
    int ncells = nx * ny * nz;
    
    solver->phi = (double*) calloc(ncells, sizeof(double));
    solver->rho = NULL;  // Will point to external density array
    solver->gx = (double*) calloc(ncells, sizeof(double));
    solver->gy = (double*) calloc(ncells, sizeof(double));
    solver->gz = (double*) calloc(ncells, sizeof(double));
    
    if (!solver->phi || !solver->gx || !solver->gy || !solver->gz) {
        fprintf(stderr, "Error: Failed to allocate gravity solver arrays.\n");
        gravity_solver_destroy(solver);
        return NULL;
    }
    
    solver->last_iterations = 0;
    solver->last_residual = 0.0;
    
    // Print boundary condition info
    printf("Gravity solver: Using multipole expansion boundary conditions\n");
    printf("  - Monopole + quadrupole terms for accurate isolated boundaries\n");
    printf("  - Eliminates artifacts from phi=0 assumption\n");
    
    // Check if we're in 2D mode and inform user
    if (is_2d_mode(solver)) {
        if (solver->nz == 1) {
            printf("Gravity solver: Initialized in TRUE 2D mode (nz=1)\n");
        } else {
            double aspect = (solver->xmax - solver->xmin) / (solver->zmax - solver->zmin);
            printf("Gravity solver: Initialized in QUASI-2D mode (aspect ratio %.1f:1)\n", aspect);
            printf("  WARNING: Extreme aspect ratio detected. Consider setting nz=1 for better performance.\n");
        }
    }
    
    return solver;
}

void gravity_solver_destroy(GravitySolver3D *solver) {
    if (!solver) return;
    
    free(solver->phi);
    free(solver->gx);
    free(solver->gy);
    free(solver->gz);
    // Note: rho is a reference, don't free it
    
    free(solver);
}

// =============================================================================
// PARAMETER SETTERS
// =============================================================================

void gravity_set_solver_params(GravitySolver3D *solver, 
                                double omega, 
                                int max_iter, 
                                double tolerance) {
    solver->omega = omega;
    solver->max_iterations = max_iter;
    solver->tolerance = tolerance;
}

void gravity_set_boundary_type(GravitySolver3D *solver, GravBCType bc_type) {
    solver->bc_type = bc_type;
}

void gravity_get_diagnostics(const GravitySolver3D *solver, 
                              int *iterations, 
                              double *residual) {
    if (iterations) *iterations = solver->last_iterations;
    if (residual) *residual = solver->last_residual;
}

// =============================================================================
// MASS MOMENT CALCULATIONS FOR MULTIPOLE EXPANSION
// =============================================================================

/**
 * @brief Compute total mass and center of mass
 */
static void compute_mass_moments(const GravitySolver3D *solver,
                                 const double *rho_field,
                                 double *M_total,
                                 double *x_cm, double *y_cm, double *z_cm) {
    const int nx = solver->nx;
    const int ny = solver->ny;
    const int nz = solver->nz;
    const double dx = solver->dx;
    const double dy = solver->dy;
    const double dz = solver->dz;
    const double dV = dx * dy * dz;
    
    double M = 0.0;
    double Mx = 0.0, My = 0.0, Mz = 0.0;
    
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                int idx = IDX3D_GRAV(i, j, k, ny, nz);
                double rho = rho_field[idx];
                
                if (rho > 0.0) {
                    double x = solver->xmin + (i + 0.5) * dx;
                    double y = solver->ymin + (j + 0.5) * dy;
                    double z = solver->zmin + (k + 0.5) * dz;
                    
                    double dm = rho * dV;
                    M += dm;
                    Mx += dm * x;
                    My += dm * y;
                    Mz += dm * z;
                }
            }
        }
    }
    
    *M_total = M;
    if (M > 0.0) {
        *x_cm = Mx / M;
        *y_cm = My / M;
        *z_cm = Mz / M;
    } else {
        *x_cm = 0.5 * (solver->xmin + solver->xmax);
        *y_cm = 0.5 * (solver->ymin + solver->ymax);
        *z_cm = 0.5 * (solver->zmin + solver->zmax);
    }
}

/**
 * @brief Compute quadrupole moment tensor (for higher-order corrections)
 */
static void compute_quadrupole_moment(const GravitySolver3D *solver,
                                      const double *rho_field,
                                      double x_cm, double y_cm, double z_cm,
                                      double Q[3][3]) {
    const int nx = solver->nx;
    const int ny = solver->ny;
    const int nz = solver->nz;
    const double dx = solver->dx;
    const double dy = solver->dy;
    const double dz = solver->dz;
    const double dV = dx * dy * dz;
    
    // Initialize quadrupole tensor
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            Q[i][j] = 0.0;
        }
    }
    
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                int idx = IDX3D_GRAV(i, j, k, ny, nz);
                double rho = rho_field[idx];
                
                if (rho > 0.0) {
                    double x = solver->xmin + (i + 0.5) * dx - x_cm;
                    double y = solver->ymin + (j + 0.5) * dy - y_cm;
                    double z = solver->zmin + (k + 0.5) * dz - z_cm;
                    
                    double dm = rho * dV;
                    double r2 = x*x + y*y + z*z;
                    
                    // Q_ij = ∫ ρ (3 r_i r_j - r² δ_ij) dV
                    Q[0][0] += dm * (3.0*x*x - r2);
                    Q[0][1] += dm * (3.0*x*y);
                    Q[0][2] += dm * (3.0*x*z);
                    Q[1][0] += dm * (3.0*y*x);
                    Q[1][1] += dm * (3.0*y*y - r2);
                    Q[1][2] += dm * (3.0*y*z);
                    Q[2][0] += dm * (3.0*z*x);
                    Q[2][1] += dm * (3.0*z*y);
                    Q[2][2] += dm * (3.0*z*z - r2);
                }
            }
        }
    }
}

// =============================================================================
// BOUNDARY CONDITIONS
// =============================================================================

static void apply_boundary_conditions(GravitySolver3D *solver) {
    const int nx = solver->nx;
    const int ny = solver->ny;
    const int nz = solver->nz;
    
    if (solver->bc_type == GRAV_BC_PERIODIC) {
        // Periodic boundaries
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                int idx_0 = IDX3D_GRAV(0, j, k, ny, nz);
                int idx_nx = IDX3D_GRAV(nx-1, j, k, ny, nz);
                solver->phi[idx_0] = solver->phi[idx_nx];
            }
        }
        
        for (int i = 0; i < nx; i++) {
            for (int k = 0; k < nz; k++) {
                int idx_0 = IDX3D_GRAV(i, 0, k, ny, nz);
                int idx_ny = IDX3D_GRAV(i, ny-1, k, ny, nz);
                solver->phi[idx_0] = solver->phi[idx_ny];
            }
        }
        
        // Z-direction periodic (safe even for nz=1)
        if (nz > 1) {
            for (int i = 0; i < nx; i++) {
                for (int j = 0; j < ny; j++) {
                    int idx_0 = IDX3D_GRAV(i, j, 0, ny, nz);
                    int idx_nz = IDX3D_GRAV(i, j, nz-1, ny, nz);
                    solver->phi[idx_0] = solver->phi[idx_nz];
                }
            }
        }
        
    } else if (solver->bc_type == GRAV_BC_ISOLATED) {
        // Multipole expansion boundary condition for isolated systems
        // This is much more accurate than phi=0 for systems like binary mergers
        
        // Compute mass moments
        double M_total, x_cm, y_cm, z_cm;
        compute_mass_moments(solver, solver->rho, &M_total, &x_cm, &y_cm, &z_cm);
        
        // Compute quadrupole moment for better accuracy
        double Q[3][3];
        compute_quadrupole_moment(solver, solver->rho, x_cm, y_cm, z_cm, Q);
        
        const double G = solver->G;
        
        // Apply multipole expansion at all boundary faces
        
        // X-boundaries (i=0 and i=nx-1)
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                // X-min boundary
                {
                    double x = solver->xmin + 0.5 * solver->dx;
                    double y = solver->ymin + (j + 0.5) * solver->dy;
                    double z = solver->zmin + (k + 0.5) * solver->dz;
                    
                    double rx = x - x_cm;
                    double ry = y - y_cm;
                    double rz = z - z_cm;
                    double r = sqrt(rx*rx + ry*ry + rz*rz);
                    
                    // Monopole term: phi = -G*M/r
                    double phi_mono = 0.0;
                    if (r > 0.0) {
                        phi_mono = -G * M_total / r;
                        
                        // Quadrupole correction: phi_quad = -G/(2*r^5) * Q_ij * r_i * r_j
                        double r5 = r * r * r * r * r;
                        double Qr = Q[0][0]*rx*rx + Q[0][1]*rx*ry + Q[0][2]*rx*rz +
                                    Q[1][0]*ry*rx + Q[1][1]*ry*ry + Q[1][2]*ry*rz +
                                    Q[2][0]*rz*rx + Q[2][1]*rz*ry + Q[2][2]*rz*rz;
                        double phi_quad = -G * Qr / (2.0 * r5);
                        
                        phi_mono += phi_quad;
                    }
                    
                    solver->phi[IDX3D_GRAV(0, j, k, ny, nz)] = phi_mono;
                }
                
                // X-max boundary
                {
                    double x = solver->xmax - 0.5 * solver->dx;
                    double y = solver->ymin + (j + 0.5) * solver->dy;
                    double z = solver->zmin + (k + 0.5) * solver->dz;
                    
                    double rx = x - x_cm;
                    double ry = y - y_cm;
                    double rz = z - z_cm;
                    double r = sqrt(rx*rx + ry*ry + rz*rz);
                    
                    double phi_mono = 0.0;
                    if (r > 0.0) {
                        phi_mono = -G * M_total / r;
                        
                        double r5 = r * r * r * r * r;
                        double Qr = Q[0][0]*rx*rx + Q[0][1]*rx*ry + Q[0][2]*rx*rz +
                                    Q[1][0]*ry*rx + Q[1][1]*ry*ry + Q[1][2]*ry*rz +
                                    Q[2][0]*rz*rx + Q[2][1]*rz*ry + Q[2][2]*rz*rz;
                        double phi_quad = -G * Qr / (2.0 * r5);
                        
                        phi_mono += phi_quad;
                    }
                    
                    solver->phi[IDX3D_GRAV(nx-1, j, k, ny, nz)] = phi_mono;
                }
            }
        }
        
        // Y-boundaries (j=0 and j=ny-1)
        for (int i = 0; i < nx; i++) {
            for (int k = 0; k < nz; k++) {
                // Y-min boundary
                {
                    double x = solver->xmin + (i + 0.5) * solver->dx;
                    double y = solver->ymin + 0.5 * solver->dy;
                    double z = solver->zmin + (k + 0.5) * solver->dz;
                    
                    double rx = x - x_cm;
                    double ry = y - y_cm;
                    double rz = z - z_cm;
                    double r = sqrt(rx*rx + ry*ry + rz*rz);
                    
                    double phi_mono = 0.0;
                    if (r > 0.0) {
                        phi_mono = -G * M_total / r;
                        
                        double r5 = r * r * r * r * r;
                        double Qr = Q[0][0]*rx*rx + Q[0][1]*rx*ry + Q[0][2]*rx*rz +
                                    Q[1][0]*ry*rx + Q[1][1]*ry*ry + Q[1][2]*ry*rz +
                                    Q[2][0]*rz*rx + Q[2][1]*rz*ry + Q[2][2]*rz*rz;
                        double phi_quad = -G * Qr / (2.0 * r5);
                        
                        phi_mono += phi_quad;
                    }
                    
                    solver->phi[IDX3D_GRAV(i, 0, k, ny, nz)] = phi_mono;
                }
                
                // Y-max boundary
                {
                    double x = solver->xmin + (i + 0.5) * solver->dx;
                    double y = solver->ymax - 0.5 * solver->dy;
                    double z = solver->zmin + (k + 0.5) * solver->dz;
                    
                    double rx = x - x_cm;
                    double ry = y - y_cm;
                    double rz = z - z_cm;
                    double r = sqrt(rx*rx + ry*ry + rz*rz);
                    
                    double phi_mono = 0.0;
                    if (r > 0.0) {
                        phi_mono = -G * M_total / r;
                        
                        double r5 = r * r * r * r * r;
                        double Qr = Q[0][0]*rx*rx + Q[0][1]*rx*ry + Q[0][2]*rx*rz +
                                    Q[1][0]*ry*rx + Q[1][1]*ry*ry + Q[1][2]*ry*rz +
                                    Q[2][0]*rz*rx + Q[2][1]*rz*ry + Q[2][2]*rz*rz;
                        double phi_quad = -G * Qr / (2.0 * r5);
                        
                        phi_mono += phi_quad;
                    }
                    
                    solver->phi[IDX3D_GRAV(i, ny-1, k, ny, nz)] = phi_mono;
                }
            }
        }
        
        // Z-boundaries (k=0 and k=nz-1, only if nz > 1)
        if (nz > 1) {
            for (int i = 0; i < nx; i++) {
                for (int j = 0; j < ny; j++) {
                    // Z-min boundary
                    {
                        double x = solver->xmin + (i + 0.5) * solver->dx;
                        double y = solver->ymin + (j + 0.5) * solver->dy;
                        double z = solver->zmin + 0.5 * solver->dz;
                        
                        double rx = x - x_cm;
                        double ry = y - y_cm;
                        double rz = z - z_cm;
                        double r = sqrt(rx*rx + ry*ry + rz*rz);
                        
                        double phi_mono = 0.0;
                        if (r > 0.0) {
                            phi_mono = -G * M_total / r;
                            
                            double r5 = r * r * r * r * r;
                            double Qr = Q[0][0]*rx*rx + Q[0][1]*rx*ry + Q[0][2]*rx*rz +
                                        Q[1][0]*ry*rx + Q[1][1]*ry*ry + Q[1][2]*ry*rz +
                                        Q[2][0]*rz*rx + Q[2][1]*rz*ry + Q[2][2]*rz*rz;
                            double phi_quad = -G * Qr / (2.0 * r5);
                            
                            phi_mono += phi_quad;
                        }
                        
                        solver->phi[IDX3D_GRAV(i, j, 0, ny, nz)] = phi_mono;
                    }
                    
                    // Z-max boundary
                    {
                        double x = solver->xmin + (i + 0.5) * solver->dx;
                        double y = solver->ymin + (j + 0.5) * solver->dy;
                        double z = solver->zmax - 0.5 * solver->dz;
                        
                        double rx = x - x_cm;
                        double ry = y - y_cm;
                        double rz = z - z_cm;
                        double r = sqrt(rx*rx + ry*ry + rz*rz);
                        
                        double phi_mono = 0.0;
                        if (r > 0.0) {
                            phi_mono = -G * M_total / r;
                            
                            double r5 = r * r * r * r * r;
                            double Qr = Q[0][0]*rx*rx + Q[0][1]*rx*ry + Q[0][2]*rx*rz +
                                        Q[1][0]*ry*rx + Q[1][1]*ry*ry + Q[1][2]*ry*rz +
                                        Q[2][0]*rz*rx + Q[2][1]*rz*ry + Q[2][2]*rz*rz;
                            double phi_quad = -G * Qr / (2.0 * r5);
                            
                            phi_mono += phi_quad;
                        }
                        
                        solver->phi[IDX3D_GRAV(i, j, nz-1, ny, nz)] = phi_mono;
                    }
                }
            }
        }
        
    } else if (solver->bc_type == GRAV_BC_DIRICHLET) {
        // Fixed φ = 0 at boundaries (already initialized to zero)
        // This is the same as isolated for now
        // Could be extended to allow non-zero boundary values
    }
}

// =============================================================================
// POISSON SOLVER USING SOR (MODIFIED FOR 2D SUPPORT)
// =============================================================================

int gravity_solve_potential(GravitySolver3D *solver, const double *rho_field) {
    const int nx = solver->nx;
    const int ny = solver->ny;
    const int nz = solver->nz;
    
    const double dx = solver->dx;
    const double dy = solver->dy;
    const double dz = solver->dz;
    
    const double dx2 = dx * dx;
    const double dy2 = dy * dy;
    const double dz2 = dz * dz;
    
    const double dx2_inv = 1.0 / dx2;
    const double dy2_inv = 1.0 / dy2;
    const double dz2_inv = 1.0 / dz2;
    
    const double four_pi_G = 4.0 * M_PI * solver->G;
    const double omega = solver->omega;
    
    // Store reference to density field
    solver->rho = (double*) rho_field;
    
    // Check if we should use 2D mode
    bool mode_2d = is_2d_mode(solver);
    bool true_2d = (nz == 1);
    
    // Pre-compute factor for the solver
    double factor;
    if (true_2d) {
        // 2D: only x and y derivatives
        factor = 2.0 * (dx2_inv + dy2_inv);
    } else {
        // 3D: all three derivatives
        factor = 2.0 * (dx2_inv + dy2_inv + dz2_inv);
    }
    
    // SOR iteration
    int iter;
    double max_residual = 0.0;
    
    // Print mode information on first call
    static bool first_call = true;
    if (first_call) {
        if (true_2d) {
            printf("Gravity solver: Using TRUE 2D Poisson equation\n");
            printf("  Solving: (∂²/∂x² + ∂²/∂y²)φ = 4πGρ\n");
        } else if (mode_2d) {
            printf("Gravity solver: Using QUASI-2D mode (WARNING: may have convergence issues)\n");
            printf("  Consider setting nz=1 and bc_z=PERIODIC for true 2D\n");
        } else {
            printf("Gravity solver: Using full 3D Poisson equation\n");
        }
        first_call = false;
    }
    
    for (iter = 0; iter < solver->max_iterations; iter++) {
        max_residual = 0.0;
        
        if (true_2d) {
            // ====================================================================
            // TRUE 2D MODE: Single z-level, no z-derivatives
            // ====================================================================
            const int k = 0;  // Single z-level
            
            // Red-Black SOR for better convergence
            // Red sweep (i+j even)
            for (int i = 1; i < nx-1; i++) {
                for (int j = 1; j < ny-1; j++) {
                    if ((i + j) % 2 == 0) {
                        int idx = IDX3D_GRAV(i, j, k, ny, nz);
                        int idx_im1 = IDX3D_GRAV(i-1, j, k, ny, nz);
                        int idx_ip1 = IDX3D_GRAV(i+1, j, k, ny, nz);
                        int idx_jm1 = IDX3D_GRAV(i, j-1, k, ny, nz);
                        int idx_jp1 = IDX3D_GRAV(i, j+1, k, ny, nz);
                        
                        double phi_old = solver->phi[idx];
                        
                        // 2D Laplacian (no z-derivative!)
                        double laplacian_phi = 
                            (solver->phi[idx_im1] + solver->phi[idx_ip1]) * dx2_inv +
                            (solver->phi[idx_jm1] + solver->phi[idx_jp1]) * dy2_inv;
                        
                        double rhs = four_pi_G * rho_field[idx];
                        
                        double phi_new = (laplacian_phi - rhs) / factor;
                        
                        // SOR update
                        solver->phi[idx] = phi_old + omega * (phi_new - phi_old);
                        
                        // Track residual
                        double residual = fabs(phi_new - phi_old);
                        max_residual = MAX(max_residual, residual);
                    }
                }
            }
            
            // Black sweep (i+j odd)
            for (int i = 1; i < nx-1; i++) {
                for (int j = 1; j < ny-1; j++) {
                    if ((i + j) % 2 == 1) {
                        int idx = IDX3D_GRAV(i, j, k, ny, nz);
                        int idx_im1 = IDX3D_GRAV(i-1, j, k, ny, nz);
                        int idx_ip1 = IDX3D_GRAV(i+1, j, k, ny, nz);
                        int idx_jm1 = IDX3D_GRAV(i, j-1, k, ny, nz);
                        int idx_jp1 = IDX3D_GRAV(i, j+1, k, ny, nz);
                        
                        double phi_old = solver->phi[idx];
                        
                        double laplacian_phi = 
                            (solver->phi[idx_im1] + solver->phi[idx_ip1]) * dx2_inv +
                            (solver->phi[idx_jm1] + solver->phi[idx_jp1]) * dy2_inv;
                        
                        double rhs = four_pi_G * rho_field[idx];
                        
                        double phi_new = (laplacian_phi - rhs) / factor;
                        
                        solver->phi[idx] = phi_old + omega * (phi_new - phi_old);
                        
                        double residual = fabs(phi_new - phi_old);
                        max_residual = MAX(max_residual, residual);
                    }
                }
            }
            
        } else {
            // ====================================================================
            // STANDARD 3D MODE: Full 3D Poisson equation
            // ====================================================================
            
            // Red-Black SOR for better convergence
            // Red sweep (i+j+k even)
            for (int i = 1; i < nx-1; i++) {
                for (int j = 1; j < ny-1; j++) {
                    for (int k = 1; k < nz-1; k++) {
                        if ((i + j + k) % 2 == 0) {
                            int idx = IDX3D_GRAV(i, j, k, ny, nz);
                            int idx_im1 = IDX3D_GRAV(i-1, j, k, ny, nz);
                            int idx_ip1 = IDX3D_GRAV(i+1, j, k, ny, nz);
                            int idx_jm1 = IDX3D_GRAV(i, j-1, k, ny, nz);
                            int idx_jp1 = IDX3D_GRAV(i, j+1, k, ny, nz);
                            int idx_km1 = IDX3D_GRAV(i, j, k-1, ny, nz);
                            int idx_kp1 = IDX3D_GRAV(i, j, k+1, ny, nz);
                            
                            double phi_old = solver->phi[idx];
                            
                            // ∇²φ = 4πGρ
                            double laplacian_phi = 
                                (solver->phi[idx_im1] + solver->phi[idx_ip1]) * dx2_inv +
                                (solver->phi[idx_jm1] + solver->phi[idx_jp1]) * dy2_inv +
                                (solver->phi[idx_km1] + solver->phi[idx_kp1]) * dz2_inv;
                            
                            double rhs = four_pi_G * rho_field[idx];
                            
                            double phi_new = (laplacian_phi - rhs) / factor;
                            
                            // SOR update
                            solver->phi[idx] = phi_old + omega * (phi_new - phi_old);
                            
                            // Track residual
                            double residual = fabs(phi_new - phi_old);
                            max_residual = MAX(max_residual, residual);
                        }
                    }
                }
            }
            
            // Black sweep (i+j+k odd)
            for (int i = 1; i < nx-1; i++) {
                for (int j = 1; j < ny-1; j++) {
                    for (int k = 1; k < nz-1; k++) {
                        if ((i + j + k) % 2 == 1) {
                            int idx = IDX3D_GRAV(i, j, k, ny, nz);
                            int idx_im1 = IDX3D_GRAV(i-1, j, k, ny, nz);
                            int idx_ip1 = IDX3D_GRAV(i+1, j, k, ny, nz);
                            int idx_jm1 = IDX3D_GRAV(i, j-1, k, ny, nz);
                            int idx_jp1 = IDX3D_GRAV(i, j+1, k, ny, nz);
                            int idx_km1 = IDX3D_GRAV(i, j, k-1, ny, nz);
                            int idx_kp1 = IDX3D_GRAV(i, j, k+1, ny, nz);
                            
                            double phi_old = solver->phi[idx];
                            
                            double laplacian_phi = 
                                (solver->phi[idx_im1] + solver->phi[idx_ip1]) * dx2_inv +
                                (solver->phi[idx_jm1] + solver->phi[idx_jp1]) * dy2_inv +
                                (solver->phi[idx_km1] + solver->phi[idx_kp1]) * dz2_inv;
                            
                            double rhs = four_pi_G * rho_field[idx];
                            
                            double phi_new = (laplacian_phi - rhs) / factor;
                            
                            solver->phi[idx] = phi_old + omega * (phi_new - phi_old);
                            
                            double residual = fabs(phi_new - phi_old);
                            max_residual = MAX(max_residual, residual);
                        }
                    }
                }
            }
        }
        
        // Apply boundary conditions
        apply_boundary_conditions(solver);
        
        // Check convergence
        if (max_residual < solver->tolerance) {
            solver->last_iterations = iter + 1;
            solver->last_residual = max_residual;
            
            if (true_2d) {
                //printf("Gravity solver converged in %d iterations (2D mode)\n", iter + 1);
            } else {
                //printf("Gravity solver converged in %d iterations (3D mode)\n", iter + 1);
            }
            //printf("  Final residual: %.3e\n", max_residual);
            
            return 0;  // Success
        }
    }
    
    // Did not converge
    solver->last_iterations = iter;
    solver->last_residual = max_residual;
    
    fprintf(stderr, "Warning: Gravity solver did not converge in %d iterations.\n", iter);
    fprintf(stderr, "         Final residual: %.3e (tolerance: %.3e)\n", 
            max_residual, solver->tolerance);
    
    if (mode_2d && !true_2d) {
        fprintf(stderr, "\n");
        fprintf(stderr, "HINT: Convergence failure may be due to extreme aspect ratio.\n");
        fprintf(stderr, "      Try setting nz=1 and bc_zmin=bc_zmax=PERIODIC for true 2D simulation.\n");
        fprintf(stderr, "      Current aspect ratio: %.1f:1\n", 
                (solver->xmax - solver->xmin) / (solver->zmax - solver->zmin));
    }
    
    return -1;  // Failed to converge
}

// =============================================================================
// COMPUTE GRAVITATIONAL ACCELERATION (MODIFIED FOR 2D SUPPORT)
// =============================================================================

void gravity_compute_acceleration(GravitySolver3D *solver) {
    const int nx = solver->nx;
    const int ny = solver->ny;
    const int nz = solver->nz;
    
    const double dx = solver->dx;
    const double dy = solver->dy;
    const double dz = solver->dz;
    
    bool true_2d = (nz == 1);
    
    if (true_2d) {
        // ====================================================================
        // TRUE 2D MODE: gz = 0 everywhere, compute only gx and gy
        // ====================================================================
        const int k = 0;
        
        // Interior points
        for (int i = 1; i < nx-1; i++) {
            for (int j = 1; j < ny-1; j++) {
                int idx = IDX3D_GRAV(i, j, k, ny, nz);
                int idx_im1 = IDX3D_GRAV(i-1, j, k, ny, nz);
                int idx_ip1 = IDX3D_GRAV(i+1, j, k, ny, nz);
                int idx_jm1 = IDX3D_GRAV(i, j-1, k, ny, nz);
                int idx_jp1 = IDX3D_GRAV(i, j+1, k, ny, nz);
                
                // g = -∇φ (2D: only x and y components)
                solver->gx[idx] = -(solver->phi[idx_ip1] - solver->phi[idx_im1]) / (2.0 * dx);
                solver->gy[idx] = -(solver->phi[idx_jp1] - solver->phi[idx_jm1]) / (2.0 * dy);
                solver->gz[idx] = 0.0;  // No z-acceleration in 2D
            }
        }
        
        // Boundary points
        for (int j = 0; j < ny; j++) {
            int idx_0 = IDX3D_GRAV(0, j, k, ny, nz);
            int idx_nx = IDX3D_GRAV(nx-1, j, k, ny, nz);
            
            solver->gx[idx_0] = 0.0;
            solver->gy[idx_0] = 0.0;
            solver->gz[idx_0] = 0.0;
            
            solver->gx[idx_nx] = 0.0;
            solver->gy[idx_nx] = 0.0;
            solver->gz[idx_nx] = 0.0;
        }
        
        for (int i = 0; i < nx; i++) {
            int idx_0 = IDX3D_GRAV(i, 0, k, ny, nz);
            int idx_ny = IDX3D_GRAV(i, ny-1, k, ny, nz);
            
            solver->gx[idx_0] = 0.0;
            solver->gy[idx_0] = 0.0;
            solver->gz[idx_0] = 0.0;
            
            solver->gx[idx_ny] = 0.0;
            solver->gy[idx_ny] = 0.0;
            solver->gz[idx_ny] = 0.0;
        }
        
    } else {
        // ====================================================================
        // STANDARD 3D MODE: Full 3D acceleration
        // ====================================================================
        
        // Compute g = -∇φ using centered differences
        // Interior points
        for (int i = 1; i < nx-1; i++) {
            for (int j = 1; j < ny-1; j++) {
                for (int k = 1; k < nz-1; k++) {
                    int idx = IDX3D_GRAV(i, j, k, ny, nz);
                    int idx_im1 = IDX3D_GRAV(i-1, j, k, ny, nz);
                    int idx_ip1 = IDX3D_GRAV(i+1, j, k, ny, nz);
                    int idx_jm1 = IDX3D_GRAV(i, j-1, k, ny, nz);
                    int idx_jp1 = IDX3D_GRAV(i, j+1, k, ny, nz);
                    int idx_km1 = IDX3D_GRAV(i, j, k-1, ny, nz);
                    int idx_kp1 = IDX3D_GRAV(i, j, k+1, ny, nz);
                    
                    // g = -∇φ
                    solver->gx[idx] = -(solver->phi[idx_ip1] - solver->phi[idx_im1]) / (2.0 * dx);
                    solver->gy[idx] = -(solver->phi[idx_jp1] - solver->phi[idx_jm1]) / (2.0 * dy);
                    solver->gz[idx] = -(solver->phi[idx_kp1] - solver->phi[idx_km1]) / (2.0 * dz);
                }
            }
        }
        
        // Boundary points - extrapolate from interior cells
        // This is more accurate than setting g=0 for systems with mass near boundaries
        
        // X-boundaries
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                // X-min boundary: copy from i=1
                int idx_0 = IDX3D_GRAV(0, j, k, ny, nz);
                int idx_1 = IDX3D_GRAV(1, j, k, ny, nz);
                solver->gx[idx_0] = solver->gx[idx_1];
                solver->gy[idx_0] = solver->gy[idx_1];
                solver->gz[idx_0] = solver->gz[idx_1];
                
                // X-max boundary: copy from i=nx-2
                int idx_nx = IDX3D_GRAV(nx-1, j, k, ny, nz);
                int idx_nx1 = IDX3D_GRAV(nx-2, j, k, ny, nz);
                solver->gx[idx_nx] = solver->gx[idx_nx1];
                solver->gy[idx_nx] = solver->gy[idx_nx1];
                solver->gz[idx_nx] = solver->gz[idx_nx1];
            }
        }
        
        // Y-boundaries
        for (int i = 0; i < nx; i++) {
            for (int k = 0; k < nz; k++) {
                // Y-min boundary: copy from j=1
                int idx_0 = IDX3D_GRAV(i, 0, k, ny, nz);
                int idx_1 = IDX3D_GRAV(i, 1, k, ny, nz);
                solver->gx[idx_0] = solver->gx[idx_1];
                solver->gy[idx_0] = solver->gy[idx_1];
                solver->gz[idx_0] = solver->gz[idx_1];
                
                // Y-max boundary: copy from j=ny-2
                int idx_ny = IDX3D_GRAV(i, ny-1, k, ny, nz);
                int idx_ny1 = IDX3D_GRAV(i, ny-2, k, ny, nz);
                solver->gx[idx_ny] = solver->gx[idx_ny1];
                solver->gy[idx_ny] = solver->gy[idx_ny1];
                solver->gz[idx_ny] = solver->gz[idx_ny1];
            }
        }
        
        // Z-boundaries
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                // Z-min boundary: copy from k=1
                int idx_0 = IDX3D_GRAV(i, j, 0, ny, nz);
                int idx_1 = IDX3D_GRAV(i, j, 1, ny, nz);
                solver->gx[idx_0] = solver->gx[idx_1];
                solver->gy[idx_0] = solver->gy[idx_1];
                solver->gz[idx_0] = solver->gz[idx_1];
                
                // Z-max boundary: copy from k=nz-2
                int idx_nz = IDX3D_GRAV(i, j, nz-1, ny, nz);
                int idx_nz1 = IDX3D_GRAV(i, j, nz-2, ny, nz);
                solver->gx[idx_nz] = solver->gx[idx_nz1];
                solver->gy[idx_nz] = solver->gy[idx_nz1];
                solver->gz[idx_nz] = solver->gz[idx_nz1];
            }
        }
    }
}

// =============================================================================
// CONVENIENCE FUNCTIONS
// =============================================================================

int gravity_solve(GravitySolver3D *solver, const double *rho_field) {
    int status = gravity_solve_potential(solver, rho_field);
    if (status == 0) {
        gravity_compute_acceleration(solver);
    }
    return status;
}

void gravity_write_potential(const GravitySolver3D *solver, const char *filename) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Error: Could not open file %s for writing potential.\n", filename);
        return;
    }
    
    fprintf(fp, "x,y,z,phi,rho,gx,gy,gz\n");
    
    for (int i = 0; i < solver->nx; i++) {
        for (int j = 0; j < solver->ny; j++) {
            for (int k = 0; k < solver->nz; k++) {
                int idx = IDX3D_GRAV(i, j, k, solver->ny, solver->nz);
                
                double x = solver->xmin + (i + 0.5) * solver->dx;
                double y = solver->ymin + (j + 0.5) * solver->dy;
                double z = solver->zmin + (k + 0.5) * solver->dz;
                
                double rho_val = solver->rho ? solver->rho[idx] : 0.0;
                
                fprintf(fp, "%.6f,%.6f,%.6f,%.6e,%.6e,%.6e,%.6e,%.6e\n",
                        x, y, z, 
                        solver->phi[idx], rho_val,
                        solver->gx[idx], solver->gy[idx], solver->gz[idx]);
            }
        }
    }
    
    fclose(fp);
    printf("Wrote gravitational potential to %s\n", filename);
}
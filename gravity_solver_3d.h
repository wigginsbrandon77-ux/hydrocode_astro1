#ifndef GRAVITY_SOLVER_3D_H
#define GRAVITY_SOLVER_3D_H

#include <stdbool.h>

// Indexing macro for 3D arrays
#define IDX3D_GRAV(i, j, k, ny, nz) (((i) * (ny) + (j)) * (nz) + (k))

// Boundary condition types for gravity
typedef enum {
    GRAV_BC_PERIODIC,
    GRAV_BC_ISOLATED,     // φ → 0 at infinity (for isolated systems)
    GRAV_BC_DIRICHLET     // Fixed value at boundary
} GravBCType;

// Gravity solver structure
typedef struct GravitySolver3D {
    int nx, ny, nz;           // Grid dimensions
    double dx, dy, dz;        // Cell spacing
    double xmin, xmax;        // Domain boundaries
    double ymin, ymax;
    double zmin, zmax;
    
    double G;                 // Gravitational constant
    
    double *phi;              // Gravitational potential
    double *rho;              // Density field (reference to hydro density)
    double *gx;               // x-component of gravitational acceleration
    double *gy;               // y-component of gravitational acceleration
    double *gz;               // z-component of gravitational acceleration
    
    // SOR solver parameters
    double omega;             // Relaxation parameter (1.0-2.0, typically ~1.5)
    int max_iterations;       // Maximum number of SOR iterations
    double tolerance;         // Convergence tolerance
    
    GravBCType bc_type;       // Boundary condition type
    
    // Diagnostics
    int last_iterations;      // Number of iterations in last solve
    double last_residual;     // Final residual in last solve
    
} GravitySolver3D;

// =============================================================================
// Public Function Prototypes
// =============================================================================

/**
 * @brief Create and initialize a 3D gravity solver
 * 
 * @param nx Number of cells in x-direction
 * @param ny Number of cells in y-direction
 * @param nz Number of cells in z-direction
 * @param xmin Minimum x coordinate
 * @param xmax Maximum x coordinate
 * @param ymin Minimum y coordinate
 * @param ymax Maximum y coordinate
 * @param zmin Minimum z coordinate
 * @param zmax Maximum z coordinate
 * @param G Gravitational constant
 * @return Pointer to new GravitySolver3D or NULL on failure
 */
GravitySolver3D* gravity_solver_create(int nx, int ny, int nz,
                                        double xmin, double xmax,
                                        double ymin, double ymax,
                                        double zmin, double zmax,
                                        double G);

/**
 * @brief Destroy gravity solver and free memory
 * 
 * @param solver Pointer to gravity solver
 */
void gravity_solver_destroy(GravitySolver3D *solver);

/**
 * @brief Solve Poisson equation for gravitational potential
 * 
 * Solves: ∇²φ = 4πGρ
 * using Successive Over-Relaxation (SOR) method
 * 
 * @param solver Pointer to gravity solver
 * @param rho_field Density field (nx × ny × nz array)
 * @return 0 on success, -1 on failure
 */
int gravity_solve_potential(GravitySolver3D *solver, const double *rho_field);

/**
 * @brief Compute gravitational acceleration from potential
 * 
 * Computes: g = -∇φ
 * using centered finite differences
 * 
 * @param solver Pointer to gravity solver
 */
void gravity_compute_acceleration(GravitySolver3D *solver);

/**
 * @brief Set SOR solver parameters
 * 
 * @param solver Pointer to gravity solver
 * @param omega Relaxation parameter (1.0-2.0, default 1.5)
 * @param max_iter Maximum iterations (default 10000)
 * @param tolerance Convergence tolerance (default 1e-6)
 */
void gravity_set_solver_params(GravitySolver3D *solver, 
                                double omega, 
                                int max_iter, 
                                double tolerance);

/**
 * @brief Set boundary condition type
 * 
 * @param solver Pointer to gravity solver
 * @param bc_type Boundary condition type
 */
void gravity_set_boundary_type(GravitySolver3D *solver, GravBCType bc_type);

/**
 * @brief Get diagnostic information from last solve
 * 
 * @param solver Pointer to gravity solver
 * @param iterations Output: number of iterations used
 * @param residual Output: final residual achieved
 */
void gravity_get_diagnostics(const GravitySolver3D *solver, 
                              int *iterations, 
                              double *residual);

/**
 * @brief Complete gravity solve: potential + acceleration
 * 
 * Convenience function that calls both gravity_solve_potential
 * and gravity_compute_acceleration
 * 
 * @param solver Pointer to gravity solver
 * @param rho_field Density field
 * @return 0 on success, -1 on failure
 */
int gravity_solve(GravitySolver3D *solver, const double *rho_field);

/**
 * @brief Write gravitational potential to file (for debugging)
 * 
 * @param solver Pointer to gravity solver
 * @param filename Output filename
 */
void gravity_write_potential(const GravitySolver3D *solver, const char *filename);

#endif // GRAVITY_SOLVER_3D_H

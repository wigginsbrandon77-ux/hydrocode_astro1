#ifndef INPUT_DECK_H
#define INPUT_DECK_H

#include <stdio.h>
#include <stdbool.h>

// Problem types (extended from 2D)
typedef enum {
    KELVIN_HELMHOLTZ,
    BLAST_WAVE,
    RAYLEIGH_TAYLOR,
    RIEMANN_2D,
    STALLED_SHOCK,
    STELLAR_CONVECTION,
    STELLAR_PROFILE,
    GALAXY_DISK,
    UNIFORM_DISK_TEST, 
    WD_MERGER,
    SUPERNOVA,
    KELVIN_HELMHOLTZ_3D,
    BLAST_WAVE_3D,
    TURBULENCE_3D,
    SEDOV_BLAST_3D,
    STELLAR_PROFILE_3D,
    BINARY_STAR_3D,
    STELLAR_OCTANT_EXPLOSION_3D
} ProblemType;

// Boundary condition types
typedef enum {
    BC_PERIODIC,
    BC_REFLECTIVE,
    BC_OUTFLOW
} BCType;

// Structure to hold all input parameters
typedef struct {
    // Grid parameters
    int nx, ny, nz;
    double xmin, xmax;
    double ymin, ymax;
    double zmin, zmax;
    
    // Time parameters
    double t_end;
    double output_dt;
    double cfl;
    
    // Physics parameters
    double gamma;
    double G_const;  // Gravitational constant
    
    // Problem setup
    ProblemType problem;
    
    // Boundary conditions
    BCType bc_xmin, bc_xmax;
    BCType bc_ymin, bc_ymax;
    BCType bc_zmin, bc_zmax;
    
    // Additional physics (for specific problems)
    double G_M_PNS;      // PNS gravitational parameter
    double L_nu;         // Neutrino luminosity
    double g_strength;   // Gravity strength
    double L_heat;       // Heating luminosity
    
    // Output control
    char output_dir[256];
    
    // Restart parameters
    bool restart;
    char restart_file[256];
    int restart_step;
    double restart_time;
    
    // Binary merger parameters (for BINARY_STAR_3D problem)
    char star1_profile[256];      // Stellar profile file for star 1
    char star2_profile[256];      // Stellar profile file for star 2
    double star1_x, star1_y, star1_z;  // Position of star 1
    double star2_x, star2_y, star2_z;  // Position of star 2
    double star1_vx, star1_vy, star1_vz;  // Velocity of star 1
    double star2_vx, star2_vy, star2_vz;  // Velocity of star 2
    double binary_separation;     // Initial separation (if positions not specified)
    double binary_orbital_velocity;  // Orbital velocity (if velocities not specified)
    bool binary_auto_orbit;       // Automatically calculate orbital parameters
    
    // Stellar octant explosion parameters (for STELLAR_OCTANT_EXPLOSION_3D)
    char stellar_profile[256];    // Stellar profile file for octant setup
    double explosion_energy;      // Energy to deposit in code units
    double explosion_radius;      // Radius below which to deposit energy
    char explosion_type[32];      // "thermal" or "kinetic"
    
    // Perturbation parameters (to seed turbulence/asymmetries)
    bool enable_perturbations;         // Master switch for perturbations
    int perturbation_seed;             // Random seed for reproducibility
    
    // Random velocity perturbations (seeds convection/turbulence)
    double perturb_velocity_amplitude; // Amplitude of velocity perturbations (v/c_s)
    double perturb_velocity_kmin;      // Minimum wavenumber (2π/L_max)
    double perturb_velocity_kmax;      // Maximum wavenumber (2π/L_min)
    int perturb_velocity_modes;        // Number of Fourier modes to sum
    
    // Random density perturbations (composition variations/shells)
    double perturb_density_amplitude;  // Fractional density perturbation (δρ/ρ)
    double perturb_density_kmin;       // Minimum wavenumber
    double perturb_density_kmax;       // Maximum wavenumber
    int perturb_density_modes;         // Number of Fourier modes to sum
    
    // Spherical harmonic mode perturbations (large-scale asymmetries)
    bool enable_sph_harmonic_mode;     // Enable single spherical harmonic mode
    int sph_harmonic_l;                // Degree l (0, 1, 2, ...)
    int sph_harmonic_m;                // Order m (-l <= m <= l)
    double sph_harmonic_amp_rho;       // Density amplitude (fractional)
    double sph_harmonic_amp_vel;       // Velocity amplitude (v/c_s)
    
    // Torus parameters (gas ring around star)
    bool enable_torus;                 // Enable torus/ring of gas
    double torus_major_radius;         // Major radius R (distance from origin to tube center)
    double torus_minor_radius;         // Minor radius r (tube/ring thickness)
    double torus_density;              // Density of torus material
    double torus_pressure;             // Pressure of torus material
    double torus_rotation_velocity;    // Azimuthal rotation velocity of torus
    
} InputDeck;

// Function prototypes
InputDeck* input_deck_create(const char *filename);
void input_deck_destroy(InputDeck *deck);
void input_deck_print(const InputDeck *deck);
ProblemType parse_problem_type(const char *str);
BCType parse_bc_type(const char *str);

#endif // INPUT_DECK_H
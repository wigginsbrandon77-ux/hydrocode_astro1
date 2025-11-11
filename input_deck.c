#include "input_deck.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

// Helper function to trim whitespace from a string
static char* trim_whitespace(char *str) {
    char *end;
    // Trim leading space
    while(isspace((unsigned char)*str)) str++;
    if(*str == 0) return str;
    // Trim trailing space
    end = str + strlen(str) - 1;
    while(end > str && isspace((unsigned char)*end)) end--;
    end[1] = '\0';
    return str;
}

// Parse problem type from string
ProblemType parse_problem_type(const char *str) {
    if (strcmp(str, "KELVIN_HELMHOLTZ") == 0) return KELVIN_HELMHOLTZ;
    if (strcmp(str, "BLAST_WAVE") == 0) return BLAST_WAVE;
    if (strcmp(str, "RAYLEIGH_TAYLOR") == 0) return RAYLEIGH_TAYLOR;
    if (strcmp(str, "RIEMANN_2D") == 0) return RIEMANN_2D;
    if (strcmp(str, "STALLED_SHOCK") == 0) return STALLED_SHOCK;
    if (strcmp(str, "STELLAR_CONVECTION") == 0) return STELLAR_CONVECTION;
    if (strcmp(str, "STELLAR_PROFILE") == 0) return STELLAR_PROFILE;
    if (strcmp(str, "GALAXY_DISK") == 0) return GALAXY_DISK;
    if (strcmp(str, "UNIFORM_DISK_TEST") == 0) return UNIFORM_DISK_TEST;
    if (strcmp(str, "WD_MERGER") == 0) return WD_MERGER;
    if (strcmp(str, "SUPERNOVA") == 0) return SUPERNOVA;
    if (strcmp(str, "KELVIN_HELMHOLTZ_3D") == 0) return KELVIN_HELMHOLTZ_3D;
    if (strcmp(str, "BLAST_WAVE_3D") == 0) return BLAST_WAVE_3D;
    if (strcmp(str, "TURBULENCE_3D") == 0) return TURBULENCE_3D;
    if (strcmp(str, "SEDOV_BLAST_3D") == 0) return SEDOV_BLAST_3D;
    if (strcmp(str, "STELLAR_PROFILE_3D") == 0) return STELLAR_PROFILE_3D;
    if (strcmp(str, "BINARY_STAR_3D") == 0) return BINARY_STAR_3D;
    if (strcmp(str, "STELLAR_OCTANT_EXPLOSION_3D") == 0) return STELLAR_OCTANT_EXPLOSION_3D;
    
    fprintf(stderr, "Warning: Unknown problem type '%s', defaulting to KELVIN_HELMHOLTZ\n", str);
    return KELVIN_HELMHOLTZ;
}

// Parse boundary condition type from string
BCType parse_bc_type(const char *str) {
    if (strcmp(str, "PERIODIC") == 0) return BC_PERIODIC;
    if (strcmp(str, "REFLECTIVE") == 0) return BC_REFLECTIVE;
    if (strcmp(str, "OUTFLOW") == 0) return BC_OUTFLOW;
    
    fprintf(stderr, "Warning: Unknown BC type '%s', defaulting to OUTFLOW\n", str);
    return BC_OUTFLOW;
}

// Create and initialize input deck with default values
InputDeck* input_deck_create(const char *filename) {
    InputDeck *deck = (InputDeck*) malloc(sizeof(InputDeck));
    if (!deck) {
        fprintf(stderr, "Error: Failed to allocate memory for InputDeck\n");
        return NULL;
    }
    
    // Set default values
    deck->nx = 128;
    deck->ny = 128;
    deck->nz = 128;
    
    deck->xmin = -1.0;
    deck->xmax = 1.0;
    deck->ymin = -1.0;
    deck->ymax = 1.0;
    deck->zmin = -1.0;
    deck->zmax = 1.0;
    
    deck->t_end = 1.0;
    deck->output_dt = 0.1;
    deck->cfl = 0.4;
    
    deck->gamma = 5.0 / 3.0;
    deck->G_const = 1.0;
    
    deck->problem = KELVIN_HELMHOLTZ_3D;
    
    deck->bc_xmin = BC_PERIODIC;
    deck->bc_xmax = BC_PERIODIC;
    deck->bc_ymin = BC_PERIODIC;
    deck->bc_ymax = BC_PERIODIC;
    deck->bc_zmin = BC_PERIODIC;
    deck->bc_zmax = BC_PERIODIC;
    
    deck->G_M_PNS = 0.0;
    deck->L_nu = 0.0;
    deck->g_strength = 0.0;
    deck->L_heat = 0.0;
    
    strcpy(deck->output_dir, "output");
    
    deck->restart = false;
    strcpy(deck->restart_file, "");
    deck->restart_step = 0;
    deck->restart_time = 0.0;
    
    // Binary merger parameters
    strcpy(deck->star1_profile, "stellar_profile_1.txt");
    strcpy(deck->star2_profile, "stellar_profile_2.txt");
    deck->star1_x = -1.5;
    deck->star1_y = 0.0;
    deck->star1_z = 0.0;
    deck->star2_x = 1.5;
    deck->star2_y = 0.0;
    deck->star2_z = 0.0;
    deck->star1_vx = 0.0;
    deck->star1_vy = 0.5;
    deck->star1_vz = 0.0;
    deck->star2_vx = 0.0;
    deck->star2_vy = -0.5;
    deck->star2_vz = 0.0;
    deck->binary_separation = 3.0;
    deck->binary_orbital_velocity = 0.0;  // Will be calculated if auto_orbit is true
    deck->binary_auto_orbit = true;
    
    // Stellar octant explosion parameters
    strcpy(deck->stellar_profile, "stellar_profile.txt");
    deck->explosion_energy = 1.0;  // Energy in code units
    deck->explosion_radius = 0.2;  // Deposit energy below this radius
    strcpy(deck->explosion_type, "thermal");  // or "kinetic"
    
    // Perturbation parameters (default: no perturbations)
    deck->enable_perturbations = false;
    deck->perturbation_seed = 42;
    
    deck->perturb_velocity_amplitude = 0.01;  // 1% of sound speed
    deck->perturb_velocity_kmin = 1.0;        // Large-scale modes
    deck->perturb_velocity_kmax = 10.0;       // Small-scale modes
    deck->perturb_velocity_modes = 20;        // Number of modes to sum
    
    deck->perturb_density_amplitude = 0.01;   // 1% density variations
    deck->perturb_density_kmin = 1.0;
    deck->perturb_density_kmax = 10.0;
    deck->perturb_density_modes = 20;
    
    deck->enable_sph_harmonic_mode = false;
    deck->sph_harmonic_l = 1;                 // Dipole mode
    deck->sph_harmonic_m = 0;
    deck->sph_harmonic_amp_rho = 0.05;        // 5% density amplitude
    deck->sph_harmonic_amp_vel = 0.05;        // 5% velocity amplitude
    
    // Torus parameters (default: no torus)
    deck->enable_torus = false;
    deck->torus_major_radius = 1.2;           // Distance from origin to tube center
    deck->torus_minor_radius = 0.2;           // Tube/ring thickness
    deck->torus_density = 1.0;                // Density of torus material
    deck->torus_pressure = 0.1;               // Pressure of torus material
    deck->torus_rotation_velocity = 0.0;      // Azimuthal rotation (0 = static)
    
    // Now read from file if provided
    if (filename == NULL) {
        printf("No input file provided, using defaults.\n");
        return deck;
    }
    
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Warning: Could not open input file '%s', using defaults.\n", filename);
        return deck;
    }
    
    printf("Reading input deck from: %s\n", filename);
    
    char line[512];
    int line_num = 0;
    
    while (fgets(line, sizeof(line), fp)) {
        line_num++;
        
        // Skip empty lines and comments
        char *trimmed = trim_whitespace(line);
        if (strlen(trimmed) == 0 || trimmed[0] == '#') continue;
        
        // Parse key-value pairs
        char key[128], value[256];
        if (sscanf(trimmed, "%127[^=]=%255[^\n]", key, value) == 2) {
            // --- ADD THIS BLOCK ---
            // Find and strip comments from the value string
            char *comment_start = strchr(value, '#');
            if (comment_start != NULL) {
                *comment_start = '\0'; // Terminate the string at the comment
            }
            // --- END ADDED BLOCK ---
            char *key_clean = trim_whitespace(key);
            char *value_clean = trim_whitespace(value);
            
            // Grid parameters
            if (strcmp(key_clean, "nx") == 0) deck->nx = atoi(value_clean);
            else if (strcmp(key_clean, "ny") == 0) deck->ny = atoi(value_clean);
            else if (strcmp(key_clean, "nz") == 0) deck->nz = atoi(value_clean);
            
            else if (strcmp(key_clean, "xmin") == 0) deck->xmin = atof(value_clean);
            else if (strcmp(key_clean, "xmax") == 0) deck->xmax = atof(value_clean);
            else if (strcmp(key_clean, "ymin") == 0) deck->ymin = atof(value_clean);
            else if (strcmp(key_clean, "ymax") == 0) deck->ymax = atof(value_clean);
            else if (strcmp(key_clean, "zmin") == 0) deck->zmin = atof(value_clean);
            else if (strcmp(key_clean, "zmax") == 0) deck->zmax = atof(value_clean);
            
            // Time parameters
            else if (strcmp(key_clean, "t_end") == 0) deck->t_end = atof(value_clean);
            else if (strcmp(key_clean, "output_dt") == 0) deck->output_dt = atof(value_clean);
            else if (strcmp(key_clean, "cfl") == 0) deck->cfl = atof(value_clean);
            
            // Physics parameters
            else if (strcmp(key_clean, "gamma") == 0) deck->gamma = atof(value_clean);
            else if (strcmp(key_clean, "G_const") == 0) deck->G_const = atof(value_clean);
            
            // Problem type
            else if (strcmp(key_clean, "problem") == 0) {
                deck->problem = parse_problem_type(value_clean);
            }
            
            // Boundary conditions
            else if (strcmp(key_clean, "bc_xmin") == 0) deck->bc_xmin = parse_bc_type(value_clean);
            else if (strcmp(key_clean, "bc_xmax") == 0) deck->bc_xmax = parse_bc_type(value_clean);
            else if (strcmp(key_clean, "bc_ymin") == 0) deck->bc_ymin = parse_bc_type(value_clean);
            else if (strcmp(key_clean, "bc_ymax") == 0) deck->bc_ymax = parse_bc_type(value_clean);
            else if (strcmp(key_clean, "bc_zmin") == 0) deck->bc_zmin = parse_bc_type(value_clean);
            else if (strcmp(key_clean, "bc_zmax") == 0) deck->bc_zmax = parse_bc_type(value_clean);
            
            // Additional physics
            else if (strcmp(key_clean, "G_M_PNS") == 0) deck->G_M_PNS = atof(value_clean);
            else if (strcmp(key_clean, "L_nu") == 0) deck->L_nu = atof(value_clean);
            else if (strcmp(key_clean, "g_strength") == 0) deck->g_strength = atof(value_clean);
            else if (strcmp(key_clean, "L_heat") == 0) deck->L_heat = atof(value_clean);
            
            // Output control
            else if (strcmp(key_clean, "output_dir") == 0) {
                strncpy(deck->output_dir, value_clean, sizeof(deck->output_dir) - 1);
            }
            
            // Restart parameters
            else if (strcmp(key_clean, "restart") == 0) {
                deck->restart = (strcmp(value_clean, "true") == 0 || strcmp(value_clean, "1") == 0);
            }
            else if (strcmp(key_clean, "restart_file") == 0) {
                strncpy(deck->restart_file, value_clean, sizeof(deck->restart_file) - 1);
            }
            else if (strcmp(key_clean, "restart_step") == 0) deck->restart_step = atoi(value_clean);
            else if (strcmp(key_clean, "restart_time") == 0) deck->restart_time = atof(value_clean);
            
            // Binary merger parameters
            else if (strcmp(key_clean, "star1_profile") == 0) {
                strncpy(deck->star1_profile, value_clean, sizeof(deck->star1_profile) - 1);
            }
            else if (strcmp(key_clean, "star2_profile") == 0) {
                strncpy(deck->star2_profile, value_clean, sizeof(deck->star2_profile) - 1);
            }
            else if (strcmp(key_clean, "star1_x") == 0) deck->star1_x = atof(value_clean);
            else if (strcmp(key_clean, "star1_y") == 0) deck->star1_y = atof(value_clean);
            else if (strcmp(key_clean, "star1_z") == 0) deck->star1_z = atof(value_clean);
            else if (strcmp(key_clean, "star2_x") == 0) deck->star2_x = atof(value_clean);
            else if (strcmp(key_clean, "star2_y") == 0) deck->star2_y = atof(value_clean);
            else if (strcmp(key_clean, "star2_z") == 0) deck->star2_z = atof(value_clean);
            else if (strcmp(key_clean, "star1_vx") == 0) deck->star1_vx = atof(value_clean);
            else if (strcmp(key_clean, "star1_vy") == 0) deck->star1_vy = atof(value_clean);
            else if (strcmp(key_clean, "star1_vz") == 0) deck->star1_vz = atof(value_clean);
            else if (strcmp(key_clean, "star2_vx") == 0) deck->star2_vx = atof(value_clean);
            else if (strcmp(key_clean, "star2_vy") == 0) deck->star2_vy = atof(value_clean);
            else if (strcmp(key_clean, "star2_vz") == 0) deck->star2_vz = atof(value_clean);
            else if (strcmp(key_clean, "binary_separation") == 0) deck->binary_separation = atof(value_clean);
            else if (strcmp(key_clean, "binary_orbital_velocity") == 0) deck->binary_orbital_velocity = atof(value_clean);
            else if (strcmp(key_clean, "binary_auto_orbit") == 0) {
                deck->binary_auto_orbit = (strcmp(value_clean, "true") == 0 || strcmp(value_clean, "1") == 0);
            }
            
            // Stellar octant explosion parameters
            else if (strcmp(key_clean, "stellar_profile") == 0) {
                strncpy(deck->stellar_profile, value_clean, sizeof(deck->stellar_profile) - 1);
            }
            else if (strcmp(key_clean, "explosion_energy") == 0) deck->explosion_energy = atof(value_clean);
            else if (strcmp(key_clean, "explosion_radius") == 0) deck->explosion_radius = atof(value_clean);
            else if (strcmp(key_clean, "explosion_type") == 0) {
                strncpy(deck->explosion_type, value_clean, sizeof(deck->explosion_type) - 1);
            }
            
            // Perturbation parameters
            else if (strcmp(key_clean, "enable_perturbations") == 0) {
                deck->enable_perturbations = (strcmp(value_clean, "true") == 0 || strcmp(value_clean, "1") == 0);
            }
            else if (strcmp(key_clean, "perturbation_seed") == 0) deck->perturbation_seed = atoi(value_clean);
            
            else if (strcmp(key_clean, "perturb_velocity_amplitude") == 0) deck->perturb_velocity_amplitude = atof(value_clean);
            else if (strcmp(key_clean, "perturb_velocity_kmin") == 0) deck->perturb_velocity_kmin = atof(value_clean);
            else if (strcmp(key_clean, "perturb_velocity_kmax") == 0) deck->perturb_velocity_kmax = atof(value_clean);
            else if (strcmp(key_clean, "perturb_velocity_modes") == 0) deck->perturb_velocity_modes = atoi(value_clean);
            
            else if (strcmp(key_clean, "perturb_density_amplitude") == 0) deck->perturb_density_amplitude = atof(value_clean);
            else if (strcmp(key_clean, "perturb_density_kmin") == 0) deck->perturb_density_kmin = atof(value_clean);
            else if (strcmp(key_clean, "perturb_density_kmax") == 0) deck->perturb_density_kmax = atof(value_clean);
            else if (strcmp(key_clean, "perturb_density_modes") == 0) deck->perturb_density_modes = atoi(value_clean);
            
            else if (strcmp(key_clean, "enable_sph_harmonic_mode") == 0) {
                deck->enable_sph_harmonic_mode = (strcmp(value_clean, "true") == 0 || strcmp(value_clean, "1") == 0);
            }
            else if (strcmp(key_clean, "sph_harmonic_l") == 0) deck->sph_harmonic_l = atoi(value_clean);
            else if (strcmp(key_clean, "sph_harmonic_m") == 0) deck->sph_harmonic_m = atoi(value_clean);
            else if (strcmp(key_clean, "sph_harmonic_amp_rho") == 0) deck->sph_harmonic_amp_rho = atof(value_clean);
            else if (strcmp(key_clean, "sph_harmonic_amp_vel") == 0) deck->sph_harmonic_amp_vel = atof(value_clean);
            
            // Torus parameters
            else if (strcmp(key_clean, "enable_torus") == 0) {
                deck->enable_torus = (strcmp(value_clean, "true") == 0 || strcmp(value_clean, "1") == 0);
            }
            else if (strcmp(key_clean, "torus_major_radius") == 0) deck->torus_major_radius = atof(value_clean);
            else if (strcmp(key_clean, "torus_minor_radius") == 0) deck->torus_minor_radius = atof(value_clean);
            else if (strcmp(key_clean, "torus_density") == 0) deck->torus_density = atof(value_clean);
            else if (strcmp(key_clean, "torus_pressure") == 0) deck->torus_pressure = atof(value_clean);
            else if (strcmp(key_clean, "torus_rotation_velocity") == 0) deck->torus_rotation_velocity = atof(value_clean);
            
            else {
                fprintf(stderr, "Warning: Unknown parameter '%s' on line %d\n", key_clean, line_num);
            }
        }
    }
    
    fclose(fp);
    printf("Successfully loaded input deck.\n");
    
    return deck;
}

// Destroy input deck
void input_deck_destroy(InputDeck *deck) {
    if (deck) {
        free(deck);
    }
}

// Print input deck parameters
void input_deck_print(const InputDeck *deck) {
    printf("\n==============================================\n");
    printf("INPUT DECK PARAMETERS\n");
    printf("==============================================\n");
    
    printf("Grid:\n");
    printf("  nx = %d, ny = %d, nz = %d\n", deck->nx, deck->ny, deck->nz);
    printf("  x: [%.2f, %.2f]\n", deck->xmin, deck->xmax);
    printf("  y: [%.2f, %.2f]\n", deck->ymin, deck->ymax);
    printf("  z: [%.2f, %.2f]\n", deck->zmin, deck->zmax);
    
    printf("\nTime:\n");
    printf("  t_end = %.3f\n", deck->t_end);
    printf("  output_dt = %.3f\n", deck->output_dt);
    printf("  cfl = %.3f\n", deck->cfl);
    
    printf("\nPhysics:\n");
    printf("  gamma = %.3f\n", deck->gamma);
    printf("  G_const = %.3f\n", deck->G_const);
    
    printf("\nProblem: ");
    switch (deck->problem) {
        case KELVIN_HELMHOLTZ_3D: printf("KELVIN_HELMHOLTZ_3D\n"); break;
        case BLAST_WAVE_3D: printf("BLAST_WAVE_3D\n"); break;
        case TURBULENCE_3D: printf("TURBULENCE_3D\n"); break;
        case SEDOV_BLAST_3D: printf("SEDOV_BLAST_3D\n"); break;
        case STELLAR_OCTANT_EXPLOSION_3D: printf("STELLAR_OCTANT_EXPLOSION_3D\n"); break;
        default: printf("OTHER\n");
    }
    
    printf("\nBoundary Conditions:\n");
    printf("  x: %s / %s\n", 
           deck->bc_xmin == BC_PERIODIC ? "PERIODIC" : (deck->bc_xmin == BC_REFLECTIVE ? "REFLECTIVE" : "OUTFLOW"),
           deck->bc_xmax == BC_PERIODIC ? "PERIODIC" : (deck->bc_xmax == BC_REFLECTIVE ? "REFLECTIVE" : "OUTFLOW"));
    printf("  y: %s / %s\n", 
           deck->bc_ymin == BC_PERIODIC ? "PERIODIC" : (deck->bc_ymin == BC_REFLECTIVE ? "REFLECTIVE" : "OUTFLOW"),
           deck->bc_ymax == BC_PERIODIC ? "PERIODIC" : (deck->bc_ymax == BC_REFLECTIVE ? "REFLECTIVE" : "OUTFLOW"));
    printf("  z: %s / %s\n", 
           deck->bc_zmin == BC_PERIODIC ? "PERIODIC" : (deck->bc_zmin == BC_REFLECTIVE ? "REFLECTIVE" : "OUTFLOW"),
           deck->bc_zmax == BC_PERIODIC ? "PERIODIC" : (deck->bc_zmax == BC_REFLECTIVE ? "REFLECTIVE" : "OUTFLOW"));
    
    printf("\nOutput:\n");
    printf("  Directory: %s\n", deck->output_dir);
    
    if (deck->restart) {
        printf("\nRestart:\n");
        printf("  File: %s\n", deck->restart_file);
        printf("  Step: %d\n", deck->restart_step);
        printf("  Time: %.3f\n", deck->restart_time);
    }
    
    if (deck->problem == BINARY_STAR_3D) {
        printf("\nBinary Merger Parameters:\n");
        printf("  Star 1 profile: %s\n", deck->star1_profile);
        printf("  Star 2 profile: %s\n", deck->star2_profile);
        printf("  Star 1 position: (%.3f, %.3f, %.3f)\n", deck->star1_x, deck->star1_y, deck->star1_z);
        printf("  Star 2 position: (%.3f, %.3f, %.3f)\n", deck->star2_x, deck->star2_y, deck->star2_z);
        printf("  Star 1 velocity: (%.3f, %.3f, %.3f)\n", deck->star1_vx, deck->star1_vy, deck->star1_vz);
        printf("  Star 2 velocity: (%.3f, %.3f, %.3f)\n", deck->star2_vx, deck->star2_vy, deck->star2_vz);
        printf("  Auto-calculate orbit: %s\n", deck->binary_auto_orbit ? "yes" : "no");
        if (deck->binary_auto_orbit) {
            printf("  Separation: %.3f\n", deck->binary_separation);
        }
    }
    
    if (deck->problem == STELLAR_OCTANT_EXPLOSION_3D) {
        printf("\nStellar Octant Explosion Parameters:\n");
        printf("  Profile file: %s\n", deck->stellar_profile);
        printf("  Explosion energy: %.3e (code units)\n", deck->explosion_energy);
        printf("  Explosion radius: %.3f\n", deck->explosion_radius);
        printf("  Explosion type: %s\n", deck->explosion_type);
        
        if (deck->enable_torus) {
            printf("\nTorus Parameters:\n");
            printf("  Enabled: YES\n");
            printf("  Major radius R: %.3f\n", deck->torus_major_radius);
            printf("  Minor radius r: %.3f\n", deck->torus_minor_radius);
            printf("  Density: %.3e\n", deck->torus_density);
            printf("  Pressure: %.3e\n", deck->torus_pressure);
            printf("  Rotation velocity: %.3f\n", deck->torus_rotation_velocity);
        }
    }
    
    printf("==============================================\n\n");
}
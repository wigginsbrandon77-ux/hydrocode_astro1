#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "input_deck.h"
#include "hydro_solver_3d.h"

int main(int argc, char *argv[]) {
    // Check command line arguments
    if (argc < 2) {
        printf("Usage: %s <input_file>\n", argv[0]);
        printf("Example: %s input_kh_3d.txt\n", argv[0]);
        return 1;
    }
    
    const char *input_filename = argv[1];
    
    printf("\n");
    printf("==============================================\n");
    printf("  3D HYDRODYNAMICS SOLVER\n");
    printf("==============================================\n");
    printf("Reading configuration from: %s\n\n", input_filename);
    
    // Load input deck
    InputDeck *deck = input_deck_create(input_filename);
    if (!deck) {
        fprintf(stderr, "Failed to load input deck.\n");
        return 1;
    }
    
    // Print parameters
    input_deck_print(deck);
    
    // Create output directory
    char mkdir_cmd[512];
    sprintf(mkdir_cmd, "mkdir -p %s", deck->output_dir);
    system(mkdir_cmd);
    
    // Create solver from input deck
    printf("Creating 3D solver...\n");
    GodunovHLLC3D *solver = solver3d_create_from_deck(deck);
    if (!solver) {
        fprintf(stderr, "Failed to create solver.\n");
        input_deck_destroy(deck);
        return 1;
    }
    
    printf("Solver created successfully.\n");
    printf("Total cells: %d x %d x %d = %d\n", 
           solver->nx, solver->ny, solver->nz,
           solver->nx * solver->ny * solver->nz);
    printf("Cell size: dx=%.4f, dy=%.4f, dz=%.4f\n\n",
           solver->dx, solver->dy, solver->dz);
    
    // Check if this is a restart run
    if (deck->restart) {
        printf("==============================================\n");
        printf("RESTART MODE\n");
        printf("==============================================\n");
        printf("Restart file: %s\n", deck->restart_file);
        printf("Restart time: %.5f\n", deck->restart_time);
        printf("Starting output count: %d\n\n", deck->restart_step);
        
        // Load state from CSV restart file
        if (load_from_csv_restart(solver, deck->restart_file) != 0) {
            fprintf(stderr, "Failed to load restart file. Exiting.\n");
            solver3d_destroy(solver);
            input_deck_destroy(deck);
            return 1;
        }
        
        printf("Restart data loaded successfully.\n\n");
        
        // Run simulation from restart point
        printf("==============================================\n");
        printf("CONTINUING SIMULATION FROM RESTART\n");
        printf("==============================================\n\n");
        
        run_simulation_3d_restart(solver, deck->restart_time, deck->t_end, 
                                   deck->output_dt, deck->restart_step);
        
    } else {
        // Normal run from initial conditions
        printf("Setting initial conditions...\n");
        set_initial_conditions_3d_from_deck(solver, deck);
        printf("Initial conditions set.\n\n");
        
        // Run simulation
        printf("==============================================\n");
        printf("STARTING SIMULATION\n");
        printf("==============================================\n\n");
        
        run_simulation_3d(solver, deck->t_end, deck->output_dt);
    }
    
    printf("\n==============================================\n");
    printf("SIMULATION COMPLETE\n");
    printf("==============================================\n");
    printf("Output files saved to: %s/\n", deck->output_dir);
    printf("\nOutput formats:\n");
    printf("  - CSV files: output_NNNNNN.csv (for analysis/restart)\n");
    printf("  - VTK files: output_NNNNNN.vtk (for ParaView/VisIt)\n");
    printf("\nRestart capability:\n");
    printf("  - To restart from any output, set in input file:\n");
    printf("      restart = true\n");
    printf("      restart_file = output/output_NNNNNN.csv\n");
    printf("      restart_time = <time_at_that_step>\n");
    printf("      restart_step = NNNNNN\n");
    printf("\nVisualization:\n");
    printf("  - Use ParaView to open VTK files\n");
    printf("  - Use Python/matplotlib for CSV files\n");
    printf("  - Use VisIt for advanced visualization\n");
    printf("==============================================\n\n");
    
    // Clean up
    solver3d_destroy(solver);
    input_deck_destroy(deck);
    
    printf("Cleanup complete. Exiting.\n");
    
    return 0;
}
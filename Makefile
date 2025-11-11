# Makefile for 3D Hydrodynamics Solver
# =====================================

CC = gcc
CFLAGS = -O3 -Wall -Wextra -march=native -std=c11
LDFLAGS = -lm

# Source files
SOURCES_3D = main_3d.c hydro_solver_3d.c input_deck.c gravity_solver_3d.c
OBJECTS_3D = $(SOURCES_3D:.c=.o)

# Target executable
TARGET_3D = hydro3d

# Phony targets
.PHONY: all clean run_kh run_blast help

# Default target
all: $(TARGET_3D)

# Build 3D solver
$(TARGET_3D): $(OBJECTS_3D)
	@echo "Linking $(TARGET_3D)..."
	$(CC) $(OBJECTS_3D) -o $(TARGET_3D) $(LDFLAGS)
	@echo "Build complete: $(TARGET_3D)"

# Compile object files
%.o: %.c
	@echo "Compiling $<..."
	$(CC) $(CFLAGS) -c $< -o $@

# Run Kelvin-Helmholtz simulation
run_kh: $(TARGET_3D)
	@echo "Running Kelvin-Helmholtz simulation..."
	./$(TARGET_3D) input_kh_3d.txt

# Run blast wave simulation
run_blast: $(TARGET_3D)
	@echo "Running blast wave simulation..."
	./$(TARGET_3D) input_blast_3d.txt

# Clean build files
clean:
	@echo "Cleaning build files..."
	rm -f $(OBJECTS_3D) $(TARGET_3D)
	@echo "Clean complete."

# Clean output files
clean_output:
	@echo "Cleaning output files..."
	rm -rf output output_*
	@echo "Output clean complete."

# Clean everything
clean_all: clean clean_output

# Display help
help:
	@echo "3D Hydrodynamics Solver - Make targets:"
	@echo "  all          - Build the 3D solver (default)"
	@echo "  clean        - Remove build files"
	@echo "  clean_output - Remove output files"
	@echo "  clean_all    - Remove build and output files"
	@echo "  run_kh       - Build and run Kelvin-Helmholtz simulation"
	@echo "  run_blast    - Build and run blast wave simulation"
	@echo "  help         - Display this help message"
	@echo ""
	@echo "Usage:"
	@echo "  make               # Build the solver"
	@echo "  make run_kh        # Run KH simulation"
	@echo "  ./hydro3d input.txt  # Run with custom input file"

# Dependencies (auto-generated would be better, but this works for small projects)
main_3d.o: main_3d.c input_deck.h hydro_solver_3d.h
hydro_solver_3d.o: hydro_solver_3d.c hydro_solver_3d.h input_deck.h gravity_solver_3d.h
input_deck.o: input_deck.c input_deck.h
gravity_solver_3d.o: gravity_solver_3d.c gravity_solver_3d.h

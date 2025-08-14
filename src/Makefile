#===============================================================================
# MAKEFILE FOR EARTHQUAKE SIMULATION
# 
# This Makefile provides an alternative compilation method for the
# improved earthquake simulation code.
#
# Author: Makefile for improved earthquake simulation code
# Last Modified: Current compilation setup
#===============================================================================

# Compiler and flags
FC = mpif90
FFLAGS = -O2 -g -Wall -Wextra -std=f2008
LDFLAGS = 

# Target executable
TARGET = earthquake_simulation

# Build directory
BUILD_DIR = build

# Source files (in dependency order)
MODULES = physical_constants \
          simulation_parameters \
          physical_variables \
          io_parameters \
          error_handling \
          time_integration \
          physics_equations \
          mpi_utilities \
          io_handling

MAIN_PROGRAM = earthquake_simulation_improved

# Object files
OBJS = $(addprefix $(BUILD_DIR)/, $(addsuffix .o, $(MODULES))) \
       $(BUILD_DIR)/$(MAIN_PROGRAM).o

# Module files
MODS = $(addprefix $(BUILD_DIR)/, $(addsuffix .mod, $(MODULES)))

# Default target
all: $(BUILD_DIR)/$(TARGET)

# Create build directory
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Compile modules in dependency order
$(BUILD_DIR)/physical_constants.o: physical_constants.f90 $(BUILD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD_DIR)/simulation_parameters.o: simulation_parameters.f90 $(BUILD_DIR)/physical_constants.o
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD_DIR)/physical_variables.o: physical_variables.f90 $(BUILD_DIR)/physical_constants.o
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD_DIR)/io_parameters.o: io_parameters.f90 $(BUILD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD_DIR)/error_handling.o: error_handling.f90 $(BUILD_DIR)/physical_constants.o
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD_DIR)/time_integration.o: time_integration.f90 $(BUILD_DIR)/physical_constants.o
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD_DIR)/physics_equations.o: physics_equations.f90 $(BUILD_DIR)/physical_constants.o $(BUILD_DIR)/physical_variables.o $(BUILD_DIR)/simulation_parameters.o
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD_DIR)/mpi_utilities.o: mpi_utilities.f90 $(BUILD_DIR)/physical_constants.o
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD_DIR)/io_handling.o: io_handling.f90 $(BUILD_DIR)/physical_constants.o $(BUILD_DIR)/io_parameters.o
	$(FC) $(FFLAGS) -c $< -o $@

# Compile main program (depends on all modules)
$(BUILD_DIR)/$(MAIN_PROGRAM).o: $(MAIN_PROGRAM).f90 $(OBJS)
	$(FC) $(FFLAGS) -c $< -o $@

# Link executable
$(BUILD_DIR)/$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

# Clean build artifacts
clean:
	rm -rf $(BUILD_DIR)

# Clean and rebuild
rebuild: clean all

# Install (create symlink in parent directory)
install: $(BUILD_DIR)/$(TARGET)
	ln -sf $(BUILD_DIR)/$(TARGET) ../$(TARGET)

# Uninstall (remove symlink)
uninstall:
	rm -f ../$(TARGET)

# Show help
help:
	@echo "Available targets:"
	@echo "  all       - Build the executable (default)"
	@echo "  clean     - Remove build artifacts"
	@echo "  rebuild   - Clean and rebuild"
	@echo "  install   - Create symlink in parent directory"
	@echo "  uninstall - Remove symlink from parent directory"
	@echo "  help      - Show this help message"
	@echo ""
	@echo "Variables:"
	@echo "  FC        - Fortran compiler (default: mpif90)"
	@echo "  FFLAGS    - Compiler flags"
	@echo "  TARGET    - Executable name (default: earthquake_simulation)"
	@echo "  BUILD_DIR - Build directory (default: build)"

# Phony targets
.PHONY: all clean rebuild install uninstall help

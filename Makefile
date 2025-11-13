FC = gfortran
FFLAGS = -O2 -Wall -Wextra -fcheck=all

SRC_DIR = src
OBJ_DIR = obj

SOURCES = $(SRC_DIR)/kinds.f90 \
          $(SRC_DIR)/params.f90 \
          $(SRC_DIR)/grid.f90 \
          $(SRC_DIR)/fields.f90 \
          $(SRC_DIR)/spectral_transforms.f90 \
          $(SRC_DIR)/torpol.f90 \
          $(SRC_DIR)/operators.f90 \
          $(SRC_DIR)/poisson_radial.f90 \
          $(SRC_DIR)/time_integration.f90 \
          $(SRC_DIR)/main.f90

OBJECTS = $(SOURCES:$(SRC_DIR)/%.f90=$(OBJ_DIR)/%.o)

TARGET = precess_sphere

all: $(TARGET)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	mkdir -p $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(TARGET): $(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) -o $@

clean:
	rm -rf $(OBJ_DIR) $(TARGET)

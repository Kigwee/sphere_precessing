############################
# COMPILATION SETTINGS
############################

FC      = gfortran
FFLAGS  = -O3 -march=native -ffast-math -funroll-loops -fopenmp -Wall

SHTNS_DIR = external/SHTns
FFTW_DIR  = /usr/lib/x86_64-linux-gnu    # change si besoin

INCLUDES = -I$(SHTNS_DIR) -I$(FFTW_DIR)/include
LIBS     = -L$(SHTNS_DIR) -lshtns \
           -L$(FFTW_DIR)/lib -lfftw3 -lfftw3_omp \
           -llapack -lblas

SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = build

SOURCES = $(wildcard $(SRC_DIR)/*.f90)
OBJECTS = $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(SOURCES))

TARGET  = $(BIN_DIR)/dns_sphere

############################
# RULES
############################

all: $(TARGET)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	mkdir -p $(OBJ_DIR)
	$(FC) $(FFLAGS) $(INCLUDES) -c $< -o $@

$(TARGET): $(OBJECTS)
	mkdir -p $(BIN_DIR)
	$(FC) $(FFLAGS) $(OBJECTS) -o $(TARGET) $(LIBS)

clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

############################
# GPU version (optionnel)
############################
gpu:
	make LIBS="$(LIBS) -lcufft" FFLAGS="$(FFLAGS) -DSHTNS_GPU=1"



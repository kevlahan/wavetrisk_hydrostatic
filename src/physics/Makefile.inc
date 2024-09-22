SRC += init_physics_module.f90 physics_simple.f90 phys_processing.f90 physics_Held_Suarez.f90 physics_trend.f90
vpath %.f90 src/physics
all: phys_package $(PREFIX)/$(BIN_DIR)/$(TEST_CASE) #compiles phys package first


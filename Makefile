

# paralel:
# 	https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.1.tar.gz
# 	https://computation.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods/download/hypre-2.10.0b.tar.gz
# 	http://glaros.dtc.umn.edu/gkhome/metis/metis/overview

ROOT_DIR=$(realpath .)

BUILD_DIR = $(ROOT_DIR)/build
INSTALL_DIR = $(ROOT_DIR)/install

# Release
# Debug
CMAKE = cmake
CMAKE_ARGS = \
	"-DCMAKE_INSTALL_PREFIX=$(INSTALL_DIR)" \
	"-DCMAKE_BUILD_TYPE=Release" \
	"-DCXX_STANDARD=17"

all: $(BUILD_DIR)/make.stamp

$(BUILD_DIR)/CMakeCache.txt: $(ROOT_DIR)/CMakeLists.txt $(MAKEFILE_LIST)
	@echo "[CMAKE] $(BUILD_DIR)"
	@mkdir -p $(dir $@)
	cd $(dir $@); $(CMAKE) $(CMAKE_ARGS) $(dir $<)

$(BUILD_DIR)/make.stamp: $(BUILD_DIR)/CMakeCache.txt
	@echo "[BUILD] $(BUILD_DIR)"
	@mkdir -p $(INSTALL_DIR)
	$(MAKE) -C $(dir $<)
	@touch $@

clean:
	@echo "[clean] $(BUILD_DIR) $(INSTALL_DIR)"
ifeq ($(BUILD_DIR)/Makefile,$(wildcard $(BUILD_DIR)/Makefile))
	$(MAKE) -C $(BUILD_DIR) clean || echo FAIL
endif
	rm -rf $(BUILD_DIR) $(INSTALL_DIR)

# ubuntu pkgs for glvis
# sudo apt-get install libglu1-mesa libpng-dev libfontconfig1-dev

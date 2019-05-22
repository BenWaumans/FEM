ROOT_DIR=$(realpath .)

BUILD_DIR ?= $(ROOT_DIR)/build
INSTALL_DIR ?= $(ROOT_DIR)/install

GCC ?= /usr/bin/x86_64-linux-gnu-
GCC-SUFFIX ?= -6

CC = $(GCC)gcc$(GCC-SUFFIX)
CXX = $(GCC)g++$(GCC-SUFFIX)
AR = $(GCC)gcc-ar$(GCC-SUFFIX)
RANLIB = $(GCC)gcc-ranlib$(GCC-SUFFIX)

USE_CUDA ?= 0
CUDA_ARCH ?= sm_61
USE_MPI ?= 0
BUILD_MPI ?= 0
DEBUG ?= 0

CMAKE = cmake
CMAKE_ARGS = \
	"-DCMAKE_INSTALL_PREFIX=$(INSTALL_DIR)" \
	"-DCMAKE_C_COMPILER=$(CC)" \
	"-DCMAKE_C_COMPILER_AR=$(AR)" \
	"-DCMAKE_C_COMPILER_RANLIB=$(RANLIB)" \
	"-DCMAKE_CXX_COMPILER=$(CXX)" \
	"-DCMAKE_CXX_COMPILER_AR=$(AR)" \
	"-DCMAKE_CXX_COMPILER_RANLIB=$(RANLIB)" \
	"-DCXX_STANDARD=17"

ifeq ($(DEBUG), 1)
	CMAKE_ARGS += "-DCMAKE_BUILD_TYPE=Debug"
else
	CMAKE_ARGS += "-DCMAKE_BUILD_TYPE=Release"
endif
ifeq ($(USE_CUDA), 1)
	CMAKE_ARGS += "-DUSE_CUDA=ON" "-DCUDA_ARCH=$(CUDA_ARCH)"
endif
ifeq ($(USE_MPI), 1)
	CMAKE_ARGS += "-DUSE_MPI=ON"
ifeq ($(BUILD_MPI), 1)
	CMAKE_ARGS += "-DBUILD_MPI=ON"
endif
endif

all: $(INSTALL_DIR)/fem_solver

$(BUILD_DIR)/CMakeCache.txt: $(ROOT_DIR)/CMakeLists.txt $(MAKEFILE_LIST)
	@echo "[CMAKE] $(BUILD_DIR)"
	@mkdir -p $(dir $@)
	cd $(dir $@); $(CMAKE) $(CMAKE_ARGS) $(dir $<)

$(BUILD_DIR)/make.stamp: $(BUILD_DIR)/CMakeCache.txt
	@echo "[BUILD] $(BUILD_DIR)"
	@mkdir -p $(INSTALL_DIR)
	$(MAKE) -C $(dir $<)
	@touch $@

$(INSTALL_DIR)/fem_solver: $(BUILD_DIR)/make.stamp
	@echo "[INSTALL] $(INSTALL_DIR)"
	@mkdir -p $(INSTALL_DIR)
	$(MAKE) -C $(dir $<) install

clean:
	@echo "[clean] $(BUILD_DIR) $(INSTALL_DIR)"
ifeq ($(BUILD_DIR)/Makefile,$(wildcard $(BUILD_DIR)/Makefile))
	$(MAKE) -C $(BUILD_DIR) clean || echo FAIL
endif
	rm -rf $(BUILD_DIR) $(INSTALL_DIR)

# ubuntu pkgs for glvis
# sudo apt-get install libglu1-mesa-dev libpng-dev libfontconfig1-dev

# Copyright (c) 2010, Lawrence Livermore National Security, LLC. Produced at the
# Lawrence Livermore National Laboratory. LLNL-CODE-443211. All Rights reserved.
# See file COPYRIGHT for details.
#
# This file is part of the MFEM library. For more information and source code
# availability see http://mfem.org.
#
# MFEM is free software; you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License (as published by the Free
# Software Foundation) version 2.1 dated February 1999.

# Defines the following variables:
#   - METIS_FOUND
#   - METIS_LIBRARIES
#   - METIS_INCLUDE_DIRS

include(FemSolverCommands)
fem_solver_find_package(METIS METIS METIS_DIR "include;Lib" "metis.h"
  "lib" "metis;metis4;metis5"
  "Paths to headers required by METIS." "Libraries required by METIS.")

# This Makefile snippet is part of the SCHISM-ESMF interface
#
# @copyright (C) 2021 Helmholtz-Zentrum Hereon
# @copyright (C) 2020-2021 Helmholtz-Zentrum Geesthacht
# 
# @author Carsten Lemmen carsten.lemmen@hereon.de
#
# @license Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# 		http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# A NUOPC compliant Makefile snippet that defines the six variables ESMF_DEP_*
# and is assembled in a driver Makefile. 

ifneq ($(origin SCHISM_CAP_INSTALL_DIR), environment)
#$(error SCHISM_CAP_INSTALL_DIR  has to be set in environment.)
# temporary fix for NEMS system
SCHISM_CAP_INSTALL_DIR=$(NEMSDIR)/../SCHISM_INSTALL
else
SCHISM_CAP_INSTALL_DIR:= $(shell readlink --canonicalize ${SCHISM_CAP_INSTALL_DIR})
endif

#LIBS=-lhydro -lcore -lparmetis -lmetis

# 1. ESMF_DEP_FRONT - The name of the Fortran module to be used in a USE statement, or (if it ends in ".h") the name of the header file to be used in an #include statement, or (if it ends in ".so") the name of the shared object to be loaded at run-time.

ESMF_DEP_FRONT=schism_cmi_nuopc

# 2. ESMF_DEP_INCPATH - The include path to find module or header files during compilation. Must be specified as absolute path.
ESMF_DEP_INCPATH=$(SCHISM_CAP_INSTALL_DIR)

# 3. ESMF_DEP_CMPL_OBJS - Object files that need to be considered as compile dependencies. Must be specified with absolute path.
ESMF_DEP_CMPL_OBJS=

# 4. ESMF_DEP_LINK_OBJS - Object files that need to be considered as link dependencies. Must be specified with absolute path.
ESMF_DEP_LINK_OBJS=\
$(SCHISM_CAP_INSTALL_DIR)/libschism_cap.a \
$(SCHISM_CAP_INSTALL_DIR)/libparmetis.a \
$(SCHISM_CAP_INSTALL_DIR)/libmetis.a \
$(SCHISM_CAP_INSTALL_DIR)/libhydro.a \
$(SCHISM_CAP_INSTALL_DIR)/libcore.a

# 5. ESMF_DEP_SHRD_PATH - The path to find shared libraries during link-time (and during run-time unless over- ridden by LD_LIBRARY_PATH). Must be specified as absolute path.
ESMF_DEP_SHRD_PATH=$(SCHISM_CAP_INSTALL_DIR)

# 6. ESMF_DEP_SHRD_LIBS - Shared libraries that need to be specified during link-time, and must be available during run-time. Must be specified with absolute path.
ESMF_DEP_SHRD_LIBS=

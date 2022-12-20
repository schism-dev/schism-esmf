<--
# SPDX-FileCopyrightText: 2021-2022 Helmholtz-Zentrum Hereon
# SPDX-FileCopyrightText: 2018-2021 Helmholtz-Zentrum Geesthacht
# SPDX-License-Identifier: CC0-1.0
# SPDX-FileContributor: Carsten Lemmen <carsten.lemmen@hereon.de
-->
This directory contains the ESMF/NUOPC caps for the SCHISM model, as well as supporting infrastructure and interfaces.

The ESMF cap `schism_esmf_cap` relies on utilities in `schism_esmf_utils` and the basic model interface `schism_bmi`.
The NUOPC cap `schism_nuopc_cap` relies on utilities in `schism_nuop_utils`and `schism_esmf_utils`and the `schism_bmi`.

For the NUOPC model instance in `schism_nuopc_cap` there are also corresponding Makefile snippets for standardized
inclusion in NUOPC compliant model systems.


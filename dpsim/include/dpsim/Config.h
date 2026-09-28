/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#pragma once

// Version
#define DPSIM_VERSION "1.1.1"
#define DPSIM_SHORT_VERSION ""
#define DPSIM_RELEASE "1.feature/ssn_dc_hvdc_release.20260828git3853e33"

#define DPSIM_MAJOR_VERSION 1
#define DPSIM_MINOR_VERSION 1
#define DPSIM_PATCH_VERSION 1

// Features
#define WITH_RT
#define WITH_SPARSE
#define WITH_VILLAS
#define WITH_CIM
/* #undef WITH_PYBIND */
#define WITH_SUNDIALS
#define WITH_OPENMP
/* #undef WITH_CUDA */
/* #undef WITH_CUDA_SPARSE */
/* #undef WITH_MAGMA */
#define WITH_KLU
#define WITH_MNASOLVERPLUGIN
#define WITH_JSON
#define CGMES_BUILD

#define HAVE_GETOPT
#define HAVE_TIMERFD

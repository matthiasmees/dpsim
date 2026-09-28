// SPDX-FileCopyrightText: 2025 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#pragma once

#ifndef DPSIM_VILLAS_CONFIG_H
#define DPSIM_VILLAS_CONFIG_H

#define DPSIM_VILLAS_VERSION "1.2.1"

// Features
#define WITH_RT
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

#endif

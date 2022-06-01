/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief unit test utilities only needed for NVIDIA.

#ifndef NVIDIA_INDEX_BIN_COMMON_UNIT_TEST_UTILITY_H
#define NVIDIA_INDEX_BIN_COMMON_UNIT_TEST_UTILITY_H

#include <cassert>
#include <cstdio>
#include <set>
#include <stdlib.h>
#include <string>

#include "common_utility.h"
#include "windows_utilities.h"

namespace nv {
namespace index_common {

//----------------------------------------------------------------------
/// unset NVIDIA unit test environment (NVIDIA unit test specific)
inline void unset_nv_unit_test_env()
{
    // unset the disable GPU env (to enable GPU)
    char * p_unit_test_env = getenv("MI_DISABLE_GPU_FOR_UNIT_TEST");
    if (p_unit_test_env != 0)
    {
        mi::Sint32 const is_success = unsetenv("MI_DISABLE_GPU_FOR_UNIT_TEST");
        assert(is_success == 0);
        no_unused_variable_warning_please(is_success);
    }
    assert(getenv("MI_DISABLE_GPU_FOR_UNIT_TEST") == 0);
}

//----------------------------------------------------------------------
/// test utility: is this host on blacklist? (NVIDIA unit test specific)
/// The test host black list has somewhat problematic test host names.
///
/// \param[in] hostname test hostname
/// \return true this is on the abort black list
inline bool unit_test_is_host_on_abort_blacklist(std::string const & hostname)
{
    std::set< std::string > known_problem_host;
    // known_problem_host.insert("trayp");
    known_problem_host.insert("plebix"); // FX3800 is not a supported card.

    // No {Tesla, Fermi} support anymore 
    known_problem_host.insert("traya");
    known_problem_host.insert("trayb");
    known_problem_host.insert("trayc");
    known_problem_host.insert("trayd");
    known_problem_host.insert("traye");
    known_problem_host.insert("trayf");
    known_problem_host.insert("trayj");
    known_problem_host.insert("trayk");
    known_problem_host.insert("trayl");
    known_problem_host.insert("trayn");
    known_problem_host.insert("trayo");
    known_problem_host.insert("trayp");
    known_problem_host.insert("trayr");
    known_problem_host.insert("trays");
    known_problem_host.insert("trayt");
    known_problem_host.insert("trayu");

    if (known_problem_host.find(hostname) != known_problem_host.end())
    {
        return true;
    }
    return false;
}

//----------------------------------------------------------------------
/// set up and check the test environment
///
/// For instance, some internal machines can not perform the test.
///
/// \return true when the test environment is suitable for the test.
inline bool set_and_check_nv_unit_test_env()
{
    // to set up NVIDIA test env
    nv::index_common::unset_nv_unit_test_env();

    // printout host_name.
    std::string const host_name = nv::index_common::get_host_name();
    // INFO_LOG << "host_name: [" << host_name << "]";

    // check the black list of the test host
    if (nv::index_common::unit_test_is_host_on_abort_blacklist(host_name))
    {
        fprintf(stderr, "warn: host [%s] cannot perform the unit test. abort.",
                host_name.c_str());
        return false;
    }
    return true;
}

//----------------------------------------------------------------------
}} // namespace nv::index_common
#endif // NVIDIA_INDEX_BIN_COMMON_UNIT_TEST_UTILITY_H

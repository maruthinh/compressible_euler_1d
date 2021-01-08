cmake_minimum_required(VERSION 3.12)
project(src)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(src
        basic_functions.h
        boundary_conditions.cpp
        BoundaryConditions.h
        flux_convar_definitions.cpp
        flux_globalmax.cpp
        flux_kfds.cpp
        flux_kfds2nd_movers.cpp
        flux_kfds2nd_order.cpp
        flux_kfds_movers.cpp
        flux_llf.cpp
        flux_movers.cpp
        flux_movers_cp.cpp
        flux_movers_h.cpp
        flux_movers_le.cpp
        global_declarations.h
        grid_compute.cpp
        initial_conditions.cpp
        initial_conditions_oscil.cpp
        initial_conditions_two.cpp
        InitialConditions.h
        initialize_variables.cpp
        kfds2nd_movers.cpp
        lr_state.cpp
        main.cpp
        MoversDissipation.h
        read_input.cpp
        read_test_case.cpp
        read_test_case_two.cpp
        solver.cpp
        solver_explicit.cpp
        time_step.cpp
        write_solution.cpp)

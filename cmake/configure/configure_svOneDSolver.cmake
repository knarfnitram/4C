#message(STATUS "Fetch content for svOneDSolver")

#find_package(OneDSolver_lib HINTS "/home/a11bmafr/software/SimVascular/svOneDSolver/cmake-build-debug" REQUIRED)

#find_package(OneDSolver_lib HINTS

#message(${BACI_OneDSolver_ROOT})

#if(NOT OneDSolver_lib_found)
#    message(
#            FATAL_ERROR
#            "OneDSolverv2 could not be found. "
#    )
#endif()
#add_subdirectory("/home/a11bmafr/software/SimVascular/install/")
add_library(onedsolver UNKNOWN IMPORTED)
set_target_properties(
  onedsolver
  PROPERTIES IMPORTED_LOCATION
             "/home/a11bmafr/software/SimVascular/install/usr/local/lib/libOneDSolver_core.a"
             INTERFACE_INCLUDE_DIRECTORIES
             "/home/a11bmafr/software/SimVascular/install/include/onedsolver"
  )

#if(NOT TARGET OneDSolver_lib)
#fetchcontent_declare(
#        OneDSolver
#        GIT_REPOSITORY git@github.com:knarfnitram/svOneDSolver.git
#        GIT_TAG d30d825f5259a72782de9fd0318cef65fa8a2955
#)
#fetchcontent_makeavailable(OneDSolver)
#baci_add_dependency(baci_all_enabled_external_dependencies OneDSolver)
#endif ()
#find_package(Qhull REQUIRED)
#baci_add_dependency(baci_all_enabled_external_dependencies OneDSolver_lib)

message(STATUS "onedsolver include directory: ${ONEDSOLVER_INCLUDE_DIR}")
message(STATUS "onedsolver library directory: ${ONEDSOLVER_LIBRARY}")
#if(ONEDSOLVER_FOUND)
#    message(STATUS "onedsolver found")
#endif()
baci_add_dependency(baci_all_enabled_external_dependencies onedsolver)

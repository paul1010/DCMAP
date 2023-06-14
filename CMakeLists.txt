cmake_minimum_required(VERSION 3.15)

# set the project name and version
project(dcmap VERSION 1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# configure a header file to pass some of the CMake settings
# to the source code
# configure_file(TutorialConfig.h.in TutorialConfig.h)

# add the MathFunctions library
# add_subdirectory(Armadillo)
# list(APPEND EXTRA_LIBS Armadillo)
# list(APPEND EXTRA_INCLUDES "${PROJECT_SOURCE_DIR}/Armadillo/include")

# add Armadillo library
# find_package(Armadillo 11.4.3 REQUIRED)
# include_directories(${ARMADILLO_INCLUDE_DIRS})

# add openBLAS .lib and .dll
add_library(libopenblas SHARED IMPORTED)
set_property(TARGET libopenblas PROPERTY
			IMPORTED_LOCATION "${PROJECT_SOURCE_DIR}/libs/libopenblas.dll")
			set_property(TARGET libopenblas PROPERTY
			IMPORTED_IMPLIB "${PROJECT_SOURCE_DIR}/libs/libopenblas.lib")


# add the executable
add_executable(dcmap dcmap.cpp)

target_link_libraries(dcmap ${ARMADILLO_LIBRARIES}
						libopenblas)
# target_link_libraries(Tutorial PUBLIC ${EXTRA_LIBS})

# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h

target_include_directories(dcmap PUBLIC
                           "${PROJECT_BINARY_DIR}"
						   "${PROJECT_SOURCE_DIR}/Armadillo/"
						   "${PROJECT_SOURCE_DIR}"
                           )


add_custom_command(TARGET dcmap POST_BUILD        # Adds a post-build event to MyTest
    COMMAND ${CMAKE_COMMAND} -E copy_if_different  # which executes "cmake - E copy_if_different..."
        "${PROJECT_SOURCE_DIR}/libs/libopenblas.dll"      # <--this is in-file
        $<TARGET_FILE_DIR:dcmap>)                 # <--this is out-file path
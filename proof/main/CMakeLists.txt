
set(CMAKE_CXX_STANDARD 14)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# GMP
find_path(GMP_INCLUDE_DIR NAMES gmp.h)
find_library(GMP_LIBRARIES NAMES gmp libgmp)
find_library(GMPXX_LIBRARIES NAMES gmpxx libgmpxx)

# OpenMP
find_package(OpenMP)

set(SOURCE main.cpp link.hpp FR.hpp fft.hpp constants.hpp)
add_executable(main ${SOURCE}) 
add_dependencies(main client)
target_include_directories(main PUBLIC ${LIBFQFFTH} ${GMP_INCLUDE_DIR})
get_target_property(CLIENT_DIR client LOCATION)
target_link_libraries(main   ${GMP_LIBRARIES} ${GMPXX_LIBRARIES})
target_link_libraries(main ${CLIENT_DIR}/libclient.so)
target_link_libraries(main ${LIBFQFFT})
if(OpenMP_CXX_FOUND)
    target_link_libraries(main OpenMP::OpenMP_CXX)
endif()


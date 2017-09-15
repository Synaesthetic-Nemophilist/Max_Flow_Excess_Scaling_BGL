cmake_minimum_required(VERSION 3.6)
project(ZaroFlow)

set(CMAKE_CXX_STANDARD 11)



# locate Boost library
find_package(Boost 1.57.0 COMPONENTS system filesystem REQUIRED)
if(Boost_FOUND)
    message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
    message(STATUS "Boost_LIBRARIES: ${Boost_LIBRARIES}")
    message(STATUS "Boost_VERSION: ${Boost_VERSION}")
    include_directories(${Boost_INCLUDE_DIRS})
endif()

add_executable(BoostTest main.cpp)
if(Boost_FOUND)
    target_link_libraries(BoostTest ${Boost_LIBRARIES})
endif()

set(SOURCE_FILES main.cpp)
add_executable(ZaroFlow ${SOURCE_FILES})

# locate include and link LEDA libraries
set(LEDAROOOT /Users/nomaterials/Documents/Code_Libraries/LEDA-6.4-free-mavericks-llvm-5.1-mt)
include_directories(SYSTEM ${LEDAROOOT}/incl)
link_directories(${LEDAROOOT})
target_link_libraries(ZaroFlow ${LEDAROOOT}/libleda.a)
project(THM)

cmake_minimum_required(VERSION 2.6)

# Add local modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      ${PROJECT_SOURCE_DIR}/../../cmake/Modules)

if(NOT IFEM_CONFIGURED)
  find_package(IFEM REQUIRED)
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${IFEM_CXX_FLAGS}")

# Generate regtest script with correct paths
configure_file(${PROJECT_SOURCE_DIR}/../Common/scripts/regtest.sh.in
               regtest.sh)

set(THM_SRCS main_THM.C
             ThermoPoroElasticity.C ThermoPoroMaterial.C)

include_directories(${IFEM_INCLUDES} ${PROJECT_SOURCE_DIR}
                    ../Common ${PROJECT_SOURCE_DIR}/../IFEM-Elasticity
                    ${PROJECT_SOURCE_DIR}/../IFEM-PoroElasticity/PoroElastic)

add_subdirectory(${PROJECT_SOURCE_DIR}/../IFEM-PoroElasticity/PoroElastic PoroElastic)
add_subdirectory(${PROJECT_SOURCE_DIR}/../IFEM-Elasticity LinEl)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

add_executable(THM ${THM_SRCS})
target_link_libraries(THM PoroElastic Elasticity IFEMAppCommon ${IFEM_LIBRARIES})

# Installation
install(TARGETS THM DESTINATION bin COMPONENT bin)

enable_testing()

# Regression tests
if(MPI_FOUND)
  # Add parallel tests here
endif()
if (NOT MPI_FOUND OR IFEM_SERIAL_TESTS_IN_PARALLEL)
  # Add serial tests here
endif()

list(APPEND TEST_APPS THM)

if(IFEM_COMMON_APP_BUILD)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
  set(UNIT_TEST_NUMBER ${UNIT_TEST_NUMBER} PARENT_SCOPE)
else()
  add_check_target()
endif()

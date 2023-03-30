cmake_minimum_required (VERSION 2.8.11)
project (ifc_to_citygml)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)
if(COMPILER_SUPPORTS_CXX14)
    add_definitions(-std=c++14)
else()
    message(FATAL_ERROR "requires a compiler with C++14 support")
endif()

file(GLOB ifc_to_citygml_src "*.h" "*.cpp")

add_executable(ifc_to_citygml ${ifc_to_citygml_src})
target_link_libraries(ifc_to_citygml ${IFC_LIBRARIES} ${OCC_LIBRARIES} ${CGAL_LIBRARIES}
    ${MPFR_LIBRARIES} ${GMP_LIBRARIES} ${WS2_LIBRARIES} ${Boost_LIBRARIES})

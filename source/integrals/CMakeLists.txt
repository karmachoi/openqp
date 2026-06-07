
set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

file(GLOB SOURCES_mathlib CONFIGURE_DEPENDS "*.F90")

list(APPEND SOURCES ${SOURCES_mathlib})
# libintRot "rotate, don't recompute" host bridge (C++); compiled into liboqp.
list(APPEND SOURCES "${dir}/openqp_iface.cpp")
set(SOURCES ${SOURCES} PARENT_SCOPE)

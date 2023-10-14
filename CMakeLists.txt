#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(gmnana LANGUAGES CXX)

#----------------------------------------------------------------------------
# Implement user-defined items here
# The generated shared library will be named lib$(PACKAGE).so/dylib 
SET(PACKAGE gmnana)

# CXX standard
#SET(CMAKE_CXX_STANDARD 17)

# sources 
SET(src 
  Cut.cpp
  Fit.cpp
  Utilities.cpp
  SetROOTVar.cpp
  KinematicVar.cpp
  ExpConstants.cpp
  EMFFFits.cpp
)	

# headers
string(REPLACE .cpp .h headers "${src}")

# additional headers
SET(allheaders ${headers} Constants.h CodaRun.h FitFns.h SimuJob.h)

# Sources live in dedicated src directory
list(TRANSFORM src PREPEND "${CMAKE_CURRENT_LIST_DIR}/src/")

# Headers live in dedicated include directory
list(TRANSFORM allheaders PREPEND "${CMAKE_CURRENT_LIST_DIR}/include/")

# checks 
# MESSAGE(${src})
# MESSAGE(${allheaders}) 

#----------------------------------------------------------------------------
# Find ROOT 
# SET(ROOT_VER 6.22) 
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED)
# Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

MESSAGE(STATUS "Found ROOT: " $ENV{ROOTSYS})

#----------------------------------------------------------------------------
#configure_file( ${PROJECT_SOURCE_DIR}/rootlogon_build.C ${PROJECT_BINARY_DIR}/rootlogon.C )
configure_file( ${PROJECT_SOURCE_DIR}/rootlogon_install.C ${PROJECT_BINARY_DIR}/runenv/rootlogon.C )
configure_file( ${PROJECT_SOURCE_DIR}/gmn_ana.rootrc ${PROJECT_BINARY_DIR}/runenv/gmn_ana.rootrc )

#----------------------------------------------------------------------------
# The following lines are required so that the executable can dynamically locate the library dependencies independent of 
# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_LIBDIR}")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

SET(CMAKE_BUILD_TYPE RelWithDebInfo)

# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_FULL_LIBDIR}" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")
ENDIF("${isSystemDir}" STREQUAL "-1")

SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

#----------------------------------------------------------------------------
# copy environment setup script(s) for sh and csh to 
#configure_file( ${PROJECT_SOURCE_DIR}/gmn_ana.sh ${CMAKE_CURRENT_BINARY_DIR}/gmn_ana.sh )
#configure_file( ${PROJECT_SOURCE_DIR}/gmn_ana.csh ${CMAKE_CURRENT_BINARY_DIR}/gmn_ana.csh @ONLY)

#----------------------------------------------------------------------------
# Library
#add_library(${PACKAGE} SHARED ${src} ${allheaders} ${PACKAGE}Dict.cpp)
add_library(${PACKAGE} SHARED ${src} ${allheaders})

target_compile_options(${PACKAGE} PRIVATE ${${PROJECT_NAME_UC}_DIAG_FLAGS_LIST})

target_include_directories(${PACKAGE}
  PUBLIC
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  )
target_link_libraries(${PACKAGE} PUBLIC ${ROOT_LIBRARIES})

include(GNUInstallDirs)

install(TARGETS ${PACKAGE}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )
install(FILES ${allheaders} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
#install(DIRECTORY ${PROJECT_SOURCE_DIR}/DB DESTINATION ${CMAKE_INSTALL_PREFIX})
#install(DIRECTORY ${PROJECT_SOURCE_DIR}/scripts DESTINATION ${CMAKE_INSTALL_PREFIX})
#install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/gmn_ana.sh DESTINATION ${CMAKE_INSTALL_BINDIR})
#install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/gmn_ana.csh DESTINATION ${CMAKE_INSTALL_BINDIR})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/runenv/rootlogon.C DESTINATION ${CMAKE_INSTALL_PREFIX}/runenv)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/runenv/gmn_ana.rootrc DESTINATION ${CMAKE_INSTALL_PREFIX}/runenv)

#----------------------------------------------------------------------------
# ROOT dictionary
root_generate_dictionary(
 ${PACKAGE}Dict
 ${allheaders} 
 LINKDEF ${CMAKE_CURRENT_LIST_DIR}/include/${PACKAGE}_LinkDef.h  
 MODULE gmnana 
)


#MESSAGE(${CMAKE_CURRENT_BINARY_DIR})	
#MESSAGE(${CMAKE_INSTALL_BINDIR})	
#MESSAGE(${CMAKE_INSTALL_INCLUDEDIR})	
#MESSAGE(${PROJECT_BINARY_DIR})	
#MESSAGE(${CMAKE_INSTALL_LIBDIR})	

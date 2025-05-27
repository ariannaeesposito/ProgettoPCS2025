cmake_minimum_required(VERSION 3.0.0)

project(progetto2025 LANGUAGES C CXX VERSION 1.0.0)

# Set C++ Standard
################################################################################
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# IMPOSE WARNINGS ON DEBUG Aggiunge opzioni di warning in modalità Debug:
#-Wall: attiva i warning più comuni
#-Wextra: attiva warning aggiuntivi
#-pedantic-errors: rende errore qualsiasi cosa non standard secondo C++

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -pedantic-errors")

# IMPOSE CXX FLAGS FOR WINDOWS
if (WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-mbig-obj")
endif (WIN32)

# Add dependencies
#target:eseguibile,libreria ( oggetto costruibile )
###############################################################################

## Eigen3
#CONFIG: cerca file tipo Eigen3Config.cmake
#REQUIRED: se non lo trova → errore e stop
#APPEND: aggiungo librerie a una lista che raccoglie tutte le librerie da linkare
#PUBLIC:target usa la libreria e la propaga
#PRIVATE: Solo il target specificato la usa (non la propaga)

find_package(Eigen3 CONFIG REQUIRED)
list(APPEND ${PROJECT_NAME}_LINKED_LIBRARIES PUBLIC Eigen3::Eigen)

## Threads
find_package(Threads REQUIRED)
list(APPEND ${PROJECT_NAME}_LINKED_LIBRARIES PRIVATE Threads::Threads)

## GTest
find_package(GTest REQUIRED)

list(APPEND ${PROJECT_NAME}_LINKED_LIBRARIES PRIVATE GTest::gtest)
list(APPEND ${PROJECT_NAME}_LINKED_LIBRARIES PRIVATE GTest::gtest_main)
list(APPEND ${PROJECT_NAME}_LINKED_LIBRARIES PRIVATE GTest::gmock)
list(APPEND ${PROJECT_NAME}_LINKED_LIBRARIES PRIVATE GTest::gmock_main)

# Add subdirectory ##################################################################################

#utilizzo tre cmake in tutto il progetto e devo dire dove li trovo

add_subdirectory(src)
add_subdirectory(src_test)
add_subdirectory(ExportParaview)

#_sources	aggiungiamo Lista di file .cpp	
#_headers	aggiungiamo Lista di file .h o .hpp	
#_includes	(aggiungiamo percorsi) Cartelle in cui si trovano gli header(Dire al compilatore dove cercare quando uso #include) 
#{src_sources}: variabile che contiene una lista di file sorgenti provenienti dalla cartella src.

list(APPEND ${PROJECT_NAME}_sources ${src_sources} ${src_paraview_sources})
list(APPEND ${PROJECT_NAME}_headers ${src_headers} ${src_paraview_headers})
list(APPEND ${PROJECT_NAME}_includes ${src_includes} ${src_paraview_includes})

list(APPEND ${PROJECT_NAME}_TEST_sources ${${PROJECT_NAME}_sources} ${src_test_sources})
list(APPEND ${PROJECT_NAME}_TEST_headers ${${PROJECT_NAME}_headers} ${src_test_headers})
list(APPEND ${PROJECT_NAME}_TEST_includes ${${PROJECT_NAME}_includes} ${src_test_includes})

#################### Copy file ################################

# COPY : permette di copiare file da cartelle di origine a destinazioni di build senza cambiare la struttura del progetto.

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/PolygonalMesh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

############################## Add executable ##################################################

# creiamo target eseguibili: progetto2025 and progetto2025_TEST

add_executable(${PROJECT_NAME} main.cpp
      ${${PROJECT_NAME}_sources}
      ${${PROJECT_NAME}_headers})

target_include_directories(${PROJECT_NAME} PUBLIC ${${CMAKE_PROJECT_NAME}_includes})
target_link_libraries(${PROJECT_NAME} PUBLIC ${${PROJECT_NAME}_LINKED_LIBRARIES})
target_compile_options(${PROJECT_NAME} PUBLIC -fPIC)

add_executable(${PROJECT_NAME}_TEST main_test.cpp
      ${${PROJECT_NAME}_TEST_sources}
      ${${PROJECT_NAME}_TEST_headers})

target_include_directories(${PROJECT_NAME}_TEST PUBLIC ${${CMAKE_PROJECT_NAME}_TEST_includes})
target_link_libraries(${PROJECT_NAME}_TEST PUBLIC ${${PROJECT_NAME}_LINKED_LIBRARIES})
target_compile_options(${PROJECT_NAME}_TEST PUBLIC -fPIC)

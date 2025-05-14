list(APPEND src_headers ${CMAKE_CURRENT_SOURCE_DIR}/UCDUtilities.hpp)
list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/UCDUtilities.cpp)

list(APPEND src_includes ${CMAKE_CURRENT_SOURCE_DIR})

#Permettiamo al file di passare variabili al CMake “padre” usando set(... PARENT_SCOPE).
set(src_paraview_sources ${src_sources} PARENT_SCOPE)
set(src_paraview_headers ${src_headers} PARENT_SCOPE)
set(src_paraview_includes ${src_includes} PARENT_SCOPE)

# src_paraview_sources è una variabile (a cui viene assegnato il contenuto di src_sources) che viene definita nel file CMakeLists.txt dentro la cartella ExpertParaview e la utilizzo nel cmake "padre"

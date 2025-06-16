cmake_minimum_required(VERSION 3.10)

project(MiProyecto)

# Estándar C++
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Ruta SDK de macOS (ajusta si cambias de versión o usas otro entorno)
set(CMAKE_OSX_SYSROOT "/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX15.5.sdk")

# Incluir directorios
include_directories(
    ${CMAKE_SOURCE_DIR}/../include              # tus headers (como io_handler.hpp)
    /opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3
    /opt/homebrew/include
    ${Python3_INCLUDE_DIRS}
    /opt/homebrew/lib/python3.13/site-packages/numpy/_core/include
)

# Directorios de librerías
link_directories(
    /opt/homebrew/lib
)

# Buscar dependencias
find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)
find_package(Ceres REQUIRED)
find_package(autodiff REQUIRED)

# Archivos fuente (añade los que vayas creando)
set(SOURCES
    main.cpp
    io_handler.cpp
    model_data.cpp
    model_functions.cpp
    visualization.cpp
    regression.cpp
    estimations.cpp
    equations.cpp
    optimization.cpp
    evolutive.cpp
)

# Crear ejecutable
add_executable(main ${SOURCES})

# Enlazar bibliotecas
target_link_libraries(main
    ${Python3_LIBRARIES}
    boost_math_tr1
    Ceres::ceres
)

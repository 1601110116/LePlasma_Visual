cmake_minimum_required(VERSION 3.5.1)
project(LePlasma_Visual)
set(SOURCE_FILES
        include/Figure.h
        include/GridViewer.h
        include/Mathx.h
        include/Plot.h
        include/Visualize.h
        include/XYZIndicator.h
        include_cal/Buffer.h
        include_cal/Case.h
        include_cal/Cell.h
        include_cal/CSPIC.h
        include_cal/DispersionRelation.h
        include_cal/Distributer.h
        include_cal/Electron.h
        include_cal/Engine.h
        include_cal/Eular.h
        include_cal/ExternalField.h
        include_cal/Grid.h
        include_cal/HostGrid.h
        include_cal/indexCache.h
        include_cal/MPIGrid.h
        include_cal/Particle.h
        include_cal/ParticleDistributer.h
        include_cal/Macros.h
        include_cal/PlateWave.h
        include_cal/Range.h
        include_cal/RunManager.h
        include_cal/Tensor3D.h
        include_cal/UniformB.h
        include_cal/UniformBFalse.h
        include_cal/UniformE.h
        include_cal/Vector3D.h
        include_cal/Vertex.h
        LeFrame/include/AScompat.h
        LeFrame/include/DisplayObject.h
        LeFrame/include/Event.h
        LeFrame/include/EventDispatcher.h
        LeFrame/include/EventListener.h
        LeFrame/include/GLEnvironment.h
        LeFrame/include/GLEventManager.h
        LeFrame/include/Graphics.h
        LeFrame/include/KeyBoard.h
        LeFrame/include/Rect.h
        LeFrame/include/Sprite.h
        LeFrame/include/Stage.h
        LeFrame/include/TextField.h
        LeFrame/src/DisplayObject.cpp
        LeFrame/src/EventDispatcher.cpp
        LeFrame/src/EventListener.cpp
        LeFrame/src/GLEnvironment.cpp
        LeFrame/src/GLEventManager.cpp
        LeFrame/src/Graphics.cpp
        LeFrame/src/Rect.cpp
        LeFrame/src/Sprite.cpp
        LeFrame/src/Stage.cpp
        LeFrame/src/TextField.cpp
        src/Figure.cpp
        src/GridViewer.cpp
        src/Mathx.cpp
        src/Plot.cpp
        src/Visualize.cpp
        src/XYZIndicator.cpp
        src_cal/Case.cpp
        src_cal/Cell.cpp
        src_cal/CSPIC.cpp
        src_cal/DispersionRelation.cpp
        src_cal/Distributer.cpp
        src_cal/Electron.cpp
        src_cal/Engine.cpp
        src_cal/Eular.cpp
        src_cal/ExternalField.cpp
        src_cal/Grid.cpp
        src_cal/HostGrid.cpp
        src_cal/MPIGrid.cpp
        src_cal/Particle.cpp
        src_cal/ParticleDistributer.cpp
        src_cal/PlateWave.cpp
        src_cal/Range.cpp
        src_cal/RunManager.cpp
        src_cal/Tensor3D.cpp
        src_cal/UniformB.cpp
        src_cal/UniformBFalse.cpp
        src_cal/UniformE.cpp
        src_cal/Vector3D.cpp
        src_cal/Vertex.cpp
        DispersionRelation.h
        LePlasma.h
        LePlasma_Visual.cpp include_cal/EngineForSingleParticle.h src_cal/EngineForSingleParticle.cpp include_cal/CaseOfDipole.h src_cal/CaseOfDipole.cpp)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

include_directories(
        library ${OPENGL_INCLUDE_DIRS} ${GLUT_INCLUDE_DIRS}
        /usr/include/mpi/
        ./
        ../
        ./include
        ./include_cal
        ./LeFrame/include)
add_executable(LePlasma_Visual ${SOURCE_FILES})
target_link_libraries(LePlasma_Visual
        ${OPENGL_LIBRARIES}
        ${GLUT_glut_LIBRARY})
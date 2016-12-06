################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Figure.cpp \
../src/GridViewer.cpp \
../src/Mathx.cpp \
../src/Plot.cpp \
../src/Visualize.cpp \
../src/XYZIndicator.cpp 

OBJS += \
./src/Figure.o \
./src/GridViewer.o \
./src/Mathx.o \
./src/Plot.o \
./src/Visualize.o \
./src/XYZIndicator.o 

CPP_DEPS += \
./src/Figure.d \
./src/GridViewer.d \
./src/Mathx.d \
./src/Plot.d \
./src/Visualize.d \
./src/XYZIndicator.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I"/home/ylang/workspace/LePlasma_Visual" -I"/home/ylang/workspace/LePlasma_Visual/include" -I"/home/ylang/workspace/LePlasma_Visual/include_cal" -I"/home/ylang/workspace/LePlasma_Visual/LeFrame/include" -O0 -g3 -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



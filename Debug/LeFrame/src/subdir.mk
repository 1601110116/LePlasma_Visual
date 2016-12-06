################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../LeFrame/src/DisplayObject.cpp \
../LeFrame/src/EventDispatcher.cpp \
../LeFrame/src/EventListener.cpp \
../LeFrame/src/GLEnvironment.cpp \
../LeFrame/src/GLEventManager.cpp \
../LeFrame/src/Graphics.cpp \
../LeFrame/src/Rect.cpp \
../LeFrame/src/Sprite.cpp \
../LeFrame/src/Stage.cpp \
../LeFrame/src/TextField.cpp 

OBJS += \
./LeFrame/src/DisplayObject.o \
./LeFrame/src/EventDispatcher.o \
./LeFrame/src/EventListener.o \
./LeFrame/src/GLEnvironment.o \
./LeFrame/src/GLEventManager.o \
./LeFrame/src/Graphics.o \
./LeFrame/src/Rect.o \
./LeFrame/src/Sprite.o \
./LeFrame/src/Stage.o \
./LeFrame/src/TextField.o 

CPP_DEPS += \
./LeFrame/src/DisplayObject.d \
./LeFrame/src/EventDispatcher.d \
./LeFrame/src/EventListener.d \
./LeFrame/src/GLEnvironment.d \
./LeFrame/src/GLEventManager.d \
./LeFrame/src/Graphics.d \
./LeFrame/src/Rect.d \
./LeFrame/src/Sprite.d \
./LeFrame/src/Stage.d \
./LeFrame/src/TextField.d 


# Each subdirectory must supply rules for building sources it contributes
LeFrame/src/%.o: ../LeFrame/src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I"/home/ylang/workspace/LePlasma_Visual" -I"/home/ylang/workspace/LePlasma_Visual/include" -I"/home/ylang/workspace/LePlasma_Visual/include_cal" -I"/home/ylang/workspace/LePlasma_Visual/LeFrame/include" -O0 -g3 -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



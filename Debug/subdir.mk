################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../LePlasma_Visual.cpp 

OBJS += \
./LePlasma_Visual.o 

CPP_DEPS += \
./LePlasma_Visual.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I"/home/ylang/workspace/LePlasma_Visual" -I"/home/ylang/workspace/LePlasma_Visual/include" -I"/home/ylang/workspace/LePlasma_Visual/include_cal" -I"/home/ylang/workspace/LePlasma_Visual/LeFrame/include" -O0 -g3 -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



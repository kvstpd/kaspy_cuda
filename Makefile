
CUDA_PATH       := /Developer/NVIDIA/CUDA-6.5
CUDA_LIBPATH := $(CUDA_PATH)/lib


NVCC := $(CUDA_PATH)/bin/nvcc

LDFLAGS := -L$(CUDA_LIBPATH) -lcudart_static -lc++

GENCODE_FLAGS := -gencode arch=compute_11,code=compute_11

OPT_FLAGS := -O3


#C_LIBPATH :=
#F_LIBPATH :=
#SDK_LIBPATH :=

#C_INCPATH :=

#CUDA_LIBS := cudart_static.lib

F_OPTS := $(OPT_FLAGS)

# /arch:SSE4.2 /QxSSE4.2 /O3 /Qparallel /Zp16

C_OPTS := $(OPT_FLAGS)

#/Ox /favor:INTEL64 /Zp16

NV_OPTS := --machine 64 --optimize 3 --compiler-options "$(C_OPTS)"

#LINK_OPTS := /MACHINE:X64 /SUBSYSTEM:CONSOLE







all : build/kaspy

test : build/kaspy
	cd build && pwd && ./kaspy

build/kaspy : kaspy.o cycler.o KaspyCycler.o
	$(NVCC) -ccbin gfortran $(OPT_FLAGS) -o $@ $(LDFLAGS) $(GENCODE_FLAGS)  $+

kaspy.o: kaspy.for
	gfortran $(OPT_FLAGS) -o $@ -c $<

cycler.o: cycler.cpp
	$(NVCC) -ccbin gcc $(OPT_FLAGS) -o $@ $(GENCODE_FLAGS)  -c $<

KaspyCycler.o: KaspyCycler.cu
	$(NVCC) -ccbin gcc $(OPT_FLAGS) -o $@ $(GENCODE_FLAGS)  -c $<

oclean:
	rm -f *.o

clean:
	rm -f *.o
	rm -f build/kaspy


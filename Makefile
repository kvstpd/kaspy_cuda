
C_LIBPATH := c:/Program Files (x86)/Microsoft Visual Studio 10.0/VC/lib/amd64
F_LIBPATH := c:/Program Files (x86)/Intel/ComposerXE-2011/compiler/lib/intel64
SDK_LIBPATH := c:/Program Files (x86)/Microsoft SDKs/Windows/v7.0A/Lib/x64
CUDA_LIBPATH := c:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v7.5/lib/x64

C_INCPATH := c:/Program Files (x86)/Microsoft Visual Studio 10.0/VC/include

CUDA_LIBS := cudart_static.lib

F_OPTS := /Zp16
#/arch:SSE4.2 /QxSSE4.2 /O3 /Qparallel /Zp16

C_OPTS := /favor:INTEL64 /Zp16
#/Ox /favor:INTEL64 /Zp16

NV_OPTS := --machine 64 --optimize 3 -I"$(C_INCPATH)" --compiler-options "$(C_OPTS)"

LINK_OPTS := /MACHINE:X64 /SUBSYSTEM:CONSOLE


all : build/kaspy_cuda.exe

test : build/kaspy_cuda.exe
	cd build && pwd && ./kaspy_cuda.exe

build/kaspy_cuda.exe : kaspy.obj KaspyCycler.obj cycler.obj
	link $+ $(CUDA_LIBS) /OUT:$@ $(LINK_OPTS) /LIBPATH:"$(F_LIBPATH)" /LIBPATH:"$(C_LIBPATH)" /LIBPATH:"$(SDK_LIBPATH)" /LIBPATH:"$(CUDA_LIBPATH)"

KaspyCycler.obj: KaspyCycler.cu
	nvcc $< --compile $(NV_OPTS)

cycler.obj : cycler.cpp
	nvcc $< --compile $(NV_OPTS)

kaspy.obj : kaspy.for
	ifort $(F_OPTS) /c $<

clean :
	rm -f *.obj
	rm -f ./build/kaspy_cuda.exe
	rm -f *.exp
	rm -f *__genmod.mod
	rm -f *__genmod.f90
	rm -f *.lib


TARGET = spectral.so

all: $(TARGET)

spectral.so: spectral.f
	f2py -c --fcompiler=gfortran spectral.f --f77flags=-O3 -m spectral


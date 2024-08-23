TARGET = spectral.so

all: $(TARGET)

spectral.so: spectral.f
	f2py -c --fcompiler=gfortran $< --f77flags=-O3 -m spectral
	@echo "--------------------------------------------------------------"
	python doc.py
	@echo "--------------------------------------------------------------"


FC = gfortran
FCFLAGS = -O2 -Wall

TARGET = spectral.so driver_fd test_fft_diff

all: $(TARGET)

spectral.so: spectral.f
	f2py -c --fcompiler=gfortran $< --f77flags=-O3 -m spectral
	@echo "--------------------------------------------------------------"
	python doc.py
	@echo "--------------------------------------------------------------"

#fd.o: fd.f
#driver_fd.o: driver_fd.f
#spectral.o: spectral.f

%.o: %.f
	$(FC) $(FCFLAGS) -c $<

driver_fd: driver_fd.o fd.o
	$(FC) $(FCFLAGS) -o $@ $^

test_fft_diff: spectral.o test_fft_diff.o
	$(FC) $(FCFLAGS) -o $@ $^

clean:
	rm -f $(TARGET) *.o

# Build the tidal-response code.
FC     ?= gfortran
FFLAGS ?=
SOURCES = k2_Dt.f90 Forced_tide.f90 Model.f90 Pot_Strain_Func.f90 Tidal_heating_map_output.f90

Forced_tide: $(SOURCES)
	$(FC) $(FFLAGS) -o $@ $(SOURCES)
	rm -f *.mod

clean:
	rm -f Forced_tide *.mod

.PHONY: clean

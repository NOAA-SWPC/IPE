lib_LIBRARIES = libshared.a

libshared_a_SOURCES = ipe_error_module.F90

libipe_a_FCFLAGS =

.PHONY: clean-modules

clean-modules:
	-test -z "$(FC_MODEXT)" || rm -f *.$(FC_MODEXT)

clean-local: clean-modules

# dependencies

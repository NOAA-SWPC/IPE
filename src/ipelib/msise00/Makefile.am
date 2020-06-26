noinst_LIBRARIES = libmsise.a

libmsise_a_SOURCES = hwm14.f90 physics_constants.f90 physics_msis.f90 utils_constants.f90 utils_spline.f90 

libmsise_a_FCFLAGS = $(FC_MODINC) .  $(FC_MODINC) ../shared

EXTRA_DIST = hwm93.f90 nrlmsise00_driver.for nrlmsise00_sub.for test77gtd7.f testgtd7.f90 utils_spline.f90


.PHONY: clean-modules

clean-modules:
	-test -z "$(FC_MODEXT)" || rm -f *.$(FC_MODEXT)

clean-local: clean-modules

# dependencies
libmsise_a-physics_constants.$(OBJEXT): libmsise_a-utils_constants.$(OBJEXT)
libmsise_a-physics_msis.$(OBJEXT): libmsise_a-physics_constants.$(OBJEXT) \
	libmsise_a-utils_constants.$(OBJEXT) libmsise_a-utils_spline.$(OBJEXT)
libmsise_a-utils_spline.$(OBJEXT): libmsise_a-utils_constants.$(OBJEXT)

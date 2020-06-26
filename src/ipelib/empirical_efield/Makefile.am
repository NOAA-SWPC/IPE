noinst_LIBRARIES = libefield.a

libefield_a_SOURCES = apex_sunloc_empirical.f efield_ipe.f weimer2005_ipe.f


libefield_a_FFLAGS = $(FC_MODINC) . $(FC_MODINC) ../shared

EXTRA_DIST = apex_sunloc.f

.PHONY: clean-modules

clean-modules:
	-test -z "$(FC_MODEXT)" || rm -f *.$(FC_MODEXT)

clean-local: clean-modules

# dependencies
libefield_a-efield_ipe.$(OBJEXT) : libefield_a-weimer2005_ipe.$(OBJEXT)



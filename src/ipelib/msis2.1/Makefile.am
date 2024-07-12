noinst_LIBRARIES = libmsis2_1.a

libmsis2_1_a_SOURCES = hwm14.f90 msis_utils.F90 msis_constants.F90 msis_init.F90 msis_gfn.F90 msis_tfn.F90 msis_dfn.F90 msis_calc.F90 msis_gtd8d.F90

libmsis2_1_a_FCFLAGS = $(FC_MODINC) .  $(FC_MODINC) ../shared

.PHONY: clean-modules

clean-modules:
	-test -z "$(FC_MODEXT)" || rm -f *.$(FC_MODEXT)

clean-local: clean-modules

# dependencies
libmsis2_1_a-msis_utils.$(OBJEXT): libmsis2_1_a-msis_constants.$(OBJEXT)
libmsis2_1_a-msis_constants.$(OBJEXT): libmsis2_1_a-msis_constants.$(OBJEXT)
libmsis2_1_a-msis_gfn.$(OBJEXT): libmsis2_1_a-msis_constants.$(OBJEXT) \
	libmsis2_1_a-msis_init.$(OBJEXT)
libmsis2_1_a-msis_tfn.$(OBJEXT): libmsis2_1_a-msis_constants.$(OBJEXT) \
	libmsis2_1_a-msis_init.$(OBJEXT) libmsis2_1_a-msis_gfn.$(OBJEXT) \
	libmsis2_1_a-msis_utils.$(OBJEXT)
libmsis2_1_a-msis_dfn.$(OBJEXT): libmsis2_1_a-msis_constants.$(OBJEXT) \
	libmsis2_1_a-msis_init.$(OBJEXT) libmsis2_1_a-msis_gfn.$(OBJEXT) \
	libmsis2_1_a-msis_utils.$(OBJEXT) libmsis2_1_a-msis_tfn.$(OBJEXT)
libmsis2_1_a-msis_calc.$(OBJEXT): libmsis2_1_a-msis_constants.$(OBJEXT) \
	libmsis2_1_a-msis_init.$(OBJEXT) libmsis2_1_a-msis_gfn.$(OBJEXT) \
	libmsis2_1_a-msis_utils.$(OBJEXT) libmsis2_1_a-msis_tfn.$(OBJEXT) \
	libmsis2_1_a-msis_dfn.$(OBJEXT)
libmsis2_1_a-msis_gtd8d.$(OBJEXT): libmsis2_1_a-msis_constants.$(OBJEXT) \
	libmsis2_1_a-msis_init.$(OBJEXT) libmsis2_1_a-msis_calc.$(OBJEXT)

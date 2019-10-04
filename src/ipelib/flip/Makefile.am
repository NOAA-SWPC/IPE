noinst_LIBRARIES = libflip.a

libflip_a_SOURCES = CTIPE-int.f CMINOR.f ELECXS.f FLIP_GRID.f INIT-PROFILES.f KEMPRN.f MINORA.f Neut_Heating.f Photoel-Freqs.f RSDENA_EXB.f RSJACA.f RSLPSD.f RSLPST.f RSPE2B.f RSPRIM.f RSTEMC_EXB.f Rates.f

libflip_a_FFLAGS = $(FC_MODINC) .

EXTRA_DIST = ASFLIPIO.f FLIPAP.f FLIP_PC.f FLIPSUBS.f FLIPWIND.f GET_PARAMS.f IRI95.f MagDec.f RSSUBC.f

.PHONY: clean-modules

clean-modules:
	-test -z "$(FC_MODEXT)" || rm -f *.$(FC_MODEXT)

clean-local: clean-modules

# dependencies
libflip_a-CMINOR.$(OBJEXT): libflip_a-CTIPE-int.$(OBJEXT)
libflip_a-INIT-PROFILES.$(OBJEXT): libflip_a-CTIPE-int.$(OBJEXT)
libflip_a-KEMPRN.$(OBJEXT): libflip_a-CTIPE-int.$(OBJEXT)
libflip_a-MINORA.$(OBJEXT): libflip_a-CTIPE-int.$(OBJEXT)
libflip_a-Neut_Heating.$(OBJEXT): libflip_a-CTIPE-int.$(OBJEXT)
libflip_a-RSDENA_EXB.$(OBJEXT): libflip_a-CTIPE-int.$(OBJEXT)
libflip_a-RSLPSD.$(OBJEXT): libflip_a-CTIPE-int.$(OBJEXT)
libflip_a-RSLPST.$(OBJEXT): libflip_a-CTIPE-int.$(OBJEXT)
libflip_a-RSPE2B.$(OBJEXT): libflip_a-CTIPE-int.$(OBJEXT)
libflip_a-RSPRIM.$(OBJEXT): libflip_a-CTIPE-int.$(OBJEXT)
libflip_a-RSTEMC_EXB.$(OBJEXT): libflip_a-CTIPE-int.$(OBJEXT)

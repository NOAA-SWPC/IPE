SUBDIRS = ipelib .

if COUPLED
lib_LIBRARIES = libipm.a

libipm_a_LIBADD  = ipelib/*.$(OBJEXT) \
	ipelib/dynamo/*.$(OBJEXT) \
	ipelib/empirical_efield/*.$(OBJEXT) \
	ipelib/flip/*.$(OBJEXT) \
	ipelib/msise00/*.$(OBJEXT)

libipm_a_FCFLAGS = $(FC_MODINC) ipelib \
	$(FC_MODINC) ipelib/dynamo \
        $(FC_MODINC) ipelib/empirical_efield \
        $(FC_MODINC) ipelib/flip \
        $(FC_MODINC) ipelib/msise00

libipm_a_SOURCES = ipeCap.F90 ipeMethods.F90 IPE_Wrapper.F90

libipm_a-ipeCap.$(OBJEXT): libipm_a-ipeMethods.$(OBJEXT) libipm_a-IPE_Wrapper.$(OBJEXT)
libipm_a-ipeMethods.$(OBJEXT): libipm_a-IPE_Wrapper.$(OBJEXT)
libipm_a-IPE_Wrapper.$(OBJEXT): $(builddir)/ipelib/libipe.a
else
bin_PROGRAMS = ipe.x

ipe_x_SOURCES = IPE_Driver.F90

ipe_x_FCFLAGS = $(FC_MODINC) ipelib \
	$(FC_MODINC) ipelib/dynamo \
        $(FC_MODINC) ipelib/empirical_efield \
        $(FC_MODINC) ipelib/flip \
        $(FC_MODINC) ipelib/msise00

ipe_x_LDADD = $(builddir)/ipelib/libipe.a
endif

.PHONY: clean-modules

clean-modules:
	-test -z "$(FC_MODEXT)" || rm -f *.$(FC_MODEXT)

clean-local: clean-modules

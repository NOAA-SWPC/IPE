SUBDIRS = shared dynamo empirical_efield flip msise00 .

lib_LIBRARIES = libipe.a

libipe_a_SOURCES = IPE_Common_Routines.F90 IPE_Constants_Dictionary.F90 \
	IPE_Electrodynamics_Class.F90 IPE_Forcing_Class.F90 \
	IPE_Grid_Class.F90 IPE_Model_Class.F90 IPE_Model_Parameters_Class.F90 \
	IPE_MPI_Layer_Class.F90 IPE_Neutrals_Class.F90 IPE_Plasma_Class.F90 \
	IPE_Precision.F90 IPE_Time_Class.F90

libipe_a_FCFLAGS = $(FC_MODINC) . \
	$(FC_MODINC) dynamo \
	$(FC_MODINC) empirical_efield \
	$(FC_MODINC) flip \
	$(FC_MODINC) msise00 \
	$(FC_MODINC) shared

libipe_a_LIBADD = dynamo/*.$(OBJEXT) \
	empirical_efield/*.$(OBJEXT) \
	flip/*.$(OBJEXT) \
	msise00/*.$(OBJEXT) \
	shared/*.$(OBJEXT)

.PHONY: clean-modules

clean-modules:
	-test -z "$(FC_MODEXT)" || rm -f *.$(FC_MODEXT)

clean-local: clean-modules

# dependencies
libipe_a-IPE_Common_Routines.$(OBJEXT): libipe_a-IPE_Precision.$(OBJEXT)

libipe_a-IPE_Constants_Dictionary.$(OBJEXT): libipe_a-IPE_Precision.$(OBJEXT)

libipe_a-IPE_Driver.$(OBJEXT): libipe_a-IPE_Model_Class.$(OBJEXT) libipe_a-IPE_Precision.$(OBJEXT)

libipe_a-IPE_Electrodynamics_Class.$(OBJEXT): libipe_a-IPE_Common_Routines.$(OBJEXT) \
	libipe_a-IPE_Constants_Dictionary.$(OBJEXT) libipe_a-IPE_Forcing_Class.$(OBJEXT) \
	libipe_a-IPE_Grid_Class.$(OBJEXT) libipe_a-IPE_MPI_Layer_Class.$(OBJEXT) \
	libipe_a-IPE_Plasma_Class.$(OBJEXT) libipe_a-IPE_Precision.$(OBJEXT) \
	libipe_a-IPE_Time_Class.$(OBJEXT)

libipe_a-IPE_Forcing_Class.$(OBJEXT): libipe_a-IPE_Common_Routines.$(OBJEXT) \
	libipe_a-IPE_Constants_Dictionary.$(OBJEXT) libipe_a-IPE_Model_Parameters_Class.$(OBJEXT) \
	libipe_a-IPE_Precision.$(OBJEXT) libipe_a-IPE_MPI_Layer_Class.$(OBJEXT)

libipe_a-IPE_Grid_Class.$(OBJEXT): libipe_a-IPE_Common_Routines.$(OBJEXT) \
	libipe_a-IPE_Constants_Dictionary.$(OBJEXT) libipe_a-IPE_Model_Parameters_Class.$(OBJEXT) \
	libipe_a-IPE_MPI_Layer_Class.$(OBJEXT) libipe_a-IPE_Precision.$(OBJEXT)

libipe_a-IPE_Model_Class.$(OBJEXT): libipe_a-IPE_Electrodynamics_Class.$(OBJEXT) \
	libipe_a-IPE_Forcing_Class.$(OBJEXT) libipe_a-IPE_Grid_Class.$(OBJEXT) \
	libipe_a-IPE_Model_Parameters_Class.$(OBJEXT) libipe_a-IPE_MPI_Layer_Class.$(OBJEXT) \
	libipe_a-IPE_Neutrals_Class.$(OBJEXT) libipe_a-IPE_Plasma_Class.$(OBJEXT) \
	libipe_a-IPE_Precision.$(OBJEXT) libipe_a-IPE_Time_Class.$(OBJEXT)

libipe_a-IPE_Model_Parameters_Class.$(OBJEXT): libipe_a-IPE_Common_Routines.$(OBJEXT) \
	libipe_a-IPE_Constants_Dictionary.$(OBJEXT) libipe_a-IPE_MPI_Layer_Class.$(OBJEXT) \
	libipe_a-IPE_Precision.$(OBJEXT)

libipe_a-IPE_MPI_Layer_Class.$(OBJEXT): libipe_a-IPE_Precision.$(OBJEXT)

libipe_a-IPE_Neutrals_Class.$(OBJEXT): libipe_a-IPE_Constants_Dictionary.$(OBJEXT) \
	libipe_a-IPE_Forcing_Class.$(OBJEXT) libipe_a-IPE_Grid_Class.$(OBJEXT) \
	libipe_a-IPE_Precision.$(OBJEXT) libipe_a-IPE_Time_Class.$(OBJEXT)

libipe_a-IPE_Plasma_Class.$(OBJEXT): libipe_a-IPE_Common_Routines.$(OBJEXT) \
	libipe_a-IPE_Constants_Dictionary.$(OBJEXT) libipe_a-IPE_Forcing_Class.$(OBJEXT) \
	libipe_a-IPE_Grid_Class.$(OBJEXT) libipe_a-IPE_MPI_Layer_Class.$(OBJEXT) \
	libipe_a-IPE_Neutrals_Class.$(OBJEXT) libipe_a-IPE_Precision.$(OBJEXT) \
	libipe_a-IPE_Time_Class.$(OBJEXT)

libipe_a-IPE_Time_Class.$(OBJEXT): libipe_a-IPE_Precision.$(OBJEXT)

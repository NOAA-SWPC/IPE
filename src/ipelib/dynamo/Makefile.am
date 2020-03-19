noinst_LIBRARIES = libdynamo.a

libdynamo_a_SOURCES = cons_module.f dynamo_module.f heelis_module.f module_ceee.f module_clearcee.f \
	module_cnm.f module_cnmmod.f module_colath.f module_divide.f module_edges.f module_flwv32.f \
	module_highlat.f module_htrpex.f module_init_cons.f module_init_heelis.f module_magfield.f \
	module_potm.f module_readin_ascii.f module_rhspde.f module_stencil.f module_stencmd.f \
	module_stenmd.f module_stenmod.f module_sub_dynamo.f module_sub_heelis.f module_threed.f \
	module_transf.f module_weimer2005Ipe.f mudcom.f mud.f mudmod.f muh2cr.f params.f readin.f util.f

libdynamo_a_FFLAGS = $(FC_MODINC) .

EXTRA_DIST = apex_sunloc.f module_cal_monthday.f module_check_err.f module_eldyn.f module_nc_close.f \
	module_nc_create.f module_plas2dyn_fli_array.f module_sub_nclplot1d.f module_sub_nclplot.f \
	module_update_fli.f nc_module.f

.PHONY: clean-modules

clean-modules:
	-test -z "$(FC_MODEXT)" || rm -f *.$(FC_MODEXT)

clean-local: clean-modules

# dependencies
libdynamo_a-cons_module.$(OBJEXT): libdynamo_a-params.$(OBJEXT)
libdynamo_a-dynamo_module.$(OBJEXT): libdynamo_a-params.$(OBJEXT)
libdynamo_a-module_cnm.$(OBJEXT): libdynamo_a-dynamo_module.$(OBJEXT)
libdynamo_a-module_cnmmod.$(OBJEXT): libdynamo_a-dynamo_module.$(OBJEXT)
libdynamo_a-module_colath.$(OBJEXT): libdynamo_a-cons_module.$(OBJEXT) libdynamo_a-dynamo_module.$(OBJEXT) \
	libdynamo_a-heelis_module.$(OBJEXT) libdynamo_a-module_magfield.$(OBJEXT) libdynamo_a-params.$(OBJEXT)
libdynamo_a-module_divide.$(OBJEXT): libdynamo_a-dynamo_module.$(OBJEXT)
libdynamo_a-module_edges.$(OBJEXT): libdynamo_a-dynamo_module.$(OBJEXT)
libdynamo_a-module_flwv32.$(OBJEXT): libdynamo_a-cons_module.$(OBJEXT) libdynamo_a-heelis_module.$(OBJEXT)
libdynamo_a-module_highlat.$(OBJEXT): libdynamo_a-dynamo_module.$(OBJEXT) libdynamo_a-module_colath.$(OBJEXT) \
	libdynamo_a-module_sub_heelis.$(OBJEXT) libdynamo_a-module_weimer2005Ipe.$(OBJEXT) libdynamo_a-params.$(OBJEXT)
libdynamo_a-module_htrpex.$(OBJEXT): libdynamo_a-dynamo_module.$(OBJEXT)
libdynamo_a-module_init_cons.$(OBJEXT): libdynamo_a-cons_module.$(OBJEXT)
libdynamo_a-module_init_heelis.$(OBJEXT): libdynamo_a-cons_module.$(OBJEXT) libdynamo_a-heelis_module.$(OBJEXT)
libdynamo_a-module_potm.$(OBJEXT): libdynamo_a-cons_module.$(OBJEXT) libdynamo_a-dynamo_module.$(OBJEXT) \
	libdynamo_a-heelis_module.$(OBJEXT) libdynamo_a-module_flwv32.$(OBJEXT) libdynamo_a-module_magfield.$(OBJEXT) \
	libdynamo_a-params.$(OBJEXT)
libdynamo_a-module_readin_ascii.$(OBJEXT): libdynamo_a-dynamo_module.$(OBJEXT) libdynamo_a-params.$(OBJEXT) \
	libdynamo_a-readin.$(OBJEXT)
libdynamo_a-module_rhspde.$(OBJEXT): libdynamo_a-cons_module.$(OBJEXT) libdynamo_a-dynamo_module.$(OBJEXT)
libdynamo_a-module_stencil.$(OBJEXT): libdynamo_a-module_cnm.$(OBJEXT) libdynamo_a-module_htrpex.$(OBJEXT)
libdynamo_a-module_stencmd.$(OBJEXT): libdynamo_a-module_cnm.$(OBJEXT) libdynamo_a-module_cnmmod.$(OBJEXT) \
	libdynamo_a-module_htrpex.$(OBJEXT)
libdynamo_a-module_stenmd.$(OBJEXT): libdynamo_a-cons_module.$(OBJEXT) libdynamo_a-dynamo_module.$(OBJEXT)
libdynamo_a-module_stenmod.$(OBJEXT): libdynamo_a-cons_module.$(OBJEXT) libdynamo_a-dynamo_module.$(OBJEXT)
libdynamo_a-module_sub_dynamo.$(OBJEXT): libdynamo_a-cons_module.$(OBJEXT) libdynamo_a-dynamo_module.$(OBJEXT) \
	libdynamo_a-module_clearcee.$(OBJEXT) libdynamo_a-module_divide.$(OBJEXT) libdynamo_a-module_edges.$(OBJEXT) \
	libdynamo_a-module_rhspde.$(OBJEXT) libdynamo_a-module_stencil.$(OBJEXT) libdynamo_a-module_stencmd.$(OBJEXT) \
	libdynamo_a-module_stenmd.$(OBJEXT) libdynamo_a-module_stenmod.$(OBJEXT) libdynamo_a-module_threed.$(OBJEXT) \
	libdynamo_a-module_transf.$(OBJEXT) libdynamo_a-params.$(OBJEXT)
libdynamo_a-module_sub_heelis.$(OBJEXT): libdynamo_a-heelis_module.$(OBJEXT) libdynamo_a-module_colath.$(OBJEXT) \
	libdynamo_a-module_potm.$(OBJEXT)
libdynamo_a-module_threed.$(OBJEXT): libdynamo_a-cons_module.$(OBJEXT) libdynamo_a-dynamo_module.$(OBJEXT) \
	libdynamo_a-params.$(OBJEXT)
libdynamo_a-module_transf.$(OBJEXT): libdynamo_a-cons_module.$(OBJEXT) libdynamo_a-dynamo_module.$(OBJEXT)
libdynamo_a-mud.$(OBJEXT): libdynamo_a-dynamo_module.$(OBJEXT) libdynamo_a-module_ceee.$(OBJEXT)
libdynamo_a-mudmod.$(OBJEXT): libdynamo_a-dynamo_module.$(OBJEXT)
libdynamo_a-muh2cr.$(OBJEXT): libdynamo_a-dynamo_module.$(OBJEXT) libdynamo_a-module_ceee.$(OBJEXT)

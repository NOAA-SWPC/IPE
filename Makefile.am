SUBDIRS = src
dist_doc_DATA = README.md

ACLOCAL_AMFLAGS = -I m4

if COUPLED
nodist_include_HEADERS = src/ipecap.$(FC_MODEXT)
nodist_data_DATA = ipe.mk
ipe.mk:
	@echo "# ESMF self-describing build dependency makefile fragment" > $(@)
	@echo >> $(@)
	@echo "ESMF_DEP_FRONT     = ipeCap"                        >> $(@)
	@echo "ESMF_DEP_INCPATH   = $(includedir)"                 >> $(@)
	@echo "ESMF_DEP_CMPL_OBJS = "                              >> $(@)
	@echo "ESMF_DEP_LINK_OBJS = $(libdir)/libipm.a"            >> $(@)
	@echo "ESMF_DEP_SHRD_PATH = $(IPE_DEP_SHRD_PATH)"          >> $(@)
	@echo "ESMF_DEP_SHRD_LIBS = $(IPE_DEP_SHRD_LIBS)"          >> $(@)
	@echo
	@echo "Finished generating ESMF self-describing build dependency makefile fragment:" $(@)
	@echo

uninstall-local:
	-rm -f $(datarootdir)/ipe.mk

CLEANFILES = ipe.mk
endif

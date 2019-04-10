#!/bin/sh

#=====================================================================
# SystemName is the filename found under ./compilation_new/arch
#=====================================================================
SystemName = muir_xeon

#=====================================================================
# Paths and internal libraries setup - This unlikely needs to be changed
#=====================================================================
WORK_DIR = $(PWD)
ODIR=$(WORK_DIR)/bin

#=====================================================================
# No need to change following settings, unless it is necessary
#=====================================================================
include ./compilation_new/arch/$(SystemName)

ifndef EXEC_NAME
EXEC_NAME = mst2
endif
MST2_DIR_LINK = MST-2
MSTLIB = $(ODIR)/libmst.a
MPPLIB = $(ODIR)/libmpp.a
IOLIB  = $(ODIR)/iolib.a
IODIR  = $(WORK_DIR)/iolib 
SLULIB = $(ODIR)/slulib.a
ACCLIB = $(ODIR)/libMatAcc.a

ifeq "$(SuperLU)" "0"
   SLULIB =
   DSUPERLULIB =
endif
 
ifeq "$(Acceleration)" "1"
   ADDLIBS += $(ACCLIB)
endif
 
PREPROC_FFLAGS = $(FPPDEFS) $(FPPFLAGS)
PREPROC_CFLAGS = $(CPPDEFS) $(CPPFLAGS)
CFLAGS += $(PREPROC_CFLAGS) -c
FFLAGS += -c
F77FLAGS += -c
CUDA_OPT += -c

#=====================================================================

all: check-bin check-git-tag main check util gen_link

#=====================================================================
# Check bin directory to make sure it exists
#=====================================================================
check-bin:
	if [ ! -d $(ODIR) ]; then echo $(ODIR) "folder does not exist and I am creating one ..."; mkdir $(ODIR); fi

check-git-tag:
	[[ `(git tag | tail -n 1 | wc -m)` -eq 0 ]] && echo '"Develop"' > $(ODIR)/git_version || git tag | tail -n 1 | sed 's/^/\"/' | sed 's/$$/\"/' > $(ODIR)/git_version
	echo 'write(myunit,'\''("# ",t5,a,t30,": ",a)'\'') "Source code version",' > $(ODIR)/prepend
	cat $(ODIR)/prepend $(ODIR)/git_version | ( tr -d '\n'; echo; ) > src/git_version.h
	rm -f $(ODIR)/git_version $(ODIR)/prepend

gen_link:
	cd ../ ; rm -f $(MST2_DIR_LINK); ln -s $(WORK_DIR) $(MST2_DIR_LINK)
	@echo
	@echo '---------------------------------------------------------------------------------------------'
	@echo '*** WARNING Generic link LSMS_2.0 is pointing now to '$(WORK_DIR)' ***'
	@echo '---------------------------------------------------------------------------------------------'
	@echo

main:
	cd $(ODIR); make "FC=$(FC)" "FFLAGS=$(FFLAGS)" "F77FLAGS=$(F77FLAGS)" "MSTLIB=$(MSTLIB)" "PREPROC=$(PREPROC)" "ODIR=$(ODIR)" \
"XLC_I=$(XLC_I)" "ARCHV_LSMS=$(ARCHV_LSMS)" "ARCHV_OPT=$(ARCHV_OPT)" "PREPROC_FFLAGS=$(PREPROC_FFLAGS)" -f ../lib/Makefile
	if test "$(Acceleration)" = "1"; then cd $(ODIR); \
make "FC=$(FC)" "FFLAGS=$(FFLAGS)" "F77FLAGS=$(F77FLAGS)" "PREPROC=$(PREPROC)" "CC=$(CC)" "CXX=$(CXX)" "CUDA_CXX=$(CUDA_CXX)" "CUDA_OPT=$(CUDA_OPT)" "CFLAGS=$(CFLAGS)" \
"XLC_I=$(XLC_I)" "ODIR=$(ODIR)" "ARCHV_LSMS=$(ARCHV_LSMS)" "ARCHV_OPT=$(ARCHV_OPT)" "ACCLIB=$(ACCLIB)" "PREPROC_FFLAGS=$(PREPROC_FFLAGS)" -f ../Accelerator/Makefile; fi
	cd $(ODIR); make "FC=$(FC)" "FFLAGS=$(FFLAGS)" "F77FLAGS=$(F77FLAGS)" "IOLIB=$(IOLIB)" "PREPROC=$(PREPROC)" "PREPROC_FFLAGS=$(PREPROC_FFLAGS)" \
"CC=$(CC)" "CFLAGS=$(CFLAGS)" "ODIR=$(ODIR)" "XLC_I=$(XLC_I)" "ARCHV_LSMS=$(ARCHV_LSMS)" "ARCHV_OPT=$(ARCHV_OPT)" -f ../iolib/Makefile
	cd $(ODIR); make "FC=$(FC)" "FFLAGS=$(FFLAGS)" "MPPLIB=$(MPPLIB)" "PREPROC=$(PREPROC)" "PREPROC_FFLAGS=$(PREPROC_FFLAGS)" \
"ODIR=$(ODIR)" "XLC_I=$(XLC_I)" "ARCHV_LSMS=$(ARCHV_LSMS)" "ARCHV_OPT=$(ARCHV_OPT)" -f ../plib/Makefile
	if test "$(SuperLU)" = "1"; then cd $(ODIR); \
make "FC=$(FC)" "FFLAGSSLU=$(FFLAGS)" "SLULIB=$(SLULIB)" "MPI2INCLUDE_PATH=$(MPI2INCLUDE_PATH)" "SLUPATH=$(SLUPATH)" "PREPROC=$(PREPROC)" \
"PREPROC_FFLAGS=$(PREPROC_FFLAGS)" "ODIR=$(ODIR)" "ARCHV_LSMS=$(ARCHV_LSMS)" "ARCHV_OPT=$(ARCHV_OPT)" "XLC_I=$(XLC_I)" -f ../slulib/Makefile; fi
	cd $(ODIR); make "FC=$(FC)" "FFLAGS=$(FFLAGS)" "F77FLAGS=$(F77FLAGS)" "MSTLIB=$(MSTLIB)" "PREPROC=$(PREPROC)" "PREPROC_FFLAGS=$(PREPROC_FFLAGS)" \
"CC=$(CC)" "CFLAGS=$(CFLAGS)" "IOLIB=$(IOLIB)" "DSUPERLULIB=$(DSUPERLULIB)" "MPPLIB=$(MPPLIB)" "SLULIB=$(SLULIB)" "SLUPATH=$(SLUPATH)" \
"ODIR=$(ODIR)" "ADDLIBS=$(ADDLIBS)" "XLC_I=$(XLC_I)" "LD=$(LD)" "Use_FFTW=$(Use_FFTW)" "FFTW_INC=$(FFTW_INC)" "EXEC_NAME=$(EXEC_NAME)" -f ../src/Makefile
	if test -d $(HOME)/bin; then echo "copy the executable $(EXEC_NAME) to $(HOME)/bin"; else mkdir $(HOME)/bin; fi
	cp $(ODIR)/$(EXEC_NAME) $(HOME)/bin

check: main
	cd $(ODIR); make "FC=$(FC)" "FFLAGS=$(FFLAGS)" "F77FLAGS=$(F77FLAGS)" "MSTLIB=$(MSTLIB)" "IOLIB=$(IOLIB)" "MPPLIB=$(MPPLIB)" \
"PREPROC=$(PREPROC)" "ODIR=$(ODIR)" "ADDLIBS=$(ADDLIBS)" "XLC_I=$(XLC_I)" \
"PREPROC_FFLAGS=$(PREPROC_FFLAGS)" "Use_FFTW=$(Use_FFTW)" "FFTW_INC=$(FFTW_INC)" "LD=$(LD)" -f ../driver/Makefile

util: check
	cd $(ODIR); make "FC=$(FC)" "FFLAGS=$(FFLAGS)" "MPICC=$(MPICC)" "F77FLAGS=$(F77FLAGS)" "MSTLIB=$(MSTLIB)" "PREPROC=$(PREPROC)" "PREPROC_FFLAGS=$(PREPROC_FFLAGS)" \
"CC=$(CC)" "CFLAGS=$(CFLAGS)" "IOLIB=$(IOLIB)" "MPPLIB=$(MPPLIB)" "ODIR=$(ODIR)" "LD=$(LD)" "ADDLIBS=$(ADDLIBS)" "XLC_I=$(XLC_I)" -f ../util/Makefile

clear:
	cd $(ODIR); make "ODIR=$(ODIR)" "MSTLIB=$(MSTLIB)" "EXEC_NAME=$(EXEC_NAME)" clear -f ../lib/Makefile
	cd $(ODIR); make "ODIR=$(ODIR)" "IOLIB=$(IOLIB)" "EXEC_NAME=$(EXEC_NAME)" clear -f ../lib/Makefile
	cd $(ODIR); make "ODIR=$(ODIR)" "MPPLIB=$(MPPLIB)" "EXEC_NAME=$(EXEC_NAME)" clear -f ../plib/Makefile 
	cd $(ODIR); make "ODIR=$(ODIR)" "EXEC_NAME=$(EXEC_NAME)" clear -f ../src/Makefile

clear_slulib:
	cd $(ODIR); make "ODIR=$(ODIR)" "SLULIB=$(SLULIB)" clear -f ../slulib/Makefile

clear_util:
	cd $(ODIR); make "ODIR=$(ODIR)" clear -f ../util/Makefile

clean_util:
	cd $(ODIR); make "ODIR=$(ODIR)" clean -f ../util/Makefile

clear_check:
	cd $(ODIR); make "ODIR=$(ODIR)" clear -f ../driver/Makefile

clean_check:
	cd $(ODIR); make "ODIR=$(ODIR)" clean -f ../driver/Makefile

clear_src:
	cd $(ODIR); make "ODIR=$(ODIR)" "EXEC_NAME=$(EXEC_NAME)" clear -f ../src/Makefile

clean:
	rm -f $(ODIR)/*.o $(ODIR)/*.mod $(ODIR)/*.a

deepclean:
	rm -f $(ODIR)/*

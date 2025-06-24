# KLU Demo Makefile

default: all

include ../../SuiteSparse_config/SuiteSparse_config.mk

CLIB = $(LDFLAGS) -L../../lib -lklu -lbtf -lamd -lcolamd -lsuitesparseconfig \
	$(LIB_WITH_PARTITION) $(LDLIBS)

CHOLMOD_LIB = -lcholmod
CHOLMOD = $(CHOLMOD_LIB) $(LAPACK) $(BLAS)

ifneq ($(GPU_CONFIG),)
CHOLMOD += $(CUDART_LIB) $(CUBLAS_LIB) 
endif

I = -I../../include

all: library klu_simple kludemo kluldemo klu_simpler klu_path klu_full

library:
	( cd ../../SuiteSparse_config ; $(MAKE) )
	( cd ../Lib ; $(MAKE) )
	( cd ../../BTF ; $(MAKE) library )
	( cd ../../AMD ; $(MAKE) library )
	( cd ../../COLAMD ; $(MAKE) library )
	( cd ../../CHOLMOD ; $(MAKE) library )
	- ( cd ../../CAMD ; $(MAKE) )
	- ( cd ../../CCOLAMD ; $(MAKE) )
	- ( cd ../.. ; $(MAKE) metis )

purge: distclean

distclean: clean
	- $(RM) kludemo kluldemo klu_simple klu_full klu_path klu_simpler
	- $(RM) -r $(PURGE)

clean:
	- $(RM) -r $(CLEAN)

kludemo: kludemo.c Makefile
	$(CC) $(CF) $(I) kludemo.c -o kludemo $(CLIB) $(CHOLMOD)

kluldemo: kludemo.c Makefile
	$(CC) $(CF) $(I) kluldemo.c -o kluldemo $(CLIB) $(CHOLMOD)

klu_simple: klu_simple.c Makefile
	$(CC) $(CF) $(I) klu_simple.c -o klu_simple $(CLIB)
	- ./klu_simple

klu_simpler: klu_simpler.c Makefile
	$(CC) $(CF) $(I) klu_simpler.c -o klu_simpler $(CLIB)

klu_path: klu_path.c Makefile
	$(CC) $(CF) $(I) klu_path.c -o klu_path $(CLIB)

klu_full: klu_full.c Makefile 
	$(CC) $(CF) $(I) klu_full.c -o klu_full $(CLIB)

#------------------------------------------------------------------------------
# KLU Makefile
#------------------------------------------------------------------------------

SUITESPARSE ?= $(realpath $(CURDIR)/..)
export SUITESPARSE

default: all

include ../SuiteSparse_config/SuiteSparse_config.mk

demos: all

all:
	( cd Lib ; $(MAKE) )

library:
	( cd Lib ; $(MAKE) )

# compile the static libraries only
static:
	( cd Lib    ; $(MAKE) static )

clean:
	( cd Lib ; $(MAKE) clean )

distclean:
	( cd Lib ; $(MAKE) distclean )
	( cd User ; $(MAKE) distclean )

purge: distclean

cov:
	( cd Tcov ; $(MAKE) )

# create PDF documents for the original distribution
docs:
	( cd Doc    ; $(MAKE) )

# install KLU
install:
	( cd Lib ; $(MAKE) install )

# uninstall KLU
uninstall:
	( cd Lib ; $(MAKE) uninstall )

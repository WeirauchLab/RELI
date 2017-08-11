########################################################################
##                                                                    ##
##   Weirauch Lab Makefile for Regulatory Element Locus Intersection  ##
##   Analysis (RELI)                                                  ##
##                                                                    ##
##   Usage:    make              # build 'RELI' binary                ##
##             make test         # run examples                       ##
##             -or- make help                                         ##
##                                                                    ##
##   Author:   ern6xv (Kevin.Ernst@cchmc.org)                         ##
##   Date:     06 July 2015                                           ##
##                                                                    ##
########################################################################
PKGNAME=RELI
# FIXME: this needs to be updated with a pre-commit hook, or we need to pass
# the value of this variable as a preprocessor flag (-D) to g++
PKGVER=0.90

# Enable debug switches by running 'make DEBUG=1' or 'make debug'
# (what this changes: don't add -O [optimize], add -ggdb)
DEBUG=0

# Prepends $(SOURCEDIR) to all the sources and builds a binary in the top
# level of the repo
SOURCEDIR=src
SOURCES=RELI.cpp RELI_impl.cpp
INCLUDES=RELI_impl.h
DATAURL=https://tf.cchmc.org/external/RELI/data.tar.bz2

# Required (third-party) libraries
LIBS=gsl gslcblas

CC=g++
CXXFLAGS=$(CFLAGS) $(LDFLAGS) -std=c++11 -I$(SOURCEDIR)
ifeq ($(DEBUG), 1)
# enable debugging with gdb
CXXFLAGS += -ggdb -Wall -Wno-sign-compare -Wno-parentheses
else
# optimize and disable warnings
CXXFLAGS += -O3 -w
endif

# ANSI terminal colors (see 'man tput') and
# https://linuxtidbits.wordpress.com/2008/08/11/output-color-on-bash-scripts/
# Don't set these if there isn't a $TERM environment variable (e.g., 'make
# clusterbuild')
ifneq ($(strip $(TERM)),)
BOLD=$(shell tput bold)
RED=$(shell tput setaf 1)
BLUE=$(shell tput setaf 4)
UL=$(shell tput sgr 0 1)
GREEN=$(shell tput setaf 2)
RESET=$(shell tput sgr0 )
endif

# To make Eclipse happy; this is the default target for Makefile projects
all: binary

binary: $(PKGNAME)

$(PKGNAME): $(addprefix $(SOURCEDIR)/,$(SOURCES) $(INCLUDES))
	$(CC) $(CXXFLAGS) -o $(PKGNAME) $(addprefix $(SOURCEDIR)/,$(SOURCES)) \
	    $(addprefix -l,$(LIBS)) 

clean:
	-rm -f a.out a.exe *.o $(PKGNAME) $(PKGNAME).exe core.* vgcore.*

test: binary validate-data
	pushd example && ./example_run.sh
	-popd

validate-data: fetch-data
	if [ ! -f .data_validated ]; then \
		if data/validate.sh; then \
			touch .data_validated; \
		else \
			echo >&2;                                                             \
			echo "$(UL)$(RED)ACK!$(RESET)" >&2;                                   \
			echo >&2;                                                             \
			echo "Data validation failed. Try deleting the 'data' directory" >&2; \
			echo "and running 'make test' again." >&2;                            \
			echo >&2; \
		fi; \
	fi

fetch-data:
	test -x data/validate.sh || curl $(DATAURL) | tar xjf -

# Preserve the complete path of this Makefile in case we were called with
# something like 'make -f Makefile.WL'. According to §5.6.3 of the manual, the
# '-f' option is not propagated, so this trickery is necessary.
ME=$(word $(words $(MAKEFILE_LIST)),$(MAKEFILE_LIST))

debug:
	make -f $(ME) DEBUG=1

help:
	@echo
	@echo "  $(UL)$(BOLD)$(BLUE)Building $(PKGNAME) v$(PKGVER)$(RESET)"
	@echo
	@echo "  Try one of these:"
	@echo
	@echo "      $(BOLD)make binary$(RESET) - (default target) build the RELI binary"
	@echo
	@echo "      $(BOLD)make test$(RESET)   - download sample data and perform a test analysis"
	@echo
	@echo "      $(BOLD)make debug$(RESET)  - build a debuggable version of the RELI binary"
	@echo
	@echo "      $(BOLD)make help$(RESET)   - you're looking at it"
	@echo
	@echo
	@echo "  For more help, see https://tf.cchmc.org/s/reli-readme"
	@echo

uninstall:
install:
	@echo
	@echo "$(UL)$(RED)OOF!$(RESET)"
	@echo
	@echo "Not a supported target (yet). Stay tuned."
	@echo

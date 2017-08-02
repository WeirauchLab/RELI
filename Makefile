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
# FIXME: this needs to be updated with a pre-commit hook
PKGVER=0.90

# FIXME: Install data on user's system, possibly in /opt/$(PKGNAME)/share
DESTROOT=/usr/local

# Prepends $(SOURCEDIR) to all the sources and builds a binary in the top
# level of the repo
SOURCEDIR=src
SOURCES=RELI.cpp RELI_impl.cpp
INCLUDES=RELI_impl.h
DATAURL=https://tf.cchmc.org/external/RELI/data.tar.bz2

# Required (third-party) libraries
LIBS=gsl gslcblas

CXXFLAGS=-std=c++11 -O3 -ggdb -w -I$(SOURCEDIR)

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

binary: $(PKGNAME)

$(PKGNAME): $(addprefix $(SOURCEDIR)/,$(SOURCES) $(INCLUDES))
	g++ $(CXXFLAGS) -o $(PKGNAME) $(addprefix $(SOURCEDIR)/,$(SOURCES)) \
	    $(addprefix -l,$(LIBS)) 

clean:
	-rm -f a.out a.exe *.o $(PKGNAME) $(PKGNAME).exe

test: fetch-data
	pushd example && example_run.sh
	-popd

fetch-data: .data_validated
	-rm -f .data_validated
	curl $(DATAURL) | tar xjf -
	data/validate.sh && touch .data_validated

help:
	@echo
	@echo "  $(UL)$(BOLD)$(BLUE)Building $(PKGNAME) v$(PKGVER)$(RESET)"
	@echo
	@echo "  Try one of these:"
	@echo
	@echo "      $(BOLD)make binary$(RESET) - to build the RELI binary"
	@echo "                                   (this is the default target)"
	@echo
	@echo "      $(BOLD)make test$(RESET)   - to download sample data and"
	@echo "                                   perform a test analysis"
	@echo
	@echo "      $(BOLD)make help$(RESET)   - you're looking at it"
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

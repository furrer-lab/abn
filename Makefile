# Determine package name and version from DESCRIPTION file
PKG_VERSION=$(shell grep -i ^version DESCRIPTION | cut -d : -d \  -f 2)
PKG_NAME=$(shell grep -i ^package DESCRIPTION | cut -d : -d \  -f 2)

# Name of built package
PKG_TAR=$(PKG_NAME)_$(PKG_VERSION).tar.gz

all: install configure clean

# Install package
install:
	cd .. && R CMD INSTALL $(PKG_NAME)

configure: configure.ac
	autoconf ./configure.ac > ./configure
	chmod +x ./configure

clean:
	./cleanup
	rm -rf docs/
	rm -rf ./revdep

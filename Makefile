# Makefile for simple installation of python project
#
# Authors:
# Geir Arne Hjelle <geir.arne.hjelle@kartverket.no>.
#
# $Revision: 10610 $
# $Date: 2016-05-30 17:14:53 +0200 (Mon, 30 May 2016) $
# $LastChangedBy: hjegei $

# Programs and directories
F2PY = f2py3.5
F2PYEXTENSION = .cpython-35m-x86_64-linux-gnu.so

EXTDIR = $(CURDIR)
SOFADIR = $(CURDIR)/sofa/src

# Define phony targets (targets that are not files)
.PHONY: develop install external sofa

# Install in developer mode (no need to reinstall after changing source)
develop:
	python3 setup.py develop --user

install:
	python3 setup.py install --user

# External libraries
external:	sofa

# SOFA
sofa:	$(EXTDIR)/sofa$(F2PYEXTENSION)

$(EXTDIR)/sofa$(F2PYEXTENSION):	$(shell find $(SOFADIR) -type f)
	( cd $(EXTDIR) && $(F2PY) -c $(SOFADIR)/sofa.pyf $(SOFADIR)/*.for )


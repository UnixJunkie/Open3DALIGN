# Makefile.am
# 
# is part of
# 
# Open3DALIGN
# -----------
# 
# An open-source software aimed at unsupervised molecular alignment
# 
# Copyright (C) 2010-2012 Paolo Tosco, Thomas Balle
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
# 
# For further information, please contact:
# 
# Paolo Tosco, PhD
# Dipartimento di Scienza e Tecnologia del Farmaco
# Universita' degli Studi di Torino
# Via Pietro Giuria, 9
# 10125 Torino (Italy)
# Phone:  +39 011 670 7680
# Mobile: +39 348 553 7206
# Fax:    +39 011 670 7687
# E-mail: paolo.tosco@unito.it


lib_LTLIBRARIES = libo3a.la
libo3a_la_SOURCES = \
align.c \
compare.c \
conf.c \
filter.c \
qmd.c \
superpose_conf.c \
tinker.c \
include/align.h \
include/basis_set.h \
include/error_messages.h \
include/extern.h \
include/ff_parm.h \
include/nice_windows.h \
include/o3header.h \
include/proc_env.h \
include/prog_exe_info.h \
include/safe_rint.h

bin_PROGRAMS = \
open3dalign
open3dalign_SOURCES = \
main.c \
parse_o3a_input.c \
include/align.h \
include/basis_set.h \
include/error_messages.h \
include/extern.h \
include/ff_parm.h \
include/keywords.h \
include/nice_windows.h \
include/o3header.h \
include/proc_env.h \
include/prog_exe_info.h \
include/safe_rint.h

open3dalign_LDADD = \
-lo3q \
-lo3a \
$(WINDOWS_RES) \
$(LIBLAPACK_ATLAS) \
$(LIBLAPACK_F77) \
$(LIBATLAS) \
$(LIBMKL) \
$(LIBSUNPERF) \
$(LIBACCELERATE) \
$(LIBEDIT)

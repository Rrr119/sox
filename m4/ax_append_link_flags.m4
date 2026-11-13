# ===========================================================================
# https://www.gnu.org/software/autoconf-archive/ax_append_link_flags.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_APPEND_LINK_FLAGS([FLAG ... [,VAR [,ACTION-IF-SUCCESS [,ACTION-IF-FAILURE]]])
#
# DESCRIPTION
#
#   AX_APPEND_LINK_FLAGS tests each linker FLAG using an AX_CHECK_LINK_FLAG
#   style test.  Each working FLAG is appended to VAR (default LDFLAGS).
#
# LICENSE
#
#   Copyright (c) 2013 Karlson2k (Evgeny Grin) <k2k@narod.ru>
#
#   Copying and distribution of this file, with or without modification,
#   are permitted in any medium without royalty provided the copyright
#   notice and this notice are preserved. This file is offered as-is,
#   without any warranty.

#serial 8

AC_DEFUN([AX_APPEND_LINK_FLAGS],
[
  m4_ifval([$2],
    [ax_append_link_flags_var="$2"],
    [ax_append_link_flags_var="LDFLAGS"])
  ax_append_link_flags_flags_success="true"
  AC_REQUIRE([AC_PROG_CC])
  for ax_append_link_flags_flag in $1
  do
    AX_CHECK_LINK_FLAG([$ax_append_link_flags_flag],
      [AS_VAR_APPEND([$ax_append_link_flags_var], [" $ax_append_link_flags_flag"])],
      [ax_append_link_flags_flags_success="false"])
  done
  AS_VAR_IF([ax_append_link_flags_flags_success], [true],
    [m4_default([$3], [:])],
    [m4_default([$4], [:])])
])

dnl Minimal fallback for environments without autoconf-archive.
m4_ifndef([AX_APPEND_LINK_FLAGS], [
AC_DEFUN([AX_APPEND_LINK_FLAGS], [
  for ax_flag in $1; do
    AC_MSG_CHECKING([whether the linker supports $ax_flag])
    ax_save_LDFLAGS="$LDFLAGS"
    LDFLAGS="$LDFLAGS $ax_flag"
    AC_LINK_IFELSE([AC_LANG_PROGRAM([[]], [[]])],
      [AC_MSG_RESULT([yes])],
      [AC_MSG_RESULT([no]); LDFLAGS="$ax_save_LDFLAGS"])
  done
])
])

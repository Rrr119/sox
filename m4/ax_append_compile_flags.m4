dnl Minimal fallback for environments without autoconf-archive.
m4_ifndef([AX_APPEND_COMPILE_FLAGS], [
AC_DEFUN([AX_APPEND_COMPILE_FLAGS], [
  for ax_flag in $1; do
    AC_MSG_CHECKING([whether $CC supports $ax_flag])
    ax_save_CFLAGS="$CFLAGS"
    CFLAGS="$CFLAGS $ax_flag"
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[]], [[]])],
      [AC_MSG_RESULT([yes])],
      [AC_MSG_RESULT([no]); CFLAGS="$ax_save_CFLAGS"])
  done
])
])

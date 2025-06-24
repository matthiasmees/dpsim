
#ifndef AMD_EXPORT_H
#define AMD_EXPORT_H

#ifdef AMD_STATIC_DEFINE
#  define AMD_EXPORT
#  define AMD_NO_EXPORT
#else
#  ifndef AMD_EXPORT
#    ifdef amd_EXPORTS
        /* We are building this library */
#      define AMD_EXPORT 
#    else
        /* We are using this library */
#      define AMD_EXPORT 
#    endif
#  endif

#  ifndef AMD_NO_EXPORT
#    define AMD_NO_EXPORT 
#  endif
#endif

#ifndef AMD_DEPRECATED
#  define AMD_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef AMD_DEPRECATED_EXPORT
#  define AMD_DEPRECATED_EXPORT AMD_EXPORT AMD_DEPRECATED
#endif

#ifndef AMD_DEPRECATED_NO_EXPORT
#  define AMD_DEPRECATED_NO_EXPORT AMD_NO_EXPORT AMD_DEPRECATED
#endif

/* NOLINTNEXTLINE(readability-avoid-unconditional-preprocessor-if) */
#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef AMD_NO_DEPRECATED
#    define AMD_NO_DEPRECATED
#  endif
#endif

#endif /* AMD_EXPORT_H */

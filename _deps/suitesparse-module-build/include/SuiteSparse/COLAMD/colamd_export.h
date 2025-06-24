
#ifndef COLAMD_EXPORT_H
#define COLAMD_EXPORT_H

#ifdef COLAMD_STATIC_DEFINE
#  define COLAMD_EXPORT
#  define COLAMD_NO_EXPORT
#else
#  ifndef COLAMD_EXPORT
#    ifdef colamd_EXPORTS
        /* We are building this library */
#      define COLAMD_EXPORT 
#    else
        /* We are using this library */
#      define COLAMD_EXPORT 
#    endif
#  endif

#  ifndef COLAMD_NO_EXPORT
#    define COLAMD_NO_EXPORT 
#  endif
#endif

#ifndef COLAMD_DEPRECATED
#  define COLAMD_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef COLAMD_DEPRECATED_EXPORT
#  define COLAMD_DEPRECATED_EXPORT COLAMD_EXPORT COLAMD_DEPRECATED
#endif

#ifndef COLAMD_DEPRECATED_NO_EXPORT
#  define COLAMD_DEPRECATED_NO_EXPORT COLAMD_NO_EXPORT COLAMD_DEPRECATED
#endif

/* NOLINTNEXTLINE(readability-avoid-unconditional-preprocessor-if) */
#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef COLAMD_NO_DEPRECATED
#    define COLAMD_NO_DEPRECATED
#  endif
#endif

#endif /* COLAMD_EXPORT_H */

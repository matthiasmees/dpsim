
#ifndef KLU_EXPORT_H
#define KLU_EXPORT_H

#ifdef KLU_STATIC_DEFINE
#  define KLU_EXPORT
#  define KLU_NO_EXPORT
#else
#  ifndef KLU_EXPORT
#    ifdef klu_EXPORTS
        /* We are building this library */
#      define KLU_EXPORT 
#    else
        /* We are using this library */
#      define KLU_EXPORT 
#    endif
#  endif

#  ifndef KLU_NO_EXPORT
#    define KLU_NO_EXPORT 
#  endif
#endif

#ifndef KLU_DEPRECATED
#  define KLU_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef KLU_DEPRECATED_EXPORT
#  define KLU_DEPRECATED_EXPORT KLU_EXPORT KLU_DEPRECATED
#endif

#ifndef KLU_DEPRECATED_NO_EXPORT
#  define KLU_DEPRECATED_NO_EXPORT KLU_NO_EXPORT KLU_DEPRECATED
#endif

/* NOLINTNEXTLINE(readability-avoid-unconditional-preprocessor-if) */
#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef KLU_NO_DEPRECATED
#    define KLU_NO_DEPRECATED
#  endif
#endif

#endif /* KLU_EXPORT_H */

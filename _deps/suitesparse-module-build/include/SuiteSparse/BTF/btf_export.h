
#ifndef BTF_EXPORT_H
#define BTF_EXPORT_H

#ifdef BTF_STATIC_DEFINE
#  define BTF_EXPORT
#  define BTF_NO_EXPORT
#else
#  ifndef BTF_EXPORT
#    ifdef btf_EXPORTS
        /* We are building this library */
#      define BTF_EXPORT 
#    else
        /* We are using this library */
#      define BTF_EXPORT 
#    endif
#  endif

#  ifndef BTF_NO_EXPORT
#    define BTF_NO_EXPORT 
#  endif
#endif

#ifndef BTF_DEPRECATED
#  define BTF_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef BTF_DEPRECATED_EXPORT
#  define BTF_DEPRECATED_EXPORT BTF_EXPORT BTF_DEPRECATED
#endif

#ifndef BTF_DEPRECATED_NO_EXPORT
#  define BTF_DEPRECATED_NO_EXPORT BTF_NO_EXPORT BTF_DEPRECATED
#endif

/* NOLINTNEXTLINE(readability-avoid-unconditional-preprocessor-if) */
#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef BTF_NO_DEPRECATED
#    define BTF_NO_DEPRECATED
#  endif
#endif

#endif /* BTF_EXPORT_H */

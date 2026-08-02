/* Native routine registration.
 *
 * Registering the entry points lets R check the number of arguments at call
 * time instead of letting a mismatch corrupt memory silently. That matters
 * here: several of these routines take 40-60 arguments, and the package has
 * already been bitten once by an argument list drifting out of step with the
 * R side. R_useDynamicSymbols(FALSE) also stops R searching this object for
 * unrelated symbols.
 */

#include <stdlib.h>
#include <R_ext/Rdynload.h>
#include <R_ext/RS.h>

/* Fortran entry points called from R. */
void F77_NAME(areq4)(void *, void *, void *, void *, void *, void *, void *);

void F77_NAME(ggfst)(void *, void *, void *, void *, void *, void *, void *,
                     void *, void *, void *, void *, void *, void *);

void F77_NAME(mcmcgld)(void *, void *, void *, void *, void *, void *, void *,
                       void *, void *, void *, void *, void *, void *, void *,
                       void *, void *, void *, void *, void *, void *, void *,
                       void *, void *, void *, void *, void *, void *, void *,
                       void *, void *, void *, void *, void *, void *, void *,
                       void *, void *, void *, void *, void *, void *, void *,
                       void *, void *, void *, void *, void *, void *, void *,
                       void *, void *, void *, void *, void *, void *, void *,
                       void *, void *, void *, void *, void *, void *, void *);

void F77_NAME(mcmchz)(void *, void *, void *, void *, void *, void *, void *,
                      void *, void *, void *, void *, void *, void *, void *,
                      void *, void *, void *, void *, void *, void *, void *,
                      void *, void *, void *, void *, void *, void *, void *,
                      void *, void *, void *, void *, void *, void *, void *,
                      void *, void *, void *, void *, void *);

void F77_NAME(postprocesschain2)(void *, void *, void *, void *, void *,
                                 void *, void *, void *, void *, void *,
                                 void *, void *, void *, void *, void *,
                                 void *, void *, void *, void *, void *,
                                 void *, void *, void *, void *, void *,
                                 void *, void *, void *, void *, void *,
                                 void *, void *, void *, void *, void *,
                                 void *, void *, void *, void *, void *,
                                 void *, void *, void *, void *, void *,
                                 void *, void *);

void F77_NAME(pppmindiv2)(void *, void *, void *, void *, void *, void *,
                          void *, void *, void *, void *, void *, void *,
                          void *, void *, void *, void *, void *, void *,
                          void *, void *, void *);

void F77_NAME(tessdyn)(void *, void *, void *, void *, void *, void *, void *,
                       void *, void *, void *, void *, void *, void *, void *,
                       void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"areq4",             (DL_FUNC) &F77_NAME(areq4),              7},
    {"ggfst",             (DL_FUNC) &F77_NAME(ggfst),             13},
    {"mcmcgld",           (DL_FUNC) &F77_NAME(mcmcgld),           63},
    {"mcmchz",            (DL_FUNC) &F77_NAME(mcmchz),            40},
    {"postprocesschain2", (DL_FUNC) &F77_NAME(postprocesschain2), 47},
    {"pppmindiv2",        (DL_FUNC) &F77_NAME(pppmindiv2),        21},
    {"tessdyn",           (DL_FUNC) &F77_NAME(tessdyn),           21},
    {NULL, NULL, 0}
};

void R_init_Geneland(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

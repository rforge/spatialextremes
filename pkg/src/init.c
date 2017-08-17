#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void bell(void *, void *);
extern void brownresnickstderr(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void circemb(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void computeWeightsBR(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void computeWeightsExtt(void *, void *, void *, void *, void *, void *, void *, void *);
extern void computeWeightsSC(void *, void *, void *, void *, void *, void *, void *);
extern void concProbKendall(void *, void *, void *, void *, void *, void *);
extern void condsimbrown(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void condsimextt(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void condsimschlather(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void DIC(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void direct(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void distance(void *, void *, void *, void *, void *);
extern void empiricalBootConcProb(void *, void *, void *, void *, void *);
extern void empiricalConcProb(void *, void *, void *, void *, void *, void *);
extern void extCoeffSmith(void *, void *, void *, void *);
extern void extCoeffST(void *, void *, void *, void *, void *, void *);
extern void extremaltstderr(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void fitbrcovariance(void *, void *, void *, void *, void *, void *, void *);
extern void fitcovariance(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void fitcovmat2d(void *, void *, void *, void *, void *, void *, void *, void *);
extern void fitcovmat3d(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void fitgcovariance(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void fiticovariance(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void fittcovariance(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void geomgaussstderr(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void getStartingPartitionBR(void *, void *, void *, void *, void *, void *);
extern void getStartingPartitionExtt(void *, void *, void *, void *, void *);
extern void getStartingPartitionSC(void *, void *, void *, void *);
extern void gev(void *, void *, void *, void *, void *, void *);
extern void gevlik(void *, void *, void *, void *, void *, void *);
extern void gibbsForPartBR(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gibbsForPartExtt(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gibbsForPartSC(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gpdlik(void *, void *, void *, void *, void *, void *);
extern void latentgev(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void listAllPartOfASet(void *, void *, void *, void *);
extern void lmadogram(void *, void *, void *, void *, void *, void *);
extern void madogram(void *, void *, void *, void *);
extern void maxLinDsgnMat(void *, void *, void *, void *, void *, void *, void *, void *);
extern void maxLinear(void *, void *, void *, void *, void *, void *, void *);
extern void rbrowndirect(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void rbrownexact(void *, void *, void *, void *, void *, void *, void *, void *);
extern void rcondMaxLin(void *, void *, void *, void *, void *, void *);
extern void rextremaltcirc(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void rextremaltdirect(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void rextremalttbm(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void rgeomcirc(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void rgeomdirect(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void rgeomtbm(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void rotation(void *, void *, void *, void *, void *, void *);
extern void rschlathercirc(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void rschlatherdirect(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void rschlathertbm(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void rsmith1d(void *, void *, void *, void *, void *, void *, void *);
extern void rsmith2d(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void schlatherindstderr(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void schlatherstderr(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void skriging(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void smithstderr(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void smithstderr3d(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void spatgevstderr(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void tbm(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void vandercorput(void *, void *);
extern void variogram(void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"bell",                     (DL_FUNC) &bell,                      2},
    {"brownresnickstderr",       (DL_FUNC) &brownresnickstderr,       29},
    {"circemb",                  (DL_FUNC) &circemb,                  10},
    {"computeWeightsBR",         (DL_FUNC) &computeWeightsBR,         11},
    {"computeWeightsExtt",       (DL_FUNC) &computeWeightsExtt,        8},
    {"computeWeightsSC",         (DL_FUNC) &computeWeightsSC,          7},
    {"concProbKendall",          (DL_FUNC) &concProbKendall,           6},
    {"condsimbrown",             (DL_FUNC) &condsimbrown,             17},
    {"condsimextt",              (DL_FUNC) &condsimextt,              11},
    {"condsimschlather",         (DL_FUNC) &condsimschlather,         10},
    {"DIC",                      (DL_FUNC) &DIC,                      13},
    {"direct",                   (DL_FUNC) &direct,                   11},
    {"distance",                 (DL_FUNC) &distance,                  5},
    {"empiricalBootConcProb",    (DL_FUNC) &empiricalBootConcProb,     5},
    {"empiricalConcProb",        (DL_FUNC) &empiricalConcProb,         6},
    {"extCoeffSmith",            (DL_FUNC) &extCoeffSmith,             4},
    {"extCoeffST",               (DL_FUNC) &extCoeffST,                6},
    {"extremaltstderr",          (DL_FUNC) &extremaltstderr,          33},
    {"fitbrcovariance",          (DL_FUNC) &fitbrcovariance,           7},
    {"fitcovariance",            (DL_FUNC) &fitcovariance,            11},
    {"fitcovmat2d",              (DL_FUNC) &fitcovmat2d,               8},
    {"fitcovmat3d",              (DL_FUNC) &fitcovmat3d,              11},
    {"fitgcovariance",           (DL_FUNC) &fitgcovariance,           13},
    {"fiticovariance",           (DL_FUNC) &fiticovariance,           12},
    {"fittcovariance",           (DL_FUNC) &fittcovariance,           12},
    {"geomgaussstderr",          (DL_FUNC) &geomgaussstderr,          33},
    {"getStartingPartitionBR",   (DL_FUNC) &getStartingPartitionBR,    6},
    {"getStartingPartitionExtt", (DL_FUNC) &getStartingPartitionExtt,  5},
    {"getStartingPartitionSC",   (DL_FUNC) &getStartingPartitionSC,    4},
    {"gev",                      (DL_FUNC) &gev,                       6},
    {"gevlik",                   (DL_FUNC) &gevlik,                    6},
    {"gibbsForPartBR",           (DL_FUNC) &gibbsForPartBR,           13},
    {"gibbsForPartExtt",         (DL_FUNC) &gibbsForPartExtt,         10},
    {"gibbsForPartSC",           (DL_FUNC) &gibbsForPartSC,            9},
    {"gpdlik",                   (DL_FUNC) &gpdlik,                    6},
    {"latentgev",                (DL_FUNC) &latentgev,                29},
    {"listAllPartOfASet",        (DL_FUNC) &listAllPartOfASet,         4},
    {"lmadogram",                (DL_FUNC) &lmadogram,                 6},
    {"madogram",                 (DL_FUNC) &madogram,                  4},
    {"maxLinDsgnMat",            (DL_FUNC) &maxLinDsgnMat,             8},
    {"maxLinear",                (DL_FUNC) &maxLinear,                 7},
    {"rbrowndirect",             (DL_FUNC) &rbrowndirect,             15},
    {"rbrownexact",              (DL_FUNC) &rbrownexact,               8},
    {"rcondMaxLin",              (DL_FUNC) &rcondMaxLin,               6},
    {"rextremaltcirc",           (DL_FUNC) &rextremaltcirc,           11},
    {"rextremaltdirect",         (DL_FUNC) &rextremaltdirect,         12},
    {"rextremalttbm",            (DL_FUNC) &rextremalttbm,            13},
    {"rgeomcirc",                (DL_FUNC) &rgeomcirc,                11},
    {"rgeomdirect",              (DL_FUNC) &rgeomdirect,              12},
    {"rgeomtbm",                 (DL_FUNC) &rgeomtbm,                 13},
    {"rotation",                 (DL_FUNC) &rotation,                  6},
    {"rschlathercirc",           (DL_FUNC) &rschlathercirc,           10},
    {"rschlatherdirect",         (DL_FUNC) &rschlatherdirect,         11},
    {"rschlathertbm",            (DL_FUNC) &rschlathertbm,            12},
    {"rsmith1d",                 (DL_FUNC) &rsmith1d,                  7},
    {"rsmith2d",                 (DL_FUNC) &rsmith2d,                 10},
    {"schlatherindstderr",       (DL_FUNC) &schlatherindstderr,       33},
    {"schlatherstderr",          (DL_FUNC) &schlatherstderr,          32},
    {"skriging",                 (DL_FUNC) &skriging,                 13},
    {"smithstderr",              (DL_FUNC) &smithstderr,              30},
    {"smithstderr3d",            (DL_FUNC) &smithstderr3d,            33},
    {"spatgevstderr",            (DL_FUNC) &spatgevstderr,            24},
    {"tbm",                      (DL_FUNC) &tbm,                      12},
    {"vandercorput",             (DL_FUNC) &vandercorput,              2},
    {"variogram",                (DL_FUNC) &variogram,                 4},
    {NULL, NULL, 0}
};

void R_init_SpatialExtremes(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

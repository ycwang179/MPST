// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// HQblkdiag
NumericMatrix HQblkdiag(NumericMatrix C, NumericVector cnt);
RcppExport SEXP _MPST_HQblkdiag(SEXP CSEXP, SEXP cntSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type C(CSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cnt(cntSEXP);
    rcpp_result_gen = Rcpp::wrap(HQblkdiag(C, cnt));
    return rcpp_result_gen;
END_RCPP
}
// HQbary
NumericMatrix HQbary(NumericMatrix V, NumericMatrix Tr, NumericVector xx, NumericVector yy);
RcppExport SEXP _MPST_HQbary(SEXP VSEXP, SEXP TrSEXP, SEXP xxSEXP, SEXP yySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type V(VSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Tr(TrSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type yy(yySEXP);
    rcpp_result_gen = Rcpp::wrap(HQbary(V, Tr, xx, yy));
    return rcpp_result_gen;
END_RCPP
}
// HQgetInd
List HQgetInd(NumericMatrix V, NumericMatrix Tr, NumericVector xx, NumericVector yy);
RcppExport SEXP _MPST_HQgetInd(SEXP VSEXP, SEXP TrSEXP, SEXP xxSEXP, SEXP yySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type V(VSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Tr(TrSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type yy(yySEXP);
    rcpp_result_gen = Rcpp::wrap(HQgetInd(V, Tr, xx, yy));
    return rcpp_result_gen;
END_RCPP
}
// mtxrbind
NumericMatrix mtxrbind(NumericMatrix a, NumericMatrix b);
RcppExport SEXP _MPST_mtxrbind(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(mtxrbind(a, b));
    return rcpp_result_gen;
END_RCPP
}
// mtxcbind
NumericMatrix mtxcbind(NumericMatrix a, NumericMatrix b);
RcppExport SEXP _MPST_mtxcbind(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(mtxcbind(a, b));
    return rcpp_result_gen;
END_RCPP
}
// vbind
NumericVector vbind(NumericVector a, NumericVector b);
RcppExport SEXP _MPST_vbind(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(vbind(a, b));
    return rcpp_result_gen;
END_RCPP
}
// BSpline2
List BSpline2(NumericMatrix V, NumericMatrix Tr, double d, int r, NumericVector xx, NumericVector yy);
RcppExport SEXP _MPST_BSpline2(SEXP VSEXP, SEXP TrSEXP, SEXP dSEXP, SEXP rSEXP, SEXP xxSEXP, SEXP yySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type V(VSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Tr(TrSEXP);
    Rcpp::traits::input_parameter< double >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type yy(yySEXP);
    rcpp_result_gen = Rcpp::wrap(BSpline2(V, Tr, d, r, xx, yy));
    return rcpp_result_gen;
END_RCPP
}
// BSpline
List BSpline(NumericMatrix V, NumericMatrix Tr, double d, int r, NumericVector xx, NumericVector yy);
RcppExport SEXP _MPST_BSpline(SEXP VSEXP, SEXP TrSEXP, SEXP dSEXP, SEXP rSEXP, SEXP xxSEXP, SEXP yySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type V(VSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Tr(TrSEXP);
    Rcpp::traits::input_parameter< double >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type yy(yySEXP);
    rcpp_result_gen = Rcpp::wrap(BSpline(V, Tr, d, r, xx, yy));
    return rcpp_result_gen;
END_RCPP
}
// CrIndices
List CrIndices(int d, int r);
RcppExport SEXP _MPST_CrIndices(SEXP dSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(CrIndices(d, r));
    return rcpp_result_gen;
END_RCPP
}
// CrCellArrays
List CrCellArrays(int d, int r);
RcppExport SEXP _MPST_CrCellArrays(SEXP dSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(CrCellArrays(d, r));
    return rcpp_result_gen;
END_RCPP
}
// CrArrays
List CrArrays(int d, int r);
RcppExport SEXP _MPST_CrArrays(SEXP dSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(CrArrays(d, r));
    return rcpp_result_gen;
END_RCPP
}
// degree
double degree(double m);
RcppExport SEXP _MPST_degree(SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(degree(m));
    return rcpp_result_gen;
END_RCPP
}
// indices
List indices(int d);
RcppExport SEXP _MPST_indices(SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(indices(d));
    return rcpp_result_gen;
END_RCPP
}
// bary
List bary(NumericVector V1, NumericVector V2, NumericVector V3, NumericVector X, NumericVector Y);
RcppExport SEXP _MPST_bary(SEXP V1SEXP, SEXP V2SEXP, SEXP V3SEXP, SEXP XSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type V1(V1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type V2(V2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type V3(V3SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(bary(V1, V2, V3, X, Y));
    return rcpp_result_gen;
END_RCPP
}
// tcord
NumericVector tcord(NumericVector V1, NumericVector V2, NumericVector V3, NumericVector u);
RcppExport SEXP _MPST_tcord(SEXP V1SEXP, SEXP V2SEXP, SEXP V3SEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type V1(V1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type V2(V2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type V3(V3SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(tcord(V1, V2, V3, u));
    return rcpp_result_gen;
END_RCPP
}
// locate
NumericVector locate(NumericVector I1, NumericVector J1, NumericVector K1, NumericVector I, NumericVector J, NumericVector K);
RcppExport SEXP _MPST_locate(SEXP I1SEXP, SEXP J1SEXP, SEXP K1SEXP, SEXP ISEXP, SEXP JSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type I1(I1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type J1(J1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type K1(K1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type I(ISEXP);
    Rcpp::traits::input_parameter< NumericVector >::type J(JSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(locate(I1, J1, K1, I, J, K));
    return rcpp_result_gen;
END_RCPP
}
// dirder
NumericMatrix dirder(NumericMatrix Bin, double lam1, double lam2, double lam3);
RcppExport SEXP _MPST_dirder(SEXP BinSEXP, SEXP lam1SEXP, SEXP lam2SEXP, SEXP lam3SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Bin(BinSEXP);
    Rcpp::traits::input_parameter< double >::type lam1(lam1SEXP);
    Rcpp::traits::input_parameter< double >::type lam2(lam2SEXP);
    Rcpp::traits::input_parameter< double >::type lam3(lam3SEXP);
    rcpp_result_gen = Rcpp::wrap(dirder(Bin, lam1, lam2, lam3));
    return rcpp_result_gen;
END_RCPP
}
// build
NumericMatrix build(double d);
RcppExport SEXP _MPST_build(SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(build(d));
    return rcpp_result_gen;
END_RCPP
}
// triarea
double triarea(NumericVector V1, NumericVector V2, NumericVector V3);
RcppExport SEXP _MPST_triarea(SEXP V1SEXP, SEXP V2SEXP, SEXP V3SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type V1(V1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type V2(V2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type V3(V3SEXP);
    rcpp_result_gen = Rcpp::wrap(triarea(V1, V2, V3));
    return rcpp_result_gen;
END_RCPP
}
// locEng
NumericMatrix locEng(NumericVector V1, NumericVector V2, NumericVector V3, NumericMatrix Mat, double d);
RcppExport SEXP _MPST_locEng(SEXP V1SEXP, SEXP V2SEXP, SEXP V3SEXP, SEXP MatSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type V1(V1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type V2(V2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type V3(V3SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Mat(MatSEXP);
    Rcpp::traits::input_parameter< double >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(locEng(V1, V2, V3, Mat, d));
    return rcpp_result_gen;
END_RCPP
}
// energy
NumericMatrix energy(NumericMatrix V, NumericMatrix Tr, double d);
RcppExport SEXP _MPST_energy(SEXP VSEXP, SEXP TrSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type V(VSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Tr(TrSEXP);
    Rcpp::traits::input_parameter< double >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(energy(V, Tr, d));
    return rcpp_result_gen;
END_RCPP
}
// mvrbind
NumericMatrix mvrbind(NumericMatrix a, NumericVector b);
RcppExport SEXP _MPST_mvrbind(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(mvrbind(a, b));
    return rcpp_result_gen;
END_RCPP
}
// mvcbind
NumericMatrix mvcbind(NumericMatrix a, NumericVector b);
RcppExport SEXP _MPST_mvcbind(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(mvcbind(a, b));
    return rcpp_result_gen;
END_RCPP
}
// newcol
NumericMatrix newcol(Rcpp::Nullable<Rcpp::NumericMatrix> B0, Rcpp::Nullable<Rcpp::NumericVector> c0);
RcppExport SEXP _MPST_newcol(SEXP B0SEXP, SEXP c0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericMatrix> >::type B0(B0SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type c0(c0SEXP);
    rcpp_result_gen = Rcpp::wrap(newcol(B0, c0));
    return rcpp_result_gen;
END_RCPP
}
// tdata
List tdata(NumericMatrix V, NumericMatrix Tr);
RcppExport SEXP _MPST_tdata(SEXP VSEXP, SEXP TrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type V(VSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Tr(TrSEXP);
    rcpp_result_gen = Rcpp::wrap(tdata(V, Tr));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MPST_HQblkdiag", (DL_FUNC) &_MPST_HQblkdiag, 2},
    {"_MPST_HQbary", (DL_FUNC) &_MPST_HQbary, 4},
    {"_MPST_HQgetInd", (DL_FUNC) &_MPST_HQgetInd, 4},
    {"_MPST_mtxrbind", (DL_FUNC) &_MPST_mtxrbind, 2},
    {"_MPST_mtxcbind", (DL_FUNC) &_MPST_mtxcbind, 2},
    {"_MPST_vbind", (DL_FUNC) &_MPST_vbind, 2},
    {"_MPST_BSpline2", (DL_FUNC) &_MPST_BSpline2, 6},
    {"_MPST_BSpline", (DL_FUNC) &_MPST_BSpline, 6},
    {"_MPST_CrIndices", (DL_FUNC) &_MPST_CrIndices, 2},
    {"_MPST_CrCellArrays", (DL_FUNC) &_MPST_CrCellArrays, 2},
    {"_MPST_CrArrays", (DL_FUNC) &_MPST_CrArrays, 2},
    {"_MPST_degree", (DL_FUNC) &_MPST_degree, 1},
    {"_MPST_indices", (DL_FUNC) &_MPST_indices, 1},
    {"_MPST_bary", (DL_FUNC) &_MPST_bary, 5},
    {"_MPST_tcord", (DL_FUNC) &_MPST_tcord, 4},
    {"_MPST_locate", (DL_FUNC) &_MPST_locate, 6},
    {"_MPST_dirder", (DL_FUNC) &_MPST_dirder, 4},
    {"_MPST_build", (DL_FUNC) &_MPST_build, 1},
    {"_MPST_triarea", (DL_FUNC) &_MPST_triarea, 3},
    {"_MPST_locEng", (DL_FUNC) &_MPST_locEng, 5},
    {"_MPST_energy", (DL_FUNC) &_MPST_energy, 3},
    {"_MPST_mvrbind", (DL_FUNC) &_MPST_mvrbind, 2},
    {"_MPST_mvcbind", (DL_FUNC) &_MPST_mvcbind, 2},
    {"_MPST_newcol", (DL_FUNC) &_MPST_newcol, 2},
    {"_MPST_tdata", (DL_FUNC) &_MPST_tdata, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_MPST(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

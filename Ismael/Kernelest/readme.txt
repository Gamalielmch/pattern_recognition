Steve Marron's Matlab Smoothing Software
(c) J. S. Marron, 1997

for Matlab version 4.


Recommended first thing to try:

% GPANAL, General Purpose data ANALysis
    Makes a single page with family analysis, and SiZer 1
    for a quick first look at new data sets (1-d density or regression). 
    Uses: gpsz1.m, bwsjpib.m, bwrswb.m, gpkde.m, gpnpr.m


Basic High Level Methods:

% GPFAM, General Purpose FAMily approach to 1-d smoothing (kde or reg)
    Does a family of kernel smooths (Gaussian kernel), for
    either density estimation, or regression, using a binned 
    implementation.  Family is 15 (usually) logarithmically
    spaced smooths, centered at a good data driven bandwidth
    (SJPI for kde, RSW for reg).  Spread is determined by peak
    height.  Curves for bandwidths less than the binwidth, and
    greater than twice the range are not shown (or computed)
    Uses:      vec2mat.m, gplbinr.m, gpkde.m, bwsjpib.m, bwos.m, bwrswb.m

% GPSZ1, General Purpose Significant derivative Zero crossings
    Creates gray level map (function of location and bandwidth),
    which is:
        "dark"   at points where deriv sig > 0
        "light"  at points where deriv sig < 0
        "gray"   at points where deriv roughly 0
        "darker gray" where "effective sample size" < 5
    Uses:      vec2mat.m, gplbinr.m, phiinv.m


Some Smoothers:

% GPKDE, General Purpose Kernel Density Estimate (1-d, Gaussian Kernel)
    Does 1-d kernel density estimation, using binned (default) or 
    direct (either matrix, or loops for bigger data sets), 
    implementations, with the bandwidth either user specified 
    (can be vector), or data driven (SJPI, Normal Reference, 
    Silverman's ROT2, Oversmoothed).
    Used by:   gpfam.m 
    Uses:      vec2mat.m, gplbinr.m, bwsjpi.m, bwos.m, bwrot, bwsnr

% GPNPR, General Purpose NonParametric Regression (1-d, Local poly))
    Does 1-d kernel local polynomial (usually linear) regression,
    using binned (default), direct (either matrix, or loops for 
    bigger data sets), or moving window (for higher degree)
    implementations, with the bandwidth either user specified 
    (can be vector), or data driven (Ruppert Sheather and Wand ROT
    or DPI).  Kernel shape is Gaussian.
    Used by:   gpfam.m, bwrswb.m
    Uses:      vec2mat.m, gplbinr.m, bwrswb.m


An interesting Error Measure:

% GPVE1, General Purpose Visual Error Measure
    As developed in Marron and Tsybakov, JASA 1995.
    These are intended to correspond more to "visual
    impression of error", than the usual integral norms.
    For each point on [xgrid,v1] curve, finds "closest
    point on [xgrid,v2] curve", then summarizes these
    distances, according to norm choice (with optional 
    symmetrization).
    Uses:      vec2mat.m


Generic Bandwidth Selectors:

% BWSJPIB, Sheather Jones Plug In (Binned version), for 1-d k.d.e.
    Does data-based bandwidth selection for 1-d kernel density 
    estimation, using the Sheather Jones Plug In method.
    A binned implementation is used, and a grid search is done, 
    taking the largest local downcrossing.  Assumes a Gaussian kernel.
    Used by:   gpkde.m, gpfam.m 
    Uses:      vec2mat.m, gplbinr.m, bwrfph.m, bwos.m

% BWRSWB, Ruppert Sheather Wand Plug In (Binned version), for 1-d reg.
    Does data-based bandwidth selection for 1-d local linear reg.
    estimation, using the Ruppert Sheather Wand Plug In method.
    A linear binned implementation is used.
    Assumes a Gaussian kernel.  Returns twice range when n <= 5
    Used by:   gpfam.m, gpnpr.m
    Uses:      vec2mat.m, gplbinr.m, gpnpr.m

% BWOS, Terrell's OverSmoothed BandWidth, for 1-d kernel density estimation
    This is the bandwidth, which maximizes the MISE asymptotically
    optimal representation, modulo "scale", chosen here to
    b "standard deviation".   Assumes the kernel is standard Gaussian.
    Used by:   gpkde.m, bwsjpib.m

% BWROT, Silverman's Rule of Thumb BandWidths, for 1-d kde
    from Silverman's density estimation book
    Variations on the Simple normal reference:  bwsnr.m
    Assumes the kernel is standard Gaussian
    Used by:   gpkde
    Uses:      iqr, cquant

% BWSNR, Simple Normal Reference BandWidth, for 1-d kernel density estimation
    This is the MISE asy optimal bandwidth, if the data are 
    normally distributed, the sample variance as variance
    Assumes the kernel is standard Gaussian
    Used by:   gpkde


Useful Numerical Methodology:

% GPLBINR, General Purpose Linear BINneR (density and regression est.)
    Does linear binning of either density or regression 1-d data,
    to an equally spaced grid.
    Used by:   gpkde.m, gpnpr.m.

% BWRFPH, Estimates Integrated Squared Density Deriv's, for 1-d k.d.e.
    Binned implementation of a "diagonals in" estimate of 
    the "Roughness of Density derivatives", (i.e. the integral of
    their square).  The kernel is K, as opposed to K*K, which
    comes from integrating the actual kde. For K*K plug in the
    bandwidth * sqrt(2) for h.  Assumes a Gaussian kernel.
    Used by:   bwsjpib.m
    Uses:      vec2mat.m, gplbinr.m


Miscellaneous Useful Functions:

% GPROOTF, General Purpose ROOT Finder.
    Approximates roots of a discretized continuous function,
    fitting a cubic to nearest subintervals.
    Used by:   bwsjpib.m, bwrswb.m
    See also: gpminr.m

% GPMINR, General Purpose MINimizeR.
    Fits quadratic near smallest value, as approx to continuous min.
    In case of ties, finds the first one.
    See also: gprootf.m

% INTERP1S, linear interpolation, with constant extrapolation outside
    slight modification of linear version of INTERP1, which
    allows values xi values outside the range of x, and 
    returns the closest values at such points.

% IQR, Inter Quartile Range
    Uses:      vec2mat, cquant

% PHI, Standard Gaussian (normal) c.d.f.
    Inverse is phiinv.m

% PHIINV, Inverse of Standard Gaussian (normal) c.d.f.
    Inverse is phi.m

% CQUANT, Continuous empirical QUANTiles.
    Calculates "how far through the 1-d data set in 'data'",
    an element with prob vprob is, assuming the underlying 
    distribution is continuous.  When there are no ties in 
    data (expected unless there is rounding), the i-th order 
    statistic of 'data' is at prob i / (n + 1) (motivated by 
    putting mass 1 / (n + 1) between each order stat, and 
    also by the fact that for U_(i) i-th i.i.d. Unif(0,1),
    order stat, E(U_(i)) = 1 / (n+1)) and everything in 
    between is linearly interpolated (with linear extension 
    at the ends).  When there are ties in the data, points
    are combined before interpolation (so the resulting cdf is
    continuous).
    Useful for generating bootstrap data.
    For inverse of this, use cprob.m

% CPROB, Continuous empirical PROBabilities.
    Calculates "how far through the 1-d data set in 'data'",
    each element of vx is, assuming the underlying distribution is
    continuous.  When there are no ties in data (expected unless
    there is rounding), the i-th order statistic of 'data' is at
    prob i / (n + 1) (motivated by putting mass 1 / (n + 1)
    between each order stat, and also by the fact that for
    U_(i) i-th i.i.d. Unif(0,1), order stat, E(U_(i)) = 1 / (n+1))
    and everything in between is linearly interpolated (with linear 
    extension at the ends).  When there are ties in the data, points
    are combined before interpolation (so the resulting cdf is
    continuous).  Useful for computing bootstrap p-values.
    For inverse of this, use cquant.m


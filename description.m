%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%       This package contains codes for denoising matrices with heteroscedastic
%       noise, estimated by a sample covariance, using pseudo-whitening. The two
%       primary user-callable functions are as follows:
%
%   fshr_approx - applies optimal shrinkage with pseudo-whitening, using estimated
%       values; uses sample covariance directly.
%
%   fshr_approx2 - applies optimal shrinkage with pseudo-whitening, using estimated
%       values; uses out-of-sample noise matrix (not sample covariance)
%
%   fshr_exact - applies optimal shrinkage with pseudo-whitening, using oracle
%       knowledge of the signal (for purposes of comparison).
%
%   fshr_draw - draws example problem
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

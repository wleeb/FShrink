%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%       This package contains codes for denoising matrices with weighted
%       loss functions, using spectral denoising. The primary user-callable
%       codes are described below.
%
%   whtd_dense - estimates the optimal spectral denoiser for specified
%       symmetric matrix of row weights.
%
%   whtd_approx - estimates the optimal spectral denoiser for specified
%       row and column weghts.
%
%   whtd_exact - computes the exact spectral denoiser when given the
%       exact low-rank target matrix. For testing purposes only.
%   
%   whtd_local - performs localized denoising, with specified submatrices.
%
%   whtd_missing - performs optimal spectral denoising for matrix wth missing
%      entries. Data is backprojected, noise is whitened, and weighting
%      is performed after denoising.
%
%   whtd_colored - performs optimal spectral denoising for matrix wth
%      additive colored noise with rank 1 variance profile. Noise is
%      whitened, and recoloring is performed after denoising.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

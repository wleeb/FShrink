%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%        This code contains four user-callable functions for analyzing random
%        matrices of the form
%
%                      N = A^{-1/2} G B^{-1/2} / sqrt(n)              (1)
%
%        where A and B are two user-specified sequences of positive weights
%        of lengths m and n, respectively, each with an associated probability
%        distribution. Here, G is an m-by-n matrix of iid N(0,1) entries.
%
%     mpbdry_stiel - evaluates the Stieltjes transform of the LSD of N at a
%        specified value lambda
%     mpbdry_edge - evaluates the asymptotic operator norm squared of the random 
%        matrix N
%     mpbdry_thresh - evaluates the asymptotic operator norm squared of the random 
%        matrix N and the minimum detectable signal eigenvalue
%     mpbdry_sback - evaluates population parameters for spiked matrix 
%        model; from observed value to population value
%     mpbdry_sforw - evaluates population parameters for spiked matrix 
%        model; from population value to observed value
%
%        Additional routines evaluate higher-order derivatives of certain
%        intermediate functions, which are useful for theoretical purposes.
%        Among these routines are:
%
%     mpbdry_evalqder2 - evaluate first and second derivatives of Q
%     mpbdry_fmin - evaluates the minimum of F
%     mpbdry_fminder - evaluates the derivative of the minimum of F
%     mpbdry_evalfmix - evaluates second and mixed derivatives of F
%     mpbdry_evalfg3 - evaluates all derivatives up to third order of both
%        F and G
%     mpbdry_evalfg3_slow - same as mpbdry_evalfg3, but not vectorized
%
%        The methods for these codes are described in the paper
%        ``Rapid evaluation of the spectral signal detection threshold and
%            Stieltjes transform'',
%        by W. Leeb.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

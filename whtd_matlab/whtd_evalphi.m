        function [phi,phi_der,rphi] = whtd_evalphi(s,gam,sig,amu,alpha)
%
%        rescale alpha and amu
%
        alpha = alpha / amu;
        amu=1;

        [cout,cout_der,rout] = whtd_evalcout(s,gam,sig);
%
%        phi and derivative phi'
%
        phi = alpha * cout / (alpha*cout^2 + amu*(1-cout^2));
        phi_der = alpha*cout_der*(1 - (alpha-1)*cout^2) ...
            / (1 + (alpha-1)*cout^2)^2;
%
%        ratio between phi' and phi
%
        dd = (alpha-1)*cout^2;
        ff = (1 - dd) / (1 + dd);
        rphi = (cout_der / cout) * ff;

        end
%

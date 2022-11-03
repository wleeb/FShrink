        function [topt,topt_der,rtopt] = whtd_evaltopt(s,gam,sig,amu,bmu,...
           alpha,beta)
%
%        evaluate optimal singular value for rank 1 signal, and its
%        derivative, and ratio of derivative to value
%
        [rlam,cout,cinn] = whtd_pop2emp(s,gam,sig);

        [phi,phi_der,rphi] = whtd_evalphi(s,gam,sig,amu,alpha);
        [psi,psi_der,rpsi] = whtd_evalpsi(s,gam,sig,bmu,beta);

        topt = s*phi*psi;
        topt_der = phi*psi + s*phi_der*psi + s*phi*psi_der;
%
%        formula for ratio of derivative to topt
%
        rtopt = phi_der / phi + psi_der / psi + 1/s;

        end
%

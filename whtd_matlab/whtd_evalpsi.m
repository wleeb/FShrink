        function [psi,psi_der,rpsi] = whtd_evalpsi(s,gam,sig,bmu,beta)
%
%        rescale beta and bmu
%
        beta = beta / bmu;
        bmu = 1;

        [cinn,cinn_der,rinn] = whtd_evalcinn(s,gam,sig);

%
%        psi and derivative psi'
%
        psi = beta * cinn / (beta*cinn^2 + bmu*(1-cinn^2));
        psi_der = beta*cinn_der*(1 - (beta-1)*cinn^2) ...
            / (1 + (beta-1)*cinn^2)^2;
%
%        ratio between psi' and psi
%
        dd = (beta-1)*cinn^2;
        ff = (1 - dd) / (1 + dd);
        rpsi = (cinn_der / cinn) * ff;

        end
%

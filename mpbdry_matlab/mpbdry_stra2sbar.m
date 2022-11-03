        function [sbar,sbar_der] = mpbdry_stra2sbar(stra,stra_der,z,gam)
%
%        given s(z) and s'(z), computes sbar(z) and sbar'(z)
%
        sbar = gam*stra - (1-gam)/z;
        sbar_der = gam*stra_der + (1-gam) / z^2;

        end
%

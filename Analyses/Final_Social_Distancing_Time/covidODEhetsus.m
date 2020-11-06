function dsdt = covidODEhetsus(t,s)

    % @ copyright
    % Authors:
    %   Ricardo Aguas
    %   Rodrigo M Corder
    %   Jessica G King
    %   Guilherme Goncalves
    %   Marcelo U Ferreira
    %   M Gabriela M Gomes
    %
    % This work is protected under the @Attribution-NonCommercial 4.0 International intellectual property license.
    % You are free to:
    %   Share - copy and redistribute the material in any medium or format
    %   Adapt - remix, transform, and build upon the material Under the following terms:
    %   Attribution - You must give appropriate credit to the authors, and indicate if any changes were made. You may do so in any reasonable manner, but not in any way that suggests the licensor endorses you or your use.
    %   NonCommercial - You may not use the material for commercial purposes.
    %   ShareAlike - If you remix, transform, or build upon the material, you must distribute your contributions under the same license as the original.

    global lam de ga rh N xvec k1 drampup drampdown dmax p inidist data R0 i

    dsdt = zeros(3*N,1);
    fsd=x;
    if(t<inidist+1)
        prox = 1;
    elseif(t>=inidist+1 && t<=inidist + drampup)
        prox = 1-(t-inidist)*(1-p)/drampup;
    elseif(t>inidist+drampup && t<=inidist+drampup+dmax)
        prox = p;
    elseif (t>inidist+drampup+dmax && t<=inidist+drampup+dmax+(1-fsd)*drampdown)
        prox = 1 - (inidist+drampup+dmax+drampdown-t)*(1-p)/drampdown;
    else
        prox = (1-(1-p)*fsd);
    end    
    
    behetsus = R0(i)/(rh/de+1/ga)/k1;
    lam = behetsus*prox/data(i).pop*(rh*sum(s(N+1:2*N))+sum(s(2*N+1:3*N)));
    

    dsdt(1:N)       = -lam*xvec'.*s(1:N);
    dsdt(N+1:2*N)   =  lam*xvec'.*s(1:N) - de*s(N+1:2*N);
    dsdt(2*N+1:3*N) =  de*s(N+1:2*N)     - ga*s(2*N+1:3*N);

end

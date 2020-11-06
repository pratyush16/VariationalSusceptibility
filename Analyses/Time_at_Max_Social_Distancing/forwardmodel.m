function LogL = forwardmodel(m)

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

global data de p model inidist var R0 N xvec qvec initinf initinfe k1...
    k2 i dmax drampup drampdown aux_pop

LogL = 0;
m    = exp(m);

%parameters
p       = m(1);
inidist = m(2);

for j=1:length(data)
    R0(j)   = m(2+j);
    vvar(j) = m(2+length(data)+j);
end

drampup   = 21;
dmax      = 30;
drampdown = 120;

for i=1:length(data)
    
    % Heterogeneous distribution
    var=vvar(i);
    Susceptibility_dist;
    k1 = xvec*qvec;         % mean susceptibility
    k2 = xvec.^2*qvec;      % mean connectivity
    
    % Initial conditions
    if (model==1)
        inf = data(i).pop/aux_pop;
        s0 = zeros(3,1);
        initinf = (inf/data(i).rep);
        s0(1) = data(i).pop-(inf/data(i).rep);
        s0(2) = initinf/(1-exp(-de));
        s0(3) = initinf;
        xvec = 1;
        N = 1;
    else
        inf = data(i).pop/aux_pop;
        s0= zeros(3*N,1);
        initinfe = (inf/data(i).rep)/(1-exp(-de))*qvec;
        initinf  = (inf/data(i).rep)*qvec;
        s0(1:N) = (data(i).pop)*qvec - initinfe - initinf;
        s0(N+1:2*N) = initinfe;
        s0(2*N+1:3*N) = initinf;  
    end 
    
    % Solving ODE's
    if(model == 1)
        soln1 = ode45(@covidODEhomo,data(i).tspan,s0);
        aux_soln1 = deval(soln1,data(i).tspan);
        incidence = de*aux_soln1(2,:)';
    elseif(model == 2)
        soln1 = ode45(@covidODEhetsus,data(i).tspan,s0);
        aux_soln1 = deval(soln1,data(i).tspan);
        incidence = de*sum(aux_soln1(N+1:2*N,:))';     
    elseif(model == 3)
        soln1 = ode45(@covidODEhetcon,data(i).tspan,s0);
        aux_soln1 = deval(soln1,data(i).tspan);
        incidence = de*sum(aux_soln1(N+1:2*N,:))';
    elseif(model == 4)
        soln1 = ode45(@covidODEhetdyn,data(i).tspan,s0);
        aux_soln1 = deval(soln1,data(i).tspan);
        incidence = de*sum(aux_soln1(N+1:2*N,:))';        
    end
    
    data(i).incidence = incidence*data(i).rep;
    
    Poi_lambda = data(i).incidence;
    Poi_y      = data(i).disease(data(i).tspan);
    
    % Computing likelihood
    data(i).LogL = sum(-Poi_lambda+Poi_y.*log(Poi_lambda));
    LogL = LogL + data(i).LogL;
    
end
end
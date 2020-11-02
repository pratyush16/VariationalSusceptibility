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

%% Parameter estimation

clc
clear all

global model data prefix ga de rh ndays aux_pop

% Aux variable to define the intial of epidemic in each country (1 case per aux_pop population)
aux_pop  = 5e6;
    
% selecting a country to be analyzed
hh = input('Choose a country: \n\n 1.Belgium \n 2.Portugal \n 3.Spain \n 4.England \n 5.Belgium - stratified by regions \n 6.Portugal - stratified by regions \n 7.Spain - stratified by regions \n 8.England - stratified by regions \n\n Country:'); 
switch hh
    case 1
        Belgium_single;
    case 2
        Portugal_single;
    case 3
        Spain_single;
    case 4
        England_single;        
    case 5
        Belgium;
    case 6
        Portugal;
    case 7
        Spain;          
    case 8
        England;      
end  

% selecting a type of model
model = input('\nChoose a model: \n\n 1.Homogeneous \n 2.Heterogeneous - susceptibility \n 3.Heterogeneous - connectivity \n 4.Heterogeneous - connectivity with dynamic CV \n\n Model:');  
switch model
    case 1
        prefix = 'Homo';
    case 2
        prefix = 'Hetsus';
    case 3
        prefix = 'Hetcon';
    case 4
        prefix = 'Hetdyn';        
otherwise
        disp('\n\n Invalid model!!')
end  

fprintf(['\n\n********************************************'])
fprintf(['\n' data(1).country ' - ' prefix '\n'])

%% Define epidemiological parameters

ga    = 1/4;    % rate of recovery/death
de    = 1/4;    % rate of progression from exposed to infectious
rh    = 0.5;    % transmissibility of exposed individuals
ndays = 366;    % simulation time length

%% Forward model

%Likelihood
logLike_P=@(m) forwardmodel(m);   

%Parameter optimization - PESTO options
opt.MSoptions                       = PestoOptions();

opt.name                            = {'p','inidist'};
for j=1:length(data)
    opt.name = [opt.name,{['R0_' num2str(j)]}];
end 
for j=1:length(data)
    opt.name = [opt.name,{['CV_' num2str(j)]}];
end 

opt.number  = length(opt.name);
opt.min     = log([0.1;min([data.lag]);1*ones(length(data),1);0.0025*ones(length(data),1)]);
opt.max     = log([0.7;min([data.lag])+40;12*ones(length(data),1);49*ones(length(data),1)]);     

opt.MSoptions.comp_type             = 'sequential';
opt.MSoptions.n_starts              = 50;
opt.MSoptions.localOptimizer        = 'fmincon';
opt.MSoptions.mode                  = 'text';
opt.MSoptions.localOptimizerOptions = optimset('Algorithm','interior-point',...
                                    'Hessian','off',...  
                                    'GradObj','off',... 
                                    'MaxIter',1000,...           
                                    'MaxFunEvals',3000,...
                                    'display', 'iter');                                  

Opt_parameters = getMultiStarts(opt,logLike_P,opt.MSoptions);    

%% MCMC

fprintf('\nMCMC\n')
optionsPesto                      	= PestoOptions();
optionsPesto.mode                   = 'text';
optionsPesto.MCMC.nIterations      	= 1e4;
optionsPesto.MCMC.mode             	= 'silent';
optionsPesto.MCMC.samplingAlgorithm	= 'PT';
optionsPesto.MCMC.theta0           	= Opt_parameters.MS.par(:,1);   % Starting MCMC chain in the best value from the optimization
optionsPesto.MCMC.sigma0          	= 0.5 * inv(squeeze(Opt_parameters.MS.hessian(:,:,1)));
parameters                        	= getParameterSamples(Opt_parameters, logLike_P, optionsPesto);
alpha                             	= [0.95]; % 95% credible interval
parameters                         	= getParameterConfidenceIntervals(parameters, alpha, optionsPesto);

%% Storing final results

Final_results.country      = data(1).country;
Final_results.model        = prefix;
Final_results.LogL         = Opt_parameters.MS.logPost(1);
Final_results.d            = exp([parameters.CI.S(1,1) median(parameters.S.par(1,:)) parameters.CI.S(1,2)]);
Final_results.p            = 1 - exp([parameters.CI.S(1,1) median(parameters.S.par(1,:)) parameters.CI.S(1,2)]);
Final_results.inidist      = exp([parameters.CI.S(2,1) median(parameters.S.par(2,:)) parameters.CI.S(2,2)]);

% Whole country analysis
if (hh<5)
    Final_results.R0_1     = exp([parameters.CI.S(3,1) median(parameters.S.par(3,:)) parameters.CI.S(3,2)]);
    Final_results.R0_2     = [];
    Final_results.R0_3     = [];
    Final_results.R0_4     = [];
    
    if (model>1)     
        Final_results.CV_1 = (exp([parameters.CI.S(4,1) median(parameters.S.par(4,:)) parameters.CI.S(4,2)])).^(1/2);
    else
        Final_results.CV_1 = [];                
    end
    Final_results.CV_2 = [];
    Final_results.CV_3 = [];
    Final_results.CV_4 = [];   
    
    Final_results.lag_1    = data(1).lag; 
    Final_results.lag_2    = []; 
    Final_results.lag_3    = []; 
    Final_results.lag_4    = []; 
end

% Belgium and Portugal (stratified by regions)
if (hh==5 || hh==6)
    Final_results.R0_1     = exp([parameters.CI.S(3,1) median(parameters.S.par(3,:)) parameters.CI.S(3,2)]);
    Final_results.R0_2     = exp([parameters.CI.S(4,1) median(parameters.S.par(4,:)) parameters.CI.S(4,2)]);
    Final_results.R0_3     = [];
    Final_results.R0_4     = [];
    
    if (model>1)
        Final_results.CV_1 = (exp([parameters.CI.S(5,1) median(parameters.S.par(5,:)) parameters.CI.S(5,2)])).^(1/2);
        Final_results.CV_2 = (exp([parameters.CI.S(6,1) median(parameters.S.par(6,:)) parameters.CI.S(6,2)])).^(1/2);
    else
        Final_results.CV_1 = [];
        Final_results.CV_2 = [];              
    end
    Final_results.CV_3 = [];
    Final_results.CV_4 = [];    
    
    Final_results.lag_1    = data(1).lag;
    Final_results.lag_2    = data(2).lag; 
    Final_results.lag_3    = []; 
    Final_results.lag_4    = []; 
end

% Spain (stratified by regions)
if (hh==7)
    Final_results.R0_1     = exp([parameters.CI.S(3,1) median(parameters.S.par(3,:)) parameters.CI.S(3,2)]);
    Final_results.R0_2     = exp([parameters.CI.S(4,1) median(parameters.S.par(4,:)) parameters.CI.S(4,2)]);
    Final_results.R0_3     = exp([parameters.CI.S(5,1) median(parameters.S.par(5,:)) parameters.CI.S(5,2)]);
    Final_results.R0_4     = [];
    
    if (model>1)           
        Final_results.CV_1 = (exp([parameters.CI.S(6,1) median(parameters.S.par(6,:)) parameters.CI.S(6,2)])).^(1/2);
        Final_results.CV_2 = (exp([parameters.CI.S(7,1) median(parameters.S.par(7,:)) parameters.CI.S(7,2)])).^(1/2);
        Final_results.CV_3 = (exp([parameters.CI.S(8,1) median(parameters.S.par(8,:)) parameters.CI.S(8,2)])).^(1/2);
    else
        Final_results.CV_1 = [];
        Final_results.CV_2 = [];
        Final_results.CV_3 = [];               
    end
    Final_results.CV_4 = [];
    
    Final_results.lag_1    = data(1).lag;
    Final_results.lag_2    = data(2).lag; 
    Final_results.lag_3    = data(3).lag; 
    Final_results.lag_4    = []; 
end

% England (stratified by regions)
if (hh==8)
    Final_results.R0_1     = exp([parameters.CI.S(3,1) median(parameters.S.par(3,:)) parameters.CI.S(3,2)]);
    Final_results.R0_2     = exp([parameters.CI.S(4,1) median(parameters.S.par(4,:)) parameters.CI.S(4,2)]);
    Final_results.R0_3     = exp([parameters.CI.S(5,1) median(parameters.S.par(5,:)) parameters.CI.S(5,2)]);
    Final_results.R0_4     = exp([parameters.CI.S(6,1) median(parameters.S.par(6,:)) parameters.CI.S(6,2)]);
    
    if (model>1)   
        Final_results.CV_1 = (exp([parameters.CI.S(7,1) median(parameters.S.par(7,:)) parameters.CI.S(7,2)])).^(1/2);
        Final_results.CV_2 = (exp([parameters.CI.S(8,1) median(parameters.S.par(8,:)) parameters.CI.S(8,2)])).^(1/2);
        Final_results.CV_3 = (exp([parameters.CI.S(9,1) median(parameters.S.par(9,:)) parameters.CI.S(9,2)])).^(1/2);
        Final_results.CV_4 = (exp([parameters.CI.S(10,1) median(parameters.S.par(10,:)) parameters.CI.S(10,2)])).^(1/2);
    else
        Final_results.CV_1 = [];
        Final_results.CV_2 = [];
        Final_results.CV_3 = [];
        Final_results.CV_4 = [];                
    end
    
    Final_results.lag_1    = data(1).lag;
    Final_results.lag_2    = data(2).lag; 
    Final_results.lag_3    = data(3).lag;
    Final_results.lag_4    = data(4).lag;
end

Final_results

%% Saving final results

if ~exist(Final_results.country,'dir')
    mkdir(Final_results.country)
end

save([Final_results.country '/' data(1).country '_' prefix '.mat'])
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

global data aux_pop

C = readtable('Spain_incidence.csv');
S = readtable('Spain_serology.csv','ReadRowNames',true);

data(1).country  = 'Spain';
data(1).pop      = S{4,20};
data(1).disease  = (C{:,20});
data(1).lag      = find(C{:,20}>=data(1).pop/aux_pop,1);
data(1).tspan    = (data(1).lag:(numel(data(1).disease)));
data(1).serology = S{1:3,20}/100;
data(1).rep      = sum(data(1).disease(1:45))/(data(1).pop*data(1).serology(1));
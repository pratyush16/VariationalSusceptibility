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

C = readtable('Belgium_incidence.csv');
S = readtable('Belgium_serology.csv','ReadRowNames',true);

data(1).country  = 'Flanders';
data(1).pop      = S{1,2};
data(1).disease  = (C{:,2});
data(1).lag      = find(data(1).disease>=data(1).pop/aux_pop,1);
data(1).tspan    = (data(1).lag:(numel(data(1).disease)));
data(1).rep      = 0.06;

data(2).country  = 'Rest';
data(2).pop      = S{1,3};
data(2).disease  = (C{:,3});
data(2).lag      = find(data(2).disease>=data(2).pop/aux_pop,1);
data(2).tspan    = (data(2).lag:(numel(data(2).disease)));
data(2).rep      = 0.06;



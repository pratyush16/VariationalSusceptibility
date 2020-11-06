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

C = readtable('England_incidence.csv');
S = readtable('England_serology.csv','ReadRowNames',true);

data(1).country  = 'London';
data(1).pop      = S{1,2};
data(1).disease  = C{:,2};
data(1).lag      = find(data(1).disease>=data(1).pop/aux_pop,1);
data(1).tspan    = (data(1).lag:(numel(data(1).disease)));
data(1).rep      = 0.024;

data(2).country  = 'NorthWest';
data(2).pop      = S{1,3};
data(2).disease  = C{:,3};
data(2).lag      = find(data(2).disease>=data(2).pop/aux_pop,1);
data(2).tspan    = (data(2).lag:(numel(data(2).disease)));
data(2).rep      = 0.024;
% 
data(3).country  = 'SouthEast';
data(3).pop      = S{1,4};
data(3).disease  = C{:,4};
data(3).lag      = find(data(3).disease>=data(3).pop/aux_pop,1);
data(3).tspan    = (data(3).lag:(numel(data(3).disease)));
data(3).rep      = 0.024;

data(4).country  = 'Rest';
data(4).pop      = S{1,5};
data(4).disease  = C{:,5};
data(4).lag      = find(data(4).disease>=data(4).pop/aux_pop,1);
data(4).tspan    = (data(4).lag:(numel(data(4).disease)));
data(4).rep      = 0.024;
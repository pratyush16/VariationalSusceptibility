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

data(1).country  = 'Madrid';
data(1).pop      = S{4,18};
data(1).disease  = (C{:,18});
data(1).lag      = find(data(1).disease>=data(1).pop/aux_pop,1);
data(1).tspan    = (data(1).lag:(numel(data(1).disease)));
data(1).serology = S{1:3,18}/100;
data(1).rep      = sum(data(1).disease(1:45))/(data(1).pop*data(1).serology(1));

data(2).country  = 'Catalunya';
data(2).pop      = S{4,19};
data(2).disease  = C{:,19};
data(2).lag      = find(data(2).disease>=data(2).pop/aux_pop,1);
data(2).tspan    = (data(2).lag:(numel(data(2).disease)));
data(2).serology = S{1:3,19}/100;
data(2).rep      = sum(data(2).disease(1:45))/(data(2).pop*data(2).serology(1));

data(3).country  = 'Rest';
data(3).pop      = sum(S{4,[1:17]});
data(3).disease  = (sum(C{:,[1:17]},2));
data(3).lag      = find(sum(C{:,[1:17]},2)>=data(3).pop/aux_pop,1);
data(3).tspan    = (data(3).lag:(numel(data(3).disease)));
data(3).serology = sum(S{1:3,[1:17]}.*S{4,[1:17]}/sum(S{4,[1:17]}),2)/100;
data(3).rep      = sum(data(3).disease(1:45))/(data(3).pop*data(3).serology(1));
function [vZ_load, vZ_St, vZ_C, A_F] = middleEarStapesCochleaInfant(varargin)
% 
% Computes the mechanical load impedance of stapes and cochlea according to
% [1, 2, 3, 4, 5]
% Optionally the function returns the impedances of the single elements, 
% vZ_St (stapes, mechanical) and vZ_C (cochlea, acoustical), and the area 
% of the stapes footplate A_F.
% 
% SYNTAX:
%   [vZ_load, vZ_St, vZ_C, A_F] = middleEarStapesCochleaInfant(vOmega)
%   [vZ_load, ___] = middleEarStapesCochleaInfant(vOmega, Name, Value)
%   [vZ_load, ___] = middleEarStapesCochleaInfant(vOmega, stStapesCochlea)
%   stStapesCochlea = middleEarStapesCochleaInfant
%   stStapesCochlea = middleEarStapesCochleaInfant(Name, Value)
%   stStapesCochlea = middleEarStapesCochleaInfant(stStapesCochlea)
% 
% INPUTS    - with no inputs the function returns a struct containing the
%           default parameters
%           - vOmega       frequency in radians
%               - with no other inputs the model will be computed using the 
%               default parameter values
%           - Different values can be used with the name/value syntax
%           and/or with a struct "stStapesCochlea". 
%           The names or fieldnames are as follows:
%             'w_St'	% stapes mechanical losses, in Ns/m
%             'm_St'	% stapes mechanical mass, in kg
%             'n_St'	% stapes mechanical compliance, in m/N
%             'A_F'     % footplate area, in m^2
%             'w_C'     % cochlea mechanical losses, in Ns/m
%             'm_C'     % cochlea mechanical mass, in kg
%             'n_C'     % cochlea mechanical compliance, in m/N
% 
% OUTPUTS   - vZ_load   Nx1 vector of the mechanical load impedance, with
%                       N=length(vOmega)
%           - vZ_St     Nx1 vector of the stapes mechanical impedance
%           - vZ_C      Nx1 vector of the cochlea acoustical impedance
%           - A_F       area of the stapes footplate
% 
% 
%       Z_load = Z_St + Z_C*(A_F^2)
%   with
%       Z_St = w_St + jw m_St + 1/(jw n_St) 
%   and
%       Z_C = w_C/(A_F^2) + jw m_C/(A_F^2) + (A_F^2)/(jw n_C) 
% 
% REFERENCES:
% [1] Sankowsky-Rothe et al. (2022), "Parametric model of young infants' 
%     eardrum and ear canal impedances supporting immittance measurement 
%     results. Part I: Development of the model", Acta Acustica
% [2] Sankowsky-Rothe et al. (2022), "Parametric model of young infants' 
%     eardrum and ear canal impedances supporting immittance measurement 
%     results. Part II: Prediction of eardrum and ear canal impedances for 
%     common middle ear pathologies", Acta Acustica
% [3] Hudde & Engel (1998). "Measuring and Modeling Basic Properties of the
%     Human Middle Ear and Ear Canal. Part I: Model Structure and Measuring
%     Techniques" ACUSTICA/acta acustica, vol. 84, no. 4, pp. 720-738
% [4] Hudde & Engel (1998). "Measuring and Modeling Basic Properties of the
%     Human Middle Ear and Ear Canal. Part II: Ear Canal, Middle Ear 
%     Cavities, Eardrum, and Ossicles" ACUSTICA/acta acustica, vol. 84, 
%     no. 5, pp. 894-913 
% [5] Hudde & Engel (1998). "Measuring and Modeling Basic Properties of the
%     Human Middle Ear and Ear Canal. Part III: Eardrum Impedances, 
%     Transfer Functions and Model Calculations" ACUSTICA/acta acustica, 
%     vol. 84, no. 6, pp. 1091-1108
% 
% AUTHOR:   Tobias Sankowsky-Rothe
% DATE:     2022
% LICENSE:  see EOF

% CHANGELOG:


p = inputParser;

% --- default values
def_w_St = 18 * 10^(-3);	% stapes mechanical losses, in Ns/m
def_m_St = 3 * 10^(-6);		% stapes mechanical mass, in kg
def_n_St = 1.2 * 10^(-3);	% stapes mechanical compliance, in m/N
def_A_F = 3 * 10^(-6);		% footplate area, in m^2
def_w_C = 70 * 10^(-3);		% cochlea mechanical losses, in Ns/m
def_m_C = 10 * 10^(-6);		% cochlea mechanical mass, in kg
def_n_C = 11 * 10^(-3);		% cochlea mechanical compliance, in m/N

def_vOmega = 'returnParameters';

addOptional(p, 'vOmega', def_vOmega);

addParameter(p, 'w_St', def_w_St);
addParameter(p, 'm_St', def_m_St);
addParameter(p, 'n_St', def_n_St);
addParameter(p, 'A_F', def_A_F);
addParameter(p, 'w_C', def_w_C);
addParameter(p, 'm_C', def_m_C);
addParameter(p ,'n_C', def_n_C);

% --- parse optional inputs
parse(p,varargin{:})

% --- check vOmega or ruturn a list of the parameter values
if ischar(p.Results.vOmega)
    if strcmpi(p.Results.vOmega,def_vOmega)
        vZ_load = rmfield(p.Results,'vOmega');
        return
    else
        error(['The value of ''vOmega'' is invalid. It must be either a \n',...
            'Nx1 frequency vector in radians, or \n',...
            'the input ''%s'' in order to get a struct of the parameter ',...
            'values.'], def_vOmega)
    end
elseif ~isnumeric(p.Results.vOmega)
    error(['The value of ''vOmega'' is invalid. It must be a \n',...
        'Nx1 frequency vector in radians, or \n.'])
elseif min(size(p.Results.vOmega))>1
    error(['The value of ''vOmega'' is invalid. It must be a \n',...
        'Nx1 frequency vector in radians.'])
else
    vOmega = p.Results.vOmega;
end

vZ_St = p.Results.w_St+vOmega*1i*p.Results.m_St+(1)./(vOmega*1i*p.Results.n_St);
% --- compute acoustical cochlea parameters
A_F = p.Results.A_F;
w_Cac = p.Results.w_C / (A_F^2);    % cochlea acoustical losses, in Ns/m^5
m_Cac = p.Results.m_C / (A_F^2);    % cochlea acoustical mass, in kg/m^4
n_Cac = p.Results.n_C * (A_F^2);    % cochlea acoustical compliance, in m^5/N
vZ_C = w_Cac + vOmega * 1i * m_Cac + 1 ./ (vOmega * 1i * n_Cac);

vZ_load = vZ_St + vZ_C * (A_F^2);


%--------------------Licence ---------------------------------------------
% Copyright (C) 2022  Tobias Sankowsky-Rothe
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.

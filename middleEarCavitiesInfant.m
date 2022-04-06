function [vZ_cav] = middleEarCavitiesInfant(varargin)
% 
% Computes the acoustic impedance of the middle ear cavities Z_cav for 
% young infants' ears. Alternatively, the function returns the
% default/current parameter values used to model the acoustic impedance of
% the middle ear cavities.
% 
% SYNTAX:
%   vZ_cav = middleEarCavitiesInfant(vOmega)
%   vZ_cav = middleEarCavitiesInfant(vOmega, Name, Value)
%   vZ_cav = middleEarCavitiesInfant(vOmega, stCav)
%   stCav = middleEarCavitiesInfant
%   stCav = middleEarCavitiesInfant(Name, Value)
%   stCav = middleEarCavitiesInfant(stCav)
% 
% DESCRIPTION:
% - middleEarCavitiesInfant(vOmega,___) returns a vector of the acoustic
%   impedance of the middle ear cavities for the frequencies in radians
%   given in vOmega
% - middleEarCavitiesInfant or middleEarCavitiesInfant(___) in which vOmega
%   is ommited returns a struct containing the complete set of parameters 
% 
% INPUTS:   - vOmega        frequency in radians
%               - with no other inputs the model will be computed using the 
%               default parameter values
%           - parameter 'AgeGroup' with one of the following values can be 
%           used to select specific default values: 
%               - 'infant' - default acc. to [1, 2], will be used if 
%               AgeGroup is not defined
%               - 'adult' - default for adult middle ear cavities acc. to
%               [3, 4, 5]
%           - Different values can be used with the name/value syntax
%           and/or with a struct "stCav". 
%           The names or fieldnames are as follows:
%             'W_t'     % tympanic cavity acoustical losses, in 
%                               Ns/m^5
%             'V_t'     % tympanic cavity volume, in m^3
% 
%             'W_ada'	% aditus ad antrum acoustical losses, in 
%                               Ns/m^5
%             'M_ada'	% aditus ad antrum acoustical mass, in 
%                               Ns^2/m^5
%             'V_ant'	% antrum cavity volume, in m^3
%             'Q_mac'	% mastoid air cells quality factor
%             'f_mac'	% mastoid air cells characteristic 
%                               frequency, in Hz
%             'V_mac'	% mastoid air cells volume, in m^3
%             'stProp_t'	% optional struct with fields c and rho 
%                               allowing to model compliance of tympanic
%                               cavity with different mass density and 
%                               speed of sound
%             'stProp_ant'	% optional struct with fields c and rho 
%                               allowing to model compliance of antrum
%                               with different mass density and speed of 
%                               sound
%             'x_A'     % used to model the pathological condition of fluid 
%                       in the middle ear, x_A is the factor of the eardrum
%                       area which is free from fluid, acc. to [2]
%                       (e.g. for the healthy middle ear x_A=1, however, 
%                       for a pathological middle ear x_A=0.2 can be used) 
%                       Effect: reduces the volume of air in the tympanic
%                       cavity
%                       V_t_effective = stCav.x_A * stCav.V_t
%             'deltaP0' % used to model the pathological condition of a
%                       pressure difference of the static air between the
%                       ear canal and the middle ear in Pa, acc. to [2]
%                       (e.g. for the healthy middle ear deltaP0=0,
%                       however, for a pathological middle ear
%                       deltaP0=-2500 can be used)
%                       Effect: reduces the mass density of air in the
%                       tympanic cavity and in the antrum
%                       stProp_t.rho_effective = stCav.stProp_t.rho * (p0 + stCav.deltaP0)/p0;
%                       stProp_ant.rho_effective = stCav.stProp_ant.rho * (p0 + stCav.deltaP0)/p0;
% 
% OUTPUTS:	- vZ_cav    Nx1 vector of the cavities' acoustic impedance, 
%                       with N=length(vOmega)
% 
% REQUIRED FUNCTIONS:
%           - propertiesOfAir.m
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
p.PartialMatching = false;


% check if AgeGroup is specified as name/value argument
ixAG = find(strcmpi(varargin,'AgeGroup'));
if ~isempty(ixAG)
    def_AgeGroup = lower(varargin{ixAG+1});
else
    def_AgeGroup = 'infant';
end
% check if AgeGroup is specified in input struct
vIdc = find(cellfun(@isstruct, varargin));
for iArg = 1:length(vIdc)
    ixArg = vIdc(iArg);
    vIdcAG = strcmp(fieldnames(varargin{ixArg}), 'AgeGroup');
    if any(vIdcAG)
        % check if it was already defined differently using name/value arg.
        if ~isempty(ixAG) && ~strcmpi(def_AgeGroup, varargin{ixArg}.AgeGroup)
            error('Inconsistent definition: AgeGroup was defined differently in ''AgeGroup'' and in input-struct')
        end
        def_AgeGroup = lower(varargin{ixArg}.AgeGroup);
    end
end

% check if other parameters are defined twice
vIdc = find(cellfun(@isstruct, varargin));
for iArg = 1:length(vIdc)   
    % we have a struct
    ixArg = vIdc(iArg);
    cFn = fieldnames(varargin{ixArg});
    for ii = 1:length(cFn)
        ixName = find(strcmp(varargin, cFn{ii}));
        if ~isempty(ixName) 
            % we have the same parameter given as N/V and as field in struct
            if isstruct(varargin{ixName +1}) || isstruct(varargin{ixArg}.(cFn{ii}))
                % we have at least one struct
                warning('Name/value definition of ''%s'' will be used, definition in struct will be ignored', ...
                    varargin{ixName})
                disp(varargin{ixName +1})
                
                varargin{ixArg}.(cFn{ii}) = varargin{ixName +1};
                
            elseif varargin{ixName +1} ~= varargin{ixArg}.(cFn{ii})
                % we have different values, throw warning and use
                % name/value definition instead of value given in struct
                szWarn = sprintf(['Inconsistent definition:\n ''%s'' was defined differently by using name/value ', ...
                    'and struct input\n the name/value input will be used\n'], varargin{ixName});
                if isnumeric(varargin{ixName +1})
                    szWarn = sprintf('%s %s = %g', szWarn, varargin{ixName}, ...
                        varargin{ixName +1});
                elseif ischar(varargin{ixName +1})
                    szWarn = sprintf('%s %s = ''%s''', szWarn, varargin{ixName}, ...
                        varargin{ixName +1});
                else
                    % -
                end
                warning(szWarn)
                
                varargin{ixArg}.(cFn{ii}) = varargin{ixName +1};
            end
        end
    end
end
                    
    
switch def_AgeGroup
    case 'infant' % --- default values for infants
        def_W_t = 2.0*10^6;	% tympanic cavity acoustical losses, in Ns/m^5
        def_V_t = 0.33e-6;  % tympanic cavity volume, in m^3
        def_W_ada = 0;      % aditus ad antrum acoustical losses, in Ns/m^5
        def_M_ada = 0;      % aditus ad antrum acoustical mass, in Ns^2/m^5
        def_V_ant = 0.5e-6; % antrum cavity volume, in m^3
        def_Q_mac = 0;		% mastoid air cells quality factor
        def_f_mac = 0;      % mastoid air cells characteristic frequency, in Hz
        def_V_mac = 0;      % mastoid air cells volume, in m^3
        
        def_stProp_ant = ''; % properties of air in the antrum given as 
                            % struct with fields c: speed of sound in m/s,
                            % and rho: mass density in kg/m^3
        def_stProp_t = '';  % properties of air in the tympanic cavity given as 
                            % struct with fields c: speed of sound in m/s,
                            % and rho: mass density in kg/m^3
                            
        % additional values modeling pathological conditions
        def_x_A = 1;        % factor reducing the eardrum area and the 
                            % tympanic cavity: 
                            % V_tcav = sqrt(x_A) * def_V_tcav
        def_deltaP0 = 0;    % pressure difference of the static air between
                            % ear canal and middle ear
    case 'adult'    % --- default values for adults
        def_W_t = 2.0*10^6;	% tympanic cavity acoustical losses, in Ns/m^5
        def_V_t = 0.5e-6;   % tympanic cavity volume, in m^3
        def_W_ada = 1.7*10^6; % aditus ad antrum acoustical losses, in Ns/m^5
        def_M_ada = 880;	% aditus ad antrum acoustical mass, in Ns^2/m^5
        def_V_ant = 0.8e-6; % antrum cavity volume, in m^3
        def_Q_mac = 0.4;	% mastoid air cells quality factor
        def_f_mac = 3500;	% mastoid air cells characteristic frequency, in Hz
        def_V_mac = 8e-6;   % mastoid air cells volume, in m^3
        
        def_stProp_ant = ''; % properties of air in the antrum given as 
                            % struct with fields c: speed of sound in m/s,
                            % and rho: mass density in kg/m^3
        def_stProp_t = '';  % properties of air in the tympanic cavity given as 
                            % struct with fields c: speed of sound in m/s,
                            % and rho: mass density in kg/m^3
                            
        % additional values modeling pathological conditions
        def_x_A = 1;        % factor reducing the eardrum area and the 
                            % tympanic cavity: 
                            % V_tcav = sqrt(x_A) * def_V_tcav
        def_deltaP0 = 0;    % pressure difference of the static air between
                            % ear canal and middle ear
    otherwise
        error('AgeGroup ''%s'' is not defined in function %s',...
            varargin{ixAG+1},mfilename)
end


def_vOmega = 'returnParameters';
def_dT = 10;    % difference to 300 K (26.85°C) d_T=10 => 36,85°C

addOptional(p,'vOmega',def_vOmega);
addOptional(p,'dT',def_dT,@isnumeric);

addParameter(p,'AgeGroup',def_AgeGroup)

addParameter(p,'W_t',def_W_t, @(x)isnumeric(x)&&isscalar(x));
addParameter(p,'V_t',def_V_t, @(x)isnumeric(x)&&isscalar(x));
addParameter(p,'W_ada',def_W_ada, @(x)isnumeric(x)&&isscalar(x));
addParameter(p,'M_ada',def_M_ada, @(x)isnumeric(x)&&isscalar(x));
addParameter(p,'V_ant',def_V_ant, @(x)isnumeric(x)&&isscalar(x));
addParameter(p,'Q_mac',def_Q_mac, @(x)isnumeric(x)&&isscalar(x));
addParameter(p,'f_mac',def_f_mac, @(x)isnumeric(x)&&isscalar(x));
addParameter(p,'V_mac',def_V_mac, @(x)isnumeric(x)&&isscalar(x));


% --- optional properties: compliances of t and ant allowing different
% mass density and speed of sound
addOptional(p, 'stProp_ant', def_stProp_ant, @validateComplianceProperties);
addOptional(p, 'stProp_t', def_stProp_t, @validateComplianceProperties);


% --- parameters introducing pathological conditions
% -- x_A used to model fluid in the middle ear
validationFcn = @(x) isnumeric(x) && isscalar(x) && (x >= 0) && (x<=1);
addOptional(p, 'x_A', def_x_A, validationFcn);
% -- deltaP0 used to model pressure difference of static air
addOptional(p, 'deltaP0', def_deltaP0, @(x)isnumeric(x)&&isscalar(x));
% ---

% --- parse inputs
parse(p,varargin{:})


% if any value of the model-parameters differs from it's default value, and
% AgeGroup was not explicitly defined
% set AgeGroup to value 'none'
cModelParams = {'W_t';'V_t';'W_ada';'M_ada';'V_ant';'Q_mac';'f_mac';...
    'V_mac'};
vbIsDefault = false(size(cModelParams));
for ii=1:length(cModelParams)
    if any(strcmp(cModelParams{ii},p.UsingDefaults))
        vbIsDefault(ii) = true;
    end
end
if any(~vbIsDefault) && any(strcmp('AgeGroup',p.UsingDefaults))
    stP = p.Results;
    stP.AgeGroup = 'none';
    parse(p, stP)
end

% --- if vOmega is not part of the input, return a list of the parameter values
if ischar(p.Results.vOmega)
    if strcmpi(p.Results.vOmega,def_vOmega)
        vZ_cav = rmfield(p.Results,'vOmega');
        return
    else
        error(['The value of ''vOmega'' is invalid. It must be either a \n',...
            'Nx1 frequency vector in radians, or \n',...
            'the input ''%s'' in order to get a struct of the parameter ',...
            'values.'],def_vOmega)
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

% --- struct of parameters
stCav = rmfield(p.Results,'vOmega');

% --- properties of air
stAir = propertiesOfAir(stCav.dT);

% --- consideration of pathological conditions
% -- deltaP0 used to model pressure difference of static air
if stCav.deltaP0 ~= 0   % do we have this pathological condition
    p0 = stAir.Ks/stAir.gmma; % static air-pressure in the ear canal
    if isstruct(stCav.stProp_t)
        % density of air in the tympanic cavity is defined already for the 
        % healthy condition
        stCav.stProp_t.rho = stCav.stProp_t.rho * (p0 + stCav.deltaP0)/p0;
    else
        % properties of air in the tympanic cavity are not yet defined
        stCav.stProp_t.c = stAir.c;
        stCav.stProp_t.rho = stAir.rho * (p0 + stCav.deltaP0)/p0;
    end
    if isstruct(stCav.stProp_ant)
        % density of air in the antrum is defined already for the healthy
        % condition
        stCav.stProp_ant.rho = stCav.stProp_ant.rho * (p0 + stCav.deltaP0)/p0;
    else
        % properties of air in the tympanic cavity are not yet defined
        stCav.stProp_ant.c = stAir.c;
        stCav.stProp_ant.rho = stAir.rho * (p0 + stCav.deltaP0)/p0;
    end
end
% - compute compliances
% tympanic cavity acoustical compliance, in m^5/N
if isstruct(stCav.stProp_t)
    N_t = stCav.V_t / ...
        (stCav.stProp_t.rho * stCav.stProp_t.c^2);
else
    N_t = stCav.V_t / (stAir.rho * stAir.c^2);
end
% antrum acoustical compliance, in m^5/N
if isstruct(stCav.stProp_ant)
    N_ant = stCav.V_ant / ...
        (stCav.stProp_ant.rho * stCav.stProp_ant.c^2);
else
    N_ant = stCav.V_ant / (stAir.rho * stAir.c^2);
end
% mastoid air cells acoustical compliance, in m^5/N
N_mac = stCav.V_mac / (stAir.rho * stAir.c^2);

if stCav.V_mac == 0
    vY_mac = zeros(size(vOmega));
else
    vY_mac = (vOmega * 1i * N_mac) ./ ...
        (1 - (vOmega ./ (stCav.f_mac * 2 * pi)).^2 ...
        + vOmega * 1i ./ (stCav.Q_mac * stCav.f_mac * 2 * pi));
end
% -- x_A used to model fluid in the middle ear
stCav.V_t = stCav.x_A * stCav.V_t;
% ---

% acoustical admittance of the middle ear cavities
vY_cav = 1 ./ (stCav.W_t + 1 ./ (vOmega * 1i * N_t)) ...
    + 1 ./ (stCav.W_ada + vOmega * 1i * stCav.M_ada + 1 ./...
    (vOmega * 1i * N_ant + vY_mac));
% acoustical impedance of the middle ear cavities
vZ_cav = 1 ./ vY_cav;


end

function bValid = validateComplianceProperties(stProp)
% validates stProp_... for compliances to be a struct with fields 'rho' and
% 'c'
bValid = true;
if isempty(stProp)
    return
end
if ~isstruct( stProp )
    bValid = false;
else
    cNames = {'rho','c'};
    cFn = fieldnames(stProp);
    for ii=1:length(cNames)
        if ~any(strcmp(cNames{ii},cFn))
            bValid = false;
        end
    end
end
if ~bValid
    szMsg = 'stProp.. has to be a struct with fields rho and c.';
    error(szMsg)
end
end

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

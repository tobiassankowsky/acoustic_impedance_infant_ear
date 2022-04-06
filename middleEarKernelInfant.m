function [mK, vY_ac, vA_D, vZ_dmi, vY_mi] = middleEarKernelInfant(varargin)
% 
% Computes the transfer matrix of the kernel element of the middle ear 
% model for young infants' ears, as well as the single acousto-mechanical 
% elements of kernel namely the eardrum acoustic shunt admittance, the 
% effective eardrum area vA_D, the mechanical impedance of eardrum, malleus
% and incus, and the mechanical admittance of the incudomalleal joint.
% 
% SYNTAX:
%   [mK, vY_ac, vA_D, vZ_dmi, vY_mi] = middleEarKernelInfant(vOmega)
%   [mK, vY_ac, vA_D, vZ_dmi, vY_mi] = middleEarKernelInfant(vOmega, Name, Value)
%   [mK, vY_ac, vA_D, vZ_dmi, vY_mi] = middleEarKernelInfant(vOmega, stKernel)
%   stKernel = middleEarKernelInfant
%   stKernel = middleEarKernelInfant(Name, Value)
%   stKernel = middleEarKernelInfant(stKernel)
% 
% DESCRIPTION:
% - middleEarKernelInfant(vOmega,___) returns the transfer matrix of the 
%   middle ear model for the frequencies in radians given in vOmega
% - middleEarKernelInfant or middleEarKernelInfant(___) in which vOmega
%   is ommited returns a struct containing the complete set of parameters 
% 
% INPUTS:   - with no inputs the function returns a struct containing the
%           default parameters
%           - vOmega       frequency in radians
%               - with no other inputs the model will be computed using the 
%               default parameter values
%           - name 'ModelStructure' can be used to compute the kernel
%           model either acc. to [1, 2] (default), or acc. to [3, 4, 5]
%           with one of the following values:
%               - 'Sankowsky-Rothe_etal_2022' or 'S' acc. [1, 2]
%               - 'Hudde_Engel_1998' or 'H' acc. to [3, 4, 5]
%           - name 'AgeGroup' with one of the following values can be 
%           used to select specific default values: 
%               - 'infant' - default acc. to [1, 2], will be used if 
%               AgeGroup is not defined
%               - 'adult' - should be used for the model acc. to [3, 4, 5]
%           - Different values can be used with the name/value syntax
%           and/or with a struct "stKernel". 
%           The names or fieldnames are as follows:
%             'a_ed' 	% physical radius of the eardrum, m
%             'T0_ed'   % tension of the eardrum, Pa
%             'h_ed'    % thickness of the eardrum, m
%             'rho_ed'  % mass density of the eardrum, kg/m^3
%             'N_ac'	% ear drum acoustical compliance, in m^5/N
%             'M_ac0'	% ear drum acoustical mass, in Ns^2/m^5
%             'f_Ym'	% characteristic frequency for ear drum 
%                       % acoustical mass, in Hz
%             'f_Yph'	% characteristic frequency for the phase of
%                       %  the ear drum acoustical admittance, in Hz
%             's_Yph'	% slope for the phase of the ear drum 
%                       % acoustical admittance
%             'W_ac'	% ear drum acoustical losses, in Ns/m^5
% 
%             'A_0',	% ear drum effective area at low frequencies, 
%                       % in m^2. 
%                       % The value is defined as 0.6333 * A_eardrum with
%                       % A_eardrum = 60e-6 m^2
%             'A_inf'	% ear drum effective area at high frequencies, 
%                       % in m^2
%             'Q_A'     % ear drum area function quality factor
%             'f_A'     % ear drum area factor characteristic 
%                       % frequency, in Hz
%             'f_Aph'	% characteristic frequency for the phase of the ear
%                       % drum area function, in Hz
%             's_Aph'	% slope for the phase of the ear drum area function
% 
%             'n_oss'	% ossicle chain mechanical compliance,in m/N
%             'm_oss'	% ossicle chain mechanical mass, in kg
%             'n_cpl'	% drum-malleaus coupling mechanical compliance, 
%                       % in m/N
%             'w_cpl'	% drum-malleaus coupling mechanical % losses, 
%                       % in Ns/m
%             'w_free'	% ear drum mechanical losses, in Ns/m
%             'm_free'	% ear drum mechanical mass, in kg
%             'n_mi'	% malleus-incus joint mechanical  compliance, 
%                       % in m/N
%             'w_mi'	% malleus-incus joint mechanical  losses, in Ns/m
%             
%             'x_A'     % used to model the pathological condition of fluid 
%                       in the middle ear, x_A is the factor of the eardrum
%                       area which is free from fluid, acc. to [2]
%                       (e.g. for the healthy middle ear x_A=1, however, 
%                       for a pathological middle ear x_A=0.2 can be used) 
%             'deltaP0' % used to model the pathological condition of a
%                       pressure difference of the static air between the
%                       ear canal and the middle ear in Pa, acc. to [2]
%                       (e.g. for the healthy middle ear deltaP0=0,
%                       however, for a pathological middle ear
%                       deltaP0=-2500 can be used)
% 
% OUTPUTS:  - mK    [Nx4] transfer matrix of the kernel with 
%                   N=length(vOmega) and the 4 columns K_11,K_12,K_21,K_22
%           - vY_ac     ear drum acoustical shunt admittance
%           - vA_D      effective drum area
%           - vZ_dmi    ear drum + malleus + suspensions mechanical 
%                       impedance + mechanical impedance of the incus
%                       Z_dm = Z_dmi*2/3
%                       Z_i  = Z_dmi/3
%           - vY_mi     mechanical admittance of the incudomalleal joint
% 
% 
%------------------------------- K ------------------------------------------
% K represents the transfer characteristics from the ear drum to the incus.
%
%  / p_delta \       / F_I \
% |           | = K |       |
%  \ q_delta /       \ v_I /
%
% k_11=(1+Y_mi*Z_dm)/A_D
% k_12=(1+Y_mi*Z_i*Z_dm/Z_dmi)*Z_dmi/A_D
% k_21=(1+Y_mi*Z_dm)*Y_ac/A_D+A_D*Y_mi
% k_22=(1+Y_mi*Z_i*Z_dm/Z_dmi)*Y_ac*Z_dmi/A_D + (1+Y_mi*Z_i)*A_D
%
% A_D  ... effective drum area
% Y_ac ... ear drum acoustical shunt admittance
% Z_dm ... ear drum + malleus + suspensions mechanical impedance
% Y_mi ... mechanical admittance of the incudomalleal joint
% Z_i  ... mechanical impedance of the incus
%
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

% --- model structure
% check if ModelStructure is specified as name/value argument
cValidModelStructures = {'Sankowsky-Rothe_etal_2022', 'Hudde_Engel_1998'};
ixMS = find(strcmpi(varargin,'ModelStructure'));
if ~isempty(ixMS)
    def_ModelStructure = validatestring(varargin{ixMS+1}, cValidModelStructures, ...
        mfilename, 'value of ''ModelStructure''');
else
    def_ModelStructure = 'Sankowsky-Rothe_etal_2022';
end
% check if ModelStructure is specified in input struct
vIdc = find(cellfun(@isstruct, varargin));
for iArg = 1:length(vIdc)
    ixArg = vIdc(iArg);
    vIdcAG = strcmp(fieldnames(varargin{ixArg}), 'ModelStructure');
    if any(vIdcAG)
        % check if it was already defined differently using name/value arg.
        if ~isempty(ixMS) && ~strcmpi(def_ModelStructure, ...
                validatestring(varargin{ixArg}.ModelStructure, cValidModelStructures, ...
                mfilename, 'field ''ModelStructure'''))
            error('Inconsistent definition: ModelStructure was defined differently in ''ModelStructure'' and in input-struct')
        end
        def_ModelStructure = validatestring(varargin{ixArg}.ModelStructure, ...
            cValidModelStructures, mfilename, 'field ''ModelStructure''');
    end
end

% --- define default values independent of ModelStructure
def_W_ac=4*10^7;		% ear drum acoustical losses, in Ns/m^5

def_A_0=38.0*10^(-6);   % ear drum effective area at low frequencies, in m^2
def_Q_A=1.3;		% ear drum area function quality factor
def_f_A=2200;		% ear drum area factor characteristic frequency, in Hz

def_n_oss=3.0*10^(-3);	% ossicle chain mechanical compliance, in m/N
def_m_oss=0.007*10^(-3);% ossicle chain mechanical mass, in kg
def_n_cpl=0.5*10^(-3);	% drum-malleaus coupling mechanical compliance, in m/N
def_w_cpl=0.08;         % drum-malleaus coupling mechanical losses, in Ns/m
def_w_free=0.02;		% ear drum mechanical losses, in Ns/m
def_m_free=0.012*10^(-3);	% ear drum mechanical mass, in kg
def_n_mi=0.04*10^(-3);	% malleus-incus joint mechanical compliance, in m/N
def_w_mi=1;             % malleus-incus joint mechanical losses, in Ns/m

% --- define default values depending on ModelStructure
switch def_ModelStructure
    case cValidModelStructures{1}
        % -- physiological based definition of N_ac and M_ac acc. to [1]
        def_a_ed = sqrt(60e-6 / pi);
        def_h_ed = 0.08e-3;    % thickness of eardrum in m (Kuypers et al. 2006)
        def_rho_ed = 1100;    % density of eardrum (Motallebzadeh et al. 2017a)
        % determine tension of drum from Y_ac(160 Hz) of Hudde and Engel model
        imY_acHE160Hz = 4.9363e-09;
        def_T0_ed = pi * def_a_ed^4 *2*pi*160/ (8 * imY_acHE160Hz * def_h_ed);
        def_N_ac = [];%pi * def_a_ed^4 / (8 * def_T0_ed * def_h_ed);
        def_M_ac0 = [];%4 * def_rho_ed * def_h_ed / (3 * pi * def_a_ed^2);
        
        % -- parameters omitted in [1]
        % - Y_ac parameters 
        def_s_Yph = 0;  % drop phase modifications
        def_f_Yph = inf;% drop phase modifications 
        def_f_Ym = inf; % drop magnitude modifications
        % - A_D parameters
        def_A_inf = 0;  % don't use min value for equiv. ED-area
        def_s_Aph = 0;  % drop phase modifications
        def_f_Aph = inf;% drop phase modifications
        
    case cValidModelStructures{2}
        % -- N_ac and M_ac acc. to [3, 4, 5]
        def_N_ac=5*10^(-12);	% ear drum acoustical compliance, in m^5/N
        def_M_ac0=2400;		% ear drum acoustical mass, in Ns^2/m^5
        % --
        % - Y_ac parameters 
        def_s_Yph=1.4*pi;	% slope for the phase of the ear drum 
                            % acoustical admittance
        def_f_Yph=8000;		% characteristic frequency for the phase of the
                            %  ear drum acoustical admittance, in Hz
        def_f_Ym=1900;		% characteristic frequency for ear drum 
                            % acoustical mass, in Hz
        % - A_D parameters
        def_A_inf=2.0*10^(-6);	% ear drum effective area at high 
                                % frequencies, in m^2
        def_s_Aph=-1.2*pi;	% slope for the phase of the ear drum area 
                            % function
        def_f_Aph=1500;		% characteristic frequency for the phase of the
                            % ear drum area function, in Hz
        
        def_a_ed = [];
        def_h_ed = [];
        def_rho_ed = [];
        def_T0_ed = [];
end

% additional values modeling pathological conditions
def_x_A = 1;        % factor reducing the eardrum area 
def_deltaP0 = 0;    % pressure difference of the static air between ear 
                    % canal and middle ear
                    
% --- age group
% check if AgeGroup is specified as name/value argument
cValidAgeGroups = {'infant', 'adult', 'none'};
ixAG = find(strcmpi(varargin,'AgeGroup'));
if ~isempty(ixAG)
    def_AgeGroup = validatestring(varargin{ixAG+1}, cValidAgeGroups, ...
        mfilename, 'value of ''AgeGroup''');
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
        if ~isempty(ixAG) && ~strcmpi(def_AgeGroup, ...
                validatestring(varargin{ixArg}.AgeGroup, cValidAgeGroups, ...
                mfilename, 'field ''AgeGroup'''))
            error('Inconsistent definition: AgeGroup was defined differently in ''AgeGroup'' and in input-struct')
        end
        def_AgeGroup = validatestring(varargin{ixArg}.AgeGroup, ...
                cValidAgeGroups, mfilename, 'field ''AgeGroup''');
    end
end

switch def_AgeGroup
    case cValidAgeGroups{1}
        thickness_factor_infant = 1.4;  % pars flaccida (PF) is 3 times
                                % as large as in adults, assuming PF is 
                                % 1/5 of the drum, the mean thickness 
                                % increased by 0.4
        def_N_ac = def_N_ac / thickness_factor_infant;
        def_M_ac0 = def_M_ac0 * thickness_factor_infant;
        def_W_ac = def_W_ac * 0.6875;
        def_f_A = def_f_A * 2;
        def_Q_A = 1;
        def_m_free = def_m_free * thickness_factor_infant;
        
        def_h_ed = def_h_ed * thickness_factor_infant;
        
    otherwise
        % nothing to do if AgeGroup is either 'adult' or 'none'
end



def_vOmega = 'returnParameters';

addOptional(p,'vOmega',def_vOmega);
addParameter(p,'ModelStructure',def_ModelStructure)
addParameter(p,'AgeGroup',def_AgeGroup)

addParameter(p, 'a_ed', def_a_ed)
addParameter(p, 'h_ed', def_h_ed)
addParameter(p, 'rho_ed', def_rho_ed)
addParameter(p, 'T0_ed', def_T0_ed)

addParameter(p,'W_ac',def_W_ac);
addParameter(p,'N_ac',def_N_ac);
addParameter(p,'M_ac0',def_M_ac0);
addParameter(p,'f_Ym',def_f_Ym);
addParameter(p,'f_Yph',def_f_Yph);
addParameter(p,'s_Yph',def_s_Yph);

addParameter(p,'A_0',def_A_0);
addParameter(p,'A_inf',def_A_inf);
addParameter(p,'Q_A',def_Q_A);
addParameter(p,'f_A',def_f_A);
addParameter(p,'f_Aph',def_f_Aph);
addParameter(p,'s_Aph',def_s_Aph);

addParameter(p,'n_oss',def_n_oss);
addParameter(p,'m_oss',def_m_oss);
addParameter(p,'n_cpl',def_n_cpl);
addParameter(p,'w_cpl',def_w_cpl);
addParameter(p,'w_free',def_w_free);
addParameter(p,'m_free',def_m_free);
addParameter(p,'n_mi',def_n_mi);
addParameter(p,'w_mi',def_w_mi);

% --- parameters introducing pathological conditions
% -- x_A used to model fluid in the middle ear
validationFcn = @(x) isnumeric(x) && isscalar(x) && (x >= 0) && (x<=1);
addOptional(p, 'x_A', def_x_A, validationFcn);
% -- deltaP0 used to model pressure difference of static air
addOptional(p, 'deltaP0', def_deltaP0, @(x)isnumeric(x)&&isscalar(x));
% ---


% --- parse inputs and verify some settings
parse(p,varargin{:})
stP = p.Results;
% -- verify correct name of ModelStructure
stP.ModelStructure = validatestring(stP.ModelStructure, ...
    cValidModelStructures, mfilename, 'field ''ModelStructure''');

% --- verify that pathological conditions are not used with the wrong model
if ~strcmp(stP.ModelStructure, cValidModelStructures{1})
    if ~(stP.x_A == def_x_A && stP.deltaP0 == def_deltaP0)
        error('ModelStructure ''%s'' can''t be used with pathological conditions', stP.ModelStructure)
    end
end

% -- verify correct name of ModelStructure
stP.AgeGroup = validatestring(stP.AgeGroup, ...
    cValidAgeGroups, mfilename, 'field ''AgeGroup''');

% if any value of the age dependent model-parameters differs from its 
% default value, and AgeGroup was not explicitly defined
% set AgeGroup to value 'none'
cModelParams = {'W_ac'; 'N_ac'; 'M_ac0'; 'f_A'; 'Q_A'; 'm_free'};
vbIsDefault = false(size(cModelParams));
for ii=1:length(cModelParams)
    if any(strcmp(cModelParams{ii},p.UsingDefaults))
        vbIsDefault(ii) = true;
    end
end
if any(~vbIsDefault) && any(strcmp('AgeGroup',p.UsingDefaults))
    stP = p.Results;
    stP.AgeGroup = 'none';
%     parse(p, stP)
end
% -- parse corrections
parse(p, stP)


% --- struct of parameters
stKern = rmfield(p.Results,'vOmega');

switch def_ModelStructure
    case cValidModelStructures{1}
        % -- physiological based definition of N_ac and M_ac acc. to [1]
        stKern.N_ac = pi * stKern.a_ed^4 / (8 * stKern.T0_ed * stKern.h_ed);
        stKern.M_ac0 = 4 * stKern.rho_ed * stKern.h_ed / ...
            (3 * pi * stKern.a_ed^2);
end



% --- check vOmega or ruturn a list of the parameter values
if ischar(p.Results.vOmega)
    if strcmpi(p.Results.vOmega,def_vOmega)
        mK = stKern;
        return
    else
        error(['The value of ''vOmega'' is invalid. It must be either a \n',...
            'Nx1 frequency vector in radians, or \n',...
            'the input ''%s'' in order to get a struct of the parameter ',...
            'values.'],def_vOmega)
    end
elseif ~isnumeric(p.Results.vOmega)
    error(['The value of ''vOmega'' is invalid. It must be a \n',...
        'Nx1 frequency vector in radians.'])
elseif min(size(p.Results.vOmega))>1
    error(['The value of ''vOmega'' is invalid. It must be a \n',...
        'Nx1 frequency vector in radians.'])
else
    vOmega = p.Results.vOmega;
    vOmega = vOmega(:);
end


% --- consideration of pathological conditions
% -- deltaP0 used to model pressure difference of static air acc. to [2]
if stKern.deltaP0 ~= 0   % do we have this pathological condition
    A_ed_phys = stKern.a_ed^2 * pi;
    A_0_re_A_phys = stKern.A_0 / A_ed_phys;
    actual_radius_drum = sqrt(stKern.A_0 /A_0_re_A_phys / pi); % effective area / (effective/physical)
    E = 1e8;    % Young's modulus
    s11 = 1/E;
    nu = 0.5; 	% Poisson's ratio
    
    Sprime = sqrt(1 ./ (6 * stKern.a_ed * ...
        (stKern.deltaP0 * stKern.a_ed^2 / (4 * stKern.T0_ed * stKern.h_ed))^2) * ...
        ((stKern.a_ed^2 + 4 * (stKern.deltaP0 * stKern.a_ed^2 / ...
        (4 * stKern.T0_ed * stKern.h_ed))^2)^(3/2) - stKern.a_ed^3)) -1;
    deltaT0 = Sprime / (s11*(1-nu));
    T0Prime = stKern.T0_ed + deltaT0;
    h_drumPrime = stKern.h_ed/(1+Sprime)^2;
    % acoustical compliance
    stKern.N_ac = pi * actual_radius_drum^4 / ...
        (8 * T0Prime * h_drumPrime);
    
end
% -- x_A used to model fluid in the middle ear acc. to [2]
if stKern.x_A ~= def_x_A
    % - decrease ED-area
    stKern.A_0 = stKern.A_0 * stKern.x_A;
    % - acoustical compliance and resistance
    stKern.N_ac = stKern.N_ac * stKern.x_A^2;
    stKern.W_ac = stKern.W_ac * stKern.x_A;
    % - decrease free vibrating ED-mass
    stKern.m_free = stKern.m_free * stKern.x_A;
end




% --- compute model elements
% - Y_ac
vM_ac = stKern.M_ac0 * (1 + sqrt(vOmega ./ (2*pi * stKern.f_Ym)));
vPhi_Y = stKern.s_Yph * log10(1 + vOmega ./ (2*pi * stKern.f_Yph));
vY_ac = 1 ./ (stKern.W_ac + vOmega * 1i .* vM_ac + 1 ./ (vOmega * 1i *...
    stKern.N_ac)) .* exp(vPhi_Y * 1i);

% - A_D
vA_D = (stKern.A_0 + stKern.A_inf) ./ (1-(vOmega ./ (2*pi * stKern.f_A)) .^2 ...
    + 1i * vOmega ./ (2*pi * stKern.f_A * stKern.Q_A)) - stKern.A_inf;
% phase modification acc. to [3, 4, 5], only relevant if ModelStructure 
% 'Hudde_Engel_1998' is used
A_D_f_Aph = (stKern.A_0 + stKern.A_inf) ./ ...
    (1 - (stKern.f_Aph / stKern.f_A)^2 + 1i * stKern.f_Aph /...
    (stKern.f_A * stKern.Q_A)) - stKern.A_inf;
vPhi_A = ones(size(vOmega)) * -1i;
vbIdc = vOmega > stKern.f_Aph*2*pi;
vPhi_A(vbIdc) = stKern.s_Aph * log10(vOmega(vbIdc) ./ (2*pi * stKern.f_Aph)) + ...
            atan2(imag(A_D_f_Aph), real(A_D_f_Aph));
vA_D(vbIdc) = abs(vA_D(vbIdc)) .* exp(vPhi_A(vbIdc)*1i);

% - Z_dmi
vZ_dmi = 1 ./ (vOmega * 1i * stKern.n_oss) + vOmega * 1i * stKern.m_oss ...
    + 1 ./ (1 ./ (vOmega * 1i * stKern.m_free + stKern.w_free) ...
    + 1 ./ (1 ./ (vOmega * 1i * stKern.n_cpl) + stKern.w_cpl));
vZ_dm = vZ_dmi * 2 ./ 3;
vZ_i = vZ_dmi ./ 3;

% - Y_mi
vY_mi = 1 ./ (1 ./ (vOmega * 1i * stKern.n_mi) + stKern.w_mi);

% - transfer matrix of the kernel
vK_11 = (1 + vY_mi .* vZ_dm) ./ vA_D;
vK_12 = (1 + vY_mi .* vZ_i .* vZ_dm ./ vZ_dmi) .* vZ_dmi ./ vA_D;
vK_21 = (1 + vY_mi .* vZ_dm) .* vY_ac ./ vA_D + vA_D .* vY_mi;
vK_22 = (1 + vY_mi .* vZ_i .* vZ_dm ./ vZ_dmi) .* vY_ac .* vZ_dmi ./ vA_D + ...
    (1 + vY_mi .* vZ_i) .* vA_D;
mK = [vK_11, vK_12, vK_21, vK_22];


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

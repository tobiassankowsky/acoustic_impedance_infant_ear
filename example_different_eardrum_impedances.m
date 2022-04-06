% example_different_eardrum_impedances
% 
% Example computing the acoustic input impedance at the eardrum (Z_D) for:
% - the healthy infant ear [1, 2],
% - the infant ear with pathological middle ear condition of fluid in the 
%   middle ear [1, 2],
% - the infant ear with pathological middle ear condition of negative
%   static air pressure difference between middle ear and ear canal (Δp0)
%   [1, 2],
% - and the adult ear according to [3, 4, 5].
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

clear

% -- frequency 
vF = logspace(2,4,400).';   % frequency vector from 100 Hz to 10 kHz
vOmega = 2*pi*vF;           % frequency vector in radians

% --- healthy infant' middle ear 
% - model parameters
stCav = middleEarCavitiesInfant;  % get default cavity parameters
stKernel = middleEarKernelInfant; % get default kernel parameters
stSC = middleEarStapesCochleaInfant;  % get default stapes/cochlea parameters
% - acoustic input impedance at the eardrum
vZ_drum = eardrumImpedanceModelInfant(vOmega, stCav, stKernel, stSC);

% --- fluid in the infant' middle ear 
x_A = 0.2;  % parameter used to model the pathological condition of fluid 
            % in the middle ear
% - model parameters
stCav_fluid = middleEarCavitiesInfant('x_A', x_A);  % get cavity parameters
stKernel_fluid = middleEarKernelInfant('x_A', x_A); % get kernel parameters
stSC_fluid = stSC;  % get stapes/cochlea parameters
% - acoustic input impedance at the eardrum
vZ_drum_fluid = eardrumImpedanceModelInfant(vOmega, stCav_fluid, ...
    stKernel_fluid, stSC_fluid);

% --- static air pressure difference between the infant' middle ear and the
% ear canal
deltaP0 = -2500;    % parameter used to model the pathological condition of 
                    % the pressure difference
% - model parameters of the 
stCav_pressure = middleEarCavitiesInfant('deltaP0', deltaP0);  % get cavity parameters
stKernel_pressure = middleEarKernelInfant('deltaP0', deltaP0); % get kernel parameters
stSC_pressure = stSC;  % get stapes/cochlea parameters
% - acoustic input impedance at the eardrum
vZ_drum_pressure = eardrumImpedanceModelInfant(vOmega, stCav_pressure, ...
    stKernel_pressure, stSC_pressure);

% --- middle ear model of the adult ear acc. to Hudde and Engel (1998)
AgeGroup_HE = 'adult';
ModelStructure_HE = 'Hudde_Engel_1998';
% - model parameters
stCav_HE = middleEarCavitiesInfant('AgeGroup', AgeGroup_HE);  % get cavity parameters
stKernel_HE = middleEarKernelInfant('AgeGroup', AgeGroup_HE, ...
    'ModelStructure', ModelStructure_HE); % get kernel parameters
stSC_HE = stSC;  % get stapes/cochlea parameters
% - acoustic input impedance at the eardrum
vZ_drum_HE = eardrumImpedanceModelInfant(vOmega, stCav_HE, stKernel_HE, stSC_HE);


% --- plotting
figure, 
subplot(2, 1, 1)
semilogx(vF, db(vZ_drum), 'DisplayName', 'infant, healthy')
hold on
semilogx(vF, db(vZ_drum_fluid), 'DisplayName', ...
    sprintf('infant, fluid: x_A=%g', x_A))
semilogx(vF, db(vZ_drum_pressure), 'DisplayName', ...
    sprintf('infant, \\Deltap_0=%gPa', deltaP0))
semilogx(vF, db(vZ_drum_HE), 'DisplayName', 'adult, Hudde and Engel 1998')
ylabel('Z_D in dB')
legend('show', 'location', 'southwest')
subplot(2, 1, 2)
semilogx(vF, 180/pi*angle(vZ_drum))
hold on
semilogx(vF, 180/pi*angle(vZ_drum_fluid), 'DisplayName', ...
    sprintf('infant, fluid: x_A=%g', x_A))
semilogx(vF, 180/pi*angle(vZ_drum_pressure), 'DisplayName', ...
    sprintf('infant, \\Deltap_0=%gPa', deltaP0))
semilogx(vF, 180/pi*angle(vZ_drum_HE), 'DisplayName', 'adult, Hudde and Engel 1998')
ylim([-180 180])
ylabel('Z_D in degrees')
xlabel('Frequency in Hz')


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


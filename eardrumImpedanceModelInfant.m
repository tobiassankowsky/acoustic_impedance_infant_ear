function vZ_drum = eardrumImpedanceModelInfant(vOmega, stCav, stKernel, stSC)
% 
% Computes the acoustic input impedance at the eardrum according to [1, 2,
% 3, 4, 5]
% 
% SYNTAX:
%   vZ_drum = eardrumImpedanceModelInfant(vOmega, stCav, stKernel, stSC)
% 
% INPUTS:   - vOmega    frequency in radians
%           - stCav     Struct containing parameters of the middle ear 
%                       cavities as used in function 
%                       'middleEarCavitiesInfant'. For more information 
%                       type help middleEarCavitiesInfant or call the
%                       function without input arguments.
%           - stKernel  Struct containing parameters of the kernel part of
%                       the middle ear model as used in function 
%                       'middleEarKernelInfant'. For more information 
%                       type help middleEarKernelInfant or call the
%                       function without input arguments.
%           - stSC      Struct containing parameters of the stapes and
%                       cochlea as used in function 
%                       'middleEarStapesCochleaInfant'. For more information 
%                       type help middleEarStapesCochleaInfant or call the
%                       function without input arguments.
% 
% OUTPUTS:  vZ_drum     Vector of the complex frequency dependent
%                       acoustic input impedance at the eardrum.
% 
% EXAMPLE:
% % -----------------------------
% stCav = middleEarCavitiesInfant;  % get default cavity parameters
% stKernel = middleEarKernelInfant; % get default kernel parameters
% stSC = middleEarStapesCochleaInfant;  % get default stapes/cochlea parameters
% vF = logspace(2,4,400).'; % frequency vector from 100 Hz to 10 kHz
% vOmega = 2*pi*vF;     % frequency vector in radians
% % - acoustic input impedance at the eardrum
% vZ_drum = eardrumImpedanceModelInfant(vOmega, stCav, stKernel, stSC);
% % - plot 
% figure, 
% subplot(2, 1, 1), semilogx(vF, db(vZ_drum))
% ylabel('Z_D in dB')
% subplot(2, 1, 2), semilogx(vF, 180/pi*angle(vZ_drum))
% ylim([-180 180])
% ylabel('Z_D in degrees')
% xlabel('Frequency in Hz')
% % -----------------------------
% 
% 
% REQUIRED FUNCTIONS:
%           - middleEarCavitiesInfant.m
%               - propertiesOfAir.m
%           - middleEarKernelInfant.m
%           - middleEarStapesCochleaInfant.m
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

% cavities
vZ_cav = middleEarCavitiesInfant(vOmega,stCav);
% kernel element
[~,vY_ac,vA_D,vZ_dmi,vY_mi]= middleEarKernelInfant(vOmega,stKernel);
% stapes and cochlea
vZ_Load = middleEarStapesCochleaInfant(vOmega,stSC);

% compute impedance in mechanical domain
vz_mech = 2/3.*vZ_dmi + 1./(vY_mi + 1./(vZ_dmi./3 + vZ_Load));
% compute impedance in acoustical domain
vZ_mech = vz_mech./(vA_D.^2);
% compute imedance at the eardrum
vZ_drum = vZ_cav + 1./(vY_ac + 1./vZ_mech);


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


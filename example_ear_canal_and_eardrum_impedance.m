% example_ear_canal_and_eardrum_impedance
% 
% Example computing the acoustic input impedance at the eardrum (Z_D) and
% at hte ear canal entrance (Z_ec) for the healthy infant ear [1, 2].
% 
% REFERENCES:
% [1] Sankowsky-Rothe et al. (2022), "Parametric model of young infants' 
%     eardrum and ear canal impedances supporting immittance measurement 
%     results. Part I: Development of the model", Acta Acustica
% [2] Sankowsky-Rothe et al. (2022), "Parametric model of young infants' 
%     eardrum and ear canal impedances supporting immittance measurement 
%     results. Part II: Prediction of eardrum and ear canal impedances for 
%     common middle ear pathologies", Acta Acustica
% 
% AUTHOR:   Tobias Sankowsky-Rothe
% DATE:     2022
% LICENSE:  see EOF

clear
close all

% -- frequency 
vF = logspace(2,4,400).';   % frequency vector from 100 Hz to 10 kHz
vOmega = 2*pi*vF;           % frequency vector in radians
nF = length(vF);

% -- ear canal parameters
stEC.r_ec = 1.7e-3; % radii of ear canal segments/slices in m
stEC.l_ec = 14e-3;  % length of ear canal segments/slices in m
stEC.dT_ec = 10;    % temperature difference to 26,85°C (300°K) in °
stEC.g_ec = 10;     % additional damping factor for the real part of the 
                    % propagation
stEC.zeta = 0.6;    % damping ratio of the ear canals soft tissue
stEC.dst = 3e-3;    % thickness of the ear canals soft tissue in m
stEC.Est = 90e3;    % Young's modulus of the ear canals soft tissue in Pa
stEC.rhoSt = 1100;  % mass density of the ear canals soft tissue in kg/m^3


% -- eardrum impedance
% - model parameters
stCav = middleEarCavitiesInfant;  % get default cavity parameters
stKernel = middleEarKernelInfant; % get default kernel parameters
stSC = middleEarStapesCochleaInfant;  % get default stapes/cochlea parameters
% - acoustic input impedance at the eardrum
vZ_drum = eardrumImpedanceModelInfant(vOmega, stCav, stKernel, stSC);

% -- ear canal model
% compliant wall
mZWall = compliantEcWallImpedance...
    (vF, stEC.l_ec, stEC.r_ec, stEC.zeta, ...
    stEC.dst, stEC.Est, stEC.rhoSt);
vAc = ones(length(vF), 1);
vBc = zeros(length(vF), 1);
vCc = zeros(length(vF), 1);
vDc = ones(length(vF), 1);

for iSlice = 1:length(stEC.r_ec) % each ec-slice
    % transfer matrix of canal
    [vA1,vB1,vC1,vD1] = ...
        chainMatTubeBenade(vF, stEC.r_ec(iSlice), ...
        stEC.l_ec(iSlice), stEC.dT_ec, stEC.g_ec);
    % transfer matrix of [canal] [compliant wall]
    [vA1,vB1,vC1,vD1] = mult2p(vA1, vB1, vC1, vD1,...
        ones(nF,1), zeros(nF,1), 1./mZWall(:,iSlice), ones(nF,1));
    
    [vAc,vBc,vCc,vDc] = mult2p(vAc, vBc, vCc, vDc, vA1, vB1, vC1, vD1);
end

% -- Zec: acoustic input impedance of the ear canal
vZ_ec = (vAc .* vZ_drum + vBc) ./ (vCc .* vZ_drum + vDc);

% --- plotting
figure, 
subplot(2, 1, 1)
semilogx(vF, db(vZ_drum), 'DisplayName', 'Z_{D}')
hold on
semilogx(vF, db(vZ_ec), 'DisplayName', 'Z_{ec}')
ylabel('Z in dB')
legend('show', 'location', 'southwest')
subplot(2, 1, 2)
semilogx(vF, 180/pi*angle(vZ_drum), 'DisplayName', 'Z_{D}')
hold on
semilogx(vF, 180/pi*angle(vZ_ec), 'DisplayName', 'Z_{ec}')
ylim([-180 180])
ylabel('Z in degrees')
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

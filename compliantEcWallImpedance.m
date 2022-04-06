function [mZWall, mZRst, mZComp, mZMass] = compliantEcWallImpedance...
    (vF, vL, vR, zeta, dst, Est, rhoSt, alphaWeight)
% 
% [mZWall, mZRst, mZComp, mZMass] = compliantEcWallImpedance...
%     (vF, vL, vR, zeta, dst, Est, rhoSt)
% Computes the acoustic impedance of the compliant wall of an infants ear canal
% 
% INPUTS:
%   vF      frequency vector in Hz
%   vL      vector of ECs slices-/segment-lengths from inside (tip) to
%           outside (entrance)
%   vR      vector of ECs slices-/segment-radii from inside (tip) to
%           outside (entrance)
%   zeta    damping ratio of the ear canals soft tissue
%   dst     thickness of the ear canals soft tissue
%   Est     Young's modulus of the ear canals soft tissue
%   rhoSt   mass density of the ear canals soft tissue
%   alphaWeight     (optional) used for weighting of damping constants
%                   alpha and beta, using alphaWeigth = alpha/(omega zeta).
%                   The default value is alphaWeight = 1.5.
% 
% OUTPUTS:
%   mZWall  matrix NxM of the EC-wall-impedance at N frequencies and M
%           segments
%   mZRst   matrix NxM of the resistance part of the EC-wall-impedance at 
%           N frequencies and M segments
%   mZComp  matrix NxM of the compliance part of the EC-wall-impedance at 
%           N frequencies and M segments
%   mZMass  matrix NxM of the mass part of the EC-wall-impedance at 
%           N frequencies and M segments
% 
% AUTHOR:   Tobias Sankowsky-Rothe
% DATE:     2020-09-03
% LICENSE:  see EOF

% CHANGELOG:
%   2021-07-08: - introduces optional input alphaWeight
%   2021-07-20: - changes default value of alphaWeight to 1

% |--- 2021-07-08 ---
if nargin < 8
    % |--- 2021-07-20 ---
%     alphaWeight = 1.5;
    % --- 2021-07-20 ---
    alphaWeight = 1;
    % --- 2021-07-20 ---|
end
% --- 2021-07-08 ---|

vF = vF(:);
vL = vL(:);
vR = vR(:);
nF = length(vF);
nSlices = length(vL);

mZRst = zeros(nF, nSlices);
mZComp = complex(zeros(nF, nSlices));
mZMass = mZComp;
mZWall = mZComp;

% |--- 2021-07-08 ---
vAlpha = (2*pi*vF)*alphaWeight*zeta;
% --- 2021-07-08 ---
% vAlpha = (2*pi*vF)*1.5*zeta;
% --- 2021-07-08 ---|
vBeta = 2 * zeta./(2*pi*vF) - vAlpha ./(2*pi*vF).^2;
% areas of soft tissue (surface of ear canal walls)
vAst = 2*pi.*vR.*vL;
% compliances of soft tissue
vCst = dst * vAst./Est;
% o-------------o
%       |
%     1/(jwC)
%       |
% o-------------o
%   1   0
%  1/Z  1

% masses of soft tissue
vMst = rhoSt*dst./vAst;

% resistance of the soft tissue
% zeta = 1/2*(alpha/w + beta*w)
% vBetaTimesOmega = 2*zeta;
% C = alpha M + beta K
% the stiffness of the spring K is given by Est*vAst/dst
%     vK = Est*vAst/dst;
%     vRst = vAlpha*vMst + vBeta*vK;
for iSlice = 1:nSlices
    mZRst(:,iSlice) = vAlpha*vMst(iSlice) + vBeta/vCst(iSlice);
    
    mZComp(:,iSlice) = 1./(2i*pi*vF*vCst(iSlice));
    mZMass(:,iSlice) = 2i*pi*vF*vMst(iSlice);
end
mZWall = mZRst + mZComp + mZMass;


%--------------------Licence ---------------------------------------------
% Copyright (C) 2020  Tobias Sankowsky-Rothe
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

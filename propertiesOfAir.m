function [stAir]=propertiesOfAir(dT)
% 
% computes temperature dependet acoustic properties of air
%   valid at 1 atm (101.325 kPa), at -10 \leq dT \leq 10
% 
% INPUT
%           dT  difference to 300 K (26.85°C)
% 
% OUTPUT
%           Struct with the following fields
%           - rho   mass density, in kg/m^3
%           - eta   shear viscosity coeff,in kg/(ms)
%           - c     speed of sound, in m/s
%           - gmma  heat capacity ratio or adiabatic index or ratio of specific heats
%           - nu    square root of the Prandtl number
%           - p0    isentropic bulk modulus [Pa]
% 
%           static pressure is approximately given by 
%           p0/gmma = rho*c^2/gmma
% 
% References
%   Benade (1968) 
%   "On the Propagation of Sound Waves in a Cylindrical Conduit"
%   JASA, vol. 44, no. 2, pp. 616-623 
% 
% AUTHOR:   Tobias Sankowsky-Rothe
% DATE:     2017
% LICENSE:  see EOF

% Changes:
%   - 2021-04-13 Fix - corrects the fieldname of the isentropic bulk
%   modulus


% -- properties of air
stAir.rho = 1.1769 * (1 - 0.00335*dT);          % mass density, in kg/m^3
stAir.eta = 1.846*10^(-5) * (1 + 0.0025*dT);    % shear viscosity coeff,in kg/(m s)
stAir.c = 347.23 * (1 + 0.00166*dT);            % speed of sound, in m/s
stAir.gmma = 1.4017 * (1 - 0.00002*dT);         % heat capacity ratio or adiabatic index or ratio of specific heats
stAir.nu = 0.8410 * (1 - 0.0002*dT);            % square root of the Prandtl number
% |--- 2021-04-13 ---
% stAir.p0 = stAir.rho*stAir.c^2;           % isentropic bulk modulus [Pa]
% --- 2021-04-13 ---
stAir.Ks = stAir.rho * stAir.c^2;           % isentropic bulk modulus [Pa]
% --- 2021-04-13 ---|

%--------------------Licence ---------------------------------------------
% Copyright (C) 2017  Tobias Sankowsky-Rothe
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

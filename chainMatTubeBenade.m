function [A,B,C,D] = chainMatTubeBenade(freq, r, l, dT, real_propagation_fact)
% Calculation of the chain parameter of a cylindrical tube.
% 
% [A,B,C,D] = chainMatTubeBenade(freq, r, l, dT, real_propagation_fact)
% Calculates the 2-port chain matrix of a cylindrical tube with given
% radius and length according to [1]
% 
% |p1|   |A  B|      |p0|
% |  | = |    |  *   |  |
% |q1|   |C	 D|      |q0|
% 
% q0                   q1
% ->                   ->
% -----------------------
% |                     |
% |p0                   |p1
% V                     V
% -----------------------
% 
% input:
%       freq    = frequency vector
%       r       = radius of the tube
%       l       = length of the tube
%       dt      = temperature difference to T0 = 26,85°C (300°K)
%   optional
%       real_propagation_fact   = additional damping factor for the real 
%                   part of the propagation constant acc. to [2]. The
%                   characteristic impedance Zw is than approximated by 
%                   Zw_approx = (rho c)/A * gamma / (j omega/c)
% 
%	References:     [1] Benade, A.H.; On the propagation of sound waves in a 
%                   cylindrical conduit. JASA 44(2):616-623, 1968.
%                   [2] Hudde & Engel (1998). 
%                   Measuring and Modeling Basic Properties of the Human 
%                   Middle Ear and Ear Canal. Part II: Ear Canal, Middle 
%                   Ear Cavities, Eardrum, and Ossicles" 
%                   ACUSTICA/acta acustica, vol. 84, no. 5, pp. 894-913
%                   
% 
%   Alternatives:
%       An alternative modeling of additional damping was proposed in  
%       Vogl, S. & Blau, M., JASA 145(2):917-930, 2019.
%       The method is implemented in "tube_benade_lossmult.m". In this
%       method the damping factor "lossmult" is multiplied with the real
%       part of the series impedance Z of the tube. The characteristic
%       impedance of the tube is given by Zw = sqrt(Z/Y) and the 
%       propagation constant is given by gamma = sqrt(Z*Y), with the shunt 
%       admittance Y of the tube. A comparable damping can achieved by
%       choosing "lossmult = sqrt(2)*real_propagation_fact".
% 
% Author:   Tobias Sankowsky-Rothe
% Date  :   2015
% License: see EOF
% modified:
%   2015-09-17 ts: replace zeros in frequncy vector with eps, to avoid
%           devision by zero


if nargin<5, real_propagation_fact = 1; end

freq = freq(:);
% |--- 2015-09-17 ts ---
freq(freq==0) = eps;
% --- 2015-09-17 ts ---|
dm = 2*r;

% dT = 10;                            % T0 = 26,85°C (300°K)
rho = 1.1769 * (1 - 0.00335 * dT);          % mass density, in kg/m^2
eta = 1.846 * 10^(-5) * (1 + 0.0025 * dT);  % shear viscosity coeff,in kg/(ms)
c = 347.23 * (1 + 0.00166 * dT); 	    % speed of sound, in m/s
gmma = 1.4017 * (1 - 0.00002 * dT);
nu = 0.8410 * (1 - 0.0002 * dT);            % square root of the Prandtl number

vOmega = 2*pi * freq;

% viscous and isothermal conditions acc. to Benade 1968
r_nu = sqrt(vOmega * rho / eta) * r;
r_t = nu * r_nu;
F_nu = besselj(1, r_nu * sqrt(-1i)) ./ besselj(0, r_nu * sqrt(-1i)) ./ ...
    (r_nu * sqrt(-1i)) * 2;
F_t = besselj(1, r_t * sqrt(-1i)) ./ besselj(0, r_t * sqrt(-1i)) ./ ...
    (r_t * sqrt(-1i)) * 2;
R0 = (rho * c) / (pi * r^2);
Z = (vOmega * 1i / c * R0) ./ (-F_nu + 1);
Y = (vOmega * 1i / c / R0) .* (F_t * (gmma - 1) + 1);

G = sqrt(Z .* Y);
if real_propagation_fact ~= 1
    G = real_propagation_fact * real(G) + 1i * imag(G);
    % Zw = rho*c/(pi*r^2) .* ((1i*vOmega) ./ (G*c)); % as given, but assumed to be wrong
    Zw = rho * c / (pi * r^2) .* ((G * c) ./ (1i * vOmega));
else
    Zw = sqrt(Z ./ Y);
end


A = cosh(G * l);
B = sinh(G * l) .* Zw;
C = sinh(G * l) ./ Zw;
D = A;

%--------------------Licence ---------------------------------------------
% Copyright (C) 2015  Tobias Sankowsky-Rothe
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

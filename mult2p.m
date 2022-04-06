function [A,B,C,D]=mult2p(A1,B1,C1,D1,A2,B2,C2,D2)
% two-port parameter multiplication
% 
% input:
%       either parameters of the two two-port as vektors in the order
%           A11,A12,A21,A22,B11,B12,B21,B22
%       or two Nx4-matrices each with the two-port-parameters
% 
% output:
%       if nargout == 1, the function returns the resulting two-port as a
%           Nx4-matrix
%       if nargout > 1 up to 4, the function returns the two-port
%           parameters as vectors in the order C11,C12,C21,C22
% 
% Author :	Tobias Sankowsky
% Date:     2008
% License: see EOF

if nargin == 2
    if size(A1,2)~=4 || size(B1,2)~=4
        error(['if there are only two arguments (A and B) to mult2p,\n',...
            'each argument has to be a Nx4 matrix where the columns\n',...
            'corresponds to the matrix parameters [11,12,21,22]\n'])
    end
    A = A1(:,1).*B1(:,1) + B1(:,3).*A1(:,2);
    B = A1(:,1).*B1(:,2) + A1(:,2).*B1(:,4);
    C = A1(:,3).*B1(:,1) + A1(:,4).*B1(:,3);
    D = A1(:,3).*B1(:,2) + A1(:,4).*B1(:,4);
else
    A = A1.*A2 + C2.*B1;
    B = A1.*B2 + B1.*D2;
    C = C1.*A2 + D1.*C2;
    D = C1.*B2 + D1.*D2;
end

if nargout==1
    A = [A(:),B(:),C(:),D(:)];
end

%--------------------Licence ---------------------------------------------
% This file is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License, 
% or (at your option) any later version.
% See the GNU General Public License for more details:
% http://www.gnu.org/licenses/gpl

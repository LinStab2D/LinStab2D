% Cylindrical Coordinate System 
%           +r
%           |
%           |
%           |______+z (axial)
%          /
%         /
%        /
%       /
%      theta (anti clock-wise)


% Continuity
% vr/r + dvr/dr + i m/r vtheta + dvz/dz = 0
A0_11 = 1./R; % vr/r
A0_12 = 1i.*m./R; % i m/r vtheta
A0_13 = Z;
A0_14 = Z;

Ar_11 = I; % dvr/dr
Ar_12 = Z;
Ar_13 = Z;
Ar_14 = Z;

Az_11 = Z;
Az_12 = Z;
Az_13 = I; % dvz/dz
Az_14 = Z;

Arr_11 = Z;
Arr_12 = Z;
Arr_13 = Z;
Arr_14 = Z;

Azz_11 = Z;
Azz_12 = Z;
Azz_13 = Z;
Azz_14 = Z;

% Radial momentum
%      1            2            3            4              5           6             7               8               9            10       11         12        13        14
% -i omega vr = -vr dVr/dr - Vr dvr/dr - imVtheta/r vr - Vz dvr/dz - dVr/dz vz + 2 Vtheta vtheta/r - dp/dr + 1/Re (1/r dvr/dr + d2vr/dr2 - m2/r2 vr + d2vrdz2 - vr/r2 - 2im/r2 vtheta)
A0_21 = - dVrdr - 1i.*m.*Vtheta./R - 1./Re.*(m.^2./R.^2 + 1./R.^2); % -vr dVr/dr - imVtheta/r vr + 1/Re(- m2/r2 vr - vr/r2 ) [2, 4, 11 and 13]
A0_22 = 2.*Vtheta./R - 1./Re.*(2.*1i.*m./R.^2); % 2 Vtheta vtheta/r - 1/Re 2im/r2 vtheta [7 and 14]
A0_23 = -dVrdz; % - dVr/dz vz [6]
A0_24 = Z;

Ar_21 = -Vr + 1./(Re.*R); %  - Vr dvr/dr + 1/Re 1/r dvr/dr [3 and 9]
Ar_22 = Z;
Ar_23 = Z;
Ar_24 = -I; % - dp/dr [8]

Az_21 = -Vz; % - Vz dvr/dz [5]
Az_22 = Z;
Az_23 = Z;
Az_24 = Z;

Arr_21 = 1./Re.*I; % 1/Re d2vr/dr2 [10]
Arr_22 = Z;
Arr_23 = Z;
Arr_24 = Z;

Azz_21 = 1./Re.*I; % 1/Re d2vrdz2 [12] 
Azz_22 = Z;
Azz_23 = Z;
Azz_24 = Z;

% Azimuthal momentum
%         1               2                3                4                 5               6             7               8           9                    10               11             12             13          14          15
% -i omega vtheta = -vr dVtheta/dr - Vr dvtheta/dr - imVtheta/r vtheta - Vz dvtheta/dz - dVtheta/dz vz  - Vtheta vr/r - Vr vtheta/r - im/r p + 1/Re (1/r dvtheta/dr + d2vtheta/dr2 - m2/r2 vtheta + d2vthetadz2 - vtheta/r2 + 2im/r2 vr)
A0_31 = -dVthetadr - Vtheta./R + 1./Re.*(2.*1i.*m./R.^2); % [2, 7 and 15]
A0_32 = -Vtheta./R.*1i.*m - Vr./R + 1./Re.*(-m.^2./R.^2 - 1./R.^2); % [4, 8, 12 and 14]
A0_33 = -dVthetadz; % [6]
A0_34 = -1i.*m./R; % [9]

Ar_31 = Z;
Ar_32 = -Vr + 1./Re.*1./R; % [3 and 10]
Ar_33 = Z;
Ar_34 = Z;

Az_31 = Z;
Az_32 = -Vz; % [5]
Az_33 = Z;
Az_34 = Z;

Arr_31 = Z;
Arr_32 = 1./Re.*I; % [11]
Arr_33 = Z;
Arr_34 = Z;

Azz_31 = Z;
Azz_32 = 1./Re.*I; % [13]
Azz_33 = Z;
Azz_34 = Z;

% Axial momentum
%      1            2            3            4              5           6           7               8           9         10         11
% -i omega vz = -vr dVz/dr - Vr dvz/dr - imVtheta/r vz - Vz dvz/dz - dVz/dz vz  - dp/dz + 1/Re (1/r dvz/dr + d2vz/dr2 - m2/r2 vz + d2vzdz2)
A0_41 = -dVzdr; % [2]
A0_42 = Z;
A0_43 = -Vtheta./R.*1i.*m - dVzdz + 1./Re.*(-m.^2./R.^2) ; % [4, 6 and 10]
A0_44 = Z;

Ar_41 = Z;
Ar_42 = Z;
Ar_43 = -Vr + 1./Re.*1./R; % [3 and 8]
Ar_44 = Z;

Az_41 = Z;
Az_42 = Z;
Az_43 = -Vz; % [5]
Az_44 = -I; % [7]

Arr_41 = Z;
Arr_42 = Z;
Arr_43 = 1./Re.*I; % [9]
Arr_44 = Z;

Azz_41 = Z;
Azz_42 = Z;
Azz_43 = 1./Re.*I; % [11]
Azz_44 = Z;
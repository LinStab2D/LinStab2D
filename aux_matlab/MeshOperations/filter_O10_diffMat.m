function [M,D] = filter_O10_diffMat(nx,alpha_f)

beta_f  = 1.0;
gamma_f = alpha_f;
f_a     = 0.5 * (  1.0  -       2.0     *alpha_f)/512.0;
f_b     = 0.5 * ( -5.0  +       10.0    *alpha_f)/256.0;
f_c     = 0.5 * ( 45.0  -       90.0    *alpha_f)/512.0;
f_d     = 0.5 * (-15.0  +       30.0    *alpha_f)/64.0;
f_e     = 0.5 * (105.0  +       302.0   *alpha_f)/256.0;
f_f     =       (193.0  +       126.0   *alpha_f)/256.0;                % central point
f_g     = f_e;
f_h     = f_d;
f_i     = f_c;
f_j     = f_b;
f_k     = f_a;

f_dw5_b = ( -1.0        +   2.0         *alpha_f)/1024.0;
f_dw5_c = (  5.0        -  10.0         *alpha_f)/512.0;
f_dw5_d = (-45.0        +  90.0         *alpha_f)/1024.0;
f_dw5_e = ( 15.0        +  98.0         *alpha_f)/128.0;
f_dw5_f = (407.0        + 210.0         *alpha_f)/512.0;                % central point
f_dw5_g = ( 63.0        + 130.0         *alpha_f)/256.0;
f_dw5_h = (-105.0       + 210.0         *alpha_f)/512.0;
f_dw5_i = ( 15.0        -  30.0         *alpha_f)/128.0;
f_dw5_j = (-45.0        +  90.0         *alpha_f)/1024.0;
f_dw5_k = (  5.0        -  10.0         *alpha_f)/512.0;
f_dw5_l = ( -1.0        +   2.0         *alpha_f)/1024.0;

f_dw4_c = (  1.0        -   2.0         *alpha_f)/1024.0;
f_dw4_d = ( -5.0        +  10.0         *alpha_f)/512.0;
f_dw4_e = ( 45.0        + 934.0         *alpha_f)/1024.0;
f_dw4_f = (113.0        +  30.0         *alpha_f)/128.0;                % central point
f_dw4_g = (105.0        + 302.0         *alpha_f)/512.0;
f_dw4_h = (-63.0        + 126.0         *alpha_f)/256.0;
f_dw4_i = (105.0        - 210.0         *alpha_f)/512.0;
f_dw4_j = (-15.0        +  30.0         *alpha_f)/128.0;
f_dw4_k = ( 45.0        -  90.0         *alpha_f)/1024.0;
f_dw4_l = ( -5.0        +  10.0         *alpha_f)/512.0;
f_dw4_m = (  1.0        -   2.0         *alpha_f)/1024.0;

f_dw3_d = ( -1.0        +   2.0         *alpha_f)/1024.0;
f_dw3_e = (  5.0        + 502.0         *alpha_f)/512.0;
f_dw3_f = (979.0        +  90.0         *alpha_f)/1024.0;               % central point
f_dw3_g = ( 15.0        +  98.0         *alpha_f)/128.0;
f_dw3_h = (-105.0       + 210.0         *alpha_f)/512.0;
f_dw3_i = ( 63.0        - 126.0         *alpha_f)/256.0;
f_dw3_j = (-105.0       + 210.0         *alpha_f)/512.0;
f_dw3_k = (  15.0       -  30.0         *alpha_f)/128.0;
f_dw3_l = (-45.0        +  90.0         *alpha_f)/1024.0;
f_dw3_m = (  5.0        -  10.0         *alpha_f)/512.0;
f_dw3_n = ( -1.0        +   2.0         *alpha_f)/1024.0;

f_dw2_e = (  1.0        +1022.0         *alpha_f)/1024.0;
f_dw2_f = (507.0        +  10.0         *alpha_f)/512.0;                % central point
f_dw2_g = ( 45.0        + 934.0         *alpha_f)/1024.0;
f_dw2_h = (-15.0        +  30.0         *alpha_f)/128.0;
f_dw2_i = (105.0        - 210.0         *alpha_f)/512.0;
f_dw2_j = (-63.0        + 126.0         *alpha_f)/256.0;
f_dw2_k = (105.0        - 210.0         *alpha_f)/512.0;
f_dw2_l = (-15.0        +  30.0         *alpha_f)/128.0;
f_dw2_m = ( 45.0        -  90.0         *alpha_f)/1024.0;
f_dw2_n = ( -5.0        +  10.0         *alpha_f)/512.0;
f_dw2_o = (  1.0        -   2.0         *alpha_f)/1024.0;

f_dw1_f = (1023.0       +   1.0         *alpha_f)/1024.0;               % central point
f_dw1_g = (  5.0        + 507.0         *alpha_f)/512.0;
f_dw1_h = (-45.0        +  45.0         *alpha_f)/1024.0;
f_dw1_i = ( 15.0        -  15.0         *alpha_f)/128.0;
f_dw1_j = (-105.0       + 105.0         *alpha_f)/512.0;
f_dw1_k = ( 63.0        -  63.0         *alpha_f)/256.0;
f_dw1_l = (-105.0       + 105.0         *alpha_f)/512.0;
f_dw1_m = ( 15.0        -  15.0         *alpha_f)/128.0;
f_dw1_n = (-45.0        +  45.0         *alpha_f)/1024.0;
f_dw1_o = (  5.0        -   5.0         *alpha_f)/512.0;
f_dw1_p = ( -1.0        +   1.0         *alpha_f)/1024.0;

diag_D      = kron([f_a f_b f_c f_d f_e f_f f_g f_h f_i f_j f_k], ones(nx,1));
D           = spdiags(diag_D,-5:5,nx,nx);
D(1,1:11)   = [f_dw1_f f_dw1_g f_dw1_h f_dw1_i f_dw1_j f_dw1_k f_dw1_l f_dw1_m f_dw1_n f_dw1_o f_dw1_p]; 
D(2,1:11)   = [f_dw2_e f_dw2_f f_dw2_g f_dw2_h f_dw2_i f_dw2_j f_dw2_k f_dw2_l f_dw2_m f_dw2_n f_dw2_o]; 
D(3,1:11)   = [f_dw3_d f_dw3_e f_dw3_f f_dw3_g f_dw3_h f_dw3_i f_dw3_j f_dw3_k f_dw3_l f_dw3_m f_dw3_n];
D(4,1:11)   = [f_dw4_c f_dw4_d f_dw4_e f_dw4_f f_dw4_g f_dw4_h f_dw4_i f_dw4_j f_dw4_k f_dw4_l f_dw4_m];
D(5,1:11)   = [f_dw5_b f_dw5_c f_dw5_d f_dw5_e f_dw5_f f_dw5_g f_dw5_h f_dw5_i f_dw5_j f_dw5_k f_dw5_l];

D(nx  ,end-10:end)   = fliplr([f_dw1_f f_dw1_g f_dw1_h f_dw1_i f_dw1_j f_dw1_k f_dw1_l f_dw1_m f_dw1_n f_dw1_o f_dw1_p]); 
D(nx-1,end-10:end)   = fliplr([f_dw2_e f_dw2_f f_dw2_g f_dw2_h f_dw2_i f_dw2_j f_dw2_k f_dw2_l f_dw2_m f_dw2_n f_dw2_o]); 
D(nx-2,end-10:end)   = fliplr([f_dw3_d f_dw3_e f_dw3_f f_dw3_g f_dw3_h f_dw3_i f_dw3_j f_dw3_k f_dw3_l f_dw3_m f_dw3_n]);
D(nx-3,end-10:end)   = fliplr([f_dw4_c f_dw4_d f_dw4_e f_dw4_f f_dw4_g f_dw4_h f_dw4_i f_dw4_j f_dw4_k f_dw4_l f_dw4_m]);
D(nx-4,end-10:end)   = fliplr([f_dw5_b f_dw5_c f_dw5_d f_dw5_e f_dw5_f f_dw5_g f_dw5_h f_dw5_i f_dw5_j f_dw5_k f_dw5_l]);

diag_M      = kron([gamma_f beta_f gamma_f], ones(nx,1));
M           = spdiags(diag_M,-1:1,nx,nx);

% % test
% x       = linspace(0,3*pi,nx);
% f       = sin(x) + 0.3*rand(1,nx);%0.1*(-1).^(1:nx);
% figure, plot(x,f,x,M\(D*f'))

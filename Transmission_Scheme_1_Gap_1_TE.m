%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename:    Transmission_Scheme_1_Gap_1_TE.m
% Author:      Shakeel Ahmed
% Date:        July 10, 2025
% MATLAB:      R2020b or later
%
% Purpose:
%   Calculate and plot the electric-field distribution through a 1D Si/SiO2 stack
%   with a central defect layer, showing the localized interface mode.
% COMPUTE_EFIELD_INTERFACE  Compute E-field distribution for a defect interface mode in a 1D Topological photonic crystal.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
c = 3e8;     
pi2 = 2*pi;
eps_in = 1^2;
miu_in = 1; 
n_in = sqrt(eps_in * miu_in);

eps_out = 1^2;
miu_out = 1; 
n_out = sqrt(eps_out * miu_out);

eps_a = 3.48^2;
miu_a = 1; 
na = sqrt(eps_a * miu_a);

eps_b = 1.45^2;
miu_b = 1; 
nb = sqrt(eps_b * miu_b);

eps_c = eps_a;
miu_c = 1; 
nc = sqrt(eps_c * miu_c);

eps_d = eps_b;
miu_d = 1; 
nd = sqrt(eps_d * miu_d);

eps_dd = 1.33^2;
miu_dd = 1;
ndd = sqrt(eps_dd * miu_dd);

ha = 110e-9;  
hb = 270e-9;   
hc = 0.5*ha;   
hd = 0.5*hb;   
hdd = 190e-9;   

f_start = 140e12;
f_end = 180e12;
num_freq = 10000;  
freq_THz = linspace(f_start/1e12, f_end/1e12, num_freq);
w_angular = 2*pi*freq_THz*1e12;  

angle_deg = [0, 15, 30, 45, 60, 75, 89];
angle_rad = angle_deg * pi/180; 
num_angles = length(angle_deg);

T_matrix = zeros(num_freq, num_angles);

for i_angle = 1:num_angles
    a = angle_rad(i_angle);

    for i_freq = 1:num_freq
        w = w_angular(i_freq);

        pa = sqrt(eps_a/miu_a) * sqrt(1 - sin(a)^2/(eps_a*miu_a));
        pb = sqrt(eps_b/miu_b) * sqrt(1 - sin(a)^2/(eps_b*miu_b));
        pc = sqrt(eps_c/miu_c) * sqrt(1 - sin(a)^2/(eps_c*miu_c));
        pd = sqrt(eps_d/miu_d) * sqrt(1 - sin(a)^2/(eps_d*miu_d));
        pdd = sqrt(eps_dd/miu_dd) * sqrt(1 - sin(a)^2/(eps_dd*miu_dd));

        % % % % Wave impedance (TM polarization)
        % pa = sqrt(miu_a/eps_a) * sqrt(1 - sin(a)^2/(eps_a*miu_a));
        % pb = sqrt(miu_b/eps_b) * sqrt(1 - sin(a)^2/(eps_b*miu_b));
        % pc = sqrt(miu_c/eps_c) * sqrt(1 - sin(a)^2/(eps_c*miu_c));
        % pd = sqrt(miu_d/eps_d) * sqrt(1 - sin(a)^2/(eps_d*miu_d));
        % pdd = sqrt(miu_dd/eps_dd) * sqrt(1 - sin(a)^2/(eps_dd*miu_dd));
        ka = (w/c) * sqrt(eps_a*miu_a) * sqrt(1 - sin(a)^2/(eps_a*miu_a));
        kb = (w/c) * sqrt(eps_b*miu_b) * sqrt(1 - sin(a)^2/(eps_b*miu_b));
        kc = (w/c) * sqrt(eps_c*miu_c) * sqrt(1 - sin(a)^2/(eps_c*miu_c));
        kd = (w/c) * sqrt(eps_d*miu_d) * sqrt(1 - sin(a)^2/(eps_d*miu_d));
        kdd = (w/c) * sqrt(eps_dd*miu_dd) * sqrt(1 - sin(a)^2/(eps_dd*miu_dd));

        ma = [cos(ka*ha), -1i/pa*sin(ka*ha); -1i*pa*sin(ka*ha), cos(ka*ha)];
        mb = [cos(kb*hb), -1i/pb*sin(kb*hb); -1i*pb*sin(kb*hb), cos(kb*hb)];
        mc = [cos(kc*hc), -1i/pc*sin(kc*hc); -1i*pc*sin(kc*hc), cos(kc*hc)];
        md = [cos(kd*hd), -1i/pd*sin(kd*hd); -1i*pd*sin(kd*hd), cos(kd*hd)];
        mdd = [cos(kdd*hdd), -1i/pdd*sin(kdd*hdd); -1i*pdd*sin(kdd*hdd), cos(kdd*hdd)];

        M = (md*ma*md)^10 *mdd* (mc*mb*mc)^10;

        t = 2*cos(a) / (M(1,1)*cos(a) + M(1,2) + (M(2,1) + M(2,2)*cos(a)));
        T = abs(t)^2;

        T_matrix(i_freq, i_angle) = T;
    end
end

%%
figure;
set(gca, 'FontSize', 12, 'FontWeight', 'normal');
set(gcf, 'Position', [100, 100, 350, 300]);
hold on;

colors = jet(num_angles);

T_dB = 10*log10(T_matrix); 
all_dB = T_dB(:);           
all_dB = all_dB(isfinite(all_dB));
min_dB = min(all_dB);
max_dB = max(all_dB);
zmin = floor(min_dB/10)*10;
zmax = ceil(max_dB/10)*10;
x = freq_THz;
for i = 1:num_angles
    y = angle_deg(i) * ones(size(freq_THz));
    z = 10*log10(T_matrix(:, i));
    
    plot3(x, y, z, 'LineWidth', 2, 'Color', colors(i, :));
end

xh = xlabel('f (THz)', 'Rotation', -15, 'FontSize', 12, 'FontWeight','bold');
yh = ylabel('\theta (°)', 'Rotation', 50, 'FontSize', 12, 'FontWeight','bold');
zlabel('T', 'FontSize', 12, 'FontWeight','bold');
grid on;

xlim([f_start/1e12, f_end/1e12]);
ylim([0 89]);
zlim([zmin, zmax]);
set(gca, 'ZTick', [zmin, zmax]); 

labels = arrayfun(@(v) sprintf('10^{%d}', v/10), [zmin, zmax], 'UniformOutput', false);
set(gca, 'ZTickLabel', labels);

view(30, 75);  
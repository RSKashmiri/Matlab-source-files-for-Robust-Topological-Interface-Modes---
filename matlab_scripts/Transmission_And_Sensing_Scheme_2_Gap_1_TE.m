%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename:    Transmission_And_Sensing_Scheme_2_Gap_1_TE.m
% Author:      Shakeel Ahmed
% Date:        July 10, 2025
% MATLAB:      R2020b or later
%
% Purpose:
%   Calculate and plot the electric-field distribution through a 1D Si/SiO2 stack
%   with a central defect layer, showing the localized interface mode.
% COMPUTE_EFIELD_INTERFACE  Compute E-field distribution for a defect interface mode in a 1D Topological photonic crystal.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;
c    = 3e8;     
pi2  = 2*pi;
eps_in  = 1^2; miu_in  = 1; n_in  = sqrt(eps_in*miu_in);
eps_out = 1^2; miu_out = 1; n_out = sqrt(eps_out*miu_out);

eps_a = 3.48^2; miu_a = 1;  
eps_b = 1.45^2; miu_b = 1;   
eps_c = eps_a;   miu_c = 1; 
eps_d = eps_b;   miu_d = 1;  

refractive_index = [ ...
    1.3333; 1.3412; 1.3425; 1.3531; 1.3951; ...
    1.4121; 1.4320; 1.4412; 1.4470; 1.4591; 1.4833];
num_samples = numel(refractive_index);

ha  = 110e-9;    
hb  = 270e-9;   
hc  = 0.5*ha;   
hd  = 0.5*hb;   
hdd = 380e-9;   

%angle_deg  = [0, 15, 30, 45, 60, 75, 89];
angle_deg  = 0;
angle_rad  = angle_deg * pi/180;
num_angles = numel(angle_deg);

% freq_ranges = [ ...
%    171, 178.1;    % 0°
%    172.8, 180.5;    % 15°
%    178.5, 187.6;    % 30°
%    187.5, 199.3;    % 45°
%    199, 214.5;    % 60°
%    209.6, 229.1;    % 75°
%    214, 235.5];   % 89°
freq_ranges = [171, 178.1];
num_freq = 100000;  
freq_vectors = cell(num_angles,1);
for ia = 1:num_angles
    f0 = freq_ranges(ia,1);
    f1 = freq_ranges(ia,2);
    freq_vectors{ia} = linspace(f0, f1, num_freq);
end

T_matrix = zeros(num_freq, num_angles, num_samples);
f_r_list  = zeros(num_samples, num_angles);
FWHM_list = zeros(num_samples, num_angles);
Q_list    = zeros(num_samples, num_angles);
S_list    = zeros(num_samples, num_angles);
for sample_idx = 1:num_samples
    eps_dd = refractive_index(sample_idx)^2;
    miu_dd = 1;
    for ia = 1:num_angles
        a      = angle_rad(ia);
        sin2a  = sin(a)^2;
        f_vec  = freq_vectors{ia};
        w_ang  = 2*pi * f_vec * 1e12; 

        for jf = 1:num_freq
            w   = w_ang(jf);
            k0  = w / c;

            pa  = sqrt(eps_a/miu_a)   * sqrt(1 - sin2a/(eps_a/miu_a));
            pb  = sqrt(eps_b/miu_b)   * sqrt(1 - sin2a/(eps_b/miu_b));
            pc  = sqrt(eps_c/miu_c)   * sqrt(1 - sin2a/(eps_c/miu_c));
            pd  = sqrt(eps_d/miu_d)   * sqrt(1 - sin2a/(eps_d/miu_d));
            pdd = sqrt(eps_dd/miu_dd) * sqrt(1 - sin2a/(eps_dd/miu_dd));

            ka  = k0*sqrt(eps_a*miu_a)*sqrt(1 - sin2a/(eps_a*miu_a));
            kb  = k0*sqrt(eps_b*miu_b)*sqrt(1 - sin2a/(eps_b*miu_b));
            kc  = k0*sqrt(eps_c*miu_c)*sqrt(1 - sin2a/(eps_c*miu_c));
            kd  = k0*sqrt(eps_d*miu_d)*sqrt(1 - sin2a/(eps_d*miu_d));
            kdd = k0*sqrt(eps_dd*miu_dd)*sqrt(1 - sin2a/(eps_dd*miu_dd));

            ma  = [cos(ka*ha),   -1i/pa*sin(ka*ha);
                   -1i*pa*sin(ka*ha), cos(ka*ha)];
            mb  = [cos(kb*hb),   -1i/pb*sin(kb*hb);
                   -1i*pb*sin(kb*hb), cos(kb*hb)];
            mc  = [cos(kc*hc),   -1i/pc*sin(kc*hc);
                   -1i*pc*sin(kc*hc), cos(kc*hc)];
            md  = [cos(kd*hd),   -1i/pd*sin(kd*hd);
                   -1i*pd*sin(kd*hd), cos(kd*hd)];
            mdd = [cos(kdd*hdd), -1i/pdd*sin(kdd*hdd);
                   -1i*pdd*sin(kdd*hdd), cos(kdd*hdd)];

            M = (md*ma*md)^9 * md *ma * mdd *mb *mc * (mc*mb*mc)^9;

            t = 2*cos(a) / ( M(1,1)*cos(a) + M(1,2) + M(2,1) + M(2,2)*cos(a) );
            T_matrix(jf,ia,sample_idx) = abs(t)^2;
        end
    end
    fprintf('Completed sample %d/%d\n', sample_idx, num_samples);
end

sample_colors = cool(num_samples); 
figure('Position',[300,300,370,90],'Color','w');
ax = axes('Position', [0.15 0.2 0.75 0.7]); 
hold(ax,'on'); 
box(ax,'on'); 
grid(ax,'off');

for s = 1:num_samples
    plot(ax, freq_vectors{1}, squeeze(T_matrix(:,1,s)), ...
         'Color', sample_colors(s,:), 'LineWidth', 1.5);
end

xlim(ax, freq_ranges(1,:));
ylim(ax, [0 1]);
xlabel(ax, 'f(THz)');

for ia = 1:num_angles
    fv = freq_vectors{ia};
    for s = 1:num_samples
        sp = squeeze(T_matrix(:,ia,s));
        [pv, loc]       = max(sp);
        f_r_list(s,ia)  = fv(loc);
        thr = pv/2;
        idx = find(sp >= thr);
        if isempty(idx)
            FWHM_list(s,ia) = NaN;
            Q_list   (s,ia) = NaN;
        else
            L = idx(1);
            R = idx(end);
            FWHM_list(s,ia) = fv(R) - fv(L);
            Q_list(s,ia)    = f_r_list(s,ia) / FWHM_list(s,ia);
        end
    end

    S_list(1,ia) = 0;
    for s = 2:num_samples
        df = f_r_list(s,ia) - f_r_list(1,ia);
        dn = refractive_index(s)   - refractive_index(1);
        S_list(s,ia) = abs(df/dn);
    end
end

for ia = 1:num_angles
    T = table( ...
      refractive_index, ...
      f_r_list(:,ia), ...
      FWHM_list(:,ia), ...
      Q_list(:,ia), ...
      S_list(:,ia), ...
      'VariableNames',{ ...
        'RefractiveIndex', ...
        'ResonantFreq_THz', ...
        'FWHM_THz', ...
        'Q_Factor', ...
        'Sensitivity_THz_per_RIU'} );
    fprintf('\n=== Angle = %d° ===\n', angle_deg(ia));
    disp(T);
end

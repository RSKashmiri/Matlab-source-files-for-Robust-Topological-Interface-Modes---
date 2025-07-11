

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename:    Code1_Dispersion().m
% Author:      Shakeel Ahhmed
% Date:        July 10, 2025
% MATLAB:      R2020b or later
%
% Purpose:
%   Compute and plot the photonic-crystal band structure and identify band gaps
%   for a silicon (Si) / silicon dioxide (SiO2) multilayer (1D stack).
%
% Usage:
%   >> [result_table, gap_ranges] = compute_dispersion();
%
% Inputs:
%   None  % all parameters are defined within the function
%
% Outputs:
%   result_table - table listing even/odd gaps and odd gap centers (THz)
%   gap_ranges   - gap start/end frequencies (THz)
%
% Dependencies:
%   None (core MATLAB functions only)
%
% Notes:
%   • Layer thicknesses: h_a (Si) and h_b (SiO2) defined below.
%   • Frequency vector wq is in THz; converted inside to angular units.
%   • To modify layer properties, edit parameters in the 'Material parameters' section.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [result_table, gap_ranges] = Code1_Dispersion()

c = 3e8;               
points = 8000;         



e_a = 3.48^2;  miu_a = 1;  h_a = 110e-9;

e_b = 1.45^2;  miu_b = 1;  h_b = 270e-9;

a = h_a + h_b;        


n_a = sqrt(e_a);
n_b = sqrt(e_b);
p_a = sqrt(e_a/miu_a);
p_b = sqrt(e_b/miu_b);


wq = linspace(0, 800, points);

temp_vals = zeros(1, points);


for i = 1:points
    omega = wq(i)*1e12*2*pi;    
    ka = omega/c * n_a;
    kb = omega/c * n_b;
    m1 = [cos(ka*h_a), -1i*(1/p_a)*sin(ka*h_a); -1i*p_a*sin(ka*h_a), cos(ka*h_a)];
    m2 = [cos(kb*h_b), -1i*(1/p_b)*sin(kb*h_b); -1i*p_b*sin(kb*h_b), cos(kb*h_b)];
    m3 = [cos(ka*h_a/2), -1i*(1/p_a)*sin(ka*h_a/2); -1i*p_a*sin(ka*h_a/2), cos(ka*h_a/2)];
    m4 = [cos(kb*h_b/2), -1i*(1/p_b)*sin(kb*h_b/2); -1i*p_b*sin(kb*h_b/2), cos(kb*h_b/2)];
    M1 = m4 * m1 * m4; % P  Total transfer matrix
    %M1 = m3 * m2 * m3;  % Q Total transfer matrix             

    temp = cos(ka*h_a).*cos(kb*h_b) - 0.5*(p_a/p_b + p_b/p_a).*sin(ka*h_a).*sin(kb*h_b);
    temp_vals(i) = real(temp);
end


ingap = abs(temp_vals) > 1;
trans = diff([0, ingap, 0]);
starts = trans == 1;
ends   = find(trans == -1) - 1;
gap_ranges = [wq(starts)', wq(ends)'];


gap_ranges = sortrows(gap_ranges, 1);
even = gap_ranges(mod((1:size(gap_ranges,1))',2)==0, :);
odd  = gap_ranges(mod((1:size(gap_ranges,1))',2)==1, :);
odd_centers = mean(odd,2);


maxr = max(size(even,1), size(odd,1));
if size(even,1) < maxr; even(end+1:maxr,:) = nan; 
end
if size(odd,1)  < maxr; odd(end+1:maxr,:)   = nan; odd_centers(end+1:maxr) = nan; 
end
result_table = table(even, odd, odd_centers, ...
    'VariableNames', {'Even_Gaps_THz','Odd_Gaps_THz','Odd_Center_THz'});


disp('Band Gap Table:');
disp(result_table);


figure('Position',[100,100,400,600]); 
hold on;
K = acos(temp_vals)/pi;
plot(K, wq, '--', 'LineWidth',1.2);
plot(-K, wq, '--', 'LineWidth',1.2);
for idx = 1:size(gap_ranges,1)
    y = gap_ranges(idx,:);
    patch([-1,1,1,-1], [y(1),y(1),y(2),y(2)], 'k', 'FaceAlpha',0.1, 'EdgeColor','none');
end
xlabel('K (\pi/a)'); ylabel('Frequency (THz)');
axis([-1 1 0 800]); box on; grid on;

end

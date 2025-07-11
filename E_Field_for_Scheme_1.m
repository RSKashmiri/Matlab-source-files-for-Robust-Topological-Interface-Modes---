
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename:    E_Field_for_Scheme_1.m
% Author:      Shakeel Ahmed
% Date:        July 10, 2025
% MATLAB:      R2020b or later
%
% Purpose:
%   Calculate and plot the electric-field distribution through a 1D Si/SiO2 stack
%   with a central defect layer, showing the localized interface mode.
% COMPUTE_EFIELD_INTERFACE  Compute E-field distribution for a defect interface mode in a 1D Topological photonic crystal.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
c=3e8; 
e0=8.854e-12; 
miu0=4*pi*1e-7; 
ea=3.48^2;
miua=1;  
eb=1.45^2;
miub=1; 
ec=ea;
miuc=1;  
ed=eb;
miud=1;   
em=1.33^2;
mium=1;
e_in=1^2;
miu_in=1; 
n_in=e_in^(1/2)*miu_in^(1/2);
e_out=1^2;
miu_out=1; 
n_out=e_out^(1/2)*miu_out^(1/2);
   
na=ea^(1/2)*miua^(1/2);
nb=eb^(1/2)*miub^(1/2);
nc=ec^(1/2)*miuc^(1/2);
nd=ed^(1/2)*miud^(1/2);
nm=em^(1/2)*mium^(1/2);

H0=1e-8; 
T0=H0/c;
f0=1/T0;
pi2=2*pi;
w0=pi2*f0;
ha=110e-9;
hb=270e-9;
hc=0.5*ha;
hd=0.5*hb;
hm=190e-9;

point=1000; 
a1=0*pi/180; 
f=147.624e12; % 147.959 THz
%f=600.357e12;
w=pi2*f;


ka=(w/c)*(ea^(1/2))*(miua^(1/2))*((1-(((sin(a1))^2)/(ea*miua)))^(1/2));
kb=(w/c)*(eb^(1/2))*(miub^(1/2))*((1-(((sin(a1))^2)/(eb*miub)))^(1/2));
kc=(w/c)*(ec^(1/2))*(miuc^(1/2))*((1-(((sin(a1))^2)/(ec*miuc)))^(1/2));
kd=(w/c)*(ed^(1/2))*(miud^(1/2))*((1-(((sin(a1))^2)/(ed*miud)))^(1/2));
km=(w/c)*(em^(1/2))*(mium^(1/2))*((1-(((sin(a1))^2)/(em*mium)))^(1/2));

p_in=(miu0*miu_in)^(1/2)/(e0*e_in)^(1/2);
p_out=(miu0*miu_out)^(1/2)/(e0*e_out)^(1/2);
pa=ea^(1/2)/miua^(1/2)*((1-((sin(a1))^2/(ea*miua)))^(1/2));
pb=eb^(1/2)/miub^(1/2)*((1-((sin(a1))^2/(eb*miub)))^(1/2));
pc=ec^(1/2)/miuc^(1/2)*((1-((sin(a1))^2/(ec*miuc)))^(1/2));
pd=ed^(1/2)/miud^(1/2)*((1-((sin(a1))^2/(ed*miud)))^(1/2));
pm=em^(1/2)/mium^(1/2)*((1-((sin(a1))^2/(em*mium)))^(1/2));

% TM wave impedance
% pa=miua^(1/2)/ea^(1/2)*((1-((sin(a1))^2/(ea*miua)))^(1/2));
% pb=miub^(1/2)/eb^(1/2)*((1-((sin(a1))^2/(eb*miub)))^(1/2));
% ... (similar for other layers)

ma=[cos(ka*ha),-1i*(1/pa)*sin(ka*ha);-1i*pa*sin(ka*ha),cos(ka*ha)];
mb=[cos(kb*hb),-1i*(1/pb)*sin(kb*hb);-1i*pb*sin(kb*hb),cos(kb*hb)];
mc=[cos(kc*hc),-1i*(1/pc)*sin(kc*hc);-1i*pc*sin(kc*hc),cos(kc*hc)];
md=[cos(kd*hd),-1i*(1/pd)*sin(kd*hd);-1i*pd*sin(kd*hd),cos(kd*hd)];
mm=[cos(km*hm),-1i*(1/pm)*sin(km*hm);-1i*pm*sin(km*hm),cos(km*hm)];


m=10; 
M = (md*ma*md)^m*mm*(mc*mb*mc)^m;

t=2*cos(a1)/((M(1)+M(3)*cos(a1))*cos(a1)+(M(2)+M(4)*cos(a1)));
r=((M(1)+M(3)*cos(a1))*cos(a1)-(M(2)+M(4)*cos(a1)))/((M(1)+M(3)*cos(a1))*cos(a1)+(M(2)+M(4)*cos(a1)));
T=(abs(t))^2; 
R=(abs(r))^2; 
veri=abs(t)^2+abs(r)^2;

H1=hd+ha+hd; 
H2=hm;       
H3=hc+hb+hc; 
L=H1*m+H2*1+H3*m; 
L1=H1*m;     
L2=H1*m+H2*1;
U0V0=[(1+r);(1-r)];
M1=md*ma*md;        
M2=(md*ma*md)^m;    
M3=(md*ma*md)^m*mm; 

is=1;
l=L/point; 

for d=0:l:L 
    if d<=L1 
        N=fix(d/H1); 
        d1=d-(N*H1);  
        if d1<=hd
            bd0=kd*d1;
            ba0=0;
            bd01=0;
        elseif d1<=hd+ha
            bd0=kd*hd;
            ba0=ka*(d1-hd);
            bd01=0;
        else
            bd0=kd*hd;
            ba0=ka*ha;
            bd01=kd*(d1-hd-ha);
        end
        md0=[cos(bd0),-1i*(1/pd)*sin(bd0);-1i*pd*sin(bd0),cos(bd0)];
        ma0=[cos(ba0),-1i*(1/pa)*sin(ba0);-1i*pa*sin(ba0),cos(ba0)];
        md01=[cos(bd01),-1i*(1/pd)*sin(bd01);-1i*pd*sin(bd01),cos(bd01)];
        m0=md0*ma0*md01; 
        M0=M1^N;          

    elseif d>L1 && d<=L2 
        N=fix((d-L1)/H2);
        d2=d-L1-(N*H2);
        bm0=km*d2;
        mm0=[cos(bm0),-1i*(1/pm)*sin(bm0);-1i*pm*sin(bm0),cos(bm0)];
        m0=mm0;   
        M0=M2;    

    else 
        N=fix((d-L2)/H3);
        d3=d-L2-(N*H3);
        if d3<=hc
            bc0=kc*d3;
            bb0=0;
            bc01=0;
        elseif d3<=hc+hb
            bc0=kc*hc;
            bb0=kb*(d3-hc);
            bc01=0;
        else
            bc0=kc*hc;
            bb0=kb*hb;
            bc01=kc*(d3-hc-hb);
        end
        mc0=[cos(bc0),-1i*(1/pc)*sin(bc0);-1i*pc*sin(bc0),cos(bc0)];
        mb0=[cos(bb0),-1i*(1/pb)*sin(bb0);-1i*pb*sin(bb0),cos(bb0)];
        mc01=[cos(bc01),-1i*(1/pc)*sin(bc01);-1i*pc*sin(bc01),cos(bc01)];
        m0=mc0*mb0*mc01; 
        M0=M3*(mc*mb*mc)^N; 
    end
    M100=M0*m0;
    UV=pinv(M100)*U0V0; 
    at(is,1)=d;
    at(is,2)=abs(UV(1)); 
    is=is+1;
end

figure('Position',[100 100 1000 200]);
hold on;
xx = at(:,1);
yy = at(:,2);
maxE = max(yy)*1.1;

for mm_i = 0:(m-1)
    x0 = mm_i*H1;
    patch([x0,      x0+hd,    x0+hd,    x0],      [0,0,maxE,maxE], 'c','FaceAlpha',0.1,'EdgeColor','none');
    patch([x0+hd,   x0+hd+ha, x0+hd+ha, x0+hd],   [0,0,maxE,maxE], 'y','FaceAlpha',0.1,'EdgeColor','none');
    patch([x0+hd+ha,x0+H1,    x0+H1,    x0+hd+ha], [0,0,maxE,maxE], 'c','FaceAlpha',0.1,'EdgeColor','none');
end

patch([L1,   L2,   L2,   L1], [0,0,maxE,maxE], 'm','FaceAlpha',0.2,'EdgeColor','none');
for mm_i = 0:(m-1)
    x0 = L2 + mm_i*H3;
    patch([x0,       x0+hc,    x0+hc,    x0],       [0,0,maxE,maxE], 'g','FaceAlpha',0.1,'EdgeColor','none');
    patch([x0+hc,    x0+hc+hb, x0+hc+hb, x0+hc],    [0,0,maxE,maxE], [0.5 0.5 0.5],'FaceAlpha',0.1,'EdgeColor','none');
    patch([x0+hc+hb, x0+H3,    x0+H3,    x0+hc+hb], [0,0,maxE,maxE], 'g','FaceAlpha',0.1,'EdgeColor','none');
end
plot(xx, yy, 'r', 'LineWidth',1.5);
xline(L1, '--m', 'LineWidth',1);
xline(L2, '--g', 'LineWidth',1);
xlabel('Position (m)');
ylabel('|E| (a.u.)');
xlim([0, L]);
ylim([0, maxE]);
box on;
hold off;



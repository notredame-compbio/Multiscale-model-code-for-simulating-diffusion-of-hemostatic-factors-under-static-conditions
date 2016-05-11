clear all;
close all;% solution of convective-diffucion equation for a spherical thrombus
% number of nodes in the cap 25, in the outer region 7*25=170; total 200

% First, impose flow field
mu = 4*1e-3; %liquid viscosity
rho=1000; % density
U = 0.01; % velocity m/s
b = 25e-6; 
a = 20e-6; 
%k = 1e-12:1e-12:1e-9;
fi = 3e-3:1e-5:16e-3;
fi0 = 0.0028;
kfi = 3e-16.*(fi-fi0).^(-1.8);

k=700e-12;
kpermeab = k;

alfa= a./sqrt(k);
beta = b./sqrt(k);
delta = beta-alfa;


J = -6*alfa+(3*alfa + 3*beta +alfa.^3 +2*beta.^3).*cosh(delta) + 3*(alfa.^2 - 1).*sinh(delta);
H = (1./J).*(3.*(alfa.^3 + 2.*beta.^3).*cosh(alfa) + 9.*alfa.*(cosh(alfa) - alfa.*sinh(alfa) - cosh(beta) + beta.*sinh(beta) ));
G = (3.*alfa - (alfa.*cosh(beta) - sinh(alfa)).*H )./(alfa.*sinh(beta) -cosh(alfa));
F = (G.*cosh(alfa) + H.*sinh(alfa))./(3.*alfa);

B0 = 3*(alfa.^4 +2*alfa.*beta.^3 + 3*alfa.^2).*cosh(alfa) + 9*alfa.^2.*(cosh(beta) - beta.*sinh(beta) -alfa.*sinh(alfa))+...
    3*cosh(delta).*( (alfa.^3 + 2*beta.^3 +3.*alfa).*(alfa.*beta.*sinh(beta) - alfa.*cosh(beta) - beta.*cosh(alfa)) + 3.*alfa.^2.*beta.*sinh(alfa) )+...
    3.*sinh(delta).*((alfa.^3 + 2.*beta.^3 + 3.*alfa).*cosh(alfa) + 3.*alfa.^2.*(alfa.*beta.*sinh(beta) - alfa.*cosh(beta) - sinh(alfa)) );
B = B0./(2.*(alfa.*sinh(beta) - cosh(alfa)).*J);
E1 = 2.*B + 2.*(beta.^3).*F;
D = 0;
C1 = -1;
A = beta.^3 - beta.^2.*B + E1 + beta.^3.*F + (cosh(beta) - beta.*sinh(beta)).*G + (sinh(beta) - beta.*cosh(beta)).*H;


radiusN_vel=50;%25
i = 1:1:radiusN_vel;
delta_r_INSIDE = (b-a)/(radiusN_vel-1);
r_INSIDE = a+(i-1)*delta_r_INSIDE;
thetaN = 200; %50
theta_INSIDE = 0:pi/(thetaN-1):pi;

f_rad = 7;
s = 1:1:f_rad*radiusN_vel;%def 7
delta_r_OUTSIDE = (b-a)/(radiusN_vel-1);
r_OUTSIDE = b+(s)*delta_r_OUTSIDE;
theta_OUTSIDE = theta_INSIDE;


for i = 1:1:radiusN_vel
    for j = 1:1:thetaN
  %INSIDE a<=r<=b
        ksi_INSIDE(i) = r_INSIDE(i)./sqrt(k);
        F_INSIDE(i) = (E1./ksi_INSIDE(i) + F.*ksi_INSIDE(i).^2 + G.*(cosh(ksi_INSIDE(i))./ksi_INSIDE(i) - sinh(ksi_INSIDE(i))) + H.*(sinh(ksi_INSIDE(i))./ksi_INSIDE(i) - cosh(ksi_INSIDE(i))));
psiINSIDE(i,j) = -0.5.*k.*U.*F_INSIDE(i).*(sin(theta_INSIDE(j))).^2;

F_INSIDE_theta(i) = (-E1.*sqrt(k)./r_INSIDE(i).^2 + 2.*F.*r_INSIDE(i)./k + G.*(sinh(r_INSIDE(i)./sqrt(k))./r_INSIDE(i) - sqrt(k).*cosh(r_INSIDE(i)./sqrt(k))./r_INSIDE(i).^2 -...
    1./sqrt(k).*cosh(r_INSIDE(i)./sqrt(k))) + H.*(cosh(r_INSIDE(i)./sqrt(k))./r_INSIDE(i) - sqrt(k).*sinh(r_INSIDE(i)./sqrt(k))./r_INSIDE(i).^2 - sinh(r_INSIDE(i)./sqrt(k))./sqrt(k))  );

        V_theta_INSIDE(i,j) = -0.5*k.*U./(r_INSIDE(i)).*F_INSIDE_theta(i).*sin(theta_INSIDE(j));
        V_R_INSIDE(i,j) = k.*U.*F_INSIDE(i).*cos(theta_INSIDE(j))./r_INSIDE(i).^2;
    end;
end; 
for i = 1:1:f_rad*radiusN_vel
    for j = 1:1:thetaN        
        %OUTSIDE r>b       
        ksi_OUTSIDE(i) = r_OUTSIDE(i)./sqrt(k);
        F_OUTSIDE(i) = (A./ksi_OUTSIDE(i) + B.*ksi_OUTSIDE(i) +C1.*ksi_OUTSIDE(i).^2 +D.*ksi_OUTSIDE(i).^4);
        F_OUTSIDE_theta(i) = -A./r_OUTSIDE(i).^2*sqrt(k) + B./sqrt(k) + 2*C1.*r_OUTSIDE(i)./k + 4*D.*r_OUTSIDE(i).^3/k.^2;
        psi_OUTSIDE(i,j) = -0.5.*k.*U.*F_OUTSIDE(i).*(sin(theta_OUTSIDE(j))).^2;
        
        V_R_OUTSIDE(i,j) = k.*U.*F_OUTSIDE(i).*cos(theta_OUTSIDE(j))./r_OUTSIDE(i).^2;

V_theta_OUTSIDE(i,j) = -0.5*k.*U./(r_OUTSIDE(i)).*F_OUTSIDE_theta(i).*sin(theta_OUTSIDE(j));
        
    end;
end;

clear A B C1 D F G;





%Now solve conv-diffusion, with protein generated at r = a_real;



a_real = 2e-5;%core radius in m
b_real = b/a*a_real; % thrombus radius in m
U0_real = 0.01; %velosity at inf, m/s
D = 1e-10; % diffusion coefficient, m^2/s
a_n = 1;
b_n = b/a*a_n;
U0 = 1;
Pe = U0_real*a_real/D;

thetamax = pi;
radiusmax = (f_rad+1)*(b_n-a_n);
timeN = 1500;
%thetaN = 50;

radiusN = (f_rad+1)*radiusN_vel;
delta_theta = thetamax/(thetaN-1);
j = 1:1:thetaN;
theta = (j-1)*delta_theta;

i = 1:1:radiusN;
delta_radius = radiusmax/(radiusN-1);
radius = a_n+(i-1)*delta_radius;
delta_time = 0.5*min(delta_radius,delta_theta);    %delta_radius
timemax = delta_time*timeN;
k = 1:1:timeN;
time = k*delta_time;
cNew = zeros(radiusN, thetaN);
cOld = zeros(radiusN, thetaN);
A = zeros(radiusN, thetaN);
B = zeros(radiusN, thetaN);
D = zeros(radiusN, thetaN);
F = zeros(radiusN, thetaN);
G = zeros(radiusN, thetaN);

u_theta = zeros(radiusN, thetaN);
u_radial = zeros(radiusN, thetaN);

  
    for i = 1:1:radiusN_vel
        for j=1:1:thetaN
            u_radial(i,j) = V_R_INSIDE(i,j)./U0_real;
            u_theta(i,j) = V_theta_INSIDE(i,j)./U0_real;
        end;
    end;
    for i = radiusN_vel+1:1:(f_rad+1)*radiusN_vel
        for j=1:1:thetaN
            u_radial(i,j) = V_R_OUTSIDE(i-radiusN_vel,j)./U0_real;
            u_theta(i,j) = V_theta_OUTSIDE(i-radiusN_vel,j)./U0_real;
        end;
    end;

% Stability condition
% for diffusion -  diffusion time step < 0.5*Pe*delta_radius^2
Kurad = 1.*delta_time/delta_radius;
Kutheta = 1.*delta_time/delta_theta;
%initial conditions:

%protein initialy
% distributed in a spherical layer near the core
for zz=1:20 %10
    cOld(zz,:) = 1;%exp(-(zz-1).^2./15);%radiusN_vel, /20
    cNew(zz,:) = cOld(zz,:);
end;

cOld(1,:) = 1; % BC
cNew(1,:) = 1;
cOld(radiusN,:) = 0; %BC on the outer boundary
cNew(radiusN,:) = 0;


cOld(:,1) = cOld(:,2);
cNew(:,1) = cNew(:,2);
cOld(:,thetaN) = cOld(:,thetaN-1);
cNew(:,thetaN) = cNew(:,thetaN-1);
for i=2:1:radiusN-1
 radius_L(i) = 0.5*(radius(i)-radius(i-1));
 radius_R(i) = 0.5*(radius(i+1)-radius(i));
end;

for j=2:1:thetaN-1
 theta_L(j) = 0.5*(theta(j)-theta(j-1));
 theta_R(j) = 0.5*(theta(j+1)-theta(j));
end;



  
    

for i=1:1:radiusN
    C(i) = (1./Pe)*(1./radius(i).^2);
end;

h = waitbar(0, 'Progress...');
for k=1:1:timeN
    for i=2:1:radiusN-1
        for j=2:1:thetaN-1
            
           
           

%------------------donor cell end ----
%general:
%radial direction:
u_radial_R = 0.5*(u_radial(i+1,j)+u_radial(i,j));
u_radial_L = 0.5*(u_radial(i,j)+u_radial(i-1,j));
if u_radial_L>=0
    TETA_L = 1;
    if i>2
     H_L = (cOld(i-1,j) - cOld(i-2,j))./(cOld(i,j)-cOld(i-1,j));
    else
        H_L=0;
    end; 
end;

if u_radial_L<0
    TETA_L = -1;
     H_L = (cOld(i+1,j) - cOld(i,j))./(cOld(i,j)-cOld(i-1,j)); 
end;

if u_radial_R>=0
    TETA_R = 1;
     H_R = (cOld(i,j) - cOld(i-1,j))./(cOld(i+1,j)-cOld(i,j));     
end;

if u_radial_R<0
    TETA_R = -1;
    if i<radiusN-1
     H_R = (cOld(i+2,j) - cOld(i+1,j))./(cOld(i+1,j)-cOld(i,j));
    else
        H_R=0;
    end;
end;

fi_rad_L = minmod(1,H_L);
fi_rad_R = minmod(1,H_R);
AL = 0.5*u_radial_L*( (1+TETA_L)*cOld(i-1,j) + (1-TETA_L)*cOld(i,j) ) +...
    0.5*abs(u_radial_L)*(1 - abs(u_radial_L*Kurad) )*fi_rad_L*(cOld(i,j) - cOld(i-1,j));
AR = 0.5*u_radial_R*( (1+TETA_R)*cOld(i,j) + (1-TETA_R)*cOld(i+1,j) ) +...
    0.5*abs(u_radial_R)*(1 - abs(u_radial_R*Kurad) )*fi_rad_R*(cOld(i+1,j) - cOld(i,j));
A(i,j) = AR-AL;


     


%general
%theta:
u_theta_R = 0.5*(u_theta(i,j+1)+u_theta(i,j));
u_theta_L = 0.5*(u_theta(i,j)+u_theta(i,j-1));
if u_theta_L>=0
    TETA_L = 1;
    if j>2
     H_L = (cOld(i,j-1) - cOld(i,j-2))./(cOld(i,j)-cOld(i,j-1));
    else
        H_L=0;
    end; 
end;

if u_theta_L<0
    TETA_L = -1;
     H_L = (cOld(i,j+1) - cOld(i,j))./(cOld(i,j)-cOld(i,j-1)); 
end;

if u_theta_R>=0
    TETA_R = 1;
     H_R = (cOld(i,j) - cOld(i,j-1))./(cOld(i,j+1)-cOld(i,j));     
end;

if u_theta_R<0
    TETA_R = -1;
    if j<thetaN-1
     H_R = (cOld(i,j+2) - cOld(i,j+1))./(cOld(i,j+1)-cOld(i,j));
    else
        H_R=0;
    end;
end;

fi_theta_L = minmod(1,H_L);
fi_theta_R = minmod(1,H_R);
BL = 0.5*u_theta_L*( (1+TETA_L)*cOld(i,j-1) + (1-TETA_L)*cOld(i,j) ) +...
    0.5*abs(u_theta_L)*(1 - abs(u_theta_L*Kutheta) )*fi_theta_L*(cOld(i,j) - cOld(i,j-1));
BR = 0.5*u_theta_R*( (1+TETA_R)*cOld(i,j) + (1-TETA_R)*cOld(i,j+1) ) +...
    0.5*abs(u_theta_R)*(1 - abs(u_theta_R*Kutheta) )*fi_theta_R*(cOld(i,j+1) - cOld(i,j));
B(i,j) = (BR-BL)/radius(i);

            G(i,j) = cOld(i,j).*(2*u_radial(i,j)/radius(i) + u_theta(i,j)/radius(i)/tan(theta(j)) );
      
            D(i,j) = (Kurad/delta_radius)*(radius_R(i).^2.*(cOld(i+1,j) - cOld(i,j)) - radius_L(i).^2.*(cOld(i,j) - cOld(i-1,j)) );
            F(i,j) = (Kutheta/delta_theta)*(1./sin(theta(j))).*(sin(theta_R(j)).*(cOld(i,j+1) - cOld(i,j)) - sin(theta_L(j)).*(cOld(i,j) - cOld(i,j-1)) );
            cNew(i,j) = cOld(i,j) - Kurad.*A(i,j) - Kutheta.*B(i,j) - G(i,j)*delta_time + C(i)*(D(i,j) + F(i,j));
        end;
    end;
%----general end


% symmetry on axis    
cNew(:,1) = cNew(:,2);
cNew(:,thetaN) = cNew(:,thetaN-1);

%Bc near the core



%BC on the outer boundary:
cNew(radiusN,:) = cNew(radiusN-1,:);

cOld = cNew;
waitbar(k/timeN);

%if rem(k,1000)==0 
%    filename = [ 'myDataFileCold700_a20_b25_coo' num2str(k) '.mat' ];
%        save(filename,'i','j','kpermeab','cOld');
end;
        
%end;

%save('myDataFileCold3_a20_b25_1','i','j','kpermeab','cOld');
%from polar radius,theta to x and y
for i=1:1:radiusN
        for j=1:1:thetaN
            x(i,j)=radius(i)*cos(theta(j));
            y(i,j) = radius(i)*sin(theta(j));     
        end;
end;
    
    for i=1:1:radiusN
        for j=1:1:thetaN
            ux(i,j)=u_radial(i,j)*cos(theta(j)) - u_theta(i,j)*sin(theta(j));
            uy(i,j) = u_radial(i,j)*sin(theta(j)) + u_theta(i,j)*cos(theta(j));     
        end;
    end;
    
figure;plot(cOld(:,:)); 
figure;
%subplot(3,2,6);
contour(x,y,cOld,200);
colorbar;
h_bar = findobj(gcf,'Tag','Colorbar');
hpos=get(h_bar,'Position');
hpos(3)=hpos(3)/2; % Halve the thickness
set(h_bar,'Position',hpos,'FontSize',12);

hold on;
%subplot(3,2,6);
quiver(x,y,ux,uy);
hold on;
ta = 0:0.01:pi;
xa = a_n*cos(ta); ya = a_n*sin(ta);%surface of the core
xb = b_n*cos(ta); yb = b_n*sin(ta);%surface of the clot
plot(xa,ya,'k-',xb,yb,'k-','Linewidth',2);
xlabel('r/a','FontSize',12);
ylabel('r/a','FontSize',12);

% % we set the units of the measures used through the file
% %
% % [ inches | centimeters | normalized | points | {pixels} | characters ]
% set(gcf, 'Units', 'centimeters');
% % we set the position and dimension of the figure ON THE SCREEN
% %
% % NOTE: measurement units refer to the previous settings!
% afFigurePosition = [1 1 20 5.5]; % [pos_x pos_y width_x width_y]
% set(gcf, 'Position', afFigurePosition); % [left bottom width height]
% % we link the dimension of the figure ON THE PAPER in such a way that
% % it is equal to the dimension on the screen
% %
% % ATTENTION: if PaperPositionMode is not ’auto’ the saved file
% % could have different dimensions from the one shown on the screen!
% set(gcf, 'PaperPositionMode', 'auto');
% % in order to make matlab to do not "cut" latex-interpreted axes labels
% set(gca, 'Units','normalized','Position',[0.15 0.2 0.75 0.7]);

figure;
hold on;
h1=subplot(3,2,6);contour(x,y,cOld,900);
h_bar_d= colorbar;
hold on;
plot(xa,ya,'k-',xb,yb,'k-','Linewidth',2);
xlabel('r/a','FontSize',12);
ylabel('r/a', 'FontSize',12);
axis([-2 2 0. 2]);
text(1.7, 1.7, '\fontsize{12}\bf F');
text(-0.7, 0.5, '\fontsize{12}\bf TC');
arrow([1.8 1.2],[1.1 1.2],14,'BaseAngle', 60, 'Width',2);
set([h1,h2,h3],'clim',[0 1]);

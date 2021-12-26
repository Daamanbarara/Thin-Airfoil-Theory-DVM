%% Characteristics of Thin Airfoil based on Discrete Vortex Method 
% using Lumped-Vortex element method derived from Katz & Plotkin (Ch-11)
% Author: Daamanjyot Barara
clc
clear all
%% Input from the user
airfoil = input('4-digit NACA: ','s');
M = input('Number of panels: ');
alpha = input('Angle of attack: '); % Input is in degrees, must be convereted to radians in the program
Geo = input('Choice of Discretization- a) Uniform b) Full Cosine: ','s');

tic

% Properties of the NACA 4 digit airfoil
c = str2double(airfoil(1))* 0.01; % Max camber
p = str2double(airfoil(2))* 0.1; % Max camber position
t = str2double(airfoil(3:4)); % Thickness

%% Geometry Discretization
N = M + 1;        % Number of data points along the camber line
x_vec=zeros(N,1); % Vector for x-coordinate of the data points
z_vec=zeros(N,1); % Vector for z-coordinate of the data points

for i=1:N
    if (Geo =='a')
        x = (i-1)/M;
        x_vec(i,1) = x;
    elseif (Geo =='b')
        x = 0.5.*(1-cos(((i-1)/(N-1))*pi));
        x_vec(i,1) = x;
    else
        print('Invalid choice for discretization');
    end
    % Distribution of mean camber line using mean camber line equations as
    % used in standard NACA 4 digit airfoil calculation
    if (0 <= x && x < p)
        z = c/(p^2)*((2*p*x)-(x^2));
        z_vec(i,1)=z;
    else
        z = c/((1-p)^2)*(1-(2*p)+(2*p*x)-(x^2));
        z_vec(i,1)=z;
    end
end

%% Computation of geometrical data for each panel

% Computation of Deltas by taking difference between subsequent coordinate
% points along the x-coordinate and z-coordinate.
delta_x=diff(x_vec); % Difference between x-coordinates of two consecutive discretized points along mean camberline
delta_z=diff(z_vec); % Difference between x-coordinates of two consecutive discretized points along mean camberline

% Computation of the panel length 
panel_length=sqrt((delta_x.^2)+(delta_z.^2));
    
% Computation of the normal vector
norm_x_vec=delta_x./(panel_length);
norm_z_vec=-delta_z./(panel_length);
norm_vec=[norm_z_vec,norm_x_vec];

for i=1:M
    % Computation of the coordinates of the vortex points located at quarter chord point of each panel
    x_vort(i,1)=x_vec(i,1)+ 0.25*delta_x(i,1);  % x-coordinate of vortex point
    z_vort(i,1)=z_vec(i,1) + 0.25*delta_z(i,1); % z-coordinate of vortex point
    % Computation of the coordinates of the collocation points located at thrid quarter chord 
    % point of each panel
    x_collp(i,1)=x_vec(i,1) + 0.75*delta_x(i,1); % x-coordinate of collocation point
    z_collp(i,1)=z_vec(i,1) + 0.75*delta_z(i,1); % z-coordinate of collocation point
end

for i=1:M 
    for j=1:M
        % Computation of induced velocities at each collocation point due to each of the vortex points
        r(i,j)=(((x_collp(i)-x_vort(j)).^2) + ((z_collp(i)-z_vort(j)).^2));
        u(i,j)= ((z_collp(i) - z_vort(j)))/(2*pi*r(i,j));
        w(i,j)= -((x_collp(i,1) - x_vort(j)))/(2*pi*r(i,j));
        % Computation of Influence Coefficient
        A(i,j)=  u(i,j)*norm_vec(i,1) +  w(i,j)*norm_vec(i,2);
    end
    % Establish RHS
    RHS(i,1)= -((cos(deg2rad(alpha))*norm_vec(i,1)) + (sin(deg2rad(alpha))*norm_vec(i,2))); 
end

% Computation of Circulation on solving the system of equations
Gamma=linsolve(A,RHS); 

%% Calculation of Loads
Cl= sum(Gamma)*2; % Computation of Lift Coefficient 
Cp=(2*(Gamma))./(panel_length); % Computation of the Pressure Distribution at each panel

toc

%% Plot
plot(x_vec(1:M,1),Cp,'*-');
xlabel('$$\frac{x}{c}$$','Interpreter','latex')
ylabel('$$\Delta C_p$$','Interpreter','latex')
xlim([0 1])

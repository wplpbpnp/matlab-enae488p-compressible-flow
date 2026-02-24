% Benjamin Stutzke
% ENAE 488P
% Homework 3

clear; close all;
%% Problem 1
theta_0 = pi/2;
theta = [pi/2 0];
M1 = 5:0.1:50;
C_d = M1;

% Calculating C_d for gamma 1.4
gamma = 1.4;
for i=1:length(M1)
    a1 = ((gamma + 1)^2).*(M1(i).^2);
    a2 = (4.*gamma.*(M1(i).^2))-(2*(gamma-1));
    a3 = (a1./a2).^(gamma/(gamma-1));
    a4 = (2.*gamma.*(M1(i).^2) - gamma+1)./(gamma+1);
    a5 = (a3 .* a4)-1;
    C_pstar = (2./(gamma.*(M1(i).^2))).*(a5);
    
    C_p = C_pstar .* (sin(theta).^2) ./ (sin(theta_0)^2);

    % Components of C_p that produce drag
    C_p_drag = C_p.*sin(theta);

    C_d(i) = sum(C_p_drag)./2;
end

figure(1);
plot(M1, C_d);
title("C_d vs M at \gamma = 1.4");
xlabel("Mach Number, M");
ylabel("Coefficient of Drag, C_d");

gamma = 1.05;
for i=1:length(M1)
    a1 = ((gamma + 1)^2).*(M1(i).^2);
    a2 = (4.*gamma.*(M1(i).^2))-(2*(gamma-1));
    a3 = (a1./a2).^(gamma/(gamma-1));
    a4 = (2.*gamma.*(M1(i).^2) - gamma+1)./(gamma+1);
    a5 = (a3 .* a4)-1;
    C_pstar = (2./(gamma.*(M1(i).^2))).*(a5);
    
    C_p = C_pstar .* (sin(theta).^2) ./ (sin(theta_0)^2);
    C_p_drag = C_p.*sin(theta);
    C_d(i) = sum(C_p_drag)./2;
end

figure(2);
plot(M1, C_d);
title("C_d vs M at \gamma = 1.05");
xlabel("Mach Number, M");
ylabel("Coefficient of Drag, C_d");

%% Problem 2
x = 0:0.01:2;
y_upper = 0.1.*x.*(2-x);
y_lower = -0.1.*x.*(2-x);

figure(1);
hold on;
plot(x, y_upper);
plot(x, y_lower);
title("Body equations");

y_u_dx = 0.2-0.2.*x;
y_l_dx = -0.2 + 0.2.*x;

theta_upper = atan2(y_u_dx, 1);
theta_lower = atan2(y_l_dx, 1);

figure(2);
hold on;
plot(x, rad2deg(theta_upper));
plot(x, rad2deg(theta_lower));
title("theta vs x");

% Newtonian Slender Body
Cp_N_upper = 2.*(theta_upper).^2;
Cp_N_lower = 2.*(theta_lower).^2;

% Need to cut off the values past x = 1, there C_p - 0
half = find(x == 1);
for j=half:length(x)
    Cp_N_upper(j) = 0;
    Cp_N_lower(j) = 0;
end

figure(3);
hold on;
plot(x, Cp_N_upper);
plot(x, Cp_N_lower);
title("C_p Newtonian");

% Tangent-wedge
gamma = 1.4;
M1 = 20;
Cp_TW_upper = (theta_upper.^2).*(((gamma+1)/2) + sqrt(((gamma+1)/2)^2 + (4./(M1.*theta_upper).^2)));
Cp_TW_lower = (theta_lower.^2).*(((gamma+1)/2) + sqrt(((gamma+1)/2)^2 + (4./(M1.*theta_lower).^2)));

for j=half:length(x)
    Cp_TW_upper(j) = 0;
    Cp_TW_lower(j) = 0;
end

figure(4);
hold on;
plot(x, Cp_TW_upper);
plot(x, Cp_TW_lower);
title("C_p Tangent Wedge");

% Shock-expansion
syms beta

RHS = (2*cot(beta)*(((M1^2)*(sin(beta)^2) - 1)/((M1^2)*(gamma+cos(2*beta)) + 2)));
LHS = tan(theta_upper(1));
beta_upper = vpasolve(RHS == LHS, beta, [0 pi/2]);

c1 = 2 + (gamma-1)*(M1^2)*sin(beta_upper)^2;
c2 = 2*gamma*(M1^2)*(sin(beta_upper)^2) - (gamma-1);
c3 = sin(beta_upper - theta_upper(1))^2;
M2 = sqrt(c1./(c2.*c3));

C_p2 = (4/(gamma+1)).*((sin(beta_upper).^2)-(1/(M1^2)));

c4 = sqrt((gamma+1)/(gamma-1));
c5 = atan(sqrt((M2^2 - 1)*(gamma-1)/(gamma+1)));
c6 = atan(sqrt(M2^2 - 1));
nu2 = c4*c5 - c6;

dtheta = theta_upper - theta_upper(1);

nu3 = nu2 - abs(dtheta);
syms M
M3 = nu3;
for i=1:length(nu3)
    M3(i) = vpasolve(sqrt((gamma+1)/(gamma-1))*atan(sqrt((M^2 - 1)*(gamma-1)/(gamma+1))) - atan(sqrt(M^2 - 1)) == nu3(i), M, [0 Inf]);
end

figure;
plot(x, M3);
title("M past shock");
%%
c7 = 1+(((gamma-1)/2)*M2^2);
c8 = 1+(((gamma-1)/2).*M3.^2);
p3_p2 = (c7./c8).^(gamma/(gamma-1));
C_p3 = (2/(gamma*(M2^2))).*((p3_p2) - 1);

figure;
plot(x, C_p3);
title("C_p past shock");

%% Problem 3
alpha = theta_upper;
for i=0:length(alpha)-1
    alpha(i+1) = 8*i/length(alpha);
end
alpha = deg2rad(alpha);

new_theta_upper = theta_upper;
Cp_Nn_l = theta_upper;
Cp_Nn_d = theta_upper;
Cp_TWn_l = theta_upper;
Cp_TWn_d = theta_upper;

for i=1:length(alpha)
    new_theta_upper = theta_upper - alpha(i);
    new_theta_lower = theta_lower - alpha(i);
    
    % Newtonian Slender Body
    Cp_Nn_upper = 2.*(new_theta_upper.^2);
    Cp_Nn_lower = 2.*(new_theta_lower.^2);

    % Cut off values past highest point of airfoil (Cp = 0)
    cutoff = find(Cp_Nn_upper == min(Cp_Nn_upper));
    for j=cutoff:length(x)
        Cp_Nn_upper(j) = 0;
    end

    cutoff = find(Cp_Nn_lower == min(Cp_Nn_lower));
    for j=cutoff:length(x)
        Cp_Nn_lower(j) = 0;
    end

    % Find the components of Cp
    Cp_Nn_upper_l = Cp_Nn_upper.*cos(new_theta_upper);
    Cp_Nn_upper_d = Cp_Nn_upper.*sin(new_theta_upper);

    Cp_Nn_lower_l = Cp_Nn_lower.*cos(new_theta_lower);
    Cp_Nn_lower_d = Cp_Nn_lower.*sin(new_theta_lower);

    Cp_Nn_upper_l_tot = sum(Cp_Nn_upper_l);
    Cp_Nn_lower_l_tot = sum(Cp_Nn_lower_l);
    Cp_Nn_l(i) = -Cp_Nn_upper_l_tot + Cp_Nn_lower_l_tot;

    Cp_Nn_upper_d_tot = sum(Cp_Nn_upper_d);
    Cp_Nn_lower_d_tot = sum(Cp_Nn_lower_d);
    Cp_Nn_d(i) = abs(Cp_Nn_upper_d_tot) + abs(Cp_Nn_lower_d_tot);

    % Tangent-wedge
    gamma = 1.4;
    M1 = 20;
    Cp_TWn_upper = (new_theta_upper.^2).*((gamma+1/2) + sqrt((gamma+1/2)^2 + (4./(M1.*new_theta_upper).^2)));
    Cp_TWn_lower = (new_theta_lower.^2).*((gamma+1/2) + sqrt((gamma+1/2)^2 + (4./(M1.*new_theta_lower).^2)));

    cutoff = find(Cp_TWn_upper == min(Cp_TWn_upper));
    for j=cutoff:length(x)
        Cp_TWn_upper(j) = 0;
    end

    cutoff = find(Cp_TWn_lower == min(Cp_TWn_lower));
    for j=cutoff:length(x)
        Cp_TWn_lower(j) = 0;
    end

    Cp_TWn_upper_l = Cp_TWn_upper.*cos(new_theta_upper);
    Cp_TWn_upper_d = Cp_TWn_upper.*sin(new_theta_upper);

    Cp_TWn_lower_l = Cp_TWn_lower.*cos(new_theta_lower);
    Cp_TWn_lower_d = Cp_TWn_lower.*sin(new_theta_lower);

    Cp_TWn_upper_l_tot = sum(Cp_TWn_upper_l);
    Cp_TWn_lower_l_tot = sum(Cp_TWn_lower_l);
    Cp_TWn_l(i) = -Cp_TWn_upper_l_tot + Cp_TWn_lower_l_tot;

    Cp_TWn_upper_d_tot = sum(Cp_TWn_upper_d);
    Cp_TWn_lower_d_tot = sum(Cp_TWn_lower_d);
    Cp_TWn_d(i) = abs(Cp_TWn_upper_d_tot) + abs(Cp_TWn_lower_d_tot);
end

figure;
hold on;
plot(alpha, Cp_Nn_l);
plot(alpha, Cp_TWn_l);
title("C_l");

figure;
hold on;
plot(alpha, Cp_Nn_d);
plot(alpha, Cp_TWn_d);
title("C_d");
% Benjamin Stutzke
% ENAE 488P
% Homework 4

%% Question 1
% Generating incident shock curve
% Using beta to generate theta because vpasolve doesn't handle well
M = 5;
beta_inc = 0.05:0.05:pi/2;
p2_p1_inc = zeros([1 length(beta_inc)]);
theta_inc = zeros([1 length(beta_inc)]);

for k=1:length(beta_inc)
    theta_inc(k) = tbm(M, beta_inc(k), false);
    p2_p1_inc(k) = p2p1(beta_inc(k), M);
end

% Clean up arrays by removing negative theta values
index = find(theta_inc > 0, 1, 'first');
theta_inc = theta_inc(index:end);
p2_p1_inc = p2_p1_inc(index:end);

% Insert theta endpoints
theta_inc = [0, theta_inc];
theta_inc(end + 1) = 0;

beta1 = tbm(M, 0, true);
beta2 = pi/2;

p1 = p2p1(beta1, M);
p2 = p2p1(beta2, M);

p2_p1_inc = [p1, p2_p1_inc];
p2_p1_inc(end + 1) = p2;

% Generating reflected shock curves and plotting
b = [25, 30, 35, 40];
b = deg2rad(b);
figureStartIndex = 1;

beta_refl = 0.1:0.1:pi/2;
p2_p1_refl = zeros([1 length(beta_refl)]);
theta_refl = zeros([1 length(beta_refl)]);

p2_p1_curves = zeros([length(b) length(beta_refl)]);
theta_curves = zeros([length(b) length(beta_refl)]);

p2_p1_refl_curves = zeros([length(b) length(beta_refl)]);
theta_refl_curves = zeros([length(b) length(beta_refl)]);

tdisp = b;
pdisp = b;

for i = 1:length(b)

    shock_theta = tbm(M, b(i), false);
    M2 = m2(b(i), shock_theta, M);

    for k=1:length(beta_refl)
        theta_refl(k) = tbm(M2, beta_refl(k), false);
        p2_p1_refl(k) = p2p1(beta_refl(k), M2);
    end

    % Clean up arrays by removing negative theta values
    index = find(theta_refl > 0, 1, 'first');
    theta_refl = theta_refl(index:end);
    p2_p1_refl = p2_p1_refl(index:end);
    
    % Insert theta endpoints
    theta_refl = [0, theta_refl];
    theta_refl(end + 1) = 0;
    
    beta1 = tbm(M2, 0, true);
    beta2 = pi/2;
    
    p1 = p2p1(beta1, M2);
    p2 = p2p1(beta2, M2);
    
    p2_p1_refl = [p1, p2_p1_refl];
    p2_p1_refl(end + 1) = p2;

    % Allocate to external arrays
    for h=1:length(theta_refl)
        theta_curves(i, h) = theta_refl(h);
        p2_p1_curves(i, h) = p2_p1_refl(h);
    end

    % Create reflection curve
    theta_disp = tbm(M, b(i), false);
    p2_p1_disp = p2p1(b(i), M);

    tdisp(i) = theta_disp;
    pdisp(i) = p2_p1_disp;

    theta_refl_curves(i, :) = -theta_curves(i, :) + theta_disp;
    p2_p1_refl_curves(i, :) = p2_p1_curves(i, :) - 1 + p2_p1_disp;
end

figure(1);
hold on;
plot(rad2deg(theta_inc), p2_p1_inc);
plot(rad2deg(theta_refl_curves(1, :)), p2_p1_refl_curves(1, :));
title("Reflected Shock for beta = 25 degrees");
xlabel("\theta (degrees)");
ylabel("p/p_1");
xline(0);

figure(2);
hold on;
plot(rad2deg(theta_inc), p2_p1_inc);
plot(rad2deg(theta_refl_curves(2, 1:14)), p2_p1_refl_curves(2, 1:14));
title("Reflected Shock for beta = 30 degrees");
xlabel("\theta (degrees)");
ylabel("p/p_1");
xline(0);

figure(3);
hold on;
plot(rad2deg(theta_inc), p2_p1_inc);
plot(rad2deg(theta_refl_curves(3, 1:14)), p2_p1_refl_curves(3, 1:14));
title("Reflected Shock for beta = 35 degrees");
xlabel("\theta (degrees)");
ylabel("p/p_1");
xline(0);

figure(4);
hold on;
plot(rad2deg(theta_inc), p2_p1_inc);
plot(rad2deg(theta_refl_curves(4, 1:13)), p2_p1_refl_curves(4, 1:13));
title("Reflected Shock for beta = 40 degrees");
xlabel("\theta (degrees)");
ylabel("p/p_1");
xline(0);

%% Question 2
n = 2/3;
Pr = 0.75;
Me = 2;
g = 1.4;
eta = 0:0.1:10;

% Me = 2
y0 = [0 0 0.460 1 0.288];

[X, Y] = ode45(@equation, eta, y0, [], n, Pr, Me, g);

eta_out = X;
u_ue = Y(:,2);
h_he = Y(:,4);

figure(5);
hold on;
plot(eta_out, u_ue);
plot(eta_out, h_he);
title("u/u_e and T/T_e for M = 2");
xlabel("sqrt(Re_x)*y/x");
legend("u/u_e", "T/T_e")

% Me = 4
Me = 4;
y0 = [0 0 0.439 1 1.097];

[X, Y] = ode45(@equation, eta, y0, [], n, Pr, Me, g);

eta_out = X;
u_ue = Y(:,2);
h_he = Y(:,4);

figure(6);
hold on;
plot(eta_out, u_ue);
plot(eta_out, h_he);
title("u/u_e and T/T_e for M = 4");
xlabel("sqrt(Re_x)*y/x");
legend("u/u_e", "T/T_e")

% Me = 6
Me = 6;
y0 = [0 0 0.415 1 2.330];

[X, Y] = ode45(@equation, eta, y0, [], n, Pr, Me, g);

eta_out = X;
u_ue = Y(:,2);
h_he = Y(:,4);

figure(7);
hold on;
plot(eta_out, u_ue);
plot(eta_out, h_he);
title("u/u_e and T/T_e for M = 6");
xlabel("sqrt(Re_x)*y/x");
legend("u/u_e", "T/T_e")
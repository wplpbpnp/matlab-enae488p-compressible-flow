%% Problem 1
close all; clear;
syms R delta Rc beta y

x = "R + delta - Rc*((cot(beta)^2)*(sqrt(1+ ((y^2)*(tan(beta)^2))/(Rc^2)) - 1))";
displayFormula(x)
x = str2sym(x);

M = 5.95;
R = 43.95;
beta = asin(1/M);
delta = (0.143*exp(3.24/(M^2))) * R;
Rc = (1.143*exp(.54/(M-1)^1.2)) * R;

x = subs(x)

Y = -500:0.1:500;
X = double(subs(x, y, Y));

figure(10);
imshow('world_cup.tif');
hold on;
plot(-X + 737.2, Y + 378.8);


%% Problem 2
clear; 
syms d_R Rc_R M gamma beta y theta

x = "1 + d_R - Rc_R*((cot(beta)^2)*(sqrt(1+ ((y^2)*(tan(beta)^2))/(Rc_R^2)) - 1))";
%displayFormula(x)
x = str2sym(x);

tbm = "tan(theta) == 2*cot(beta)*(((M^2)*(sin(beta)^2) - 1)/((M^2)*(gamma+cos(2*beta)) + 2))";
%displayFormula(tbm)
tbm = str2sym(tbm);

tbm_sonic = "asin(sqrt(((gamma+1)/(4*gamma))*(1+sqrt(1-((2*(-gamma+3))/((gamma+1)*M^2))+((gamma+9)/((gamma+1)*M^4)))) - ((3-gamma)/(4*gamma*M^2))))";
%displayFormula(tbm_sonic)
tbm_sonic = str2sym(tbm_sonic);

tbm_max = "asin(sqrt(((gamma+1)/(4*gamma))*(1+sqrt(1+((8*(gamma-1))/((gamma+1)*M^2))+((16)/((gamma+1)*M^4)))) - (1/(gamma*M^2))))";
%displayFormula(tbm_max)
tbm_max = str2sym(tbm_max);

theta = deg2rad(10);
g = 1.4;
tbm = subs(tbm);

M1 = [3 5 10];
beta_far = M1;
Y = -50:0.1:50;

for k=1:length(M1)
    figure(k+1);
    hold on;

    % Solving Theta-Beta-Mach to get relevant values of beta
    tbM_ = subs(tbm, [M, gamma], [M1(k), g]);
    beta_sonic = double(subs(tbm_sonic, [M, gamma], [M1(k), g]));
    beta_max = double(subs(tbm_max, [M, gamma], [M1(k), g]));
    beta_far = vpasolve(tbM_, beta, [0 pi/2]);

    % Solving for shock curve
    d_R_ = 0.386*exp(4.67/(M1(k)^2));
    Rc_R_ = 1.386*exp(1.8/(M1(k)-1)^.75);
    x_ = subs(x, [beta, d_R, Rc_R], [beta_far, d_R_, Rc_R_]);
    X = double(subs(x_, y, Y));
    plot(X, Y);

    % Finding sonic point by finding smallest difference between beta_sonic
    % and shock angle
    angle = abs(atan(Y./X));
    angle = angle(1:find(angle == max(angle)));
    diff = abs(beta_sonic-abs(angle));
    plot(X(find(diff == min(diff))), Y(find(diff == min(diff))), 'r-o');

    % Same process as above for max deflection angle
    diff = abs(beta_max-abs(angle));
    plot(X(find(diff == min(diff))), Y(find(diff == min(diff))), 'b+');
   
    % Mach number independence holds for Msin(beta) >= 3
    ind = M1(k)*sin(angle);
    for j = 1:length(ind)
        if ind(j) >= 3
            plot(X(j), Y(j), 'k+');
            break;
        end
    end

    title("M"+M1(k));
    legend("shock", "sonic point", "max deflection", "independence invalid");
    xlabel("x/R");
    ylabel("y/R");
end

%% Problem 3
clear;
gamma = 1.4;
x_d = [2 5 10];
r_rs = 0:0.01:1;
eta0_eta = 1./r_rs;

%syms v
%for i=2:length(eta0_eta)
%    eq = eta0_eta(i)^2 == v*(gamma+1-(gamma*v))*(((2*gamma*v)-gamma-1)/(gamma-1))^(-(gamma-1)/gamma)
%    v_h(i) = vpasolve(eq, v, [0, inf]);
%end

function [p, rho, T, M2] = oshock(M, beta, theta, gamma)
    % Outputs post-shock conditions for given M1 and beta

    Mn = M*sin(beta);
    p = 1 + ((2*gamma)/(gamma+1))*((Mn^2) - 1);
    rho = (gamma+1)/(gamma - 1 + (2/(Mn^2)));
    
    temp1 = (2*gamma*(Mn^2) - gamma + 1);
    temp2 = (((gamma-1)*(Mn^2)) + 2);
    temp3 = ((gamma+1)^2)*(Mn^2);
    T = temp1*temp2/temp3;

    temp4 = (((gamma - 1)*(Mn^2)) + 2)/((2*gamma*(Mn^2)) - (gamma - 1));
    temp5 = temp4/(sin(beta-theta)^2);
    M2 = sqrt(temp5);
end
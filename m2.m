function [M2] = m2(beta, theta, M)
    c1 = 2 + (0.4)*M^2*sin(beta)^2;
    c2 = 2*1.4*M^2*sin(beta)^2 - .4;
    c3 = 1/(sin(beta - theta));
    M2 = c3*sqrt(c1/c2);
end
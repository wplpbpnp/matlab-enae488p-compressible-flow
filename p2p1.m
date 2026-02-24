function [pratio] = p2p1(beta, M)
    C = 2*1.4/(2.4);
    pratio = 1 + C*(M^2*sin(beta)^2 - 1);
end
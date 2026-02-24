function dF = equation(eta, y, n, Pr, Me, g)
% F, U, S, H, Q
F = y(1);
U = y(2);
S = y(3);
H = y(4);
Q = y(5);

dF = zeros(5, 1);
dF(1) = U;
dF(2) = S;
dF(3) = (-F*S - S*(n-1)*(H^(n-2))*Q)/H^(n-1);
dF(4) = Q;
dF(5) = ((-F*Q) - (g-1)*(Me^2)*(H^(n-1))*S*S - (Q*Q*(n-1)*(H^(n-2))/Pr))*(Pr/(H^(n-1)));
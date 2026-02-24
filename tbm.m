function [output] = tbm(m, angle, theta_flag)
    % Solves Theta-Beta-Mach equation for Theta given M and Beta.
    % With theta_flag == true, uses angle for Theta and 
    % solves numerically for Beta.
    % Assume gamma = 1.4
    g = 1.4;
    s1 = sym("M");
    s2 = sym("theta");
    s3 = sym("beta");
    s4 = sym("gamma");

    eq = tan(s2) == 2*cot(s3)*(((s1^2)*(sin(s3)^2) - 1)/((s1^2)*(s4+cos(2*s3)) + 2));

    if theta_flag == true
        eq = subs(eq, [s1, s4, s2], [m, g, angle]);
        res = vpasolve(eq, s3, [0 inf]);
        while res > pi/2
            res = res - (2*pi);
        end
        if res < 0
            res = double(res + pi);
            if res < 0
                res = -1*res;
            end
            output = res;
        else
            output = res;
        end
    else
        eq = subs(eq, [s1, s4, s3], [m, g, angle]);
        output = double(solve(eq, s2));
    end
end
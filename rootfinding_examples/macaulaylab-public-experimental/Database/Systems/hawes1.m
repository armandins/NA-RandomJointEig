function [system] = hawes1(a,b)
    eqs = cell(4,1);
    eqs{1} = [1 2 0 0; -2 1 1 0; 1 0 2 0; a^2*b 1 0 0; -a^2 0 0 1; a^2 3 0 0; -2*a*b^2 2 0 0; 4*a*b 1 0 1; -4*a*b 4 0 0; -2*a 0 0 2; 4*a 3 0 1; -2*a 6 0 0; b^3 3 0 0; -3*b^2 2 0 1; 3*b^2 5 0 0; 3*b 1 0 2; -6*b 4 0 1; 3*b 7 0 0; -1 0 0 3; 3 3 0 2; -3 6 0 1; 1 9 0 0];
    eqs{2} = [4 3 0 0; -4 1 1 0; a^2*b 0 0 0; 3*a^2 2 0 0; -4*a*b^2 1 0 0; 4*a*b 0 0 1; -16*a*b 3 0 0; 12*a 2 0 1; -12*a 5 0 0; 3*b^3 2 0 0; -6*b^2 1 0 1; 15*b^2 4 0 0; 3*b 0 0 2; -24*b 3 0 1; 21*b 6 0 0; 9 2 0 2; -18 5 0 1; 9 8 0 0];
    eqs{3} = [12 2 0 0; -4 0 1 0; 6*a^2 1 0 0; -4*a*b^2 0 0; -48*a*b 2 0 0; 24*a 1 0 1; -60*a 4 0 0; 6*b^3 1 0 0; -6*b^2 0 0 1; 60*b^2 3 0 0; -72*b 2 0 1; 126*b 5 0 0; 18 1 0 2; -90 4 0 1; 72 7 0 0];
    eqs{4} = [6*a^2 0 0 0; -96*a*b 1 0 0; 24*a 0 0 1; -240*a 3 0 0; 6*b^3 0 0 0; 180*b^2 2 0 0; -144*b 1 0 1; 630*b 4 0 0; 18 0 0 2; -360 3 0 1; 504 6 0 0; 24 1 0 0];

    system = systemstruct(eqs);
end
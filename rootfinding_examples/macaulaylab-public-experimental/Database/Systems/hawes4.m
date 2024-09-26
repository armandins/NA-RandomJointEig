function [system] = hawes4(a,b,c)
    eqs = cell(5,1);
    eqs{1} = [3 0 0 0 2 0; 1 0 2 0 0 0; b 0 0 0 0 0];
    eqs{2} = [3 0 0 0 0 2; 1 0 0 2 0 0; b 0 0 0 0 0];
    eqs{3} = [1 1 0 0 0 0; 2 0 1 0 1 0; 3*a 0 2 0 0 0; 5 0 4 0 0 0; 2*c 0 1 0 0 0];
    eqs{4} = [1 1 0 0 0 0; 2 0 0 1 0 1; 3*a 0 0 2 0 0; 5 0 0 4 0 0; 2*c 0 0 1 0 0];
    eqs{5} = [1 1 1 0 0 0; 1 0 0 0 3 0; 1 0 2 0 1 0; a 0 3 0 0 0; 1 0 5 0 0 0; b 0 0 0 1 0; c 0 2 0 0 0; -1 1 0 1 0 0; -1 0 0 0 0 3; -1 0 0 2 0 1; -a 0 0 3 0 0; -1 0 0 5 0 0; -b 0 0 0 0 1; -c 0 0 2 0 0];
    eqs{6} = [162*a^2 0 1 2 1 1; 18*a*b 0 1 1 1 0; 54*a*c 0 2 2 1 1; -54*a 1 1 0 1 1; -9*a 0 3 2 1 1; -18*a 0 2 2 0 1; 810*a 0 1 4 1 1; -18*a 0 1 3 1 0; 54*a 0 1 2 0 1; 54*a 0 1 1 1 2; 54*a 0 0 2 2 1; 6*b*c 0 2 1 1 0; -b 0 3 1 1; -2*b 0 2 1 0 0; 6*b 0 1 1 0 0; 6*b 0 0 1 2 0; -18*c 1 2 0 1 1; 270*c 0 2 4 1 1; -6*c 0 2 3 1 0; 18*c 0 2 1 1 2; 3 1 3 0 1 1; 6 1 2 0 0 1; -18 1 1 0 0 1; -18 1 1 0 0 2 1; -45 0 3 4 1 1; 1 0 3 3 1 0; -3 0 3 1 1 2; -90 0 2 4 0 1; 2 0 2 3 0 0; -6 0 2 1 0 2; 270 0 1 4 0 1; -6 0 1 3 0 0; 18 0 1 1 0 2; 270 0 0 4 2 1; -6 0 0 3 2 0; -18 0 0 1 2 2; -162*a^2 0 2 1 1 1; -18*a*b 0 1 1 0 1; -54*a*c 0 2 2 1 1; 54*a 1 0 1 1 1; -810*a 0 4 1 1 1; 18*a 0 3 1 0 1; -540*a 0 2 3 1 1; 18*a 0 2 2 1 0; 54*a 0 2 0 1 2; -54*a 0 1 1 2 1; -6*b*c 0 1 2 0 1; -60*b 0 1 3 0 1; 2*b 0 3 2 0 1; -6*b 0 1 0 0 2; 18*c 1 0 2 1 1; -270*c 0 4 2 1 1; 6*c 0 3 2 0 1; -18*c 0 1 2 2 1; 180 1 0 3 1 1; -6 1 0 2 1 0; 18 1 0 0 1 2; -2700 0 4 3 1 1; 90 0 4 2 1 0; -270 0 4 0 1 2; 60 0 3 3 0 1; -2 0 3 2 0 0; 6 0 3 0 0 2; -180 0 1 3 2 1; 6 0 1 2 2 0; -18 0 1 0 2 2];

    system = systemstruct(eqs);
end
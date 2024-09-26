function [system] = pb601()
    eqs = cell(3,1);
    eqs{1} = [-62500000000000 2 0 1; 538354663.40508 0 6 0; 503135199.444 0 5 0; 89525836.2 0 4 0; 5775860.4 0 3 0; 107358 0 2 0; 617 0 1 0; 1 0 0 0];
    eqs{2} = [62500000000000 2 1 0; 187500000000000 2 0 1; 20250000 1 1 0; -503135199.444 0 5 0; -179051672.4 0 4 0; -17327581.2 0 3 0; -429432 0 2 0; -3085 0 1 0; -6 0 0 0];
    eqs{3} = [-5.55555555555555e+15 2 0 0; 1.11111111111111e+16 1 0 1; 1800000000 0 1 0; 1 0 0 0];

    system = systemstruct(eqs);
end
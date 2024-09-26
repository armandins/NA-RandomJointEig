function [system] = kotsireas7(m)
    eqs = cell(6,1);
    eqs{1} = [-(m^2+m) 1 0 0 1 0 0; (m^2+m) 1 0 0 0 1 0; (m^3-m) 1 0 0 0 0 0; (m^2+m) 0 1 0 1 0 0; -(m^2+m) 0 1 0 0 1 0; (m^3-m) 0 1 0 0 0 0; 2*(m^2+m) 0 0 1 0 0 0; -2*(m^3+m^2) 0 0 0 0 0 0];
    eqs{2} = [m 1 0 0 1 0 0; -m 1 0 0 0 1 0; (m^2+m) 1 0 0 0 0 0; m 0 1 0 1 0 0; -m 0 1 0 0 1 0; -(m^2+m) 0 1 0 0 0 0; -2*m 0 0 1 1 0 0; 2*m 0 0 1 0 1 0];
    eqs{3} = [1 0 0 0 2 0 0; -2 0 0 0 1 1 0; -2 0 0 0 1 0 0; 1 0 0 0 0 2 0; -2 0 0 0 0 1 0; 1 0 0 0 0 0 1; 1 0 0 0 0 0 0];
    eqs{4} = [1 2 0 0 3 0 0; -1 0 0 0 0 0 0];
    eqs{5} = [1 0 2 0 0 3 0; -1 0 0 0 0 0 0];
    eqs{6} = [1 0 0 2 0 0 3; -1 0 0 0 0 0 0];

    system = systemstruct(eqs);
end
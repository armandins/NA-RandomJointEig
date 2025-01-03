function [system] = hairer2()
    eqs = cell(11,1);
    eqs{1} = [1 0 0 0 1 0 0 0 0 0 0 0 0 0; 1 0 0 0 0 1 0 0 0 0 0 0 0 0; 1 0 0 0 0 1 0 0 0 0 0 0 0 0; 1 0 0 0 0 0 0 1 0 0 0 0 0 0; -1 0 0 0 0 0 0 0 0 0 0 0 0 0];
    eqs{2} = [1 0 0 1 1 0 0 0 0 0 0 0 0 0; 1 0 1 0 0 1 0 0 0 0 0 0 0 0; 1 1 0 0 0 0 1 0 0 0 0 0 0 0; -1/2 0 0 0 0 0 0 0 0 0 0 0 0 0];
    eqs{3} = [1 0 0 2 1 0 0 0 0 0 0 0 0 0; 1 0 2 0 0 1 0 0 0 0 0 0 0 0; 1 2 0 0 0 0 1 0 0 0 0 0 0 0; -1/3 0 0 0 0 0 0 0 0 0 0 0 0 0];
    eqs{4} = [1 1 0 0 0 1 0 0 0 1 0 0 0 0; 1 1 0 0 1 0 0 0 0 0 0 0 1 0; 1 0 1 0 1 0 0 0 0 0 0 0 0 1; -1/6 0 0 0 0 0 0 0 0 0 0 0 0 0];
    eqs{5} = [1 0 0 3 1 0 0 0 0 0 0 0 0 0; 1 0 3 0 0 1 0 0 0 0 0 0 0 0; 1 3 0 0 0 0 1 0 0 0 0 0 0 0; -1/4 0 0 0 0 0 0 0 0 0 0 0 0 0];
    eqs{6} = [1 1 1 0 0 1 0 0 0 1 0 0 0 0; 1 1 0 1 1 0 0 0 0 0 0 0 1 0; 1 0 1 1 1 0 0 0 0 0 0 0 0 1; -1/8 0 0 0 0 0 0 0 0 0 0 0 0 0];
    eqs{7} = [1 2 0 0 0 1 0 0 0 1 0 0 0 0; 1 2 0 0 1 0 0 0 0 0 0 0 1 0; 1 0 2 0 1 0 0 0 0 0 0 0 0 1; -1/2 0 0 0 0 0 0 0 0 0 0 0 0 0];
    eqs{8} = [1 1 0 0 1 0 0 0 0 1 0 0 0 1; -1/24 0 0 0 0 0 0 0 0 0 0 0 0 0];
    eqs{9} = [1 1 0 0 0 0 0 0 0 0 0 0 0 0; -1 0 0 0 0 0 0 0 1 0 0 0 0 0];
    eqs{10} = [1 0 1 0 0 0 0 0 0 0 0 0 0 0; -1 0 0 0 0 0 0 0 0 1 0 0 0 0; -1 0 0 0 0 0 0 0 0 0 1 0 0 0];
    eqs{11} = [1 0 0 1 0 0 0 0 0 0 0 0 0 0; -1 0 0 0 0 0 0 0 0 0 1 0 0 0; -1 0 0 0 0 0 0 0 0 0 0 0 1 0; -1 0 0 0 0 0 0 0 0 0 0 0 0 1];

    system = systemstruct(eqs);
end
function [system] = ku10()
    eqs = cell(10,1);
    eqs{1} = [5 1 0 1 0 0 0 0 0 0 0; 5 1 0 0 0 0 0 0 0 0 0; 3 0 0 1 0 0 0 0 0 0 0; 55 0 0 0 0 0 0 0 0 0 0];
    eqs{2} = [7 0 0 1 1 0 0 0 0 0 0; 9 0 0 1 0 0 0 0 0 0 0; 9 0 0 0 1 0 0 0 0 0 0; 19 0 0 0 0 0 0 0 0 0 0];
    eqs{3} = [3 0 0 0 1 1 0 0 0 0 0; 6 0 0 0 1 0 0 0 0 0 0; 5 0 0 0 0 1 0 0 0 0 0; -4 0 0 0 0 0 0 0 0 0 0];
    eqs{4} = [6 0 0 0 0 1 1 0 0 0 0; 6 0 0 0 0 1 0 0 0 0 0; 7 0 0 0 0 0 1 0 0 0 0; 118 0 0 0 0 0 0 0 0 0 0];
    eqs{5} = [1 0 0 0 0 0 1 1 0 0 0; 3 0 0 0 0 0 1 0 0 0 0; 9 0 0 0 0 0 0 1 0 0 0; 27 0 0 0 0 0 0 0 0 0 0];
    eqs{6} = [6 0 0 0 0 0 0 1 1 0 0; 7 0 0 0 0 0 0 1 0 0 0; 1 0 0 0 0 0 0 0 1 0 0; 72 0 0 0 0 0 0 0 0 0 0];
    eqs{7} = [9 0 0 0 0 0 0 0 1 1 0; 7 0 0 0 0 0 0 0 1 0 0; 1 0 0 0 0 0 0 0 0 1 0; 35 0 0 0 0 0 0 0 0 0 0];
    eqs{8} = [4 0 0 0 0 0 0 0 0 1 1; 4 0 0 0 0 0 0 0 0 1 0; 6 0 0 0 0 0 0 0 0 0 1; 16 0 0 0 0 0 0 0 0 0 0];
    eqs{9} = [8 0 1 0 0 0 0 0 0 0 1; 3 0 1 0 0 0 0 0 0 0 0; 4 0 0 0 0 0 0 0 0 0 1; -51 0 0 0 0 0 0 0 0 0 0];
    eqs{10} = [3 1 1 0 0 0 0 0 0 0 0; -6 1 0 0 0 0 0 0 0 0 0; 1 0 1 0 0 0 0 0 0 0 0; 5 0 0 0 0 0 0 0 0 0 0];

    system = systemstruct(eqs);
end
function [system] = sparse5()
    eqs = cell(5,1);
    eqs{1} = [1 2 2 2 2 2; 3 2 0 0 0 0; 1 1 1 1 1 1; 1 0 2 0 0 0; 1 0 0 2 0 0; 1 0 0 0 2 0; 1 0 0 0 0 2; 5 0 0 0 0 0];
    eqs{2} = [1 2 2 2 2 2; 1 2 0 0 0 0; 1 1 1 1 1 1; 3 0 2 0 0 0; 1 0 0 2 0 0; 1 0 0 0 2 0; 1 0 0 0 0 2; 5 0 0 0 0 0];
    eqs{3} = [1 2 2 2 2 2; 1 2 0 0 0 0; 1 1 1 1 1 1; 1 0 2 0 0 0; 3 0 0 2 0 0; 1 0 0 0 2 0; 1 0 0 0 0 2; 5 0 0 0 0 0];
    eqs{4} = [1 2 2 2 2 2; 1 2 0 0 0 0; 1 1 1 1 1 1; 1 0 2 0 0 0; 1 0 0 2 0 0; 3 0 0 0 2 0; 1 0 0 0 0 2; 5 0 0 0 0 0];
    eqs{5} = [1 2 2 2 2 2; 1 2 0 0 0 0; 1 1 1 1 1 1; 1 0 2 0 0 0; 1 0 0 2 0 0; 1 0 0 0 2 0; 3 0 0 0 0 2; 5 0 0 0 0 0];

    system = systemstruct(eqs);
end
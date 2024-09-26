function [system] = caeqsrasse() 
    eqs = cell(4,1); 
    eqs{1} = [2 1 1 1 0; -2 0 1 0 0; 1 0 0 2 1; -1 0 0 0 1]; 
    eqs{2} = [4 1 2 1 0; 2 1 0 3 0; -10 1 0 1 0; -1 0 3 0 1; 4 0 2 0 0; ...
        4 0 1 2 1; 4 0 1 0 1; -10 0 0 2 0; 2 0 0 0 0]; 
    eqs{3} = [1 2 1 0 0; 2 1 0 1 1; -1 0 1 0 0; -2 0 0 0 1]; 
    eqs{4} = [2 3 0 1 0; 4 2 1 0 1; -10 2 0 0 0; 4 1 0 1 2; -10 1 0 1 0; ...
        -1 0 1 0 3; 4 0 1 0 1; 4 0 0 0 2; 2 0 0 0 0]; 

    system = systemstruct(eqs); 
end
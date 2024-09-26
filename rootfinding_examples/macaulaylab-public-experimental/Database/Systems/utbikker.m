function [system] = utbikker() 
    eqs = cell(4,1); 
    eqs{1} = [1 2 0 0 0; -3 1 1 0 0; 2 1 0 1 0; -2 1 0 0 0; 1 0 2 0 0; ...
        1 0 1 1 0; -2 0 1 0 1; -3 0 1 0 0; 1 0 0 2 0; -4 0 0 1 1; ...
        -2 0 0 1 0; 3 0 0 0 2; 3 0 0 0 1; -2 0 0 0 0]; 
    eqs{2} = [-3 2 0 0 0; -1 1 1 0 0; 1 1 0 1 0; -5 1 0 0 1; 2 1 0 0 0; ...
        2 0 2 0 0; -1 0 1 1 0; -1 0 1 0 1; -5 0 1 0 0; 1 0 0 2 0; ...
        -1 0 0 1 1; 1 0 0 1 0; -6 0 0 0 2; 5 0 0 0 1; 5 0 0 0 0]; 
    eqs{3} = [2 3 0 0 0; 1 2 1 0 0; -3 2 0 0 0; 1 1 1 1 0; -3 1 1 0 0; ...
        -5 1 0 0 2; -2 1 0 0 1; -2 1 0 0 0; 1 0 3 0 0; -1 0 2 0 1; ...
        1 0 2 0 0; 1 0 1 1 1; -3 0 1 1 0; -5 0 1 0 2; 2 0 1 0 1; ...
        -1 0 1 0 0; 1 0 0 3 0; -5 0 0 2 1; -1 0 0 2 0; 7 0 0 1 2; ...
        1 0 0 1 0; -3 0 0 0 3; 2 0 0 0 2; 11 0 0 0 1; -3 0 0 0 0]; 
    eqs{4} = [4 3 0 0 0; 11 2 1 0 0; 5 2 0 0 0; 6 1 1 1 0; -1 1 1 0 1; ...
        -7 1 0 2 0; 2 1 0 1 0; -1 1 0 0 1; -10 1 0 0 0; -1 0 3 0 0; ...
        6 0 2 1 0; -1 0 2 0 1; 3 0 2 0 0; -12 0 1 2 0; -4 0 1 1 1; ...
        2 0 1 1 0; 5 0 1 0 2; -35 0 1 0 0; 6 0 0 3 0; 6 0 0 2 1; ...
        2 0 0 2 0; 4 0 0 1 2; -14 0 0 1 0; 15 0 0 0 3; -1 0 0 0 2; ...
        4 0 0 0 1; -15 0 0 0 0]; 

    system = systemstruct(eqs); 
end
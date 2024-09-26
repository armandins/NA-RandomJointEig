function [system] = quadgrid() 
    eqs = cell(5,1); 
    eqs{1} = [1 0 1 0 0 0; 1 0 0 1 0 0; 1 0 0 0 1 0; 1 0 0 0 0 1; ...
        -1 0 0 0 0 0]; 
    eqs{2} = [1 1 1 0 0 0; 1 1 0 1 0 0; 1 1 0 0 1 0; 1 1 0 0 0 1; ...
        0.5 0 0 1 0 0; 1 0 0 0 1 0; 1.5 0 0 0 0 1; ...
        -0.6339745962 0 0 0 0 0]; 
    eqs{3} = [1 2 1 0 0 0; 1 2 0 1 0 0; 1 2 0 0 1 0; 1 2 0 0 0 1; ...
        1 1 0 1 0 0; 2 1 0 0 1 0; 3 1 0 0 0 1; 0.25 0 0 1 0 0; ...
        1 0 0 0 1 0; 4.5 0 0 0 0 1; -0.4019237886 0 0 0 0 0]; 
    eqs{4} = [1 3 1 0 0 0; 1 3 0 1 0 0; 1 3 0 0 1 0; 1 3 0 0 0 1; ...
        1.5 2 0 1 0 0; 3 2 0 0 1 0; 4.5 2 0 0 0 1; 0.75 1 0 1 0 0; ...
        3 1 0 0 1 0; 6.75 1 0 0 0 1; 0.125 0 0 1 0 0; 1 0 0 0 1 0; ...
        3.375 0 0 0 0 1; -0.1310915568 0 0 0 0 0]; 
    eqs{5} = [1 4 1 0 0 0; 1 4 0 1 0 0; 1 4 0 0 1 0; 1 4 0 0 0 1; ...
        2 3 0 1 0 0; 4 3 0 0 1 0; 6 3 0 0 0 1; 1.5 2 0 1 0 0; ...
        6 2 0 0 1 0; 13.5 2 0 0 0 1; 0.5 1 0 1 0 0; 4 1 0 0 1 0; ...
        13.5 1 0 0 0 1; 0.0625 0 0 1 0 0; 1 0 0 0 1 0; ...
        5.0625 0 0 0 0 1; 0.3021933285 0 0 0 0 0]; 

    system = systemstruct(eqs); 
end
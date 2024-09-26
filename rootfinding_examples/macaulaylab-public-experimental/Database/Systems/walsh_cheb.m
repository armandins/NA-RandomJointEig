function [system] = walsh_cheb() 
    eqs = cell(6,1); 
    eqs{1} = [-1 0 0 0 1 0 0; -1 0 0 0 0 1 0; 1 0 0 0 0 0 0]; 
    eqs{2} = [2 0 1 0 0 1 0; 1 0 1 0 0 0 0; -1 0 0 1 0 0 0; -12 0 0 0 1 0 0; ...
        9 0 0 0 0 0 0]; 
    eqs{3} = [-10 0 0 0 0 0 0; 1 1 0 0 0 0 0; 9 0 1 0 0 0 0; -12 0 0 1 0 0 0; ...
        -49 0 0 0 1 0 0; -0.5 0 0 0 0 1 0; -2 1 0 0 0 1 0; -0.5 0 2 0 0 1 0]; 
    eqs{4} = [-1 1 0 0 0 0 0; -10 0 1 0 0 0 0; -127 0 0 1 0 0 0; -78 0 0 0 1 0 0; ...
        -0.5 0 0 0 0 1 0; -0.5 2 0 0 0 1 0; 2 1 1 0 0 1 0]; 
    eqs{5} = [-10 1 0 0 0 0 0; -78 0 0 1 0 0 0; -0.5 0 0 0 0 1 0; -0.5 2 0 0 0 1 0]; 
    eqs{6} = [-0.024302413273002 0 0 0 0 0 0; 0.000980392156863 1 0 0 0 0 0; ...
        1 0 0 0 0 0 1;  -0.000037707390649 2 0 0 0 0 0; -0.000245098039216 0 2 0 0 0 0; ...
        -0.024302413273002 0 0 0 0 2 0; 0.000980392156863 1 0 0 0 2 0; ...
        -0.000037707390649 2 0 0 0 2 0; -0.000245098039216 0 2 0 0 2 0]; 
    
    system = systemstruct(eqs); 
end
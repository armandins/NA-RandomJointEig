function [system] = walsh() 
    eqs = cell(6,1); 
    eqs{1} = [-1 0 0 0 1 0 0; -1 0 0 0 0 1 0; 1 0 0 0 0 0 0]; 
    eqs{2} = [2 0 1 0 0 1 0; 1 0 1 0 0 0 0; -1 0 0 1 0 0 0; -12 0 0 0 1 0 0; ...
        9 0 0 0 0 0 0]; 
    eqs{3} = [-2 1 0 0 0 1 0; 1 1 0 0 0 0 0; -1 0 2 0 0 1 0; 9 0 1 0 0 0 0; ...
        -12 0 0 1 0 0 0; -49 0 0 0 1 0 0; -10 0 0 0 0 0 0]; 
    eqs{4} = [2 1 1 0 0 1 0; 9 1 0 0 0 0 0; -10 0 1 0 0 0 0; -49 0 0 1 0 0 0; ...
        -78 0 0 0 1 0 0; -1 2 0 0 0 1 0; -10 1 0 0 0 0 0; -78 0 0 1 0 0 0]; 
    eqs{5} = [-1 2 0 0 0 1 0; -10 1 0 0 0 0 0; -78 0 0 1 0 0 0]; 
    eqs{6} = [-0.000150829562594 2 0 0 0 2 0; 0.001960784313725 1 0 0 0 2 0; ...
        -0.000980392156863 0 2 0 0 2 0; -0.048039215686275 0 0 0 0 2 0; 1 0 0 0 0 0 1]; 
    
    system = systemstruct(eqs); 
end
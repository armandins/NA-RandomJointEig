function [system] = batselier9() 
    eqs = cell(3,1); 
    eqs{1} = [1 1 4; 3 3 0; -1 0 4; -3 2 0]; 
    eqs{2} = [1 2 1; -2 2 0]; 
    eqs{3} = [2 1 4; -1 3 0; -2 0 4; 1 2 0]; 

    system = systemstruct(eqs); 
end
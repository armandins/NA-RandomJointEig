function [system] = bender1() 
    eqs = cell(2,1); 
    eqs{1} = [1 0 0; 3 1 0; 2 0 1; 4 1 1]; 
    eqs{2} = [3 0 0; -2 1 0; 4 0 1; -4 1 1]; 

    system = systemstruct(eqs); 
end
function [system] = phcpack3() 
    eqs = cell(2,1); 
    eqs{1} = [1 2 0]; 
    eqs{2} = [1 1 0; -1 0 1]; 

    system = systemstruct(eqs); 
end
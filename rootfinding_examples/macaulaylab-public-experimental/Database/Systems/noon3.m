function [system] = noon3() 
    eqs = cell(3,1); 
    eqs{1} = [1 1 2 0; 1 1 0 2; -1.1 1 0 0; 1 0 0 0]; 
    eqs{2} = [1 2 1 0; 1 0 1 2; -1.1 0 1 0; 1 0 0 0]; 
    eqs{3} = [1 2 0 1; 1 0 2 1; -1.1 0 0 1; 1 0 0 0]; 

    system = systemstruct(eqs); 
end
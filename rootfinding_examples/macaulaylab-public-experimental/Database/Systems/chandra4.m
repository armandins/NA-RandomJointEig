function [system] = chandra4() 
    eqs = cell(4,1); 
    eqs{1} = [-0.25617 2 0 0 0; -0.17078 1 1 0 0; -0.128085 1 0 1 0; ...
        7.48766 1 0 0 0; -8 0 0 0 0]; 
    eqs{2} = [-0.34156 1 1 0 0; -0.25617 0 2 0 0; -0.204936 0 1 1 0; ...
        7.48766 0 1 0 0; -8 0 0 0 0]; 
    eqs{3} = [-0.384255 1 0 1 0; -0.307404 0 1 1 0; -0.25617 0 0 2 0; ...
        7.48766 0 0 1 0; -8 0 0 0 0]; 
    eqs{4} = [-0.409872 1 0 0 1; -0.34156 0 1 0 1; -0.2927657143 0 0 1 1; ...
        7.48766 0 0 0 1; -8 0 0 0 0]; 

    system = systemstruct(eqs); 
end
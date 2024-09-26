function [system] = batselier10() 
    eqs = cell(2,1); 
    eqs{1} = [2 5 2; 2 4 3; 2 4 2; 800 4 1; 400 3 2; 400 3 1; ...
        80000 3 0; 1 2 4; -2 2 3; -700 2 2; 3 1 5; 1 1 4; 900 1 3; ...
        -300 1 2; -40000 1 1; 3 0 5; 600 0 4; 700 0 3; 140000 0 2; ...
        20000 0 1; 4000000 0 0]; 
    eqs{2} = [-2 5 2; -2 4 3; -2 4 2; -800 4 1; -400 3 2; -400 3 1; ...
        -80000 3 0; -1 2 4; 2 2 3; 900 2 2; 3 1 5; -1 1 4; -700 1 3; ...
        500 1 2; 120000 1 1; -3 0 5; -600 0 4; -600 0 4; -500 0 3; ...
        -100000 0 2; 20000 0 1; 4000000 0 0]; 

    system = systemstruct(eqs); 
end
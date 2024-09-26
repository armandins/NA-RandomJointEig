function [system] = dreesen1() 
    % DREESEN1 contains a system of multivariate polynomial equations. 
    % 
    %   [system] = DREESEN1() returns the system of multivariate polynomial 
    %   equations. 

    % MacaulayLab (2022) - Christof Vermeersch. 

    eqs = cell(2,1); 
    eqs{1} = [-1 2 0; 2 1 1; 1 0 2; 5 1 0; -3 0 1; -4 0 0]; 
    eqs{2} = [1 2 0; 2 1 1; 1 0 2; -1 0 0]; 

    system = systemstruct(eqs); 
end
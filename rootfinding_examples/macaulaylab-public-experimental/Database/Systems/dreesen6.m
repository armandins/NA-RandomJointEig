function [system] = dreesen6() 
    % DREESEN6 contains a system of multivariate polynomial equations. 
    % 
    %   [system] = DREESEN6() returns the system of multivariate polynomial 
    %   equations. 

    % MacaulayLab (2022) - Christof Vermeersch. 

    eqs = cell(2,1); 
    eqs{1} = [1 2 0; -1 0 1]; 
    eqs{2} = [1 1 0; -5 0 0]; 

    system = systemstruct(eqs); 
end
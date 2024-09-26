function [system] = gaukwa2()
    eqs = cell(4,1);
    eqs{1} = [1 1 0 0 0; 1 0 1 0 0; -0.998250904334731+0.059119641363025i 0 0 0 0];
    eqs{2} = [1 1 0 1 0; 1 0 1 0 1; -0.892749639148806+0.450553084330444i 0 0 0 0];
    eqs{3} = [1 1 0 2 0; 1 0 1 0 2; 0.160088552022675+0.98710265702777i 0 0 0 0];
    eqs{4} = [1 1 0 3 0; 1 0 1 0 3; -0.725369971319578+0.688359211972815i 0 0 0 0];

    system = systemstruct(eqs);
end
function [system] = puma()
    eqs = cell(8,1);
    eqs{1} = [1 2 0 0 0 0 0 0 0; 1 0 2 0 0 0 0 0 0; -1 0 0 0 0 0 0 0 0];
    eqs{2} = [1 0 0 2 0 0 0 0 0; 1 0 0 0 2 0 0 0 0; -1 0 0 0 0 0 0 0 0];
    eqs{3} = [1 0 0 0 0 2 0 0 0; 1 0 0 0 0 0 2 0 0; -1 0 0 0 0 0 0 0 0];
    eqs{4} = [1 0 0 0 0 0 0 2 0; 1 0 0 0 0 0 0 0 2; -1 0 0 0 0 0 0 0 0];
    eqs{5} = [0.004731 1 0 1 0 0 0 0 0; -0.1238 1 0 0 0 0 0 0 0; -0.3578 0 1 1 0 0 0 0 0; -0.001637 0 1 0 0 0 0 0 0; -0.9338 0 0 0 1 0 0 0 0; 1 0 0 0 0 0 0 1 0; -0.3571 0 0 0 0 0 0 0 0];
    eqs{6} = [0.2238 1 0 1 0 0 0 0 0; 0.2638 1 0 0 0 0 0 0 0; 0.7623 0 1 1 0 0 0 0 0; -0.07745 0 1 0 0 0 0 0 0; -0.6734 0 0 0 1 0 0 0 0; -0.6022 0 0 0 0 0 0 0 0];
    eqs{7} = [0.3578 1 0 0 0 0 0 0 0; 0.004731 0 1 0 0 0 0 0 0; 1 0 0 0 0 0 1 0 1];
    eqs{8} = [-0.7623 1 0 0 0 0 0 0 0; 0.2238 0 1 0 0 0 0 0 0; 0.3461 0 0 0 0 0 0 0 0];

    system = systemstruct(eqs);
end
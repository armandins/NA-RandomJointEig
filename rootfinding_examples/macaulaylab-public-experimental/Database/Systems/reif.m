function [system] = reif()
    eqs = cell(14,1);
    M = zeros(5, 17);
    M(1,1) = 1;
    M(1,5) = 1;
    M(1,14) = 1;
    M(2,1) = 1;
    M(2,6) = 1;
    M(2,15) = 1;
    M(3,1) = 1;
    M(3,7) = 1;
    M(4,1) = -1;
    M(4,7) = 1;
    M(4,14) = 1;
    M(5,1) = -1;
    M(5,7) = 1;
    M(5,15) = 1;
    eqs{1} = M;
    M = zeros(4, 17);
    M(1,1) = 1;
    M(1,5) = 1;
    M(1,16) = 1;
    M(2,1) = 1;
    M(2,6) = 1;
    M(2,17) = 1;
    M(3,1) = -1;
    M(3,7) = 1;
    M(3,16) = 1;
    M(4,1) = -1;
    M(4,7) = 1;
    M(4,17) = 1;
    eqs{2} = M;
    M = zeros(5, 17);
    M(1,1) = 1;
    M(1,8) = 1;
    M(1,14) = 1;
    M(2,1) = 1;
    M(2,9) = 1;
    M(2,15) = 1;
    M(3,1) = 1;
    M(3,10) = 1;
    M(4,1) = -1;
    M(4,10) = 1;
    M(4,14) = 1;
    M(5,1) = -1;
    M(5,10) = 1;
    M(5,15) = 1;
    eqs{3} = M;
    M = zeros(4, 17);
    M(1,1) = 1;
    M(1,8) = 1;
    M(1,16) = 1;
    M(2,1) = 1;
    M(2,9) = 1;
    M(2,17) = 1;
    M(3,1) = -1;
    M(3,10) = 1;
    M(3,16) = 1;
    M(4,1) = -1;
    M(4,10) = 1;
    M(4,17) = 1;
    eqs{4} = M;
    M = zeros(5, 17);
    M(1,1) = 1;
    M(1,11) = 1;
    M(1,14) = 1;
    M(2,1) = 1;
    M(2,12) = 1;
    M(2,15) = 1;
    M(3,1) = 1;
    M(3,13) = 1;
    M(4,1) = -1;
    M(4,13) = 1;
    M(4,14) = 1;
    M(5,1) = -1;
    M(5,13) = 1;
    M(5,15) = 1;
    eqs{5} = M;
    M = zeros(4, 17);
    M(1,1) = 1;
    M(1,11) = 1;
    M(1,16) = 1;
    M(2,1) = 1;
    M(2,12) = 1;
    M(2,17) = 1;
    M(3,1) = -1;
    M(3,13) = 1;
    M(3,16) = 1;
    M(4,1) = -1;
    M(4,13) = 1;
    M(4,17) = 1;
    eqs{6} = M;
    M = zeros(5, 17);
    M(1,1) = 1;
    M(1,2) = 1;
    M(1,14) = 1;
    M(2,1) = 1;
    M(2,3) = 1;
    M(2,15) = 1;
    M(3,1) = 1;
    M(3,4) = 1;
    M(4,1) = -1;
    M(4,4) = 1;
    M(4,14) = 1;
    M(5,1) = -1;
    M(5,4) = 1;
    M(5,15) = 1;
    eqs{7} = M;
    M = zeros(4, 17);
    M(1,1) = 1;
    M(1,2) = 1;
    M(1,16) = 1;
    M(2,1) = 1;
    M(2,3) = 1;
    M(2,17) = 1;
    M(3,1) = -1;
    M(3,4) = 1;
    M(3,16) = 1;
    M(4,1) = -1;
    M(4,4) = 1;
    M(4,17) = 1;
    eqs{8} = M;
    M = zeros(6, 17);
    M(1,1) = 1;
    M(1,2) = 1;
    M(1,5) = 1;
    M(1,14) = 1;
    M(2,1) = 1;
    M(2,3) = 1;
    M(2,6) = 1;
    M(2,15) = 1;
    M(3,1) = 1;
    M(3,4) = 1;
    M(3,7) = 1;
    M(4,1) = -1;
    M(4,4) = 1;
    M(4,7) = 1;
    M(4,14) = 1;
    M(5,1) = -1;
    M(5,4) = 1;
    M(5,7) = 1;
    M(5,15) = 1;
    M(6,1) = -1;
    eqs{9} = M;
    M = zeros(4, 17);
    M(1,1) = 1;
    M(1,2) = 1;
    M(1,5) = 1;
    M(1,16) = 1;
    M(2,1) = 1;
    M(2,3) = 1;
    M(2,6) = 1;
    M(2,17) = 1;
    M(3,1) = -1;
    M(3,4) = 1;
    M(3,7) = 1;
    M(3,16) = 1;
    M(4,1) = -1;
    M(4,4) = 1;
    M(4,10) = 1;
    M(4,17) = 1;
    eqs{10} = M;
    M = zeros(6, 17);
    M(1,1) = 1;
    M(1,2) = 1;
    M(1,8) = 1;
    M(1,14) = 1;
    M(2,1) = 1;
    M(2,3) = 1;
    M(2,9) = 1;
    M(2,15) = 1;
    M(3,1) = 1;
    M(3,4) = 1;
    M(3,10) = 1;
    M(4,1) = -1;
    M(4,4) = 1;
    M(4,10) = 1;
    M(4,14) = 1;
    M(5,1) = -1;
    M(5,4) = 1;
    M(5,10) = 1;
    M(5,15) = 1;
    M(6,1) = -1;
    eqs{11} = M;
    M = zeros(4, 17);
    M(1,1) = 1;
    M(1,2) = 1;
    M(1,8) = 1;
    M(1,16) = 1;
    M(2,1) = 1;
    M(2,3) = 1;
    M(2,9) = 1;
    M(2,17) = 1;
    M(3,1) = -1;
    M(3,4) = 1;
    M(3,10) = 1;
    M(3,16) = 1;
    M(4,1) = -1;
    M(4,4) = 1;
    M(4,10) = 1;
    M(4,17) = 1;
    eqs{12} = M;
    M = zeros(6, 17);
    M(1,1) = 1;
    M(1,2) = 1;
    M(1,11) = 1;
    M(1,14) = 1;
    M(2,1) = 1;
    M(2,3) = 1;
    M(2,12) = 1;
    M(2,15) = 1;
    M(3,1) = 1;
    M(3,4) = 1;
    M(3,13) = 1;
    M(4,1) = -1;
    M(4,4) = 1;
    M(4,13) = 1;
    M(4,14) = 1;
    M(5,1) = -1;
    M(5,4) = 1;
    M(5,13) = 1;
    M(5,15) = 1;
    M(6,1) = -1;
    eqs{13} = M;
    M = zeros(4, 17);
    M(1,1) = 1;
    M(1,2) = 1;
    M(1,11) = 1;
    M(1,16) = 1;
    M(2,1) = 1;
    M(2,3) = 1;
    M(2,12) = 1;
    M(2,17) = 1;
    M(3,1) = -1;
    M(3,4) = 1;
    M(3,13) = 1;
    M(3,16) = 1;
    M(4,1) = -1;
    M(4,4) = 1;
    M(4,13) = 1;
    M(4,17) = 1;
    eqs{14} = M;

    system = systemstruct(eqs);
end
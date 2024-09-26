function [system] = rbpl24()
    eqs = cell(9,1);
    eqs{1} = [62500 2 0 0 0 0 0 0 0 0; 62500 0 0 0 2 0 0 0 0 0; 62500 0 0 0 0 0 0 2 0 0; -74529 0 0 0 0 0 0 0 0 0];
    eqs{2} = [625 0 2 0 0 0 0 0 0 0; -1250 0 1 0 0 0 0 0 0 0; 625 0 0 0 0 2 0 0 0 0; 625 0 0 0 0 0 0 0 2 0; -2624 0 0 0 0 0 0 0 0 0];
    eqs{3} = [12500 0 0 2 0 0 0 0 0 0; 2500 0 0 1 0 0 0 0 0 0; 12500 0 0 0 0 0 2 0 0 0; -44975 0 0 0 0 0 1 0 0 0; 12500 0 0 0 0 0 0 0 0 2; -10982 0 0 0 0 0 0 0 0 0];
    eqs{4} = [400000 1 1 0 0 0 0 0 0 0; -400000 0 1 0 0 0 0 0 0 0; 400000 0 0 0 1 1 0 0 0 0; 400000 0 0 0 0 0 0 1 1 0; 178837 0 0 0 0 0 0 0 0 0];
    eqs{5} = [1000000 1 0 1 0 0 0 0 0 0; 100000 0 0 1 0 0 0 0 0 0; 1000000 0 0 0 1 0 1 0 0 0; -1799000 0 0 0 0 0 1 0 0 0; 1000000 0 0 0 0 0 0 1 0 1; -805427 0 0 0 0 0 0 0 0 0];
    eqs{6} = [2000000 0 1 1 0 0 0 0 0 0; -2000000 0 1 0 0 0 0 0 0 0; 200000 0 0 1 0 0 0 0 0 0; 2000000 0 0 0 0 1 1 0 0 0; -3598000 0 0 0 0 0 1 0 0 0; 2000000 0 0 0 0 0 0 0 1 1; -1403 0 0 0 0 0 0 0 0 0];
    eqs{7} = [-113800000000000 1 0 0 0 1 0 0 0 1; 206888400000000 1 0 0 0 1 0 0 0 0; 113800000000000 1 0 0 0 0 1 0 1 0; -206888400000000 1 0 0 0 0 1 0 0 0; 2014260000000 1 0 0 0 0 0 0 1 0; -2014260000000 1 0 0 0 0 0 0 0 1; -362960716800000 1 0 0 0 0 0 0 0 0; 113800000000000 0 1 0 1 0 0 0 0 1; -206888400000000 0 1 0 1 0 0 0 0 0; -113800000000000 0 1 0 0 0 1 1 0 0; 206888400000000 0 1 0 0 0 1 0 0 0; -2014260000000 0 1 0 0 0 0 1 0 0; 2014260000000 0 1 0 0 0 0 0 0 1; 38025201600000 0 1 0 0 0 0 0 0 0; -113800000000000 0 0 1 1 0 0 0 1 0; 206888400000000 0 0 1 1 0 0 0 0 0; 113800000000000 0 0 1 0 1 0 1 0 0; -206888400000000 0 0 1 0 1 0 0 0 0; 2014260000000 0 0 1 0 0 0 1 0 0; -2014260000000 0 0 1 0 0 0 0 1 0; 292548849600000 0 0 1 0 0 0 0 0 0; 61907200000000 0 0 0 1 0 0 0 1 0; -61907200000000 0 0 0 1 0 0 0 0 1; 11809567440000 0 0 0 1 0 0 0 0 0; -61907200000000 0 0 0 0 1 0 1 0 0; 61907200000000 0 0 0 0 1 0 0 0 1; 1475978220000 0 0 0 0 1 0 0 0 0; 61907200000000 0 0 0 0 0 1 1 0 0; -61907200000000 0 0 0 0 0 1 0 1 0; -825269402280000 0 0 0 0 0 1 0 0 0; -1.2129826896e+15 0 0 0 0 0 0 1 0 0; -151600474800000 0 0 0 0 0 0 0 1 0; 825859951200000 0 0 0 0 0 0 0 0 1; -19295432410527 0 0 0 0 0 0 0 0 0];
    eqs{8} = [777600000000 1 0 0 0 1 0 0 0 1; 1409011200000 1 0 0 0 1 0 0 0 0; -777600000000 1 0 0 0 0 1 0 1 0; -1409011200000 1 0 0 0 0 1 0 0 0; 1065312000000 1 0 0 0 0 0 0 1 0; -1065312000000 1 0 0 0 0 0 0 0 1; 235685027200 1 0 0 0 0 0 0 0 0; -777600000000 0 1 0 1 0 0 0 0 1; -1409011200000 0 1 0 1 0 0 0 0 0; 777600000000 0 1 0 0 0 1 1 0 0; 1409011200000 0 1 0 0 0 1 0 0 0; -1065312000000 0 1 0 0 0 0 1 0 0; 1065312000000 0 1 0 0 0 0 0 0 1; 398417510400 0 1 0 0 0 0 0 0 0; 777600000000 0 0 1 1 0 0 0 1 0; 1409011200000 0 0 1 1 0 0 0 0 0; -777600000000 0 0 1 0 1 0 1 0 0; -1409011200000 0 0 1 0 1 0 0 0 0; 1065312000000 0 0 1 0 0 0 1 0 0; -1065312000000 0 0 1 0 0 0 0 1 0; 158626915200 0 0 1 0 0 0 0 0 0; 805593600000 0 0 0 1 0 0 0 1 0; -805593600000 0 0 0 1 0 0 0 0 1; -311668424000 0 0 0 1 0 0 0 0 0; -805593600000 0 0 0 0 1 0 1 0 0; 805593600000 0 0 0 0 1 0 0 0 1; -268090368000 0 0 0 0 1 0 0 0 0; 805593600000 0 0 0 0 0 1 1 0 0; -805593600000 0 0 0 0 0 1 0 1 0; 72704002800 0 0 0 0 0 1 0 0 0; 412221302400 0 0 0 0 0 0 1 0 0; 354583756800 0 0 0 0 0 0 0 1 0; 307085438400 0 0 0 0 0 0 0 0 1; 282499646407 0 0 0 0 0 0 0 0 0];
    eqs{9} = [3200 0 1 0 0 0 0 0 0 0; 1271 0 0 0 0 0 0 0 0 0];

    system = systemstruct(eqs);
end
function [system] = kruppa()
    eqs = cell(6,1);
    eqs{1} = [-6.15115567158004e-05 2 0 0 0 0; -0.00043664138510483 1 1 0 0 0; 2.26071559827849e-07 1 0 1 0 0; 2.11511465115345e-07 1 0 0 1 0; 5.76239877677818e-07 1 0 0 0 1; 0.079081046812384 1 0 0 0 0; 1.80797761045375e-05 0 2 0 0 0; -1.41387960643007e-07 0 1 1 0 0; 8.50811727332718e-07 0 1 0 1 0; -3.79317162317744e-08 0 1 0 0 1; 0.0419624438965366 0 1 0 0 0; -7.61153829922538e-10 0 0 1 1 0; 2.8786696090516e-10 0 0 1 0 1; 5.74460303784112e-07 0 0 1 0 0; -1.48478169590623e-11 0 0 0 2 0; -6.79218657484955e-10 0 0 0 1 1; -0.000175161395035901 0 0 0 1 0; 2.28267982419871e-12 0 0 0 0 2; -6.65943977391716e-05 0 0 0 0 1; -6.97692477559067 0 0 0 0 0];
    eqs{2} = [0.000130546289981528 2 0 0 0 0; 1.40622165312777e-05 1 1 0 0 0; -3.44609936855853e-07 1 0 1 0 0; -2.59013301909069e-07 1 0 0 1 0; -2.71953420104016e-08 1 0 0 0 1; -0.0214033108963866 1 0 0 0 0; -3.08742081009702e-06 0 2 0 0 0; 4.5206337405973e-08 0 1 1 0 0; 8.39472055445605e-08 0 1 0 1 0; 2.6800860494094e-10 0 1 0 0 1; -0.00346746304174839 0 1 0 0 0; -1.62584884422345e-10 0 0 2 0 0; 3.04870430133298e-10 0 0 1 1 0; -1.70262311565956e-12 0 0 1 0 1; 4.14428595359516e-05 0 0 1 0 0; 1.57602690882905e-10 0 0 0 2 0; 8.76086414050974e-11 0 0 0 1 1; -5.43945527306181e-06 0 0 0 1 0; 2.11006021466014e-06 0 0 0 0 1; 1.18267820506517 0 0 0 0 0];
    eqs{3} = [-2.73338631829837e-05 2 0 0 0 0; -0.000301424125478683 1 1 0 0 0; 1.0694256762866e-07 1 0 1 0 0; 8.88388814932693e-08 1 0 0 1 0; 4.82696454079478e-07 1 0 0 0 1; 0.0468617055588229 1 0 0 0 0; 8.8812735132687e-06 0 2 0 0 0; -7.36189802491257e-08 0 1 1 0 0; 5.46904606080249e-07 0 1 0 1 0; -2.03849283077562e-08 0 1 0 0 1; 0.030514264226846 0 1 0 0 0; -3.39695973068974e-10 0 0 1 1 0; 1.62584884422345e-10 0 0 1 0 1; -8.31437565649053e-08 0 0 1 0 0; -6.39949110828445e-12 0 0 0 2 0; -6.26570189478977e-10 0 0 0 1 1; -9.57036922257151e-05 0 0 0 1 0; 1.70262311565956e-12 0 0 0 0 2; -5.51207104648816e-05 0 0 0 0 1; -4.40144726692841 0 0 0 0 0];
    eqs{4} = [0.00020432250740847 2 0 0 0 0; 1.56849356143919e-05 1 1 0 0 0; -4.04587350244954e-07 1 0 1 0 0; -4.87305302729242e-07 1 0 0 1 0; -4.37601758133548e-08 1 0 0 0 1; -0.0336358923748597 1 0 0 0 0; -5.46658290219168e-06 0 2 0 0 0; 7.95604828075377e-08 0 1 1 0 0; 1.04006712183549e-07 0 1 0 1 0; 3.38991257849075e-10 0 1 0 0 1; -0.00404932285692824 0 1 0 0 0; -2.8786696090516e-10 0 0 2 0 0; 3.39304230142816e-10 0 0 1 1 0; -2.28267982419871e-12 0 0 1 0 1; 4.87371291653175e-05 0 0 1 0 0; 3.67088373419699e-10 0 0 0 2 0; 1.49222268850375e-10 0 0 0 1 1; 4.47962288468087e-06 0 0 0 1 0; 3.20769606792884e-06 0 0 0 0 1; 1.80765948728211 0 0 0 0 0];
    eqs{5} = [-1.58964327245376e-07 2 0 0 0 0; -6.3493815217782e-06 1 1 0 0 0; 5.48812398775846e-10 1 0 1 0 0; 3.60872410814741e-10 1 0 0 1 0; -6.09097372823032e-09 1 0 0 0 1; 0.00180088751701278 1 0 0 0 0; -1.14740252397262e-06 0 2 0 0 0; 4.40298856360943e-09 0 1 1 0 0; 1.94297453317849e-08 0 1 0 1 0; 3.11269996945191e-09 0 1 0 0 1; 0.00050930863127144 0 1 0 0 0; -1.23884787674437e-12 0 0 1 1 0; -1.25885325939368e-11 0 0 1 0 1; -4.45780648509274e-07 0 0 1 0 0; -4.60106346383124e-15 0 0 0 2 0; -8.23833445341353e-13 0 0 0 1 1; -4.47136812555249e-06 0 0 0 1 0; 4.79826571419022e-13 0 0 0 0 2; 1.53044969512683e-06 0 0 0 0 1; -0.171772909691699 0 0 0 0 0];
    eqs{6} = [-1.82442494606386e-05 2 0 0 0 0; -1.77258860204536e-05 1 1 0 0 0; 7.75702888189264e-08 1 0 1 0 0; 3.89076811646233e-08 1 0 0 1 0; -2.71833797487333e-09 1 0 0 0 1; 0.00660929284056171 1 0 0 0 0; 7.26639981391655e-07 0 2 0 0 0; -6.06889914479223e-09 0 1 1 0 0; 3.80801465844382e-08 0 1 0 1 0; 1.06282125755784e-10 0 1 0 0 1; 0.00208594002826038 0 1 0 0 0; 1.25885325939368e-11 0 0 2 0 0; -1.65497689996062e-10 0 0 1 1 0; -4.79826571419022e-13 0 0 1 0 1; -9.18833264708608e-06 0 0 1 0 0; -6.16356638207838e-13 0 0 0 2 0; 6.32558972702553e-12 0 0 0 1 1; -9.49191582964126e-06 0 0 0 1 0; 2.91764102389389e-07 0 0 0 0 1; -0.524497069977953 0 0 0 0 0];

    system = systemstruct(eqs);
end
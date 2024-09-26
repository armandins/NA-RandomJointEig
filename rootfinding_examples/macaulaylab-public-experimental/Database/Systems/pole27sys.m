function [system] = pole27sys()
    eqs = cell(14,1);
    eqs{1} = [0.0712448724229539 1 0 0 0 0 0 0 0 1 0 0 0 0 0; 0.107391200992692 1 0 0 0 0 0 0 0 0 1 0 0 0 0; 0.0606046247191059 1 0 0 0 0 0 0 0 0 0 1 0 0 0; -0.056429815411452 1 0 0 0 0 0 0 0 0 0 0 1 0 0; -0.152165362133568 1 0 0 0 0 0 0 0 0 0 0 0 1 0; -0.124372362421262 1 0 0 0 0 0 0 0 0 0 0 0 0 1; 0.0563452683679251 1 0 0 0 0 0 0 0 0 0 0 0 0 0; -0.0712448724229539 0 1 0 0 0 0 0 1 0 0 0 0 0 0; 0.143578946113947 0 1 0 0 0 0 0 0 0 1 0 0 0 0; 0.193137079981164 0 1 0 0 0 0 0 0 0 0 1 0 0 0; 0.031085886947234 0 1 0 0 0 0 0 0 0 0 0 1 0 0; -0.232220459198047 0 1 0 0 0 0 0 0 0 0 0 0 1 0; -0.321605074744818 0 1 0 0 0 0 0 0 0 0 0 0 0 1; 0.0417330238116869 0 1 0 0 0 0 0 0 0 0 0 0 0 0; -0.107391200992692 0 0 1 0 0 0 0 1 0 0 0 0 0 0; -0.143578946113947 0 0 1 0 0 0 0 0 1 0 0 0 0 0; 0.168990053868345 0 0 1 0 0 0 0 0 0 0 1 0 0 0; 0.160579755009589 0 0 1 0 0 0 0 0 0 0 0 1 0 0; -0.0433812507912106 0 0 1 0 0 0 0 0 0 0 0 0 1 0; -0.234126357908529 0 0 1 0 0 0 0 0 0 0 0 0 0 1; -0.0506455353188055 0 0 1 0 0 0 0 0 0 0 0 0 0 0; -0.0606046247191059 0 0 0 1 0 0 0 1 0 0 0 0 0 0; -0.193137079981164 0 0 0 1 0 0 0 0 1 0 0 0 0 0; -0.168990053868345 0 0 0 1 0 0 0 0 0 1 0 0 0 0; 0.179418361633886 0 0 0 1 0 0 0 0 0 0 0 1 0 0; 0.21496480257607 0 0 0 1 0 0 0 0 0 0 0 0 1 0; 0.0635857696373987 0 0 0 1 0 0 0 0 0 0 0 0 0 1; -0.11724557954495 0 0 0 1 0 0 0 0 0 0 0 0 0 0; 0.056429815411452 0 0 0 0 1 0 0 1 0 0 0 0 0 0; -0.031085886947234 0 0 0 0 1 0 0 0 1 0 0 0 0 0; -0.160579755009589 0 0 0 0 1 0 0 0 0 1 0 0 0 0; -0.179418361633886 0 0 0 0 1 0 0 0 0 0 1 0 0 0; 0.250324722121778 0 0 0 0 1 0 0 0 0 0 0 0 1 0; 0.3089954329644 0 0 0 0 1 0 0 0 0 0 0 0 0 1; -0.0576396494666057 0 0 0 0 1 0 0 0 0 0 0 0 0 0; 0.152165362133568 0 0 0 0 0 1 0 1 0 0 0 0 0 0; 0.232220459198047 0 0 0 0 0 1 0 0 1 0 0 0 0 0; 0.0433812507912106 0 0 0 0 0 1 0 0 0 1 0 0 0 0; -0.21496480257607 0 0 0 0 0 1 0 0 0 0 1 0 0 0; -0.250324722121778 0 0 0 0 0 1 0 0 0 0 0 1 0 0; 0.28149879236875 0 0 0 0 0 1 0 0 0 0 0 0 0 1; 0.0945219379834688 0 0 0 0 0 1 0 0 0 0 0 0 0 0; 0.124372362421262 0 0 0 0 0 0 1 1 0 0 0 0 0 0; 0.321605074744818 0 0 0 0 0 0 1 0 1 0 0 0 0 0; 0.234126357908529 0 0 0 0 0 0 1 0 0 1 0 0 0 0; -0.0635857696373987 0 0 0 0 0 0 1 0 0 0 1 0 0 0; -0.3089954329644 0 0 0 0 0 0 1 0 0 0 0 1 0 0; -0.28149879236875 0 0 0 0 0 0 1 0 0 0 0 0 1 0; 0.181493615509306 0 0 0 0 0 0 1 0 0 0 0 0 0 0; -0.0762045088275659 0 0 0 0 0 0 0 1 0 0 0 0 0 0; -0.010407713551175 0 0 0 0 0 0 0 0 1 0 0 0 0 0; 0.137885940063516 0 0 0 0 0 0 0 0 0 1 0 0 0 0; 0.197728766482615 0 0 0 0 0 0 0 0 0 0 1 0 0 0; 0.0414933734919406 0 0 0 0 0 0 0 0 0 0 0 1 0 0; -0.22615736379943 0 0 0 0 0 0 0 0 0 0 0 0 1 0; -0.325824498608256 0 0 0 0 0 0 0 0 0 0 0 0 0 1; 0.0364070996284862 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    eqs{2} = [0.104143143049402 1 0 0 0 0 0 0 0 1 0 0 0 0 0; 0.13852403963522 1 0 0 0 0 0 0 0 0 1 0 0 0 0; 0.000622912817688869 1 0 0 0 0 0 0 0 0 0 1 0 0 0; -0.197816126487464 1 0 0 0 0 0 0 0 0 0 0 1 0 0; -0.195973956736754 1 0 0 0 0 0 0 0 0 0 0 0 1 0; 0.0743571310980861 1 0 0 0 0 0 0 0 0 0 0 0 0 1; 0.0584052525730866 1 0 0 0 0 0 0 0 0 0 0 0 0 0; -0.104143143049402 0 1 0 0 0 0 0 1 0 0 0 0 0 0; 0.180993113756719 0 1 0 0 0 0 0 0 0 1 0 0 0 0; 0.136892030161293 0 1 0 0 0 0 0 0 0 0 1 0 0 0; -0.154310073351042 0 1 0 0 0 0 0 0 0 0 0 1 0 0; -0.312425784925551 0 1 0 0 0 0 0 0 0 0 0 0 1 0; -0.0501533334257578 0 1 0 0 0 0 0 0 0 0 0 0 0 1; 0.0199499186252916 0 1 0 0 0 0 0 0 0 0 0 0 0 0; -0.13852403963522 0 0 1 0 0 0 0 1 0 0 0 0 0 0; -0.180993113756719 0 0 1 0 0 0 0 0 1 0 0 0 0 0; 0.181001778219758 0 0 1 0 0 0 0 0 0 0 1 0 0 0; 0.138537224293507 0 0 1 0 0 0 0 0 0 0 0 1 0 0; -0.0749789658775526 0 0 1 0 0 0 0 0 0 0 0 0 1 0; -0.195937729909761 0 0 1 0 0 0 0 0 0 0 0 0 0 1; -0.0749680197464147 0 0 1 0 0 0 0 0 0 0 0 0 0 0; -0.000622912817688869 0 0 0 1 0 0 0 1 0 0 0 0 0 0; -0.136892030161293 0 0 0 1 0 0 0 0 1 0 0 0 0 0; -0.181001778219758 0 0 0 1 0 0 0 0 0 1 0 0 0 0; 0.25909847389687 0 0 0 1 0 0 0 0 0 0 0 1 0 0; 0.255731275152658 0 0 0 1 0 0 0 0 0 0 0 0 1 0; -0.0980394818925705 0 0 0 1 0 0 0 0 0 0 0 0 0 1; -0.0766520608370986 0 0 0 1 0 0 0 0 0 0 0 0 0 0; 0.197816126487464 0 0 0 0 1 0 0 1 0 0 0 0 0 0; 0.154310073351042 0 0 0 0 1 0 0 0 1 0 0 0 0 0; -0.138537224293507 0 0 0 0 1 0 0 0 0 1 0 0 0 0; -0.25909847389687 0 0 0 0 1 0 0 0 0 0 1 0 0 0; 0.303064628411147 0 0 0 0 1 0 0 0 0 0 0 0 1 0; 0.205440241922479 0 0 0 0 1 0 0 0 0 0 0 0 0 1; 0.0486455760226156 0 0 0 0 1 0 0 0 0 0 0 0 0 0; 0.195973956736754 0 0 0 0 0 1 0 1 0 0 0 0 0 0; 0.312425784925551 0 0 0 0 0 1 0 0 1 0 0 0 0 0; 0.0749789658775526 0 0 0 0 0 1 0 0 0 1 0 0 0 0; -0.255731275152658 0 0 0 0 0 1 0 0 0 0 1 0 0 0; -0.303064628411147 0 0 0 0 0 1 0 0 0 0 0 1 0 0; 0.317446077342152 0 0 0 0 0 1 0 0 0 0 0 0 0 1; 0.137672457057934 0 0 0 0 0 1 0 0 0 0 0 0 0 0; -0.0743571310980861 0 0 0 0 0 0 1 1 0 0 0 0 0 0; 0.0501533334257578 0 0 0 0 0 0 1 0 1 0 0 0 0 0; 0.195937729909761 0 0 0 0 0 0 1 0 0 1 0 0 0 0; 0.0980394818925705 0 0 0 0 0 0 1 0 0 0 1 0 0 0; -0.205440241922479 0 0 0 0 0 0 1 0 0 0 0 1 0 0; -0.317446077342152 0 0 0 0 0 0 1 0 0 0 0 0 1 0; 0.0423708819565526 0 0 0 0 0 0 1 0 0 0 0 0 0 0; -0.0689249826639022 0 0 0 0 0 0 0 1 0 0 0 0 0 0; 0.0572227116175119 0 0 0 0 0 0 0 0 1 0 0 0 0 0; 0.195900256154442 0 0 0 0 0 0 0 0 0 1 0 0 0 0; 0.0909414224394049 0 0 0 0 0 0 0 0 0 0 1 0 0 0; -0.2108193938362 0 0 0 0 0 0 0 0 0 0 0 1 0 0; -0.314452800844557 0 0 0 0 0 0 0 0 0 0 0 0 1 0; 0.00766348133208784 0 0 0 0 0 0 0 0 0 0 0 0 0 1; 0.0452949141173368 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    eqs{3} = [0.114632606545634 1 0 0 0 0 0 0 0 1 0 0 0 0 0; 0.104355127922911 1 0 0 0 0 0 0 0 0 1 0 0 0 0; -0.108498626086074 1 0 0 0 0 0 0 0 0 0 1 0 0 0; -0.226439615085336 1 0 0 0 0 0 0 0 0 0 0 1 0 0; 0.0108456873899591 1 0 0 0 0 0 0 0 0 0 0 0 1 0; 0.306846878338656 1 0 0 0 0 0 0 0 0 0 0 0 0 1; 0.0474799712724454 1 0 0 0 0 0 0 0 0 0 0 0 0 0; -0.114632606545634 0 1 0 0 0 0 0 1 0 0 0 0 0 0; 0.151696780234845 0 1 0 0 0 0 0 0 0 1 0 0 0 0; 0.0397970231052206 0 1 0 0 0 0 0 0 0 0 1 0 0 0; -0.185188921460442 0 1 0 0 0 0 0 0 0 0 0 1 0 0; -0.120404672788264 0 1 0 0 0 0 0 0 0 0 0 0 1 0; 0.183888999265512 0 1 0 0 0 0 0 0 0 0 0 0 0 1; 0.004576233612467 0 1 0 0 0 0 0 0 0 0 0 0 0 0; -0.104355127922911 0 0 1 0 0 0 0 1 0 0 0 0 0 0; -0.151696780234845 0 0 1 0 0 0 0 0 1 0 0 0 0 0; 0.179808488137732 0 0 1 0 0 0 0 0 0 0 1 0 0 0; 0.131068701915729 0 0 1 0 0 0 0 0 0 0 0 1 0 0; -0.123962119644981 0 0 1 0 0 0 0 0 0 0 0 0 1 0; -0.238657431349847 0 0 1 0 0 0 0 0 0 0 0 0 0 1; -0.0586657280706742 0 0 1 0 0 0 0 0 0 0 0 0 0 0; 0.108498626086074 0 0 0 1 0 0 0 1 0 0 0 0 0 0; -0.0397970231052206 0 0 0 1 0 0 0 0 1 0 0 0 0 0; -0.179808488137732 0 0 0 1 0 0 0 0 0 1 0 0 0 0; 0.253892561770563 0 0 0 1 0 0 0 0 0 0 0 1 0 0; 0.110196530296909 0 0 0 1 0 0 0 0 0 0 0 0 1 0; -0.280577202672743 0 0 0 1 0 0 0 0 0 0 0 0 0 1; -0.0208149901260111 0 0 0 1 0 0 0 0 0 0 0 0 0 0; 0.226439615085336 0 0 0 0 1 0 0 1 0 0 0 0 0 0; 0.185188921460442 0 0 0 0 1 0 0 0 1 0 0 0 0 0; -0.131068701915729 0 0 0 0 1 0 0 0 0 1 0 0 0 0; -0.253892561770563 0 0 0 0 1 0 0 0 0 0 1 0 0 0; 0.255362673788986 0 0 0 0 1 0 0 0 0 0 0 0 1 0; 0.132465698011361 0 0 0 0 1 0 0 0 0 0 0 0 0 1; 0.0676642041641225 0 0 0 0 1 0 0 0 0 0 0 0 0 0; -0.0108456873899591 0 0 0 0 0 1 0 1 0 0 0 0 0 0; 0.120404672788264 0 0 0 0 0 1 0 0 1 0 0 0 0 0; 0.123962119644981 0 0 0 0 0 1 0 0 0 1 0 0 0 0; -0.110196530296909 0 0 0 0 0 1 0 0 0 0 1 0 0 0; -0.255362673788986 0 0 0 0 0 1 0 0 0 0 0 1 0 0; 0.339695674349431 0 0 0 0 0 1 0 0 0 0 0 0 0 1; 0.0503036873888375 0 0 0 0 0 1 0 0 0 0 0 0 0 0; -0.306846878338656 0 0 0 0 0 0 1 1 0 0 0 0 0 0; -0.183888999265512 0 0 0 0 0 0 1 0 1 0 0 0 0 0; 0.238657431349847 0 0 0 0 0 0 1 0 0 1 0 0 0 0; 0.280577202672743 0 0 0 0 0 0 1 0 0 0 1 0 0 0; -0.132465698011361 0 0 0 0 0 0 1 0 0 0 0 1 0 0; -0.339695674349431 0 0 0 0 0 0 1 0 0 0 0 0 1 0; -0.0639158580154462 0 0 0 0 0 0 1 0 0 0 0 0 0 0; -0.0405423244682057 0 0 0 0 0 0 0 1 0 0 0 0 0 0; 0.0990982148072104 0 0 0 0 0 0 0 0 1 0 0 0 0 0; 0.143864363422955 0 0 0 0 0 0 0 0 0 1 0 0 0 0; -0.0797203920066704 0 0 0 0 0 0 0 0 0 0 1 0 0 0; -0.261249847323615 0 0 0 0 0 0 0 0 0 0 0 1 0 0; -0.0332078033264689 0 0 0 0 0 0 0 0 0 0 0 0 1 0; 0.330301006649526 0 0 0 0 0 0 0 0 0 0 0 0 0 1; 0.042664226937987 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    eqs{4} = [0.133691512146562 1 0 0 0 0 0 0 0 1 0 0 0 0 0; 0.0503388945210272 1 0 0 0 0 0 0 0 0 1 0 0 0 0; -0.206927694163866 1 0 0 0 0 0 0 0 0 0 1 0 0 0; -0.130228450502176 1 0 0 0 0 0 0 0 0 0 0 1 0 0; 0.257323211908165 1 0 0 0 0 0 0 0 0 0 0 0 1 0; 0.233882169890492 1 0 0 0 0 0 0 0 0 0 0 0 0 1; 0.0461818739467038 1 0 0 0 0 0 0 0 0 0 0 0 0 0; -0.133691512146562 0 1 0 0 0 0 0 1 0 0 0 0 0 0; 0.140842317176365 0 1 0 0 0 0 0 0 0 1 0 0 0 0; 0.00301731534119662 0 1 0 0 0 0 0 0 0 0 1 0 0 0; -0.151903372926202 0 1 0 0 0 0 0 0 0 0 0 1 0 0; -0.0103332213591955 0 1 0 0 0 0 0 0 0 0 0 0 1 0; 0.165643079002682 0 1 0 0 0 0 0 0 0 0 0 0 0 1; 0.000308329423446791 0 1 0 0 0 0 0 0 0 0 0 0 0 0; -0.0503388945210272 0 0 1 0 0 0 0 1 0 0 0 0 0 0; -0.140842317176365 0 0 1 0 0 0 0 0 1 0 0 0 0 0; 0.219131819083488 0 0 1 0 0 0 0 0 0 0 1 0 0 0; 0.0799978150607582 0 0 1 0 0 0 0 0 0 0 0 1 0 0; -0.274977519351958 0 0 1 0 0 0 0 0 0 0 0 0 1 0; -0.184022282914927 0 0 1 0 0 0 0 0 0 0 0 0 0 1; -0.0485359247695727 0 0 1 0 0 0 0 0 0 0 0 0 0 0; 0.206927694163866 0 0 0 1 0 0 0 1 0 0 0 0 0 0; -0.00301731534119662 0 0 0 1 0 0 0 0 1 0 0 0 0 0; -0.219131819083488 0 0 0 1 0 0 0 0 0 1 0 0 0 0; 0.238055165102804 0 0 0 1 0 0 0 0 0 0 0 1 0 0; 0.0101861694309639 0 0 0 1 0 0 0 0 0 0 0 0 1 0; -0.261660864551533 0 0 0 1 0 0 0 0 0 0 0 0 0 1; -0.00151952184637281 0 0 0 1 0 0 0 0 0 0 0 0 0 0; 0.130228450502176 0 0 0 0 1 0 0 1 0 0 0 0 0 0; 0.151903372926202 0 0 0 0 1 0 0 0 1 0 0 0 0 0; -0.0799978150607582 0 0 0 0 1 0 0 0 0 1 0 0 0 0; -0.238055165102804 0 0 0 0 1 0 0 0 0 0 1 0 0 0; 0.302442111530847 0 0 0 0 1 0 0 0 0 0 0 0 1 0; 0.104389940203804 0 0 0 0 1 0 0 0 0 0 0 0 0 1; 0.0521725653746122 0 0 0 0 1 0 0 0 0 0 0 0 0 0; -0.257323211908165 0 0 0 0 0 1 0 1 0 0 0 0 0 0; 0.0103332213591955 0 0 0 0 0 1 0 0 1 0 0 0 0 0; 0.274977519351958 0 0 0 0 0 1 0 0 0 1 0 0 0 0; -0.0101861694309639 0 0 0 0 0 1 0 0 0 0 1 0 0 0; -0.302442111530847 0 0 0 0 0 1 0 0 0 0 0 1 0 0; 0.336899213941113 0 0 0 0 0 1 0 0 0 0 0 0 0 1; 0.00416292579016346 0 0 0 0 0 1 0 0 0 0 0 0 0 0; -0.233882169890492 0 0 0 0 0 0 1 1 0 0 0 0 0 0; -0.165643079002682 0 0 0 0 0 0 1 0 1 0 0 0 0 0; 0.184022282914927 0 0 0 0 0 0 1 0 0 1 0 0 0 0; 0.261660864551533 0 0 0 0 0 0 1 0 0 0 1 0 0 0; -0.104389940203804 0 0 0 0 0 0 1 0 0 0 0 1 0 0; -0.336899213941113 0 0 0 0 0 0 1 0 0 0 0 0 1 0; -0.0566797017879628 0 0 0 0 0 0 1 0 0 0 0 0 0 0; -0.0171870498597946 0 0 0 0 0 0 0 1 0 0 0 0 0 0; 0.131220654885444 0 0 0 0 0 0 0 0 1 0 0 0 0 0; 0.0675148817453912 0 0 0 0 0 0 0 0 0 1 0 0 0 0; -0.202715403227657 0 0 0 0 0 0 0 0 0 0 1 0 0 0; -0.147349918388012 0 0 0 0 0 0 0 0 0 0 0 1 0 0; 0.251239007277294 0 0 0 0 0 0 0 0 0 0 0 0 1 0; 0.25085427502714 0 0 0 0 0 0 0 0 0 0 0 0 0 1; 0.0453679887295185 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    eqs{5} = [0.142420903705832 1 0 0 0 0 0 0 0 1 0 0 0 0 0; -0.0273865492488036 1 0 0 0 0 0 0 0 0 1 0 0 0 0; -0.232972818871747 1 0 0 0 0 0 0 0 0 0 1 0 0 0; 0.0724666332853743 1 0 0 0 0 0 0 0 0 0 0 1 0 0; 0.317003006996352 1 0 0 0 0 0 0 0 0 0 0 0 1 0; -0.134402243485427 1 0 0 0 0 0 0 0 0 0 0 0 0 1; 0.0479152148010873 1 0 0 0 0 0 0 0 0 0 0 0 0 0; -0.142420903705832 0 1 0 0 0 0 0 1 0 0 0 0 0 0; 0.144398806337843 0 1 0 0 0 0 0 0 0 1 0 0 0 0; -0.000423516620252108 0 1 0 0 0 0 0 0 0 0 1 0 0 0; -0.147634142603894 0 1 0 0 0 0 0 0 0 0 0 1 0 0; 0.0014740606068652 0 1 0 0 0 0 0 0 0 0 0 0 1 0; 0.152036153038626 0 1 0 0 0 0 0 0 0 0 0 0 0 1; -4.25882464391501e-05 0 1 0 0 0 0 0 0 0 0 0 0 0 0; 0.0273865492488036 0 0 1 0 0 0 0 1 0 0 0 0 0 0; -0.144398806337843 0 0 1 0 0 0 0 0 1 0 0 0 0 0; 0.236289721082867 0 0 1 0 0 0 0 0 0 0 1 0 0 0; -0.0450840112752414 0 0 1 0 0 0 0 0 0 0 0 1 0 0; -0.321688909823353 0 0 1 0 0 0 0 0 0 0 0 0 1 0; 0.107033290331583 0 0 1 0 0 0 0 0 0 0 0 0 0 1; -0.0485724588005584 0 0 1 0 0 0 0 0 0 0 0 0 0 0; 0.232972818871747 0 0 0 1 0 0 0 1 0 0 0 0 0 0; 0.000423516620252108 0 0 0 1 0 0 0 0 1 0 0 0 0 0; -0.236289721082867 0 0 0 1 0 0 0 0 0 1 0 0 0 0; 0.241716154665455 0 0 0 1 0 0 0 0 0 0 0 1 0 0; -0.00146860472861732 0 0 0 1 0 0 0 0 0 0 0 0 1 0; -0.249101162853298 0 0 0 1 0 0 0 0 0 0 0 0 0 1; 0.000212151396801636 0 0 0 1 0 0 0 0 0 0 0 0 0 0; -0.0724666332853743 0 0 0 0 1 0 0 1 0 0 0 0 0 0; 0.147634142603894 0 0 0 0 1 0 0 0 1 0 0 0 0 0; 0.0450840112752414 0 0 0 0 1 0 0 0 0 1 0 0 0 0; -0.241716154665455 0 0 0 0 1 0 0 0 0 0 1 0 0 0; 0.329356759644547 0 0 0 0 1 0 0 0 0 0 0 0 1 0; -0.0619628973206706 0 0 0 0 1 0 0 0 0 0 0 0 0 1; 0.0496474551418917 0 0 0 0 1 0 0 0 0 0 0 0 0 0; -0.317003006996352 0 0 0 0 0 1 0 1 0 0 0 0 0 0; -0.0014740606068652 0 0 0 0 0 1 0 0 1 0 0 0 0 0; 0.321688909823353 0 0 0 0 0 1 0 0 0 1 0 0 0 0; 0.00146860472861732 0 0 0 0 0 1 0 0 0 0 1 0 0 0; -0.329356759644547 0 0 0 0 0 1 0 0 0 0 0 1 0 0; 0.339795868996559 0 0 0 0 0 1 0 0 0 0 0 0 0 1; -0.000590717588518114 0 0 0 0 0 1 0 0 0 0 0 0 0 0; 0.134402243485427 0 0 0 0 0 0 1 1 0 0 0 0 0 0; -0.152036153038626 0 0 0 0 0 0 1 0 1 0 0 0 0 0; -0.107033290331583 0 0 0 0 0 0 1 0 0 1 0 0 0 0; 0.249101162853298 0 0 0 0 0 0 1 0 0 0 1 0 0 0; 0.0619628973206706 0 0 0 0 0 0 1 0 0 0 0 1 0 0; -0.339795868996559 0 0 0 0 0 0 1 0 0 0 0 0 1 0; -0.0511099198580016 0 0 0 0 0 0 1 0 0 0 0 0 0 0; 0.00918550243895438 0 0 0 0 0 0 0 1 0 0 0 0 0 0; 0.14175545754757 0 0 0 0 0 0 0 0 1 0 0 0 0 0; -0.0365716567700609 0 0 0 0 0 0 0 0 0 1 0 0 0 0; -0.2318569638525 0 0 0 0 0 0 0 0 0 0 1 0 0 0; 0.0816497735422322 0 0 0 0 0 0 0 0 0 0 0 1 0 0; 0.315426774753603 0 0 0 0 0 0 0 0 0 0 0 0 1 0; -0.143579906061253 0 0 0 0 0 0 0 0 0 0 0 0 0 1; 0.0476940829282326 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    eqs{6} = [0.117849377768783 1 0 0 0 0 0 0 0 1 0 0 0 0 0; -0.0876860647276761 1 0 0 0 0 0 0 0 0 1 0 0 0 0; -0.140609191366246 1 0 0 0 0 0 0 0 0 0 1 0 0 0; 0.205557289895931 1 0 0 0 0 0 0 0 0 0 0 1 0 0; 0.0938323198896144 1 0 0 0 0 0 0 0 0 0 0 0 1 0; -0.31768738271886 1 0 0 0 0 0 0 0 0 0 0 0 0 1; 0.045293582361765 1 0 0 0 0 0 0 0 0 0 0 0 0 0; -0.117849377768783 0 1 0 0 0 0 0 1 0 0 0 0 0 0; 0.142939316698851 0 1 0 0 0 0 0 0 0 1 0 0 0 0; -0.0215225245442627 0 1 0 0 0 0 0 0 0 0 1 0 0 0; -0.172445006761772 0 1 0 0 0 0 0 0 0 0 0 1 0 0; 0.0687288281380629 0 1 0 0 0 0 0 0 0 0 0 0 1 0; 0.191049331501888 0 1 0 0 0 0 0 0 0 0 0 0 0 1; -0.00235329316165676 0 1 0 0 0 0 0 0 0 0 0 0 0 0; 0.0876860647276761 0 0 1 0 0 0 0 1 0 0 0 0 0 0; -0.142939316698851 0 0 1 0 0 0 0 0 1 0 0 0 0 0; 0.186558534563466 0 0 1 0 0 0 0 0 0 0 1 0 0 0; -0.121012047796254 0 0 1 0 0 0 0 0 0 0 0 1 0 0; -0.164946888393417 0 0 1 0 0 0 0 0 0 0 0 0 1 0; 0.243171868227018 0 0 1 0 0 0 0 0 0 0 0 0 0 1; -0.05318553916708 0 0 1 0 0 0 0 0 0 0 0 0 0 0; 0.140609191366246 0 0 0 1 0 0 0 1 0 0 0 0 0 0; 0.0215225245442627 0 0 0 1 0 0 0 0 1 0 0 0 0 0; -0.186558534563466 0 0 0 1 0 0 0 0 0 1 0 0 0 0; 0.243289063682708 0 0 0 1 0 0 0 0 0 0 0 1 0 0; -0.0648658201248332 0 0 0 1 0 0 0 0 0 0 0 0 1 0; -0.285964399164089 0 0 0 1 0 0 0 0 0 0 0 0 0 1; 0.0110796247829899 0 0 0 1 0 0 0 0 0 0 0 0 0 0; -0.205557289895931 0 0 0 0 1 0 0 1 0 0 0 0 0 0; 0.172445006761772 0 0 0 0 1 0 0 0 1 0 0 0 0 0; 0.121012047796254 0 0 0 0 1 0 0 0 0 1 0 0 0 0; -0.243289063682708 0 0 0 0 1 0 0 0 0 0 1 0 0 0; 0.257181049755608 0 0 0 0 1 0 0 0 0 0 0 0 1 0; -0.131625811988251 0 0 0 0 1 0 0 0 0 0 0 0 0 1; 0.0621718645504763 0 0 0 0 1 0 0 0 0 0 0 0 0 0; -0.0938323198896144 0 0 0 0 0 1 0 1 0 0 0 0 0 0; -0.0687288281380629 0 0 0 0 0 1 0 0 1 0 0 0 0 0; 0.164946888393417 0 0 0 0 0 1 0 0 0 1 0 0 0 0; 0.0648658201248332 0 0 0 0 0 1 0 0 0 0 1 0 0 0; -0.257181049755608 0 0 0 0 0 1 0 0 0 0 0 1 0 0; 0.337387301227059 0 0 0 0 0 1 0 0 0 0 0 0 0 1; -0.0282885650968677 0 0 0 0 0 1 0 0 0 0 0 0 0 0; 0.31768738271886 0 0 0 0 0 0 1 1 0 0 0 0 0 0; -0.191049331501888 0 0 0 0 0 0 1 0 1 0 0 0 0 0; -0.243171868227018 0 0 0 0 0 0 1 0 0 1 0 0 0 0; 0.285964399164089 0 0 0 0 0 0 1 0 0 0 1 0 0 0; 0.131625811988251 0 0 0 0 0 0 1 0 0 0 0 1 0 0; -0.337387301227059 0 0 0 0 0 0 1 0 0 0 0 0 1 0; -0.0670830617515417 0 0 0 0 0 0 1 0 0 0 0 0 0 0; 0.0322519389538854 0 0 0 0 0 0 0 1 0 0 0 0 0 0; 0.108159456398762 0 0 0 0 0 0 0 0 1 0 0 0 0 0; -0.119594583168914 0 0 0 0 0 0 0 0 0 1 0 0 0 0; -0.123157803884071 0 0 0 0 0 0 0 0 0 0 1 0 0 0; 0.235848937781088 0 0 0 0 0 0 0 0 0 0 0 1 0 0; 0.0673080748714457 0 0 0 0 0 0 0 0 0 0 0 0 1 0; -0.343850826908064 0 0 0 0 0 0 0 0 0 0 0 0 0 1; 0.0422134389521502 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    eqs{7} = [0.111092461406876 1 0 0 0 0 0 0 0 1 0 0 0 0 0; -0.134751367758009 1 0 0 0 0 0 0 0 0 1 0 0 0 0; -0.0358679150221221 1 0 0 0 0 0 0 0 0 0 1 0 0 0; 0.229123900385785 1 0 0 0 0 0 0 0 0 0 0 1 0 0; -0.149394425689347 1 0 0 0 0 0 0 0 0 0 0 0 1 0; -0.178244285854098 1 0 0 0 0 0 0 0 0 0 0 0 0 1; 0.0558946565710762 1 0 0 0 0 0 0 0 0 0 0 0 0 0; -0.111092461406876 0 1 0 0 0 0 0 1 0 0 0 0 0 0; 0.177562334692695 0 1 0 0 0 0 0 0 0 1 0 0 0 0; -0.102374113386776 0 1 0 0 0 0 0 0 0 0 1 0 0 0; -0.18648403787105 0 1 0 0 0 0 0 0 0 0 0 1 0 0; 0.262306836773212 0 1 0 0 0 0 0 0 0 0 0 0 1 0; 0.0649687664514805 0 1 0 0 0 0 0 0 0 0 0 0 0 1; -0.0136069730368173 0 1 0 0 0 0 0 0 0 0 0 0 0 0; 0.134751367758009 0 0 1 0 0 0 0 1 0 0 0 0 0 0; -0.177562334692695 0 0 1 0 0 0 0 0 1 0 0 0 0 0; 0.181505047943072 0 0 1 0 0 0 0 0 0 0 1 0 0 0; -0.14001666108798 0 0 1 0 0 0 0 0 0 0 0 1 0 0; -0.0793877631326525 0 0 1 0 0 0 0 0 0 0 0 0 1 0; 0.206088163958722 0 0 1 0 0 0 0 0 0 0 0 0 0 1; -0.0728332722795214 0 0 1 0 0 0 0 0 0 0 0 0 0 0; 0.0358679150221221 0 0 0 1 0 0 0 1 0 0 0 0 0 0; 0.102374113386776 0 0 0 1 0 0 0 0 1 0 0 0 0 0; -0.181505047943072 0 0 0 1 0 0 0 0 0 1 0 0 0 0; 0.271351893722558 0 0 0 1 0 0 0 0 0 0 0 1 0 0; -0.222360013390142 0 0 0 1 0 0 0 0 0 0 0 0 1 0; -0.185232145045196 0 0 0 1 0 0 0 0 0 0 0 0 0 1; 0.055901359853485 0 0 0 1 0 0 0 0 0 0 0 0 0 0; -0.229123900385785 0 0 0 0 1 0 0 1 0 0 0 0 0 0; 0.18648403787105 0 0 0 0 1 0 0 0 1 0 0 0 0 0; 0.14001666108798 0 0 0 0 1 0 0 0 0 1 0 0 0 0; -0.271351893722558 0 0 0 0 1 0 0 0 0 0 1 0 0 0; 0.290218520618396 0 0 0 0 1 0 0 0 0 0 0 0 1 0; -0.165212083236475 0 0 0 0 1 0 0 0 0 0 0 0 0 1; 0.0657630448153702 0 0 0 0 1 0 0 0 0 0 0 0 0 0; 0.149394425689347 0 0 0 0 0 1 0 1 0 0 0 0 0 0; -0.262306836773212 0 0 0 0 0 1 0 0 1 0 0 0 0 0; 0.0793877631326525 0 0 0 0 0 1 0 0 0 1 0 0 0 0; 0.222360013390142 0 0 0 0 0 1 0 0 0 0 1 0 0 0; -0.290218520618396 0 0 0 0 0 1 0 0 0 0 0 1 0 0; 0.333494485353364 0 0 0 0 0 1 0 0 0 0 0 0 0 1; -0.113677782232455 0 0 0 0 0 1 0 0 0 0 0 0 0 0; 0.178244285854098 0 0 0 0 0 0 1 1 0 0 0 0 0 0; -0.0649687664514805 0 0 0 0 0 0 1 0 1 0 0 0 0 0; -0.206088163958722 0 0 0 0 0 0 1 0 0 1 0 0 0 0; 0.185232145045196 0 0 0 0 0 0 1 0 0 0 1 0 0 0; 0.165212083236475 0 0 0 0 0 0 1 0 0 0 0 1 0 0; -0.333494485353364 0 0 0 0 0 0 1 0 0 0 0 0 1 0; -0.0108561974574617 0 0 0 0 0 0 1 0 0 0 0 0 0 0; 0.0610375024133259 0 0 0 0 0 0 0 1 0 0 0 0 0 0; 0.0759824241868603 0 0 0 0 0 0 0 0 1 0 0 0 0 0; -0.189722117506034 0 0 0 0 0 0 0 0 0 1 0 0 0 0; 0.0317152848571259 0 0 0 0 0 0 0 0 0 0 1 0 0 0; 0.259170684828122 0 0 0 0 0 0 0 0 0 0 0 1 0 0; -0.246298483803003 0 0 0 0 0 0 0 0 0 0 0 0 1 0; -0.157607131519002 0 0 0 0 0 0 0 0 0 0 0 0 0 1; 0.0457055959570352 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    eqs{8} = [0.0796589673561539 1 0 0 0 0 0 0 0 1 0 0 0 0 0; -0.118214547454375 1 0 0 0 0 0 0 0 0 1 0 0 0 0; 0.0502196723763858 1 0 0 0 0 0 0 0 0 0 1 0 0 0; 0.0962288305776844 1 0 0 0 0 0 0 0 0 0 0 1 0 0; -0.183826215464085 1 0 0 0 0 0 0 0 0 0 0 0 1 0; 0.0930518287164687 1 0 0 0 0 0 0 0 0 0 0 0 0 1; 0.0565305559379184 1 0 0 0 0 0 0 0 0 0 0 0 0 0; -0.0796589673561539 0 1 0 0 0 0 0 1 0 0 0 0 0 0; 0.156892786171454 0 1 0 0 0 0 0 0 0 1 0 0 0 0; -0.180960821935206 0 1 0 0 0 0 0 0 0 0 1 0 0 0; -0.0316791495604159 0 1 0 0 0 0 0 0 0 0 0 1 0 0; 0.287656203064056 0 1 0 0 0 0 0 0 0 0 0 0 1 0; -0.268191744107457 0 1 0 0 0 0 0 0 0 0 0 0 0 1; -0.0336853681856347 0 1 0 0 0 0 0 0 0 0 0 0 0 0; 0.118214547454375 0 0 1 0 0 0 0 1 0 0 0 0 0 0; -0.156892786171454 0 0 1 0 0 0 0 0 1 0 0 0 0 0; 0.169636863253279 0 0 1 0 0 0 0 0 0 0 1 0 0 0; -0.142515944998032 0 0 1 0 0 0 0 0 0 0 0 1 0 0; -0.0648281156979872 0 0 1 0 0 0 0 0 0 0 0 0 1 0; 0.214727927849111 0 0 1 0 0 0 0 0 0 0 0 0 0 1; -0.061350730887235 0 0 1 0 0 0 0 0 0 0 0 0 0 0; -0.0502196723763858 0 0 0 1 0 0 0 1 0 0 0 0 0 0; 0.180960821935206 0 0 0 1 0 0 0 0 1 0 0 0 0 0; -0.169636863253279 0 0 0 1 0 0 0 0 0 1 0 0 0 0; 0.198630892268182 0 0 0 1 0 0 0 0 0 0 0 1 0 0; -0.236248891911516 0 0 0 1 0 0 0 0 0 0 0 0 1 0; 0.0423082798564274 0 0 0 1 0 0 0 0 0 0 0 0 0 1; 0.107183760927322 0 0 0 1 0 0 0 0 0 0 0 0 0 0; -0.0962288305776844 0 0 0 0 1 0 0 1 0 0 0 0 0 0; 0.0316791495604159 0 0 0 0 1 0 0 0 1 0 0 0 0 0; 0.142515944998032 0 0 0 0 1 0 0 0 0 1 0 0 0 0; -0.198630892268182 0 0 0 0 1 0 0 0 0 0 1 0 0 0; 0.274386708513746 0 0 0 0 1 0 0 0 0 0 0 0 1 0; -0.286973028473548 0 0 0 0 1 0 0 0 0 0 0 0 0 1; -0.0182109271553282 0 0 0 0 1 0 0 0 0 0 0 0 0 0; 0.183826215464085 0 0 0 0 0 1 0 1 0 0 0 0 0 0; -0.287656203064056 0 0 0 0 0 1 0 0 1 0 0 0 0 0; 0.0648281156979872 0 0 0 0 0 1 0 0 0 1 0 0 0 0; 0.236248891911516 0 0 0 0 0 1 0 0 0 0 1 0 0 0; -0.274386708513746 0 0 0 0 0 1 0 0 0 0 0 1 0 0; 0.282877601218352 0 0 0 0 0 1 0 0 0 0 0 0 0 1; -0.126402734837165 0 0 0 0 0 1 0 0 0 0 0 0 0 0; -0.0930518287164687 0 0 0 0 0 0 1 1 0 0 0 0 0 0; 0.268191744107457 0 0 0 0 0 0 1 0 1 0 0 0 0 0; -0.214727927849111 0 0 0 0 0 0 1 0 0 1 0 0 0 0; -0.0423082798564274 0 0 0 0 0 0 1 0 0 0 1 0 0 0; 0.286973028473548 0 0 0 0 0 0 1 0 0 0 0 1 0 0; -0.282877601218352 0 0 0 0 0 0 1 0 0 0 0 0 1 0; 0.150975385205826 0 0 0 0 0 0 1 0 0 0 0 0 0 0; 0.0741114933988293 0 0 0 0 0 0 0 1 0 0 0 0 0 0; 0.0106049119596627 0 0 0 0 0 0 0 0 1 0 0 0 0 0; -0.161704500852779 0 0 0 0 0 0 0 0 0 1 0 0 0 0; 0.175044347513942 0 0 0 0 0 0 0 0 0 0 1 0 0 0; 0.042283844136786 0 0 0 0 0 0 0 0 0 0 0 1 0 0; -0.292096325092408 0 0 0 0 0 0 0 0 0 0 0 0 1 0; 0.261902681100377 0 0 0 0 0 0 0 0 0 0 0 0 0 1; 0.0388653608430193 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    eqs{9} = [0.0686720365790188 1 0 0 0 0 0 0 0 1 0 0 0 0 0; -0.100526456833352 1 0 0 0 0 0 0 0 0 1 0 0 0 0; 0.0710675860144179 1 0 0 0 0 0 0 0 0 0 1 0 0 0; 0.00909375336074949 1 0 0 0 0 0 0 0 0 0 0 1 0 0; -0.0933332792890137 1 0 0 0 0 0 0 0 0 0 0 0 1 0; 0.123556695608219 1 0 0 0 0 0 0 0 0 0 0 0 0 1; 0.0647222554402329 1 0 0 0 0 0 0 0 0 0 0 0 0 0; -0.0686720365790188 0 1 0 0 0 0 0 1 0 0 0 0 0 0; 0.12896006925326 0 1 0 0 0 0 0 0 0 1 0 0 0 0; -0.219076103190301 0 1 0 0 0 0 0 0 0 0 1 0 0 0; 0.155671333861986 0 1 0 0 0 0 0 0 0 0 0 1 0 0; 0.0691398784252439 0 1 0 0 0 0 0 0 0 0 0 0 1 0; -0.307748264710023 0 1 0 0 0 0 0 0 0 0 0 0 0 1; -0.0696638882849485 0 1 0 0 0 0 0 0 0 0 0 0 0 0; 0.100526456833352 0 0 1 0 0 0 0 1 0 0 0 0 0 0; -0.12896006925326 0 0 1 0 0 0 0 0 1 0 0 0 0 0; 0.187238711082961 0 0 1 0 0 0 0 0 0 0 1 0 0 0; -0.244958785625722 0 0 1 0 0 0 0 0 0 0 0 1 0 0; 0.0740604095950238 0 0 1 0 0 0 0 0 0 0 0 0 1 0; 0.218472661842609 0 0 1 0 0 0 0 0 0 0 0 0 0 1; -0.0195643343667746 0 0 1 0 0 0 0 0 0 0 0 0 0 0; -0.0710675860144179 0 0 0 1 0 0 0 1 0 0 0 0 0 0; 0.219076103190301 0 0 0 1 0 0 0 0 1 0 0 0 0 0; -0.187238711082961 0 0 0 1 0 0 0 0 0 1 0 0 0 0; 0.190112462207817 0 0 0 1 0 0 0 0 0 0 0 1 0 0; -0.226198138884709 0 0 0 1 0 0 0 0 0 0 0 0 1 0; 0.0756842724166312 0 0 0 1 0 0 0 0 0 0 0 0 0 1; 0.134381556139309 0 0 0 1 0 0 0 0 0 0 0 0 0 0; -0.00909375336074949 0 0 0 0 1 0 0 1 0 0 0 0 0 0; -0.155671333861986 0 0 0 0 1 0 0 0 1 0 0 0 0 0; 0.244958785625722 0 0 0 0 1 0 0 0 0 1 0 0 0 0; -0.190112462207817 0 0 0 0 1 0 0 0 0 0 1 0 0 0; 0.220731142362197 0 0 0 0 1 0 0 0 0 0 0 0 1 0; -0.320841255436336 0 0 0 0 1 0 0 0 0 0 0 0 0 1; -0.155942747392222 0 0 0 0 1 0 0 0 0 0 0 0 0 0; 0.0933332792890137 0 0 0 0 0 1 0 1 0 0 0 0 0 0; -0.0691398784252439 0 0 0 0 0 1 0 0 1 0 0 0 0 0; -0.0740604095950238 0 0 0 0 0 1 0 0 0 1 0 0 0 0; 0.226198138884709 0 0 0 0 0 1 0 0 0 0 1 0 0 0; -0.220731142362197 0 0 0 0 0 1 0 0 0 0 0 1 0 0; 0.29386721048664 0 0 0 0 0 1 0 0 0 0 0 0 0 1; 0.0295181324174065 0 0 0 0 0 1 0 0 0 0 0 0 0 0; -0.123556695608219 0 0 0 0 0 0 1 1 0 0 0 0 0 0; 0.307748264710023 0 0 0 0 0 0 1 0 1 0 0 0 0 0; -0.218472661842609 0 0 0 0 0 0 1 0 0 1 0 0 0 0; -0.0756842724166312 0 0 0 0 0 0 1 0 0 0 1 0 0 0; 0.320841255436336 0 0 0 0 0 0 1 0 0 0 0 1 0 0; -0.29386721048664 0 0 0 0 0 0 1 0 0 0 0 0 1 0; 0.164706371379088 0 0 0 0 0 0 1 0 0 0 0 0 0 0; 0.0907376891101431 0 0 0 0 0 0 0 1 0 0 0 0 0 0; -0.0588184998584 0 0 0 0 0 0 0 0 1 0 0 0 0 0; -0.0842952033007785 0 0 0 0 0 0 0 0 0 1 0 0 0 0; 0.228599169718138 0 0 0 0 0 0 0 0 0 0 1 0 0 0; -0.213480460991175 0 0 0 0 0 0 0 0 0 0 0 1 0 0; -0.0114146799493276 0 0 0 0 0 0 0 0 0 0 0 0 1 0; 0.300805799760116 0 0 0 0 0 0 0 0 0 0 0 0 0 1; 0.036612781420333 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    eqs{10} = [0.0634507161126456 1 0 0 0 0 0 0 0 1 0 0 0 0 0; -0.0919620308733866 1 0 0 0 0 0 0 0 0 1 0 0 0 0; 0.0679575969632192 1 0 0 0 0 0 0 0 0 0 1 0 0 0; -0.00300443501367983 1 0 0 0 0 0 0 0 0 0 0 1 0 0; -0.0677597788838842 1 0 0 0 0 0 0 0 0 0 0 0 1 0; 0.10435358295116 1 0 0 0 0 0 0 0 0 0 0 0 0 1; 0.0619361353947554 1 0 0 0 0 0 0 0 0 0 0 0 0 0; -0.0634507161126456 0 1 0 0 0 0 0 1 0 0 0 0 0 0; 0.078611561030047 0 1 0 0 0 0 0 0 0 1 0 0 0 0; -0.147782392940417 0 1 0 0 0 0 0 0 0 0 1 0 0 0; 0.174135484214182 0 1 0 0 0 0 0 0 0 0 0 1 0 0; -0.131620014695342 0 1 0 0 0 0 0 0 0 0 0 0 1 0; 0.0203173052336943 0 1 0 0 0 0 0 0 0 0 0 0 0 1; -0.101545834581481 0 1 0 0 0 0 0 0 0 0 0 0 0 0; 0.0919620308733866 0 0 1 0 0 0 0 1 0 0 0 0 0 0; -0.078611561030047 0 0 1 0 0 0 0 0 1 0 0 0 0 0; 0.129992484030551 0 0 1 0 0 0 0 0 0 0 1 0 0 0; -0.248660226639642 0 0 1 0 0 0 0 0 0 0 0 1 0 0; 0.274713146139065 0 0 1 0 0 0 0 0 0 0 0 0 1 0; -0.158734515906069 0 0 1 0 0 0 0 0 0 0 0 0 0 1; 0.0704399439615642 0 0 1 0 0 0 0 0 0 0 0 0 0 0; -0.0679575969632192 0 0 0 1 0 0 0 1 0 0 0 0 0 0; 0.147782392940417 0 0 0 1 0 0 0 0 1 0 0 0 0 0; -0.129992484030551 0 0 0 1 0 0 0 0 0 1 0 0 0 0; 0.179506665255695 0 0 0 1 0 0 0 0 0 0 0 1 0 0; -0.2987875211034 0 0 0 1 0 0 0 0 0 0 0 0 1 0; 0.264809264107425 0 0 0 1 0 0 0 0 0 0 0 0 0 1; 0.0354962014035918 0 0 0 1 0 0 0 0 0 0 0 0 0 0; 0.00300443501367983 0 0 0 0 1 0 0 1 0 0 0 0 0 0; -0.174135484214182 0 0 0 0 1 0 0 0 1 0 0 0 0 0; 0.248660226639642 0 0 0 0 1 0 0 0 0 1 0 0 0 0; -0.179506665255695 0 0 0 0 1 0 0 0 0 0 1 0 0 0; 0.192193665161998 0 0 0 0 1 0 0 0 0 0 0 0 1 0; -0.287352213449315 0 0 0 0 1 0 0 0 0 0 0 0 0 1; -0.165170571878262 0 0 0 0 1 0 0 0 0 0 0 0 0 0; 0.0677597788838842 0 0 0 0 0 1 0 1 0 0 0 0 0 0; 0.131620014695342 0 0 0 0 0 1 0 0 1 0 0 0 0 0; -0.274713146139065 0 0 0 0 0 1 0 0 0 1 0 0 0 0; 0.2987875211034 0 0 0 0 0 1 0 0 0 0 1 0 0 0; -0.192193665161998 0 0 0 0 0 1 0 0 0 0 0 1 0 0; 0.194770441825301 0 0 0 0 0 1 0 0 0 0 0 0 0 1; 0.236920231474818 0 0 0 0 0 1 0 0 0 0 0 0 0 0; -0.10435358295116 0 0 0 0 0 0 1 1 0 0 0 0 0 0; -0.0203173052336943 0 0 0 0 0 0 1 0 1 0 0 0 0 0; 0.158734515906069 0 0 0 0 0 0 1 0 0 1 0 0 0 0; -0.264809264107425 0 0 0 0 0 0 1 0 0 0 1 0 0 0; 0.287352213449315 0 0 0 0 0 0 1 0 0 0 0 1 0 0; -0.194770441825301 0 0 0 0 0 0 1 0 0 0 0 0 1 0; -0.18683866418635 0 0 0 0 0 0 1 0 0 0 0 0 0 0; 0.0873429611491346 0 0 0 0 0 0 0 1 0 0 0 0 0 0; -0.113861083558443 0 0 0 0 0 0 0 0 1 0 0 0 0 0; 0.0568114937293848 0 0 0 0 0 0 0 0 0 1 0 0 0 0; 0.0814809744571268 0 0 0 0 0 0 0 0 0 0 1 0 0 0; -0.234314465101891 0 0 0 0 0 0 0 0 0 0 0 1 0 0; 0.302774891323165 0 0 0 0 0 0 0 0 0 0 0 0 1 0; -0.21522823486282 0 0 0 0 0 0 0 0 0 0 0 0 0 1; 0.0286395254233433 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    eqs{11} = [0.0826005505766948 1 0 0 0 0 0 0 0 1 0 0 0 0 0; -0.157249261206188 1 0 0 0 0 0 0 0 0 1 0 0 0 0; 0.196518587756542 1 0 0 0 0 0 0 0 0 0 1 0 0 0; -0.180741998867581 1 0 0 0 0 0 0 0 0 0 0 1 0 0; 0.104055420899098 1 0 0 0 0 0 0 0 0 0 0 0 1 0; 0.02281518442848 1 0 0 0 0 0 0 0 0 0 0 0 0 1; 0.063076046883314 1 0 0 0 0 0 0 0 0 0 0 0 0 0; -0.0826005505766948 0 1 0 0 0 0 0 1 0 0 0 0 0 0; 0.0856326835924026 0 1 0 0 0 0 0 0 0 1 0 0 0 0; -0.152840687407958 0 1 0 0 0 0 0 0 0 0 1 0 0 0; 0.184074701375133 0 1 0 0 0 0 0 0 0 0 0 1 0 0; -0.16992986298103 0 1 0 0 0 0 0 0 0 0 0 0 1 0; 0.111588008725639 0 1 0 0 0 0 0 0 0 0 0 0 0 1; -0.141860837990795 0 1 0 0 0 0 0 0 0 0 0 0 0 0; 0.157249261206188 0 0 1 0 0 0 0 1 0 0 0 0 0 0; -0.0856326835924026 0 0 1 0 0 0 0 0 1 0 0 0 0 0; 0.0872351465148534 0 0 1 0 0 0 0 0 0 0 1 0 0 0; -0.163052041458224 0 0 1 0 0 0 0 0 0 0 0 1 0 0; 0.215625687145652 0 0 1 0 0 0 0 0 0 0 0 0 1 0; -0.236086288344242 0 0 1 0 0 0 0 0 0 0 0 0 0 1; 0.20467346385818 0 0 1 0 0 0 0 0 0 0 0 0 0 0; -0.196518587756542 0 0 0 1 0 0 0 1 0 0 0 0 0 0; 0.152840687407958 0 0 0 1 0 0 0 0 1 0 0 0 0 0; -0.0872351465148534 0 0 0 1 0 0 0 0 0 1 0 0 0 0; 0.10350256682128 0 0 0 1 0 0 0 0 0 0 0 1 0 0; -0.211747676133372 0 0 0 1 0 0 0 0 0 0 0 0 1 0; 0.30770020513503 0 0 0 1 0 0 0 0 0 0 0 0 0 1; -0.220793990451125 0 0 0 1 0 0 0 0 0 0 0 0 0 0; 0.180741998867581 0 0 0 0 1 0 0 1 0 0 0 0 0 0; -0.184074701375133 0 0 0 0 1 0 0 0 1 0 0 0 0 0; 0.163052041458224 0 0 0 0 1 0 0 0 0 1 0 0 0 0; -0.10350256682128 0 0 0 0 1 0 0 0 0 0 1 0 0 0; 0.139944497867327 0 0 0 0 1 0 0 0 0 0 0 0 1 0; -0.295014232194373 0 0 0 0 1 0 0 0 0 0 0 0 0 1; 0.169847619992437 0 0 0 0 1 0 0 0 0 0 0 0 0 0; -0.104055420899098 0 0 0 0 0 1 0 1 0 0 0 0 0 0; 0.16992986298103 0 0 0 0 0 1 0 0 1 0 0 0 0 0; -0.215625687145652 0 0 0 0 0 1 0 0 0 1 0 0 0 0; 0.211747676133372 0 0 0 0 0 1 0 0 0 0 1 0 0 0; -0.139944497867327 0 0 0 0 0 1 0 0 0 0 0 1 0 0; 0.187508657883292 0 0 0 0 0 1 0 0 0 0 0 0 0 1; -0.0489450151813522 0 0 0 0 0 1 0 0 0 0 0 0 0 0; -0.02281518442848 0 0 0 0 0 0 1 1 0 0 0 0 0 0; -0.111588008725639 0 0 0 0 0 0 1 0 1 0 0 0 0 0; 0.236086288344242 0 0 0 0 0 0 1 0 0 1 0 0 0 0; -0.30770020513503 0 0 0 0 0 0 1 0 0 0 1 0 0 0; 0.295014232194373 0 0 0 0 0 0 1 0 0 0 0 1 0 0; -0.187508657883292 0 0 0 0 0 0 1 0 0 0 0 0 1 0; -0.124395195675994 0 0 0 0 0 0 1 0 0 0 0 0 0 0; 0.0865478409889395 0 0 0 0 0 0 0 1 0 0 0 0 0 0; -0.1633512279268 0 0 0 0 0 0 0 0 1 0 0 0 0 0; 0.221251988007423 0 0 0 0 0 0 0 0 0 1 0 0 0 0; -0.228491468619491 0 0 0 0 0 0 0 0 0 0 1 0 0 0; 0.164564998336673 0 0 0 0 0 0 0 0 0 0 0 1 0 0; -0.0277300574846218 0 0 0 0 0 0 0 0 0 0 0 0 1 0; -0.162039956559369 0 0 0 0 0 0 0 0 0 0 0 0 0 1; 0.0239005614858983 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    eqs{12} = [0.0550615500574426 1 0 0 0 0 0 0 0 1 0 0 0 0 0; 0.117680618616067 1 0 0 0 0 0 0 0 0 1 0 0 0 0; 0.175869560187308 1 0 0 0 0 0 0 0 0 0 1 0 0 0; 0.21740269285875 1 0 0 0 0 0 0 0 0 0 0 1 0 0; 0.231631285949858 1 0 0 0 0 0 0 0 0 0 0 0 1 0; 0.211149287455176 1 0 0 0 0 0 0 0 0 0 0 0 0 1; 0.0374824859744506 1 0 0 0 0 0 0 0 0 0 0 0 0 0; -0.0550615500574426 0 1 0 0 0 0 0 1 0 0 0 0 0 0; 0.0644689106284467 0 1 0 0 0 0 0 0 0 1 0 0 0 0; 0.127778822569236 0 1 0 0 0 0 0 0 0 0 1 0 0 0; 0.17898877502978 0 1 0 0 0 0 0 0 0 0 0 1 0 0; 0.208303279847868 0 1 0 0 0 0 0 0 0 0 0 0 1 0; 0.208483988744256 0 1 0 0 0 0 0 0 0 0 0 0 0 1; 0.0919611633862619 0 1 0 0 0 0 0 0 0 0 0 0 0 0; -0.117680618616067 0 0 1 0 0 0 0 1 0 0 0 0 0 0; -0.0644689106284467 0 0 1 0 0 0 0 0 1 0 0 0 0 0; 0.0671788557376545 0 0 1 0 0 0 0 0 0 0 1 0 0 0; 0.127998485098761 0 0 1 0 0 0 0 0 0 0 0 1 0 0; 0.173991508586389 0 0 1 0 0 0 0 0 0 0 0 0 1 0; 0.198359112908501 0 0 1 0 0 0 0 0 0 0 0 0 0 1; 0.152658098959448 0 0 1 0 0 0 0 0 0 0 0 0 0 0; -0.175869560187308 0 0 0 1 0 0 0 1 0 0 0 0 0 0; -0.127778822569236 0 0 0 1 0 0 0 0 1 0 0 0 0 0; -0.0671788557376545 0 0 0 1 0 0 0 0 0 1 0 0 0 0; 0.0671833070849081 0 0 0 1 0 0 0 0 0 0 0 1 0 0; 0.127795770664042 0 0 0 1 0 0 0 0 0 0 0 0 1 0; 0.175904602378396 0 0 0 1 0 0 0 0 0 0 0 0 0 1; 0.206745023022343 0 0 0 1 0 0 0 0 0 0 0 0 0 0; -0.21740269285875 0 0 0 0 1 0 0 1 0 0 0 0 0 0; -0.17898877502978 0 0 0 0 1 0 0 0 1 0 0 0 0 0; -0.127998485098761 0 0 0 0 1 0 0 0 0 1 0 0 0 0; -0.0671833070849081 0 0 0 0 1 0 0 0 0 0 1 0 0 0; 0.0694912118439611 0 0 0 0 1 0 0 0 0 0 0 0 1 0; 0.136785619966463 0 0 0 0 1 0 0 0 0 0 0 0 0 1; 0.241251114345689 0 0 0 0 1 0 0 0 0 0 0 0 0 0; -0.231631285949858 0 0 0 0 0 1 0 1 0 0 0 0 0 0; -0.208303279847868 0 0 0 0 0 1 0 0 1 0 0 0 0 0; -0.173991508586389 0 0 0 0 0 1 0 0 0 1 0 0 0 0; -0.127795770664042 0 0 0 0 0 1 0 0 0 0 1 0 0 0; -0.0694912118439611 0 0 0 0 0 1 0 0 0 0 0 1 0 0; 0.078245623195126 0 0 0 0 0 1 0 0 0 0 0 0 0 1; 0.245059533434823 0 0 0 0 0 1 0 0 0 0 0 0 0 0; -0.211149287455176 0 0 0 0 0 0 1 1 0 0 0 0 0 0; -0.208483988744256 0 0 0 0 0 0 1 0 1 0 0 0 0 0; -0.198359112908501 0 0 0 0 0 0 1 0 0 1 0 0 0 0; -0.175904602378396 0 0 0 0 0 0 1 0 0 0 1 0 0 0; -0.136785619966463 0 0 0 0 0 0 1 0 0 0 0 1 0 0; -0.078245623195126 0 0 0 0 0 0 1 0 0 0 0 0 1 0; 0.210728465262044 0 0 0 0 0 0 1 0 0 0 0 0 0 0; -0.0507499507280962 0 0 0 0 0 0 0 1 0 0 0 0 0 0; -0.104931527195624 0 0 0 0 0 0 0 0 1 0 0 0 0 0; -0.164844850632663 0 0 0 0 0 0 0 0 0 1 0 0 0 0; -0.217383865432538 0 0 0 0 0 0 0 0 0 0 1 0 0 0; -0.249334155143539 0 0 0 0 0 0 0 0 0 0 0 1 0 0; -0.249430743920153 0 0 0 0 0 0 0 0 0 0 0 0 1 0; -0.210231368904348 0 0 0 0 0 0 0 0 0 0 0 0 0 1; 0.0133292654056755 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    eqs{13} = [0.0613272910033801 1 0 0 0 0 0 0 0 1 0 0 0 0 0; 0.0945410920231372 1 0 0 0 0 0 0 0 0 1 0 0 0 0; 0.0793744778856648 1 0 0 0 0 0 0 0 0 0 1 0 0 0; 0.0187275995292419 1 0 0 0 0 0 0 0 0 0 0 1 0 0; -0.0613347985668647 1 0 0 0 0 0 0 0 0 0 0 0 1 0; -0.122647067909255 1 0 0 0 0 0 0 0 0 0 0 0 0 1; 0.0566441763218766 1 0 0 0 0 0 0 0 0 0 0 0 0 0; -0.0613272910033801 0 1 0 0 0 0 0 1 0 0 0 0 0 0; 0.0660158412945904 0 1 0 0 0 0 0 0 0 1 0 0 0 0; 0.119084807355863 0 1 0 0 0 0 0 0 0 0 1 0 0 0; 0.141829336338518 0 1 0 0 0 0 0 0 0 0 0 1 0 0; 0.122662091747595 0 1 0 0 0 0 0 0 0 0 0 0 1 0; 0.0613197747259972 0 1 0 0 0 0 0 0 0 0 0 0 0 1; 0.103912755419511 0 1 0 0 0 0 0 0 0 0 0 0 0 0; -0.0945410920231372 0 0 1 0 0 0 0 1 0 0 0 0 0 0; -0.0660158412945904 0 0 1 0 0 0 0 0 1 0 0 0 0 0; 0.0981363223025041 0 0 1 0 0 0 0 0 0 0 1 0 0 0; 0.19848230536285 0 0 1 0 0 0 0 0 0 0 0 1 0 0; 0.255117683753933 0 0 1 0 0 0 0 0 0 0 0 0 1 0; 0.226553099089066 0 0 1 0 0 0 0 0 0 0 0 0 0 1; 0.0992154115829214 0 0 1 0 0 0 0 0 0 0 0 0 0 0; -0.0793744778856648 0 0 0 1 0 0 0 1 0 0 0 0 0 0; -0.119084807355863 0 0 0 1 0 0 0 0 1 0 0 0 0 0; -0.0981363223025041 0 0 0 1 0 0 0 0 0 1 0 0 0 0; 0.147201299631271 0 0 0 1 0 0 0 0 0 0 0 1 0 0; 0.277858060931396 0 0 0 1 0 0 0 0 0 0 0 0 1 0; 0.317519773645382 0 0 0 1 0 0 0 0 0 0 0 0 0 1; 0.0245006726594114 0 0 0 1 0 0 0 0 0 0 0 0 0 0; -0.0187275995292419 0 0 0 0 1 0 0 1 0 0 0 0 0 0; -0.141829336338518 0 0 0 0 1 0 0 0 1 0 0 0 0 0; -0.19848230536285 0 0 0 0 1 0 0 0 0 1 0 0 0 0; -0.147201299631271 0 0 0 0 1 0 0 0 0 0 1 0 0 0; 0.179304191118717 0 0 0 0 1 0 0 0 0 0 0 0 1 0; 0.302366599377504 0 0 0 0 1 0 0 0 0 0 0 0 0 1; -0.0992668902554587 0 0 0 0 1 0 0 0 0 0 0 0 0 0; 0.0613347985668647 0 0 0 0 0 1 0 1 0 0 0 0 0 0; -0.122662091747595 0 0 0 0 0 1 0 0 1 0 0 0 0 0; -0.255117683753933 0 0 0 0 0 1 0 0 0 1 0 0 0 0; -0.277858060931396 0 0 0 0 0 1 0 0 0 0 1 0 0 0; -0.179304191118717 0 0 0 0 0 1 0 0 0 0 0 1 0 0; 0.183981873010699 0 0 0 0 0 1 0 0 0 0 0 0 0 1; -0.217220765130656 0 0 0 0 0 1 0 0 0 0 0 0 0 0; 0.122647067909255 0 0 0 0 0 0 1 1 0 0 0 0 0 0; -0.0613197747259972 0 0 0 0 0 0 1 0 1 0 0 0 0 0; -0.226553099089066 0 0 0 0 0 0 1 0 0 1 0 0 0 0; -0.317519773645382 0 0 0 0 0 0 1 0 0 0 1 0 0 0; -0.302366599377504 0 0 0 0 0 0 1 0 0 0 0 1 0 0; -0.183981873010699 0 0 0 0 0 0 1 0 0 0 0 0 1 0; -0.264450012984886 0 0 0 0 0 0 1 0 0 0 0 0 0 0; -0.0793866296075661 0 0 0 0 0 0 0 1 0 0 0 0 0 0; -0.119072652109299 0 0 0 0 0 0 0 0 1 0 0 0 0 0; -0.0981045032477805 0 0 0 0 0 0 0 0 0 1 0 0 0 0; 3.93283605357748e-05 0 0 0 0 0 0 0 0 0 0 1 0 0 0; 0.147233114329272 0 0 0 0 0 0 0 0 0 0 0 1 0 0; 0.277870209128634 0 0 0 0 0 0 0 0 0 0 0 0 1 0; 0.317507614874154 0 0 0 0 0 0 0 0 0 0 0 0 0 1; 0.0245324895357735 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    eqs{14} = [0.0713406844524503 1 0 0 0 0 0 0 0 1 0 0 0 0 0; 0.101894746713401 1 0 0 0 0 0 0 0 0 1 0 0 0 0; 0.0725919228555626 1 0 0 0 0 0 0 0 0 0 1 0 0 0; -0.00113261912308242 1 0 0 0 0 0 0 0 0 0 0 1 0 0; -0.0768863179978937 1 0 0 0 0 0 0 0 0 0 0 0 1 0; -0.109061791346399 1 0 0 0 0 0 0 0 0 0 0 0 0 1; 0.070396263166772 1 0 0 0 0 0 0 0 0 0 0 0 0 0; -0.0713406844524503 0 1 0 0 0 0 0 1 0 0 0 0 0 0; 0.121014518776515 0 1 0 0 0 0 0 0 0 1 0 0 0 0; 0.220554583835156 0 1 0 0 0 0 0 0 0 0 1 0 0 0; 0.20219488165898 0 1 0 0 0 0 0 0 0 0 0 1 0 0; 0.0286236322874945 0 1 0 0 0 0 0 0 0 0 0 0 1 0; -0.220343637084719 0 1 0 0 0 0 0 0 0 0 0 0 0 1; 0.0881961297792858 0 1 0 0 0 0 0 0 0 0 0 0 0 0; -0.101894746713401 0 0 1 0 0 0 0 1 0 0 0 0 0 0; -0.121014518776515 0 0 1 0 0 0 0 0 1 0 0 0 0 0; 0.191877565375177 0 0 1 0 0 0 0 0 0 0 1 0 0 0; 0.290712932889619 0 0 1 0 0 0 0 0 0 0 0 1 0 0; 0.171304195197489 0 0 1 0 0 0 0 0 0 0 0 0 1 0; -0.129712785418193 0 0 1 0 0 0 0 0 0 0 0 0 0 1; 0.00655660087593419 0 0 1 0 0 0 0 0 0 0 0 0 0 0; -0.0725919228555626 0 0 0 1 0 0 0 1 0 0 0 0 0 0; -0.220554583835156 0 0 0 1 0 0 0 0 1 0 0 0 0 0; -0.191877565375177 0 0 0 1 0 0 0 0 0 1 0 0 0 0; 0.209242730218869 0 0 0 1 0 0 0 0 0 0 0 1 0 0; 0.26682494738595 0 0 0 1 0 0 0 0 0 0 0 0 1 0; 0.112963728335479 0 0 0 1 0 0 0 0 0 0 0 0 0 1; -0.127891846668245 0 0 0 1 0 0 0 0 0 0 0 0 0 0; 0.00113261912308242 0 0 0 0 1 0 0 1 0 0 0 0 0 0; -0.20219488165898 0 0 0 0 1 0 0 0 1 0 0 0 0 0; -0.290712932889619 0 0 0 0 1 0 0 0 0 1 0 0 0 0; -0.209242730218869 0 0 0 0 1 0 0 0 0 0 1 0 0 0; 0.217457968262388 0 0 0 0 1 0 0 0 0 0 0 0 1 0; 0.312602851836633 0 0 0 0 1 0 0 0 0 0 0 0 0 1; -0.200918407685917 0 0 0 0 1 0 0 0 0 0 0 0 0 0; 0.0768863179978937 0 0 0 0 0 1 0 1 0 0 0 0 0 0; -0.0286236322874945 0 0 0 0 0 1 0 0 1 0 0 0 0 0; -0.171304195197489 0 0 0 0 0 1 0 0 0 1 0 0 0 0; -0.26682494738595 0 0 0 0 0 1 0 0 0 0 1 0 0 0; -0.217457968262388 0 0 0 0 0 1 0 0 0 0 0 1 0 0; 0.281230208482166 0 0 0 0 0 1 0 0 0 0 0 0 0 1; -0.123296720506736 0 0 0 0 0 1 0 0 0 0 0 0 0 0; 0.109061791346399 0 0 0 0 0 0 1 1 0 0 0 0 0 0; 0.220343637084719 0 0 0 0 0 0 1 0 1 0 0 0 0 0; 0.129712785418193 0 0 0 0 0 0 1 0 0 1 0 0 0 0; -0.112963728335479 0 0 0 0 0 0 1 0 0 0 1 0 0 0; -0.312602851836633 0 0 0 0 0 0 1 0 0 0 0 1 0 0; -0.281230208482166 0 0 0 0 0 0 1 0 0 0 0 0 1 0; 0.082597199690713 0 0 0 0 0 0 1 0 0 0 0 0 0 0; -0.0993782882388173 0 0 0 0 0 0 0 1 0 0 0 0 0 0; -0.086709518011466 0 0 0 0 0 0 0 0 1 0 0 0 0 0; 0.044728633838126 0 0 0 0 0 0 0 0 0 1 0 0 0 0; 0.219004436000958 0 0 0 0 0 0 0 0 0 0 1 0 0 0; 0.283036113868953 0 0 0 0 0 0 0 0 0 0 0 1 0 0; 0.13332284695869 0 0 0 0 0 0 0 0 0 0 0 0 1 0; -0.174384058860174 0 0 0 0 0 0 0 0 0 0 0 0 0 1; 0.0372964512210185 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

    system = systemstruct(eqs);
end
function [system] = rps10()
    eqs = cell(10,1);
    eqs{1} = [-0.207515496428277 2 0 0 2 0 0 0 0 0 0; -1.81633498321741 2 0 0 1 1 0 0 0 0 0; -0.326880264194788 2 0 0 1 0 1 0 0 0 0; 0.327228678790004 2 0 0 1 0 0 0 0 0 0; 0.0871726886752159 2 0 0 0 2 0 0 0 0 0; 0.581046929846164 2 0 0 0 1 1 0 0 0 0; 1.56146164401314 2 0 0 0 1 0 0 0 0 0; 0.120342807753062 2 0 0 0 0 2 0 0 0 0; 0.289766299873838 2 0 0 0 0 1 0 0 0 0; -0.127970368707512 2 0 0 0 0 0 0 0 0 0; -0.377029684790685 1 1 0 2 0 0 0 0 0 0; -0.91671440576214 1 1 0 1 1 0 0 0 0 0; 0.657630238424344 1 1 0 1 0 1 0 0 0 0; 0.842682927567249 1 1 0 1 0 0 0 0 0 0; 0.950415464426347 1 1 0 0 2 0 0 0 0 0; -1.15573823637833 1 1 0 0 1 1 0 0 0 0; -1.73883806758226 1 1 0 0 1 0 0 0 0 0; -0.573385779635662 1 1 0 0 0 2 0 0 0 0; -1.27789279652515 1 1 0 0 0 1 0 0 0 0; -0.485961231255263 1 1 0 0 0 0 0 0 0 0; 0.798695431832303 1 0 1 2 0 0 0 0 0 0; -0.231946468231119 1 0 1 1 1 0 0 0 0 0; 0.455139341948007 1 0 1 1 0 1 0 0 0 0; -1.13714055986675 1 0 1 1 0 0 0 0 0 0; -1.0650624231277 1 0 1 0 2 0 0 0 0 0; 1.70771409335099 1 0 1 0 1 1 0 0 0 0; -0.430912104468477 1 0 1 0 1 0 0 0 0 0; 0.266366991295394 1 0 1 0 0 2 0 0 0 0; -0.581261259115421 1 0 1 0 0 1 0 0 0 0; 0.377897769852767 1 0 1 0 0 0 0 0 0 0; -0.166889068191594 0 2 0 2 0 0 0 0 0 0; 1.02033685044889 0 2 0 1 1 0 0 0 0 0; 1.10939193639721 0 2 0 1 0 1 0 0 0 0; 0.607564575703416 0 2 0 1 0 0 0 0 0 0; -0.482067565714208 0 2 0 0 2 0 0 0 0 0; -0.113679611876378 0 2 0 0 1 1 0 0 0 0; 0.0679091571307072 0 2 0 0 1 0 0 0 0 0; 0.648956633905802 0 2 0 0 0 2 0 0 0 0; 0.908789677888625 0 2 0 0 0 1 0 0 0 0; 0.306995563707175 0 2 0 0 0 0 0 0 0 0; 0.866826144775651 0 1 1 2 0 0 0 0 0 0; 0.539670777307627 0 1 1 1 1 0 0 0 0 0; 1.85538525130694 0 1 1 1 0 1 0 0 0 0; 0.229293271620915 0 1 1 1 0 0 0 0 0 0; 0.120995290927416 0 1 1 0 2 0 0 0 0 0; -0.365479427671087 0 1 1 0 1 1 0 0 0 0; 0.908627200628342 0 1 1 0 1 0 0 0 0 0; -0.987821435703067 0 1 1 0 0 2 0 0 0 0; -0.759590462498355 0 1 1 0 0 1 0 0 0 0; -0.234045440765696 0 1 1 0 0 0 0 0 0 0; 0.374404564619872 0 0 2 2 0 0 0 0 0 0; 0.795998132768522 0 0 2 1 1 0 0 0 0 0; -0.782511672202421 0 0 2 1 0 1 0 0 0 0; -0.21948911177438 0 0 2 1 0 0 0 0 0 0; 0.394894877038992 0 0 2 0 2 0 0 0 0 0; -0.467367317969786 0 0 2 0 1 1 0 0 0 0; -0.27649317513946 0 0 2 0 1 0 0 0 0 0; -0.769299441658863 0 0 2 0 0 2 0 0 0 0; 0.508489276049675 0 0 2 0 0 1 0 0 0 0; 0.0156362617850807 0 0 2 0 0 0 0 0 0 0; -0.71530414271904 0 0 0 1 0 0 1 0 0 0; 1.22441854786319 0 0 0 1 0 0 0 1 0 0; 1.28988607045863 0 0 0 1 0 0 0 0 1 0; 0.517207385575697 0 0 0 1 0 0 0 0 0 1; -1.35287762600439 0 0 0 0 1 0 1 0 0 0; -0.0584456769855244 0 0 0 0 1 0 0 1 0 0; -0.665594849718029 0 0 0 0 1 0 0 0 1 0; 0.691709405412229 0 0 0 0 1 0 0 0 0 1; -1.70704525381214 0 0 0 0 0 1 1 0 0 0; -0.377061499535853 0 0 0 0 0 1 0 1 0 0; 0.697758704890495 0 0 0 0 0 1 0 0 1 0; -1.45796722508605 0 0 0 0 0 1 0 0 0 1; -0.194661456784744 0 0 0 0 0 0 1 0 0 0; -1.05166358226696 0 0 0 0 0 0 0 1 0 0; 0.580102254517945 0 0 0 0 0 0 0 0 1 0; -0.0429214367475854 0 0 0 0 0 0 0 0 0 1];
    eqs{2} = [-0.274164348420013 2 0 0 2 0 0 0 0 0 0; 0.145555496398316 2 0 0 1 1 0 0 0 0 0; -0.615099683261694 2 0 0 1 0 1 0 0 0 0; -0.138445244156793 2 0 0 1 0 0 0 0 0 0; 0.18078854641092 2 0 0 0 2 0 0 0 0 0; 0.755639699059209 2 0 0 0 1 1 0 0 0 0; 0.413774454738696 2 0 0 0 1 0 0 0 0 0; 0.0933758020090927 2 0 0 0 0 2 0 0 0 0; 0.402729889121092 2 0 0 0 0 1 0 0 0 0; 0.160110343036881 2 0 0 0 0 0 0 0 0 0; -0.342120127755503 1 1 0 2 0 0 0 0 0 0; -1.1674745895518 1 1 0 1 1 0 0 0 0 0; 0.099371922391061 1 1 0 1 0 1 0 0 0 0; 1.55680856443337 1 1 0 1 0 0 0 0 0 0; 0.940454186982468 1 1 0 0 2 0 0 0 0 0; -1.4126149515014 1 1 0 0 1 1 0 0 0 0; -1.57893837362116 1 1 0 0 1 0 0 0 0 0; -0.598334059226964 1 1 0 0 0 2 0 0 0 0; -0.827248467395868 1 1 0 0 0 1 0 0 0 0; -0.900546882440308 1 1 0 0 0 0 0 0 0 0; 0.374589872787203 1 0 1 2 0 0 0 0 0 0; 0.00249167758187343 1 0 1 1 1 0 0 0 0 0; 0.473243820363174 1 0 1 1 0 1 0 0 0 0; -1.74094581211543 1 0 1 1 0 0 0 0 0 0; -1.02574184475855 1 0 1 0 2 0 0 0 0 0; -0.978157634258566 1 0 1 0 1 1 0 0 0 0; 0.616379366719068 1 0 1 0 1 0 0 0 0 0; 0.651151971971351 1 0 1 0 0 2 0 0 0 0; -0.123532266650023 1 0 1 0 0 1 0 0 0 0; 0.52025861583069 1 0 1 0 0 0 0 0 0 0; -0.0754243659911413 0 2 0 2 0 0 0 0 0 0; -0.94280644898765 0 2 0 1 1 0 0 0 0 0; 0.822604277549155 0 2 0 1 0 1 0 0 0 0; 1.68638623822392 0 2 0 1 0 0 0 0 0 0; -0.578003051555137 0 2 0 0 2 0 0 0 0 0; -0.0967654551556513 0 2 0 0 1 1 0 0 0 0; 1.26831951729494 0 2 0 0 1 0 0 0 0 0; 0.653427417546279 0 2 0 0 0 2 0 0 0 0; 1.05713963692447 0 2 0 0 0 1 0 0 0 0; -0.351901583868926 0 2 0 0 0 0 0 0 0 0; 0.478256199646769 0 1 1 2 0 0 0 0 0 0; 0.529162155528347 0 1 1 1 1 0 0 0 0 0; 1.69460505803347 0 1 1 1 0 1 0 0 0 0; -0.138723560936029 0 1 1 1 0 0 0 0 0 0; 0.0925177817398932 0 1 1 0 2 0 0 0 0 0; 2.50066179951447 0 1 1 0 1 1 0 0 0 0; -0.433745742064066 0 1 1 0 1 0 0 0 0 0; -0.570773981386662 0 1 1 0 0 2 0 0 0 0; -2.57418557618624 0 1 1 0 0 1 0 0 0 0; 0.908682123022068 0 1 1 0 0 0 0 0 0 0; 0.349588714411154 0 0 2 2 0 0 0 0 0 0; 0.797250952589334 0 0 2 1 1 0 0 0 0 0; -0.207504594287461 0 0 2 1 0 1 0 0 0 0; -0.515904708485933 0 0 2 1 0 0 0 0 0 0; 0.397214505144217 0 0 2 0 2 0 0 0 0 0; -0.658874243903558 0 0 2 0 1 1 0 0 0 0; -0.206145801724319 0 0 2 0 1 0 0 0 0 0; -0.746803219555371 0 0 2 0 0 2 0 0 0 0; 1.56047400768576 0 0 2 0 0 1 0 0 0 0; -0.446456217064578 0 0 2 0 0 0 0 0 0 0; -1.0320362855812 0 0 0 1 0 0 1 0 0 0; -0.215759067016851 0 0 0 1 0 0 0 1 0 0; -0.0617112344292475 0 0 0 1 0 0 0 0 1 0; 0.432148237421863 0 0 0 1 0 0 0 0 0 1; -1.47594817030931 0 0 0 0 1 0 1 0 0 0; -0.12705583446957 0 0 0 0 1 0 0 1 0 0; -0.731495415588662 0 0 0 0 1 0 0 0 1 0; 0.66776248682605 0 0 0 0 1 0 0 0 0 1; -3.02034353373132 0 0 0 0 0 1 1 0 0 0; -0.514859363952448 0 0 0 0 0 1 0 1 0 0; 0.718934807521354 0 0 0 0 0 1 0 0 1 0; 0.516257114442282 0 0 0 0 0 1 0 0 0 1; 0.638247457896623 0 0 0 0 0 0 1 0 0 0; -0.228124386753508 0 0 0 0 0 0 0 1 0 0; 1.36677938008601 0 0 0 0 0 0 0 0 1 0; -0.902118536026858 0 0 0 0 0 0 0 0 0 1];
    eqs{3} = [0.194615627960121 2 0 0 2 0 0 0 0 0 0; -0.350892341051474 2 0 0 1 1 0 0 0 0 0; -0.0694189987244619 2 0 0 1 0 1 0 0 0 0; -0.542424299098038 2 0 0 1 0 0 0 0 0 0; -0.0505869210529748 2 0 0 0 2 0 0 0 0 0; -0.291585397036852 2 0 0 0 1 1 0 0 0 0; 0.865528115844618 2 0 0 0 1 0 0 0 0 0; -0.144028706907147 2 0 0 0 0 2 0 0 0 0; 0.41615153726461 2 0 0 0 0 1 0 0 0 0; 0.208164758092194 2 0 0 0 0 0 0 0 0 0; 1.07265142556711 1 1 0 2 0 0 0 0 0 0; -1.86373257471055 1 1 0 1 1 0 0 0 0 0; -0.561272501958868 1 1 0 1 0 1 0 0 0 0; 0.412606879197033 1 1 0 1 0 0 0 0 0 0; -0.854553109316494 1 1 0 0 2 0 0 0 0 0; -1.2521117933147 1 1 0 0 1 1 0 0 0 0; -1.42270075336129 1 1 0 0 1 0 0 0 0 0; -0.218098316250618 1 1 0 0 0 2 0 0 0 0; -1.60311321241731 1 1 0 0 0 1 0 0 0 0; -0.446247956964454 1 1 0 0 0 0 0 0 0 0; 0.347415005165549 1 0 1 2 0 0 0 0 0 0; 0.742138963310435 1 0 1 1 1 0 0 0 0 0; -0.831020434150999 1 0 1 1 0 1 0 0 0 0; -0.713545443616961 1 0 1 1 0 0 0 0 0 0; -1.48243855687384 1 0 1 0 2 0 0 0 0 0; 1.23091291787156 1 0 1 0 1 1 0 0 0 0; -0.700205382747984 1 0 1 0 1 0 0 0 0 0; 1.1350235517083 1 0 1 0 0 2 0 0 0 0; 0.123669434766218 1 0 1 0 0 1 0 0 0 0; 0.24325115365766 1 0 1 0 0 0 0 0 0 0; -0.29746151974578 0 2 0 2 0 0 0 0 0 0; -0.764391087841086 0 2 0 1 1 0 0 0 0 0; 1.48353631082628 0 2 0 1 0 1 0 0 0 0; 0.0262284934902558 0 2 0 1 0 0 0 0 0 0; -0.250458095627853 0 2 0 0 2 0 0 0 0 0; 0.38706376702247 0 2 0 0 1 1 0 0 0 0; -0.356160898625373 0 2 0 0 1 0 0 0 0 0; 0.547919615373632 0 2 0 0 0 2 0 0 0 0; 1.16527685308026 0 2 0 0 0 1 0 0 0 0; 0.157345778181886 0 2 0 0 0 0 0 0 0 0; 1.56627942536379 0 1 1 2 0 0 0 0 0 0; 0.804389089622383 0 1 1 1 1 0 0 0 0 0; 1.36508876117873 0 1 1 1 0 1 0 0 0 0; 0.385667250235707 0 1 1 1 0 0 0 0 0 0; -0.276031198589472 0 1 1 0 2 0 0 0 0 0; 2.00133869763712 0 1 1 0 1 1 0 0 0 0; 0.145126372137632 0 1 1 0 1 0 0 0 0 0; -1.29024822677432 0 1 1 0 0 2 0 0 0 0; -0.0335102717324866 0 1 1 0 0 1 0 0 0 0; -0.458742499196916 0 1 1 0 0 0 0 0 0 0; 0.102845891785658 0 0 2 2 0 0 0 0 0 0; 1.11528342889256 0 0 2 1 1 0 0 0 0 0; -1.41411731210182 0 0 2 1 0 1 0 0 0 0; 0.0664193570818274 0 0 2 1 0 0 0 0 0 0; 0.301045016680828 0 0 2 0 2 0 0 0 0 0; -0.0954783699856177 0 0 2 0 1 1 0 0 0 0; -0.444642792945758 0 0 2 0 1 0 0 0 0 0; -0.403890908466486 0 0 2 0 0 2 0 0 0 0; 0.662502386860574 0 0 2 0 0 1 0 0 0 0; -0.0487485648968095 0 0 2 0 0 0 0 0 0 0; 0.449776448525955 0 0 0 1 0 0 1 0 0 0; 0.238318786510823 0 0 0 1 0 0 0 1 0 0; 1.37658917888669 0 0 0 1 0 0 0 0 1 0; 0.0892904130224598 0 0 0 1 0 0 0 0 0 1; -0.0647244242734867 0 0 0 0 1 0 1 0 0 0; 0.0319777767304853 0 0 0 0 1 0 0 1 0 0; 0.530652390187601 0 0 0 0 1 0 0 0 1 0; 0.943551163213123 0 0 0 0 1 0 0 0 0 1; -2.24393077720544 0 0 0 0 0 1 1 0 0 0; 0.166874554065645 0 0 0 0 0 1 0 1 0 0; 0.475489118193304 0 0 0 0 0 1 0 0 1 0; -1.25272501307127 0 0 0 0 0 1 0 0 0 1; -0.31676197137727 0 0 0 0 0 0 1 0 0 0; -0.554716522369026 0 0 0 0 0 0 0 1 0 0; 0.94233779062752 0 0 0 0 0 0 0 0 1 0; 0.0967323009365533 0 0 0 0 0 0 0 0 0 1];
    eqs{4} = [-0.542231318156468 2 0 0 2 0 0 0 0 0 0; -1.55876752645388 2 0 0 1 1 0 0 0 0 0; -0.16763110456854 2 0 0 1 0 1 0 0 0 0; 0.428144429217309 2 0 0 1 0 0 0 0 0 0; 0.433190480841401 2 0 0 0 2 0 0 0 0 0; 0.629720710084467 2 0 0 0 1 1 0 0 0 0; 1.13569115672556 2 0 0 0 1 0 0 0 0 0; 0.109040837315068 2 0 0 0 0 2 0 0 0 0; 0.211985556976183 2 0 0 0 0 1 0 0 0 0; -0.0409504982462884 2 0 0 0 0 0 0 0 0 0; 0.487194647173144 1 1 0 2 0 0 0 0 0 0; -1.13175815270035 1 1 0 1 1 0 0 0 0 0; -1.32554638157209 1 1 0 1 0 1 0 0 0 0; -0.882896093611704 1 1 0 1 0 0 0 0 0 0; -0.138928570583023 1 1 0 0 2 0 0 0 0 0; -2.06777132116189 1 1 0 0 1 1 0 0 0 0; -1.22210068517825 1 1 0 0 1 0 0 0 0 0; -0.348266076590121 1 1 0 0 0 2 0 0 0 0; -0.0164770915412989 1 1 0 0 0 1 0 0 0 0; 0.0431160255118422 1 1 0 0 0 0 0 0 0 0; -0.102737729000881 1 0 1 2 0 0 0 0 0 0; 1.00204398148402 1 0 1 1 1 0 0 0 0 0; -0.966917627888521 1 0 1 1 0 1 0 0 0 0; -1.39245278810297 1 0 1 1 0 0 0 0 0 0; -0.237693851366647 1 0 1 0 2 0 0 0 0 0; 0.0203913235490343 1 0 1 0 1 1 0 0 0 0; -1.65249593965271 1 0 1 0 1 0 0 0 0 0; 0.340431580367528 1 0 1 0 0 2 0 0 0 0; 0.77528301920843 1 0 1 0 0 1 0 0 0 0; 0.762977033403646 1 0 1 0 0 0 0 0 0 0; 0.33193585698171 0 2 0 2 0 0 0 0 0 0; 0.0182908073673928 0 2 0 1 1 0 0 0 0 0; 0.807847851433961 0 2 0 1 0 1 0 0 0 0; 0.567668288627952 0 2 0 1 0 0 0 0 0 0; -0.0121837101271557 0 2 0 0 2 0 0 0 0 0; -0.178628196977515 0 2 0 0 1 1 0 0 0 0; -0.0933400214333203 0 2 0 0 1 0 0 0 0 0; -0.319752146854554 0 2 0 0 0 2 0 0 0 0; -0.189452181152727 0 2 0 0 0 1 0 0 0 0; 0.00394049919878622 0 2 0 0 0 0 0 0 0 0; 0.662487411536578 0 1 1 2 0 0 0 0 0 0; 0.798404976028356 0 1 1 1 1 0 0 0 0 0; -0.368729267937399 0 1 1 1 0 1 0 0 0 0; -0.0867133900253777 0 1 1 1 0 0 0 0 0 0; 0.243828370682188 0 1 1 0 2 0 0 0 0 0; 2.7775636697444 0 1 1 0 1 1 0 0 0 0; 1.55697251241841 0 1 1 0 1 0 0 0 0 0; -0.906315782218766 0 1 1 0 0 2 0 0 0 0; -2.0514046696465 0 1 1 0 0 1 0 0 0 0; -0.849235076014679 0 1 1 0 0 0 0 0 0 0; 0.210295461174758 0 0 2 2 0 0 0 0 0 0; 1.54047671908649 0 0 2 1 1 0 0 0 0 0; -0.640216746865421 0 0 2 1 0 1 0 0 0 0; -0.725610409598415 0 0 2 1 0 0 0 0 0 0; -0.421006770714245 0 0 2 0 2 0 0 0 0 0; -0.451092513106952 0 0 2 0 1 1 0 0 0 0; 0.766938606845301 0 0 2 0 1 0 0 0 0 0; 0.210711309539487 0 0 2 0 0 2 0 0 0 0; 0.104320284945126 0 0 2 0 0 1 0 0 0 0; -0.278417478342462 0 0 2 0 0 0 0 0 0 0; -0.270202308246846 0 0 0 1 0 0 1 0 0 0; 1.43325839652983 0 0 0 1 0 0 0 1 0 0; 0.649216154064302 0 0 0 1 0 0 0 0 1 0; 1.03862986490006 0 0 0 1 0 0 0 0 0 1; -1.80928974213754 0 0 0 0 1 0 1 0 0 0; -0.365829296961495 0 0 0 0 1 0 0 1 0 0; 0.143547014784455 0 0 0 0 1 0 0 0 1 0; -0.273987073183022 0 0 0 0 1 0 0 0 0 1; -0.126853660768583 0 0 0 0 0 1 1 0 0 0; -0.458292808629767 0 0 0 0 0 1 0 1 0 0; 1.80496860452622 0 0 0 0 0 1 0 0 1 0; -0.389922892948351 0 0 0 0 0 1 0 0 0 1; 0.315427477389965 0 0 0 0 0 0 1 0 0 0; -0.972105406031357 0 0 0 0 0 0 0 1 0 0; 1.10574800017004 0 0 0 0 0 0 0 0 1 0; 0.361964167551302 0 0 0 0 0 0 0 0 0 1];
    eqs{5} = [-0.123076583143715 2 0 0 2 0 0 0 0 0 0; 0.137696942578023 2 0 0 1 1 0 0 0 0 0; -0.587515731220698 2 0 0 1 0 1 0 0 0 0; 0.516979121184011 2 0 0 1 0 0 0 0 0 0; 0.0321246390056673 2 0 0 0 2 0 0 0 0 0; 0.801736157469905 2 0 0 0 1 1 0 0 0 0; -0.478426051816978 2 0 0 0 1 0 0 0 0 0; 0.0909519441380479 2 0 0 0 0 2 0 0 0 0; 0.600257226540674 2 0 0 0 0 1 0 0 0 0; -0.416157646089455 2 0 0 0 0 0 0 0 0 0; 1.36499155854057 1 1 0 2 0 0 0 0 0 0; -1.3462021315217 1 1 0 1 1 0 0 0 0 0; -0.267947779185724 1 1 0 1 0 1 0 0 0 0; -0.328163318667352 1 1 0 1 0 0 0 0 0 0; -0.543299069419884 1 1 0 0 2 0 0 0 0 0; -0.0611990506923743 1 1 0 0 1 1 0 0 0 0; -2.0921341984233 1 1 0 0 1 0 0 0 0 0; -0.821692489120687 1 1 0 0 0 2 0 0 0 0; 0.180205552106969 1 1 0 0 0 1 0 0 0 0; -1.23311710017938 1 1 0 0 0 0 0 0 0 0; 1.28252229495899 1 0 1 2 0 0 0 0 0 0; -0.32360789087359 1 0 1 1 1 0 0 0 0 0; 0.0821787313133683 1 0 1 1 0 1 0 0 0 0; -1.78387184821123 1 0 1 1 0 0 0 0 0 0; -1.83489085782808 1 0 1 0 2 0 0 0 0 0; -0.512032167814848 1 0 1 0 1 1 0 0 0 0; -0.932207034691206 1 0 1 0 1 0 0 0 0 0; 0.552368562869093 1 0 1 0 0 2 0 0 0 0; -0.526767653424667 1 0 1 0 0 1 0 0 0 0; 0.445174124091856 1 0 1 0 0 0 0 0 0 0; -0.0164114427593356 0 2 0 2 0 0 0 0 0 0; -0.358310068008011 0 2 0 1 1 0 0 0 0 0; 1.81617668797614 0 2 0 1 0 1 0 0 0 0; 0.10768082059655 0 2 0 1 0 0 0 0 0 0; -0.108131323761818 0 2 0 0 2 0 0 0 0 0; 0.532014075603258 0 2 0 0 1 1 0 0 0 0; -0.265247887538097 0 2 0 0 1 0 0 0 0 0; 0.124542766521154 0 2 0 0 0 2 0 0 0 0; 1.56443317441967 0 2 0 0 0 1 0 0 0 0; 0.104235944986372 0 2 0 0 0 0 0 0 0 0; 0.131183896772422 0 1 1 2 0 0 0 0 0 0; -0.233129735700999 0 1 1 1 1 0 0 0 0 0; 0.416252286227597 0 1 1 1 0 1 0 0 0 0; 0.196238595543859 0 1 1 1 0 0 0 0 0 0; -0.337948478287682 0 1 1 0 2 0 0 0 0 0; 1.39949592457995 0 1 1 0 1 1 0 0 0 0; 0.282296457934625 0 1 1 0 1 0 0 0 0 0; 0.20676458151526 0 1 1 0 0 2 0 0 0 0; -1.67615733877489 0 1 1 0 0 1 0 0 0 0; -0.0807794759847403 0 1 1 0 0 0 0 0 0 0; 0.139488025903051 0 0 2 2 0 0 0 0 0 0; 0.220613125429988 0 0 2 1 1 0 0 0 0 0; -1.22866095675544 0 0 2 1 0 1 0 0 0 0; -0.0932755727182936 0 0 2 1 0 0 0 0 0 0; 0.0760066847561506 0 0 2 0 2 0 0 0 0 0; -1.33375023307316 0 0 2 0 1 1 0 0 0 0; -0.0743800369279021 0 0 2 0 1 0 0 0 0 0; -0.215494710659201 0 0 2 0 0 2 0 0 0 0; 0.400830791369246 0 0 2 0 0 1 0 0 0 0; 0.0155848221518674 0 0 2 0 0 0 0 0 0 0; -0.531384369062268 0 0 0 1 0 0 1 0 0 0; 1.1609512634986 0 0 0 1 0 0 0 1 0 0; -0.0139734087649947 0 0 0 1 0 0 0 0 1 0; 0.124974829380608 0 0 0 1 0 0 0 0 0 1; 0.818053976282977 0 0 0 0 1 0 1 0 0 0; -1.52895122887586 0 0 0 0 1 0 0 1 0 0; -0.300803742720313 0 0 0 0 1 0 0 0 1 0; 0.0633797450007162 0 0 0 0 1 0 0 0 0 1; -2.56552119232959 0 0 0 0 0 1 1 0 0 0; -0.174245615059662 0 0 0 0 0 1 0 1 0 0; 1.54732862760874 0 0 0 0 0 1 0 0 1 0; -1.12232290624583 0 0 0 0 0 1 0 0 0 1; 0.296336878951216 0 0 0 0 0 0 1 0 0 0; -1.26020892235827 0 0 0 0 0 0 0 1 0 0; 0.103649018766031 0 0 0 0 0 0 0 0 1 0; -0.0427614568633893 0 0 0 0 0 0 0 0 0 1];
    eqs{6} = [-0.0682076866558154 2 0 0 2 0 0 0 0 0 0; 0.0457910614219832 2 0 0 1 1 0 0 0 0 0; 0.0331360143111403 2 0 0 1 0 1 0 0 0 0; -0.0884128665796767 2 0 0 1 0 0 0 0 0 0; 0.0627040723671076 2 0 0 0 2 0 0 0 0 0; 0.0406719410026242 2 0 0 0 1 1 0 0 0 0; 0.219647969289755 2 0 0 0 1 0 0 0 0 0; 0.00550361428870782 2 0 0 0 0 2 0 0 0 0; 0.0913689660756994 2 0 0 0 0 1 0 0 0 0; 0.0995239208560676 2 0 0 0 0 0 0 0 0 0; 1.30430608641081 1 1 0 2 0 0 0 0 0 0; -1.61361776738871 1 1 0 1 1 0 0 0 0 0; 0.999575794109116 1 1 0 1 0 1 0 0 0 0; -0.543744076782795 1 1 0 1 0 0 0 0 0 0; -0.191627925861076 1 1 0 0 2 0 0 0 0 0; -2.52714783213803 1 1 0 0 1 1 0 0 0 0; -0.721386997138274 1 1 0 0 1 0 0 0 0 0; -1.11267816054973 1 1 0 0 0 2 0 0 0 0; -1.58312714998811 1 1 0 0 0 1 0 0 0 0; -0.406983564612776 1 1 0 0 0 0 0 0 0 0; 1.30775741185764 1 0 1 2 0 0 0 0 0 0; -0.318087002835764 1 0 1 1 1 0 0 0 0 0; -1.87684508222007 1 0 1 1 0 1 0 0 0 0; -2.30082830180562 1 0 1 1 0 0 0 0 0 0; -1.96485022822957 1 0 1 0 2 0 0 0 0 0; 0.542624987335489 1 0 1 0 1 1 0 0 0 0; 0.676969390772226 1 0 1 0 1 0 0 0 0 0; 0.657092816371931 1 0 1 0 0 2 0 0 0 0; 1.5619440370847 1 0 1 0 0 1 0 0 0 0; 0.845110778688351 1 0 1 0 0 0 0 0 0 0; -0.00905406506209638 0 2 0 2 0 0 0 0 0 0; -0.0831370521557332 0 2 0 1 1 0 0 0 0 0; 0.000567364575571421 0 2 0 1 0 1 0 0 0 0; 0.831018754794303 0 2 0 1 0 0 0 0 0 0; -0.00906823242262303 0 2 0 0 2 0 0 0 0 0; -0.112214431521285 0 2 0 0 1 1 0 0 0 0; 0.0642888006202359 0 2 0 0 1 0 0 0 0 0; 0.0181222974847194 0 2 0 0 0 2 0 0 0 0; 1.15862032388824 0 2 0 0 0 1 0 0 0 0; 0.282266186700428 0 2 0 0 0 0 0 0 0 0; 0.0676847422464324 0 1 1 2 0 0 0 0 0 0; -0.132325340467337 0 1 1 1 1 0 0 0 0 0; 0.155743818487839 0 1 1 1 0 1 0 0 0 0; 1.75762120445086 0 1 1 1 0 0 0 0 0 0; -0.0969698396601749 0 1 1 0 2 0 0 0 0 0; 0.0184466497266884 0 1 1 0 1 1 0 0 0 0; 1.09388891868304 0 1 1 0 1 0 0 0 0 0; 0.0292850974137425 0 1 1 0 0 2 0 0 0 0; 0.537376562278144 0 1 1 0 0 1 0 0 0 0; -0.531670897845279 0 1 1 0 0 0 0 0 0 0; 0.0772617517179118 0 0 2 2 0 0 0 0 0 0; 0.03734599073375 0 0 2 1 1 0 0 0 0 0; -0.0337033788867117 0 0 2 1 0 1 0 0 0 0; 0.812356773433574 0 0 2 1 0 0 0 0 0 0; -0.0536358399444846 0 0 2 0 2 0 0 0 0 0; 0.0715424905186608 0 0 2 0 1 1 0 0 0 0; 1.03377817636989 0 0 2 0 1 0 0 0 0 0; -0.0236259117734272 0 0 2 0 0 2 0 0 0 0; -0.751007808577168 0 0 2 0 0 1 0 0 0 0; -0.879810896388137 0 0 2 0 0 0 0 0 0 0; -1.5549626616482 0 0 0 1 0 0 1 0 0 0; -0.0611151070662345 0 0 0 1 0 0 0 1 0 0; 1.16825335990701 0 0 0 1 0 0 0 0 1 0; 1.17536979107655 0 0 0 1 0 0 0 0 0 1; -1.31771494627988 0 0 0 0 1 0 1 0 0 0; -0.0415701442807332 0 0 0 0 1 0 0 1 0 0; 0.128978831838364 0 0 0 0 1 0 0 0 1 0; 1.30374626562344 0 0 0 0 1 0 0 0 0 1; -0.49898148138677 0 0 0 0 0 1 1 0 0 0; -0.00799683232729561 0 0 0 0 0 1 0 1 0 0; 1.61660744020337 0 0 0 0 0 1 0 0 1 0; -0.955052946324799 0 0 0 0 0 1 0 0 0 1; 0.498020788831642 0 0 0 0 0 0 1 0 0 0; -0.123388245737918 0 0 0 0 0 0 0 1 0 0; 0.395354545362724 0 0 0 0 0 0 0 0 1 0; -1.16343765243919 0 0 0 0 0 0 0 0 0 1];
    eqs{7} = [-0.169712144553532 2 0 0 2 0 0 0 0 0 0; -1.12984565416907 2 0 0 1 1 0 0 0 0 0; -1.31149805422042 2 0 0 1 0 1 0 0 0 0; -1.57442020054952 2 0 0 1 0 0 0 0 0 0; 0.546165736846632 2 0 0 0 2 0 0 0 0 0; 0.210400444896558 2 0 0 0 1 1 0 0 0 0; 0.342888715346696 2 0 0 0 1 0 0 0 0 0; -0.3764535922931 2 0 0 0 0 2 0 0 0 0; -0.818696931313652 2 0 0 0 0 1 0 0 0 0; -0.439462593920411 2 0 0 0 0 0 0 0 0 0; 0.176491210445206 1 1 0 2 0 0 0 0 0 0; -2.03952445704438 1 1 0 1 1 0 0 0 0 0; 0.745869398754165 1 1 0 1 0 1 0 0 0 0; -0.286166350617928 1 1 0 1 0 0 0 0 0 0; -0.40076456254394 1 1 0 0 2 0 0 0 0 0; -2.09467429476926 1 1 0 0 1 1 0 0 0 0; -2.05201729549009 1 1 0 0 1 0 0 0 0 0; 0.224273352098734 1 1 0 0 0 2 0 0 0 0; -0.0292507674785184 1 1 0 0 0 1 0 0 0 0; -0.279419440031289 1 1 0 0 0 0 0 0 0 0; 1.54872221545522 1 0 1 2 0 0 0 0 0 0; -0.793722785103128 1 0 1 1 1 0 0 0 0 0; -0.560905070395481 1 0 1 1 0 1 0 0 0 0; 0.533176901737327 1 0 1 1 0 0 0 0 0 0; -1.26069133915014 1 0 1 0 2 0 0 0 0 0; -0.774182834231516 1 0 1 0 1 1 0 0 0 0; -1.16611569860173 1 0 1 0 1 0 0 0 0 0; -0.288030876305083 1 0 1 0 0 2 0 0 0 0; 0.128489937325171 1 0 1 0 0 1 0 0 0 0; 0.483784293228967 1 0 1 0 0 0 0 0 0 0; 0.230772384999563 0 2 0 2 0 0 0 0 0 0; -0.461588020714783 0 2 0 1 1 0 0 0 0 0; 1.30944305602726 0 2 0 1 0 1 0 0 0 0; 0.73869418690517 0 2 0 1 0 0 0 0 0 0; -0.578944896323666 0 2 0 0 2 0 0 0 0 0; 0.901493574776607 0 2 0 0 1 1 0 0 0 0; 0.254843437739574 0 2 0 0 1 0 0 0 0 0; 0.348172511324102 0 2 0 0 0 2 0 0 0 0; 0.739213520073845 0 2 0 0 0 1 0 0 0 0; 0.286335138172419 0 2 0 0 0 0 0 0 0 0; 1.15491900621773 0 1 1 2 0 0 0 0 0 0; 1.39258639807434 0 1 1 1 1 0 0 0 0 0; -0.560863866588139 0 1 1 1 0 1 0 0 0 0; 0.375431323864404 0 1 1 1 0 0 0 0 0 0; -0.883987084162639 0 1 1 0 2 0 0 0 0 0; 1.01797051417409 0 1 1 0 1 1 0 0 0 0; 1.14180163111881 0 1 1 0 1 0 0 0 0 0; -0.270931922055091 0 1 1 0 0 2 0 0 0 0; -0.924672196116106 0 1 1 0 0 1 0 0 0 0; -0.278944762251348 0 1 1 0 0 0 0 0 0 0; -0.0610602404460308 0 0 2 2 0 0 0 0 0 0; 1.59143367488385 0 0 2 1 1 0 0 0 0 0; 0.0020549981931534 0 0 2 1 0 1 0 0 0 0; -0.495572723429795 0 0 2 1 0 0 0 0 0 0; 0.0327791594770331 0 0 2 0 2 0 0 0 0 0; -1.11189401967317 0 0 2 0 1 1 0 0 0 0; -0.0371732602372415 0 0 2 0 1 0 0 0 0 0; 0.0282810809689977 0 0 2 0 0 2 0 0 0 0; 0.347118065573936 0 0 2 0 0 1 0 0 0 0; 0.00840391192868436 0 0 2 0 0 0 0 0 0 0; 1.33129873707415 0 0 0 1 0 0 1 0 0 0; 1.1967487089081 0 0 0 1 0 0 0 1 0 0; 0.905704515833505 0 0 0 1 0 0 0 0 1 0; -0.0473888920688874 0 0 0 1 0 0 0 0 0 1; -0.560558892849028 0 0 0 0 1 0 1 0 0 0; -0.541678872953628 0 0 0 0 1 0 0 1 0 0; 0.790785041402257 0 0 0 0 1 0 0 0 1 0; 1.23608801426924 0 0 0 0 1 0 0 0 0 1; -0.267634654334129 0 0 0 0 0 1 1 0 0 0; 0.357321613739484 0 0 0 0 0 1 0 1 0 0; 0.253306442722846 0 0 0 0 0 1 0 0 1 0; -0.0314635718174723 0 0 0 0 0 1 0 0 0 1; 0.144723543819308 0 0 0 0 0 0 1 0 0 0; 0.344749681587688 0 0 0 0 0 0 0 1 0 0; 0.408683409015256 0 0 0 0 0 0 0 0 1 0; -0.385416258634141 0 0 0 0 0 0 0 0 0 1];
    eqs{8} = [0.128673369147191 2 0 0 2 0 0 0 0 0 0; -1.38234500115645 2 0 0 1 1 0 0 0 0 0; -0.802841955600137 2 0 0 1 0 1 0 0 0 0; -1.06554289952459 2 0 0 1 0 0 0 0 0 0; 0.067353014566104 2 0 0 0 2 0 0 0 0 0; -0.282995108376947 2 0 0 0 1 1 0 0 0 0; -0.365065432035237 2 0 0 0 1 0 0 0 0 0; -0.196026383713295 2 0 0 0 0 2 0 0 0 0; -0.513700726163696 2 0 0 0 0 1 0 0 0 0; -0.336498513669433 2 0 0 0 0 0 0 0 0 0; 1.24221646845706 1 1 0 2 0 0 0 0 0 0; -0.109218236628474 1 1 0 1 1 0 0 0 0 0; 1.49793621886232 1 1 0 1 0 1 0 0 0 0; 1.20155250366535 1 1 0 1 0 0 0 0 0 0; -0.858540836260963 1 1 0 0 2 0 0 0 0 0; -1.67180387717421 1 1 0 0 1 1 0 0 0 0; -1.8191248166962 1 1 0 0 1 0 0 0 0 0; -0.383675632196097 1 1 0 0 0 2 0 0 0 0; -0.982032732017591 1 1 0 0 0 1 0 0 0 0; -0.622071337766898 1 1 0 0 0 0 0 0 0 0; 0.988903111633841 1 0 1 2 0 0 0 0 0 0; 0.641617617224501 1 0 1 1 1 0 0 0 0 0; -1.06256691641945 1 0 1 1 0 1 0 0 0 0; -0.197381874002467 1 0 1 1 0 0 0 0 0 0; -0.400141848372886 1 0 1 0 2 0 0 0 0 0; -0.455733535175845 1 0 1 0 1 1 0 0 0 0; 0.226509921918815 1 0 1 0 1 0 0 0 0 0; -0.588761263260955 1 0 1 0 0 2 0 0 0 0; -0.700348728440972 1 0 1 0 0 1 0 0 0 0; 0.0969124161690369 1 0 1 0 0 0 0 0 0 0; -0.32698199621343 0 2 0 2 0 0 0 0 0 0; 0.94112047467914 0 2 0 1 1 0 0 0 0 0; 0.679863927282093 0 2 0 1 0 1 0 0 0 0; 0.760110382580537 0 2 0 1 0 0 0 0 0 0; -0.309076225851502 0 2 0 0 2 0 0 0 0 0; 0.228628408297259 0 2 0 0 1 1 0 0 0 0; -0.0470235306249258 0 2 0 0 1 0 0 0 0 0; 0.636058222064932 0 2 0 0 0 2 0 0 0 0; 0.926095293449879 0 2 0 0 0 1 0 0 0 0; 0.302536639890857 0 2 0 0 0 0 0 0 0 0; 0.0159980209440115 0 1 1 2 0 0 0 0 0 0; 0.605017898829773 0 1 1 1 1 0 0 0 0 0; 1.64997375301907 0 1 1 1 0 1 0 0 0 0; 0.615156982334222 0 1 1 1 0 0 0 0 0 0; 0.477967548464756 0 1 1 0 2 0 0 0 0 0; -0.0303995602496141 0 1 1 0 1 1 0 0 0 0; 0.675739415961882 0 1 1 0 1 0 0 0 0 0; -0.493965569408768 0 1 1 0 0 2 0 0 0 0; -0.943037489506003 0 1 1 0 0 1 0 0 0 0; -0.234752870814125 0 1 1 0 0 0 0 0 0 0; 0.198308627066239 0 0 2 2 0 0 0 0 0 0; 0.441224526477308 0 0 2 1 1 0 0 0 0 0; 0.122978028318044 0 0 2 1 0 1 0 0 0 0; -0.198964445570808 0 0 2 1 0 0 0 0 0 0; 0.241723211285398 0 0 2 0 2 0 0 0 0 0; 0.054366700079688 0 0 2 0 1 1 0 0 0 0; -0.207840065145108 0 0 2 0 1 0 0 0 0 0; -0.440031838351637 0 0 2 0 0 2 0 0 0 0; 0.0886821730285199 0 0 2 0 0 1 0 0 0 0; 0.0375920488021885 0 0 2 0 0 0 0 0 0 0; 0.504396962514857 0 0 0 1 0 0 1 0 0 0; 0.914185822067171 0 0 0 1 0 0 0 1 0 0; -0.234410230701992 0 0 0 1 0 0 0 0 1 0; 0.188375974173305 0 0 0 1 0 0 0 0 0 1; 0.619929027805271 0 0 0 0 1 0 1 0 0 0; -0.0447464202297084 0 0 0 0 1 0 0 1 0 0; 0.586055133336572 0 0 0 0 1 0 0 0 1 0; 0.235297469969392 0 0 0 0 1 0 0 0 0 1; -0.501076740314703 0 0 0 0 0 1 1 0 0 0; 0.215101073556974 0 0 0 0 0 1 0 1 0 0; 0.651461084826877 0 0 0 0 0 1 0 0 1 0; 0.345029693650027 0 0 0 0 0 1 0 0 0 1; -0.00363017502361266 0 0 0 0 0 0 1 0 0 0; 0.278457550638103 0 0 0 0 0 0 0 1 0 0; 0.626115216359136 0 0 0 0 0 0 0 0 1 0; -0.14143984284375 0 0 0 0 0 0 0 0 0 1];
    eqs{9} = [-0.117002115381675 2 0 0 2 0 0 0 0 0 0; -1.74847883564849 2 0 0 1 1 0 0 0 0 0; -0.0681141308664966 2 0 0 1 0 1 0 0 0 0; -1.55440416275564 2 0 0 1 0 0 0 0 0 0; 0.0524174914208927 2 0 0 0 2 0 0 0 0 0; 0.891835794098135 2 0 0 0 1 1 0 0 0 0; -0.532212191693489 2 0 0 0 1 0 0 0 0 0; 0.0645846239607825 2 0 0 0 0 2 0 0 0 0; 0.726322812157602 2 0 0 0 0 1 0 0 0 0; -0.499715382680063 2 0 0 0 0 0 0 0 0 0; 1.42921686648242 1 1 0 2 0 0 0 0 0 0; -0.74626698952358 1 1 0 1 1 0 0 0 0 0; -0.31001429706461 1 1 0 1 0 1 0 0 0 0; 2.37399503336766 1 1 0 1 0 0 0 0 0 0; -1.03490808839995 1 1 0 0 2 0 0 0 0 0; -1.50042030622453 1 1 0 0 1 1 0 0 0 0; -0.570583676220818 1 1 0 0 1 0 0 0 0 0; -0.394308778082465 1 1 0 0 0 2 0 0 0 0; -2.22106589619302 1 1 0 0 0 1 0 0 0 0; 1.13798548335481 1 1 0 0 0 0 0 0 0 0; 0.117561497872033 1 0 1 2 0 0 0 0 0 0; -0.836428661723292 1 0 1 1 1 0 0 0 0 0; -1.85264538393266 1 0 1 1 0 1 0 0 0 0; -1.25116789358166 1 0 1 1 0 0 0 0 0 0; -0.870211964862699 1 0 1 0 2 0 0 0 0 0; -1.08784599474853 1 0 1 0 1 1 0 0 0 0; 0.16266860394272 1 0 1 0 1 0 0 0 0 0; 0.752650466990666 1 0 1 0 0 2 0 0 0 0; -1.0504867050872 1 0 1 0 0 1 0 0 0 0; 0.544117767318216 1 0 1 0 0 0 0 0 0 0; -0.085618875067609 0 2 0 2 0 0 0 0 0 0; 0.902493014258047 0 2 0 1 1 0 0 0 0 0; 1.41534491202543 0 2 0 1 0 1 0 0 0 0; -0.626498203472872 0 2 0 1 0 0 0 0 0 0; -0.445613806863466 0 2 0 0 2 0 0 0 0 0; -0.436330028670723 0 2 0 0 1 1 0 0 0 0; 1.33855299774799 0 2 0 0 1 0 0 0 0 0; 0.531232681931075 0 2 0 0 0 2 0 0 0 0; 1.61088929689305 0 2 0 0 0 1 0 0 0 0; -0.64743092481944 0 2 0 0 0 0 0 0 0 0; 1.2753331275028 0 1 1 2 0 0 0 0 0 0; 0.103152118690796 0 1 1 1 1 0 0 0 0 0; 1.33218223933555 0 1 1 1 0 1 0 0 0 0; 1.01452814587935 0 1 1 1 0 0 0 0 0 0; 0.134763379898546 0 1 1 0 2 0 0 0 0 0; -0.151028039266065 0 1 1 0 1 1 0 0 0 0; 0.825593864852827 0 1 1 0 1 0 0 0 0 0; -1.41009650740134 0 1 1 0 0 2 0 0 0 0; 1.06590993597413 0 1 1 0 0 1 0 0 0 0; -0.635620936262422 0 1 1 0 0 0 0 0 0 0; 0.202620990449284 0 0 2 2 0 0 0 0 0 0; 0.845985821390447 0 0 2 1 1 0 0 0 0 0; -1.34723078115893 0 0 2 1 0 1 0 0 0 0; -0.0736031466992901 0 0 2 1 0 0 0 0 0 0; 0.393196315442573 0 0 2 0 2 0 0 0 0 0; -0.455505765427412 0 0 2 0 1 1 0 0 0 0; -0.287740373903234 0 0 2 0 1 0 0 0 0 0; -0.595817305891858 0 0 2 0 0 2 0 0 0 0; 0.5672833993063 0 0 2 0 0 1 0 0 0 0; -0.00249167421336129 0 0 2 0 0 0 0 0 0 0; 2.2545055129278 0 0 0 1 0 0 1 0 0 0; 1.15613847743635 0 0 0 1 0 0 0 1 0 0; -0.0696028370001405 0 0 0 1 0 0 0 0 1 0; 1.04770954693153 0 0 0 1 0 0 0 0 0 1; -0.518600432151271 0 0 0 0 1 0 1 0 0 0; -0.034590519682761 0 0 0 0 1 0 0 1 0 0; 0.697522840619187 0 0 0 0 1 0 0 0 1 0; 0.558176868433451 0 0 0 0 1 0 0 0 0 1; -2.90449550835695 0 0 0 0 0 1 1 0 0 0; -0.586010417443472 0 0 0 0 0 1 0 1 0 0; 1.17614457684379 0 0 0 0 0 1 0 0 1 0; 0.436057648153209 0 0 0 0 0 1 0 0 0 1; 1.14963798171286 0 0 0 0 0 0 1 0 0 0; 0.3811371457709 0 0 0 0 0 0 0 1 0 0; -0.422623054759932 0 0 0 0 0 0 0 0 1 0; -0.413250478876494 0 0 0 0 0 0 0 0 0 1];
    eqs{10} = [0.933614330874605 1 0 0 0 0 0 0 0 0 0; 1.17815802717665 0 1 0 0 0 0 0 0 0 0; 0.551650235747964 0 0 1 0 0 0 0 0 0 0; -1 0 0 0 0 0 0 0 0 0 0];

    system = systemstruct(eqs);
end
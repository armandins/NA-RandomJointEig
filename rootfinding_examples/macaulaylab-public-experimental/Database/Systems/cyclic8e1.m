function [system] = cyclic8e1()
    eqs = cell(9,1);
    eqs{1} = [1 1 0 0 0 0 0 0 0 0; 1 0 1 0 0 0 0 0 0 0; 1 0 0 1 0 0 0 0 0 0; 1 0 0 0 1 0 0 0 0 0; 1 0 0 0 0 1 0 0 0 0; 1 0 0 0 0 0 1 0 0 0; 1 0 0 0 0 0 0 1 0 0; 1 0 0 0 0 0 0 0 1 0; -0.667083740234375-0.40496826171875i 0 0 0 0 0 0 0 0 1];
    eqs{2} = [1 1 1 0 0 0 0 0 0 0; 1 1 0 0 0 0 0 0 1 0; 1 0 1 1 0 0 0 0 0 0; 1 0 0 1 1 0 0 0 0 0; 1 0 0 0 1 1 0 0 0 0; 1 0 0 0 0 1 1 0 0 0; 1 0 0 0 0 0 1 1 0 0; 1 0 0 0 0 0 0 1 1 0; -0.745269775390625-0.704345703125i 0 0 0 0 0 0 0 0 1];
    eqs{3} = [1 1 1 1 0 0 0 0 0 0; 1 1 1 0 0 0 0 0 1 0; 1 1 0 0 0 0 0 1 1 0; 1 0 1 1 1 0 0 0 0 0; 1 0 0 1 1 1 0 0 0 0; 1 0 0 0 1 1 1 0 0 0; 1 0 0 0 0 1 1 1 0 0; 1 0 0 0 0 0 1 1 1 0; 0.200958251953125+0.49725341796875i 0 0 0 0 0 0 0 0 1];
    eqs{4} = [1 1 1 1 1 0 0 0 0 0; 1 1 1 1 0 0 0 0 1 0; 1 1 1 0 0 0 0 1 1 0; 1 1 0 0 0 0 1 1 1 0; 1 0 1 1 1 1 0 0 0 0; 1 0 0 1 1 1 1 0 0 0; 1 0 0 0 1 1 1 1 0 0; 1 0 0 0 0 1 1 1 1 0; 0.405975341796875+0.8092041015625i 0 0 0 0 0 0 0 0 1];
    eqs{5} = [1 1 1 1 1 1 0 0 0 0; 1 1 1 1 1 0 0 0 1 0; 1 1 1 1 0 0 0 1 1 0; 1 1 1 0 0 0 1 1 1 0; 1 1 0 0 0 1 1 1 1 0; 1 0 1 1 1 1 1 0 0 0; 1 0 0 1 1 1 1 1 0 0; 1 0 0 0 1 1 1 1 1 0; 0.354156494140625+0.09088134765625i 0 0 0 0 0 0 0 0 1];
    eqs{6} = [1 1 1 1 1 1 1 0 0 0; 1 1 1 1 1 1 0 0 1 0; 1 1 1 1 1 0 0 1 1 0; 1 1 1 1 0 0 1 1 1 0; 1 1 1 0 0 1 1 1 1 0; 1 1 0 0 1 1 1 1 1 0; 1 0 1 1 1 1 1 1 0 0; 1 0 0 1 1 1 1 1 1 0; 0.779876708984375+0.45166015625i 0 0 0 0 0 0 0 0 1];
    eqs{7} = [1 1 1 1 1 1 1 1 0 0; 1 1 1 1 1 1 1 0 1 0; 1 1 1 1 1 1 0 1 1 0; 1 1 1 1 1 0 1 1 1 0; 1 1 1 1 0 1 1 1 1 0; 1 1 1 0 1 1 1 1 1 0; 1 1 0 1 1 1 1 1 1 0; 1 0 1 1 1 1 1 1 1 0; 0.667510986328125+0.25091552734375i 0 0 0 0 0 0 0 0 1];
    eqs{8} = [1 1 1 1 1 1 1 1 1 0; -0.748565673828125+0.0980224609375i 0 0 0 0 0 0 0 0 1; -1 0 0 0 0 0 0 0 0 0];
    eqs{9} = [-0.786651611328125+0.0364990234375i 1 0 0 0 0 0 0 0 0; -0.787689208984375-0.59979248046875i 0 1 0 0 0 0 0 0 0; -0.623687744140625-0.969482421875i 0 0 1 0 0 0 0 0 0; -0.310272216796875+0.28680419921875i 0 0 0 1 0 0 0 0 0; 0.386932373046875+0.7784423828125i 0 0 0 0 1 0 0 0 0; -0.047698974609375+0.36480712890625i 0 0 0 0 0 1 0 0 0; -0.879791259765625-0.8447265625i 0 0 0 0 0 0 1 0 0; 0.875030517578125-0.49078369140625i 0 0 0 0 0 0 0 1 0; 0.016021728515625+0.85235595703125i 0 0 0 0 0 0 0 0 1; -0.354949951171875-0.16998291015625i 0 0 0 0 0 0 0 0 0];

    system = systemstruct(eqs);
end
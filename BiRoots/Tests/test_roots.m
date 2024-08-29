function rez = test_roots(P1,P2,x,y)

if isempty(x)
    rez = [0 0 0 0 (size(P1,1)-1)*(size(P2,1)-1)];
else
    v1 = bipolyval(P1,x,y);
    v2 = bipolyval(P2,x,y);
    oc = ([abs(v1) abs(v2)]);
    
    [dx1, dy1] = bipolyder(P1,x,y);
    [dx2, dy2] = bipolyder(P2,x,y);
    for k = 1:length(x)
        obc(k,1) = norm([dx1(k) dy1(k); dx2(k) dy2(k)]);
    end
    
    reloc = [oc(:,1)./obc oc(:,2)./obc];
    reloc = max(reloc.').';
    rez = [max(max(oc)) max(reloc) max(obc) sum(reloc<1e-10) sum(reloc>=1e-10)];
    %oc
end 
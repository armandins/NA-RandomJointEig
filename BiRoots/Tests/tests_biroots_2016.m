for degr = 2:15   
    runs = 10;
for kompleks = 0:1
    
Example = [degr runs kompleks]

t1 = []; t2 = []; t3 = []; t4 = []; v1 = 0; w1 = 0; l1 = 0;
s1 = []; s2 = []; s3 = []; s4 = []; v2 = 0; w2 = 0; l2 = 0;
o1 = []; o2 = []; o3 = []; o4 = []; v3 = 0; w3 = 0; l3 = 0;
z1 = []; z2 = []; z3 = []; z4 = []; v4 = 0; w4 = 0; l4 = 0;

opts = [];
%opts.mingap = 2.5;
opts.rankeps = 1e-12;
opts.max_linerr = 2e-11;
opts.fixedrankgap = 0.5;
opts.mingap = 0.5;
opts.minimalgap = 2;

% warm up
rand('seed',0);
P1 = rot90(triu(rand(degr))); 
P2 = rot90(triu(rand(degr)));
P1 = P1 + 1i * rot90(triu(rand(degr))); 
P2 = P2 + 1i * rot90(triu(rand(degr)));
[x,y,stat] = biroots(P1,P2,6);

%test 1 : Default
if 1==1
fprintf('Biroots default:')
for k = 1:runs
    rand('seed',k);
    fprintf('.')
    P1 = rot90(triu(rand(degr))); P2 = rot90(triu(rand(degr)));
    if kompleks
        P1 = P1 + 1i * rot90(triu(rand(degr))); 
        P2 = P2 + 1i * rot90(triu(rand(degr)));
    end
    tic
    [x,y,stat,lin] = biroots(P1,P2,0,opts);
    t1(k,1) = toc;
    z1(k,1) = length(x);
    ocena = test_roots(P1,P2,x,y);
    o1(k,1) = ocena(2);
    s1(k,1) = stat;
    v1 = v1 + ocena(4);
    w1 = w1 + ocena(5);
    l1 = l1 + lin;
end
fprintf('\n')
rez1 = [sum(t1)/runs sum(s1) sum(z1) max(o1) sum(o1>1e-10) sum(log10(o1))/runs v1 w1 l1];
disp(rez1)
end

% test 2 : nn
if degr<=10
fprintf('Testing nn     :')
%disp('--------------------------------------------------');
for k = 1:runs
    rand('seed',k);
    fprintf('.')
    P1 = rot90(triu(rand(degr))); P2 = rot90(triu(rand(degr)));
    if kompleks
        P1 = P1 + 1i * rot90(triu(rand(degr))); 
        P2 = P2 + 1i * rot90(triu(rand(degr)));
    end
    tic
    [x,y,stat,lin] = biroots(P1,P2,5,opts);
    t2(k,1) = toc;
    z2(k,1) = length(x);
    s2(k,1) = stat;
    v2 = v2 + ocena(4);
    w2 = w2 + ocena(5);
    l2 = l2 + lin;
    if length(x)>0
        ocena = test_roots(P1,P2,x,y);
        o2(k,1) = ocena(2);
    else
        o2(k,1) = 1;
    end
end
fprintf('\n')
rez2 = [sum(t2)/runs sum(s2) sum(z2) max(o2) length(find(o2>1e-10)) sum(log10(o2))/runs v2 w2 l2];
disp(rez2)
end

% test 3 : Lin2
if degr<=8
fprintf('Testing Lin 2  :')
%disp('--------------------------------------------------');
for k = 1:runs
    rand('seed',k);
    fprintf('.')
    P1 = rot90(triu(rand(degr))); P2 = rot90(triu(rand(degr)));
    if kompleks
        P1 = P1 + 1i * rot90(triu(rand(degr))); 
        P2 = P2 + 1i * rot90(triu(rand(degr)));
    end
    tic
    [x,y,stat] = biroots(P1,P2,3,opts);
    t3(k,1) = toc;
    z3(k,1) = length(x);
    s3(k,1) = stat;
    v3 = v3 + ocena(4);
    w3 = w3 + ocena(5);
    if length(x)>0
        ocena = test_roots(P1,P2,x,y);
        o3(k,1) = ocena(2);
    else
        o3(k,1) = 1;
    end
end    
fprintf('\n')
rez3 = [sum(t3)/runs sum(s3) sum(z3) max(o3) length(find(o3>1e-10)) sum(log10(o3))/runs v3 w3];
disp(rez3)
end

% test 4 : LinUnif
if degr<=10
opts.mingap  = 0.5;
fprintf('Testing UnifL  :')
%disp('--------------------------------------------------');
for k=1:runs
    rand('seed',k);
    fprintf('.')
    P1 = rot90(triu(rand(degr))); P2 = rot90(triu(rand(degr)));
    if kompleks
        P1 = P1 + 1i * rot90(triu(rand(degr))); 
        P2 = P2 + 1i * rot90(triu(rand(degr)));
    end
    tic
    [x,y,stat] = biroots(P1,P2,6,opts);
    t4(k,1) = toc;
    z4(k,1) = length(x);
    v4 = v4 + ocena(4);
    w4 = w4 + ocena(5);
    ocena = test_roots(P1,P2,x,y);
    o4(k,1) = ocena(2);
    s4(k,1) = stat;
end
fprintf('\n')
rez4 = [sum(t4)/runs sum(s4) sum(z4) max(o4) length(find(o4>1e-10)) sum(log10(o4))/runs v4 w4];
disp(rez4)
end

end
end
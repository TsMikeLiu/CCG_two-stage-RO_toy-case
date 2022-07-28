clear all; warning off;
%% 初始化
% 主问题
MP =  MPParams();
SP = SPParams();

[d,dconstrains] = UncertaintySet();                                             % 不确定性参数
[MPconstrains,MPFunc] =  MPconstrainsAndFunc(MP);                  %调用下文子程序求解
[SPconstrains,SPFunc] =  SPConstrainsAndFunc(MP,SP,d);              %调用下文子程序求解

theta = sdpvar(1);
constrains = [MPconstrains;   SPconstrains;   dconstrains;  theta >= SPFunc ];

opt = sdpsettings('verbose',1,'solver','cplex');
optimize(constrains,MPFunc + theta,opt);                   % 主问题求解
value(MP.Y)
disp('*********')
value(MP.Z)
LB = value(MPFunc + theta)
value(MPFunc)
% 子问题
KKT =  KKTParams();
[SPconstrains,pi,SPFunc] =  SPKKT(MP,SP,KKT,d);
SPKKTconstrains = [SPconstrains;dconstrains];

optimize(SPKKTconstrains,-SPFunc,opt);                  % 主问题求解
UB = value(MPFunc+SPFunc) 
%% CCG
n = 1;
cons = [];
a = 1;

while abs(UB-LB) >1e-5
    disp(['迭代',num2str(n+1),'次'])
    
    [SPconstrains,SPFunc] = SPConstrainsAndFunc(MP,SP,value(d));   
    value(d)
    
    if a ==1
        cons = [ theta >= SPFunc; SPconstrains];
    else
        cons = [SPconstrains];
    end
    
    constrains = [cons;MPconstrains];
    result = optimize(constrains,MPFunc + theta,opt);             % 主问题求解
    if result.problem == 0
        a =1;
    else
        a =0;
    end
    LB = max(LB, value(MPFunc + theta));
    
    [SPconstrains,pi,SPFunc] =  SPKKT(MP,SP,KKT,d); 
    constrains = [SPconstrains;dconstrains];
    optimize(constrains,-SPFunc,opt);                       % 子问题求解
    
 
    UB = min(UB,value(SPFunc+MPFunc))   
    n = n+1;
end

function MP =  MPParams()
MP.Y = binvar(1,3);
MP.Z = sdpvar(1,3);
end

function SP = SPParams()
SP.X = sdpvar(3,3,'full');
end

function [constrains,MPFunc] =  MPconstrainsAndFunc(MP)
f = [400,414,326];
d = [18,25,20];
constrains = [];
for i = 1:length(MP.Y)
    constrains = [constrains; MP.Z(i) <= 800 * MP.Y(i);];
end
constrains = [constrains; sum(MP.Z) >=772;];
MPFunc = sum(f.*MP.Y +d.*MP.Z);
end

function [SPconstrains,SPFunc] =SPConstrainsAndFunc(MP,SP,d)      % 子问题约束和目标函数
a = [22,33,24;33,23,30;20,25,27];
SPconstrains = [SP.X >= 0;
    sum(SP.X,1) >= d;
    sum(SP.X,2) <=MP.Z';];
SPFunc = sum(sum(SP.X .* a));
end

function KKT =  KKTParams()                                   % 定义对偶变量
b = [22,33,24,33,23,30,20,25,27]';
G = [1,1,1,0,0,0,0,0,0;
    0,0,0,1,1,1,0,0,0;
    0,0,0,0,0,0,1,1,1;
    -1,0,0,-1,0,0,-1,0,0;
    0,-1,0,0,-1,0,0,-1,0;
    0,0,-1,0,0,-1,0,0,-1];
KKT.pi =  sdpvar(size(G,1),1);                           
KKT.v = binvar(size(G,1),1);
KKT.W = binvar(length(b),1);
end

function [constrains,pi,SPKKTFunc] =  SPKKT(MP,SP,KKT,d)
a = sdpvar(9,1);
k = 1;
for i = 1:size(SP.X,1)
    for j = 1:size(SP.X,1)
        a(k) = SP.X(i,j);
        k = k+1;
    end
end
M = 10000;

pi = KKT.pi;
v = KKT.v;
W = KKT.W;

b = [22,33,24,33,23,30,20,25,27]';
G = [-1,-1,-1,0,0,0,0,0,0;
    0,0,0,-1,-1,-1,0,0,0;
    0,0,0,0,0,0,-1,-1,-1;
    1,0,0,1,0,0,1,0,0;
    0,1,0,0,1,0,0,1,0;
    0,0,1,0,0,1,0,0,1];
C = [-value(MP.Z)';d'];

constrains = [G *a >= C;         a >= 0;
                   G'* pi<= b;        pi >= 0];         

Mid = G * a -C;
for i = 1:size(G,1)
    constrains = [constrains;
                        pi(i) <= M * v(i);                      
                      Mid(i) <= M * (1 - v(i))   ];      
end

Mid = b -  G' * pi ;
for i = 1:length(b)
    constrains = [constrains;
                        a(i) <= M * W(i);
                     Mid(i) <= M * (1 - W(i))  ];
end

SPKKTFunc = b'* a;
end

function [d,constrains] = UncertaintySet()
g0 = sdpvar(1);
g1 = sdpvar(1);
g2 = sdpvar(1);
constrains = [ 0  <= g0 <= 1;
    0  <= g1 <= 1;
    0  <= g2<= 1;
    g1 +  g2 + g0 <= 1.8;
    g0 +  g1 <= 1.2];
d = sdpvar(1,3);
constrains = [constrains;
    d(1) == 206 + 40 * g0;
    d(2) == 274 + 40 * g1;
    d(3) == 220 + 40 * g2;];
end
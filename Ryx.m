R_Inner = [120 110 100 85 60 45];
AA = [0.84 0.87 0.9 0.93 0.96 1];
Spur_Length = [1.1 3.2 3.2 2.7 2.9 2.9];
dist = linspace(0,16,160);
EI = zeros(6,2);
Q = zeros(160,1);
M = zeros(160,1);
d_theta = zeros(160,1);
theta = zeros(160,1);
d_delta = zeros(160,1);
delta = zeros(160,1);
tan_theta = zeros(160,1);
kk = 1;
Uinf = 7.2;
BestPower = 100000;
count = 0;



% introduction
disp("これは桁構造を考慮した最適循環分布を求めるプログラムです。")
disp("スパン長を入力[m]")
span = input(">>");

for beta = AA
    for p1 = 7:11
        for p2 = 7:11
            for p3 = 7:11
                for p4 = 6:10
               
                
                    % 初期ベクトルを定義
                    V_V = [p1 p2 p3 p4 4 3];
                    
                    % 初期桁設定
                    k=1;
                    for i = R_Inner
                        A = Enumerate_EI_Weight(i);
                        EI(k,1) = A(V_V(k),1);
                        EI(k,2) = A(V_V(k),2);
                        k = k + 1;
                    end
                    
                    %機体重量計算
                    Weight = 2 * Spur_Length * EI(:,2);% 桁重量
                    SpurWeight = Weight;
                    
                    Weight = 20 + Weight + 8 + 60; %フレーム等＋桁＋二次構造
                    Weight = Weight * 9.8; %Nに変換
                    
                    %循環分布の取得
                    Gamma = TR(Weight,span,Uinf,beta); 
                    
                    %曲げモーメントの算出
                    LL = Gamma * 1.154 * Uinf; %LL計算
                    fdx = LL * 0.1 ;
                    for i = 1:160
                        Q(i) = sum(fdx(i:160));
                    end
                    Qdx = Q * 0.1;
                    for i = 1:160
                        M(i) = sum(Qdx(i:160));
                    end
                    M = M * 1000 / 9.8;
                    
                    %局所曲率計算
                    for i = 1:160
                        p = EI(SpurPosition(i))*10^10;
                        d_theta(i) = M(i) / p * 0.1 * 1000;
                    end
                    for i = 1:160
                        theta(i) = sum(d_theta(1:i) * 180 / pi);
                    end
                    d_delta = tan(theta*pi/180) * 0.1;
                    for i = 1:160
                        delta(i) = sum(d_delta(1:i)) * 1000;
                    end
                    
    
    
                    S = Weight * 2 / (1.154 * Uinf^2 * 1.05);
                    Power = (TR2(Weight,span,Uinf,beta) + 4.5 + 1/2 * 1.154 * Uinf^2 * S * 0.009) * Uinf;
                    disp(BestPower/0.85/0.85)
                    if Power < BestPower && delta(160) < 2400
                        BestPower = Power;
                        disp("パワー")
                        disp(Power/0.85/0.85)
                        disp("たわみ量")
                        disp(delta(160))
                        disp("機体重量")
                        disp(Weight/9.8)
                        disp(V_V)
                        disp(beta)
                        count = 0;
                        disp("----------------------------")
                    end
                    %if mod(kk/20, 1) == 0
                    %    disp(kk)
                    %end
                    kk = kk + 1;
                    
                end
            end
        end
    end
end
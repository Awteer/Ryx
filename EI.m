function y=EI(PlyNumber,pp,R)   %積層構成ベクトル

    thickness = zeros(11,1);    %各層の厚みの縦ベクトル
    Inner = zeros(11,1);    %各層の内径の縦ベクトル
    Outer = zeros(11,1);    %各層の外径の縦ベクトル
    Ix = zeros(11,1);   %各層の上下断面二次M
    Iy = zeros(11,1);   %各層の前後断面二次M
    EIx = zeros(11,1);  %各層の上下剛性
    EIy = zeros(11,1);  %各層の上下剛性
    Inner(1) = R;   %最内径=マンドレル径

    %プリプレグの剛性値,それぞれ24t,40tの0,45,90度の値
    Material = [13000 1900 900 22000 1900 800];

    %%各層の厚みの計算
    thickness = PlyNumber .* 0.111;
    thickness(1) = 0.125;
    thickness(11) = 0.125;

    %%各層の内径・外径の計算
    for i=2:11
        Inner(i) = Inner(i-1) + 2*thickness(i-1);
        Outer(1) = Inner(2);
        Outer(i) = Inner(i) + 2*thickness(i);
    end

    %%断面二次モーメントの計算(上下・前後)

    Ix = pi/64 * (Outer.^4 - Inner.^4) .* (1 - cos(pp .* pi/180));
    Iy = 1/64 * (Outer.^4 - Inner.^4) .* (2*pp.*pi/180 + sin(2*pp.*pi/180)); 
    
    %弾性率計算

    for i=1:11
        if i==1 || i==11
            EIx(i) = Ix(i) * Material(3);
            EIy(i) = Iy(i) * Material(3);
        elseif i==2
            EIx(i) = Ix(i) * Material(5);
            EIy(i) = Iy(i) * Material(5);
        else
            EIx(i) = Ix(i) * Material(4);
            EIy(i) = Iy(i) * Material(4);
        end
    end
    
    EIX = 0;
    EIY = 0;
    for i=1:11
        EIX = EIX + EIx(i);
        EIY = EIY + EIy(i);
    end
    
    y = EIX;
end
    
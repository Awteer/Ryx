function y=SpurWeight(PlyNumber,pp,R)   %積層構成ベクトル

    thickness = zeros(11,1);    %各層の厚みの縦ベクトル
    weight = zeros(11,1);    %各層の1m毎の重量の縦ベクトル
    rS = zeros(11,1);
    Inner = zeros(11,1);    %各層の内径の縦ベクトル
    Outer = zeros(11,1);    %各層の外径の縦ベクトル
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

    %%各層の断面積
    rS = pi/4 * (Outer.^2 - Inner.^2) .* pp/90;

    %%各層の重量
    for i=2:10
        weight(i) = 0.001559  * 1000 * rS(i);
    end
    weight(1) = 0.001496 * 1000 * rS(1);
    weight(11) = 0.001496 * 1000 * rS(11);

    omosa = 0;
    for i=1:11
        omosa = omosa + weight(i);
    end
    
    y = omosa;
    
    end
    
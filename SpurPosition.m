function y = SpurPosition(l)
ppp = 0;

Spur_Length = [1.1 4.3 7.5 10.2 13.1 16];
Spur_Length = Spur_Length * 10;

for k = 2:6
    if l >= Spur_Length(k-1) && l <= Spur_Length(k)
        ppp = k;
    end
    if l <= 11
        ppp = 1;
    end
end

y = ppp;

end

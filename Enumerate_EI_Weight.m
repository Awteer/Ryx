function y = Enumerate_EI_Weight(R)

PlyNumber = [1 2 6 1 0 0 0 0 1 1 1];
pp = [90 90 90 50 45 40 35 30 25 20 90];
AA = zeros(2000,1);
WW = zeros(2000,1);

numbers = [1 2];
count = 1;

for i=1:9
    for j=numbers
        for k=numbers
            for l=numbers
                for m=numbers
                    for n=numbers
                        for o=numbers
                            for p=numbers
                                PlyNumber(10)=p;
                                AA(count)=round(EI(PlyNumber-1 , pp, R)/10^10, 9);
                                WW(count)=round(SpurWeight(PlyNumber-1, pp, R)/1000, 9);
                                count = count + 1;
                            end
                            PlyNumber(9)=o;
                        end
                        PlyNumber(8)=n;
                    end
                    PlyNumber(7)=m;
                end
                PlyNumber(6)=l;
            end
            PlyNumber(5)=k;
        end
        PlyNumber(4)=j;
    end
    PlyNumber(3)=i;
end

D = cat(2,AA,WW);
D = unique(D,'stable','rows');
D = sortrows(D);

y = D;

end

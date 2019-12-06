% test s^0/s^e > s_i^o/s_e_i


N = 8;
p = 2;

go = 0;
for i=1:2:N
    i
    rowSum = hilbertRowSum(N,p,i);
    go = go + rowSum;
end
so = abs(go)
ge =0;
for i = 2:2:N
    rowSum = hilbertRowSum(N,p,i);
    ge = ge + rowSum;
end
se = abs(ge)

row = 1;

gio = 0;
for j = 1:2:N
    entry = hilbertEntryVal(N,p,row,j);
    gio = gio + entry;
end
sio = abs(gio)

gie = 0;
for j = 2:2:N
    entry = hilbertEntryVal(N,p,row,j);
    gie = gie + entry;
end 
sie = abs(gie)

r = so/se

ind = sio/sie

r/ind

assert( r < ind)

splus = se + so;
sminus = se - so;

siplus = sie + sio;
siminus = sie - sio;

(splus - sminus)/ (splus + sminus)

(siplus - siminus)/(siplus + siminus)

(1-sminus/splus)/(1+sminus/splus)

(1-siminus/siplus)/(1+siminus/siplus)

splus/sminus

siplus/siminus
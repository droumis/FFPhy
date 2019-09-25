pushd /data19a/mkarlsso/alex
%clusterdayprocess('alex', '/data/mkarlsso/Ale/', 'ale', [1:8], [2 3 4 5])
%cd ../bond;
%clusterdayprocess('bond', '/data/mkarlsso/Bon/', 'bon', [3:10], [2 3 4 5])
cd ../conley
clusterdayprocess('conley', '/data/mkarlsso/Con/', 'con', [1:6], [2 3 4 5])
% check day 2
cd ../dudley
clusterdayprocess('dudley', '/data/mkarlsso/Dud/', 'dud', [1:6], [2 3 4 5])
cd ../frank
clusterdayprocess('frank', '/data/mkarlsso/Fra/', 'fra', [1:12], [2 3 4 5])
cd ../miles
clusterdayprocess('miles', '/data/mkarlsso/Mil/', 'mil', [1:5], [2 3 4 5])
cd ../nine
clusterdayprocess('nine', '/data/mkarlsso/Nin/', 'nin', [2], [2 3 4 5])
cd ../ten
clusterdayprocess('ten', '/data/mkarlsso/Ten/', 'ten', [1:7], [2 3 4 5])
popd

addpath('toolbox/');
mynorm = @(a)norm(a(:));
aa = staggered([3,3,3]);
aa.M{1} = rand(size(aa.M{1}));
aa.M{2} = rand(size(aa.M{2}));
aa.M{3} = rand(size(aa.M{3}));
old_aa = [1 1 1];
tans =[];
for ii = 1:100
%     fprintf("%d %d %d %d %d %d\n",mynorm(aa.M{1}),mynorm(aa.M{2}),mynorm(aa.M{3}),mynorm(aa.M{1})-old_aa(1),mynorm(aa.M{2})-old_aa(2),mynorm(aa.M{3})-old_aa(3));
    old_aa = [mynorm(aa.M{1}),mynorm(aa.M{2}),mynorm(aa.M{3})];
    tt = incompressible_R(aa);
    aa = aa+tt*(1/(L*4*sqrt(2)));
    tans =[tans;mynorm(div_x(aa))];
    aa.M{1}(:,:,1);
end
plot(tans)

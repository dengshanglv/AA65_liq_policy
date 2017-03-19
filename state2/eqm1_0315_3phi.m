clear;
global equ irate delta zmc zlc zmd zld yc yd q;

beta = 0.98;
gma = [0 0];
gma_new = [0 0];

pi = 0.5; % fix creditor's fraction.
omega = 30;

syms x1 x2 x3 y;
v = sqrt(x2);
u = log(x1);
udif = diff(u);
vdif = diff(v);
use = @(symfun, para)(cell2sym(arrayfun(@subs, repmat(symfun, [1 length(para)]), para, ...
    'UniformOutput',false)));

% theta = sym('theta',[1 2],'positive');
theta = [1 1];

pr_mat = [0.7 0.3; 0.3 0.7];
phi = [0.8 1.25];

q = sym('q', [1 2], 'real');
qnew = sym('qnew', [1 2], 'real');
% qnew = q;

delta = sym('delta', [1 2], 'real');
irate = sym('irate', [1 2], 'real');
zmc = sym('zmc', [1 2], 'positive');
zmd = sym('zmd', [1 2], 'positive');
zlc = sym('zlc', [1 2], 'real');
zld = sym('zld', [1 2], 'real');
mc = sym('mc', [1 2], 'positive');
md = sym('md', [1 2], 'positive');
yc = sym('yc', [1 2], 'positive');
yd = sym('yd', [1 2], 'positive');

eq1 = pi.*zmc + (1-pi).*zmd - (1+delta).*q;
eq2 = delta.*q + pi.*zlc + (1-pi).*zld;

xnewd = kron(omega.*theta, [1 1]) + kron((zmd+zld.*(1.+irate))./q, qnew./(1+gma_new));
xnewc = kron((zmc+zlc.*(1.+irate))./q, qnew./(1+gma_new));

tmp = (use(udif, xnewc)+use(udif,subs(xnewc, zmc, zmc-mc))).*kron((1+irate)./q, qnew./(1+gma_new));
tmp = reshape(tmp, [2 2])';
right = beta*[dot(tmp(1,:), pr_mat(1,:)), dot(tmp(2,:), pr_mat(2,:))];

eq3 = 2*use(udif, omega-zmc-zlc) - right;  % Zlc 

tmp2 = (use(udif, xnewc)+kron(phi.*use(vdif, yc), ones(1,...
    2)).*use(udif,subs(xnewc, zmc, zmc+mc))).*kron(1./q, qnew./(1+gma_new));
tmp2 = reshape(tmp2, [2 2])';
right2 = beta*[dot(tmp2(1,:), pr_mat(1,:)), dot(tmp2(2,:), pr_mat(2,:))];

eq4 = 2*use(udif, omega-zmc-zlc) - right2; % Zmc

tmp = (use(udif, xnewd)+use(udif,subs(xnewd, zmd, zmd-md))).*kron((1+irate)./q, qnew./(1+gma_new));
tmp = reshape(tmp, [2 2])';
right = beta*[dot(tmp(1,:), pr_mat(1,:)), dot(tmp(2,:), pr_mat(2,:))];

eq5 = 2*use(udif, -zmd-zld) - right;       % Zld

tmp2 = (use(udif, xnewd)+kron(phi.*use(vdif, yd), ones(1,...
    2)).*use(udif,subs(xnewd, zmd, zmd+md))).*kron(1./q, qnew./(1+gma_new));
tmp2 = reshape(tmp2, [2 2])';
right2 = beta*[dot(tmp2(1,:), pr_mat(1,:)), dot(tmp2(2,:), pr_mat(2,:))];

eq6 = 2*use(udif, -zmd-zld) - right2;      % Zmd

tmp3 = (use(u, subs(xnewc, zmc, zmc+mc))-use(u,xnewc)).*kron(phi, ones(1,2));
tmp3 = reshape(tmp3, [2 2])';
right3 = beta*[dot(tmp3(1,:), pr_mat(1,:)), dot(tmp3(2,:), pr_mat(2,:))];
eq7 = yc - right3;                         % yc

tmp3 = (use(u, subs(xnewd, zmd, zmd+md))-use(u,xnewd)).*kron(phi, ones(1,2));
tmp3 = reshape(tmp3, [2 2])';
right3 = beta*[dot(tmp3(1,:), pr_mat(1,:)), dot(tmp3(2,:), pr_mat(2,:))];
eq8 = yd - right3;                         % yd


tmp2 = (kron(phi.*use(vdif, yc), ones(1,2)).*use(udif,subs(xnewc, zmc, zmc+mc))...
    - use(udif,subs(xnewc, zmc, zmc-mc))).*kron([1 1], qnew./(1+gma_new));
tmp2 = reshape(tmp2, [2 2])';
right2 = beta*[dot(tmp2(1,:), pr_mat(1,:)), dot(tmp2(2,:), pr_mat(2,:))];

eqn4 = right2;                             % FOC for mc, when i==0

tmp2 = (kron(phi.*use(vdif, yd), ones(1,2)).*use(udif,subs(xnewd, zmd, zmd+md))...
    - use(udif,subs(xnewd, zmd, zmd-md))).*kron([1 1], qnew./(1+gma_new));
tmp2 = reshape(tmp2, [2 2])';
right2 = beta*[dot(tmp2(1,:), pr_mat(1,:)), dot(tmp2(2,:), pr_mat(2,:))];

eqn6 = right2;                             % FOC for md, when i==0

tmp4 = use(u, xnewc)+use(u, subs(xnewc, zmc, zmc-mc));
tmp4 = reshape(tmp4, [2 2])';
right4 = beta*[dot(tmp4(1,:), pr_mat(1,:)), dot(tmp4(2,:), pr_mat(2,:))];

welc = use(u, omega - zmc - zlc) + 0.5*(right4 + use(v,yc));

tmp5 = use(u, xnewd)+use(u, subs(xnewd, zmd, zmd-md));
tmp5 = reshape(tmp5, [2 2])';
right4 = beta*[dot(tmp5(1,:), pr_mat(1,:)), dot(tmp5(2,:), pr_mat(2,:))];

weld = use(u,  - zmd - zld) + 0.5*(right4 + use(v,yd));

lbondd = omega.*theta.*q + (zmd+zld.*(1.+irate));
lbondc = (zmc+zlc.*(1.+irate));
lbond = [subs(lbondd, zmd, zmd-md), subs(lbondc, zmc, zmc-mc)];

%% if i == 0
conc = sym('conc', [1 2], 'positive');
cond = sym('cond', [1 2], 'positive');

eqn1 = eq1+eq2;
noweq = [eqn1, eq3, eqn4, eq5, eqn6, eq7, eq8];
noweq = subs(noweq,irate,[0,0]);
noweq = subs(noweq, qnew, q);
noweq = subs(noweq, [zmc, zmd], [omega-zlc-conc, -zld-cond]);
noweq = simplify(noweq);

[zconc1, zconc2] = solve(noweq([5,6]), conc);

nconc1 = sidx(solve(noweq(3), conc(1)),1);
nconc2 = sidx(solve(noweq(4), conc(2)),1);

noweq(3) = conc(1) - nconc1;
noweq(4) = conc(2) - nconc2;
noweq(5) = conc(1) - zconc1;
noweq(6) = conc(2) - zconc2;

objf = @(x)(sum(double(subs(noweq, [conc, cond, q, yc, yd, mc, md], x)).^2));
x0 = [10,10,2,2,10,10,2,2,1,1,1,1,1,1];

s00 = [sconc1,sconc2,scond1,scond2,sq1,sq2,syc1,syc2,syd1,syd2,smc1,smc2,smd1,smd2];

mindelta1 = (pi*smc1+(1-pi)*smd1)/sq1-1;
mindelta2 = (pi*smc2+(1-pi)*smd2)/sq2-1;

welcr = simplify(subs(welc, [zmc, irate, qnew], ...
    [omega-zlc-conc, 0, 0, q]));
welcr = double(subs(welcr,[conc, cond, q, yc, yd, mc, md],s00));

weldr = simplify(subs(weld, [zmd, irate, qnew], ...
    [-zld-cond, 0, 0, q]));
weldr = double(subs(weldr,[conc, cond, q, yc, yd, mc, md],s00));
wel00 = [welcr, weldr];

%% i1 = 0, i2 > 0
conc = sym('conc', [1 2], 'positive');
cond = sym('cond', [1 2], 'positive');

eqa(1) = eq1(1) + eq2(1);
eqa([2 3]) = [eq1(2), eq2(2)];
eqa([4 5]) = eq3;
eqa(6) = eq4(2);
eqa([7 8]) = eq5;
eqa(9) = eq6(2);
eqa([10 11 12 13]) = [eq7 eq8];
eqa([14 15]) = [eqn4(1) eqn6(1)];

eqa(16:19) = [welc weld];

eqa = subs(eqa,irate(1),0);
eqa = subs(eqa, qnew, q);
eqa = subs(eqa, [zmc(1), zmd(1)], ...
    [omega-zlc(1)-conc(1), -zld(1)-cond(1)]);
eqa = simplify(eqa);
eqa = subs(eqa, mc(2), zmc(2));
eqa = subs(eqa, md(2), zmd(2));
eqa = subs(eqa, delta(2), mindelta2*0.99);
eqa = simplify(eqa);
welr2 = eqa(16:19);
eqa = eqa(1:15);

nconc1 = sidx(solve(eqa(4), conc(1)), 1);  % smaller
nzlc2 = sidx(solve(eqa(5), zlc(2)), 2);    % bigger
nzmc2 = sidx(solve(eqa(6), zmc(2)), 2);    % bigger, so positive
ncond1 = sidx(solve(eqa(7), cond(1)), 1);  % smaller
nzld2 = sidx(solve(eqa(8), zld(2)), 1);    % smaller absolute value, bigger
nzmd2 = sidx(solve(eqa(9), zmd(2)), 2);    % bigger, so positive
[nmc1, nmd1] = solve(eqa(14:15), [mc(1), md(1)]);

eqa(4) = nconc1 - conc(1);
eqa(5) = nzlc2 - zlc(2);
eqa(6) = nzmc2 - zmc(2);
eqa(7) = ncond1 - cond(1);
eqa(8) = nzld2 - zld(2);
eqa(9) = nzmd2 - zmd(2);
eqa(14:15) = [nmc1, nmd1] - [mc(1), md(1)];

[s1irate2, s1conc1, s1cond1, s1zmc2, s1zlc2, s1zmd2, s1zld2, s1q1, s1q2, s1yc1, s1yc2, ...
    s1yd1, s1yd2, s1mc1, s1md1] = ...
vpasolve(eqa, [irate(2), conc(1), cond(1), zmc(2), zlc(2), ...
    zmd(2), zld(2), q, yc, yd, mc(1), md(1)], [ 0.0039006652460165745897872645020515, ...
    479.22546708966659657638826415176, 345.84658079698541667935255926502, ...
    76.482143074030296262525194436873, 444.21694558946432881612536954185, ...
    55.216808490540469271676888877504, -401.25139338694291200123419116813, ...
    0.80186468594211136183315878028401, 0.65888452755638439372584421711419, ...
    0.13782646603130568239102037822635, 0.13662346373651220077462728739497, ...
    0.13782646603130568239102037822635, 0.13662346373651220077462728739497, ...
    76.95816214879973223169501490536, 55.539029269916441424211715011069]);

s1 = [s1irate2, s1conc1, s1cond1, s1zmc2, s1zlc2, s1zmd2, s1zld2, s1q1, s1q2, s1yc1, s1yc2, ...
    s1yd1, s1yd2, s1mc1, s1md1];

wel01 = subs(welr2, [irate(2), conc(1), cond(1), zmc(2), zlc(2), ...
    zmd(2), zld(2), q, yc, yd, mc(1), md(1)],s1);

%%  i1>0 i2=0
conc = sym('conc', [1 2], 'positive');
cond = sym('cond', [1 2], 'positive');
clear eqa;

eqa(1) = eq1(2) + eq2(2);
eqa([2 3]) = [eq1(1), eq2(1)];
eqa([4 5]) = eq3;
eqa(6) = eq4(1);
eqa([7 8]) = eq5;
eqa(9) = eq6(1);
eqa([10 11 12 13]) = [eq7 eq8];
eqa([14 15]) = [eqn4(2) eqn6(2)];

eqa(16:19) = [welc weld];

eqa = subs(eqa,irate(2),0);
eqa = subs(eqa, qnew, q);
eqa = subs(eqa, [zmc(2), zmd(2)], ...
    [omega-zlc(2)-conc(2), -zld(2)-cond(2)]);
eqa = simplify(eqa);
eqa = subs(eqa, mc(1), zmc(1));
eqa = subs(eqa, md(1), zmd(1));
eqa = simplify(eqa);
welr2 = eqa(16:19);
eqa = eqa(1:15);

nconc2 = sidx(solve(eqa(5), conc(2)), 1);  % smaller
nzlc1 = sidx(solve(eqa(4), zlc(1)), 2);    % bigger
nzmc1 = sidx(solve(eqa(6), zmc(1)), 1);    % bigger, so positive
ncond2 = sidx(solve(eqa(8), cond(2)), 1);  % smaller
nzld1 = sidx(solve(eqa(7), zld(1)), 1);    % smaller absolute value, bigger
nzmd1 = sidx(solve(eqa(9), zmd(1)), 2);    % bigger, so positive
[nmc2, nmd2] = solve(eqa(14:15), [mc(2), md(2)]);

eqa(5) = nconc2 - conc(2);
eqa(4) = nzlc1 - zlc(1);
eqa(6) = nzmc1 - zmc(1);
eqa(8) = ncond2 - cond(2);
eqa(7) = nzld1 - zld(1);
eqa(9) = nzmd1 - zmd(1);
eqa(14:15) = [nmc2, nmd2] - [mc(2), md(2)];

eqa_now = subs(eqa, delta(1), mindelta1*0.99);
[s2irate1, s2conc2, s2cond2, s2zmc1, s2zlc1, s2zmd1, s2zld1, s2q1, s2q2, s2yc1, s2yc2, ...
    s2yd1, s2yd2, s2mc2, s2md2] = ...
vpasolve(eqa_now, [irate(1), conc(2), cond(2), zmc(1), zlc(1), ...
    zmd(1), zld(1), q, yc, yd, mc(2), md(2)], [0.0038921699826169293898369469572849, ...
    479.22546708966659657638826415176, 346.08880426572194863144094449406, ...
    76.483175138853248614440850560023, 444.2160775221868985885810114003, ...
    55.178801393033537316886456029833, -400.97041761365623209728689812491, ...
    0.8047313320424420129684795270819, 0.65653060470019219058012764285834, ...
    0.13662607015679904198760894043726, 0.13782646603130568239102037822635, ...
    0.13662607015679904198760894043726, 0.13782646603130568239102037822635, ...
    76.95816214879973223169501490536, 55.577927605383623236611918035265]);

s2 = [s2irate1, s2conc2, s2cond2, s2zmc1, s2zlc1, s2zmd1, s2zld1, s2q1, s2q2, s2yc1, s2yc2, ...
    s2yd1, s2yd2, s2mc2, s2md2];

wel10 = subs(welr2, [irate(1), conc(2), cond(2), zmc(1), zlc(1), ...
    zmd(1), zld(1), q, yc, yd, mc(2), md(2)],s2);

%%  i1>0 i2>0
equ = [eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, welc, weld];
equ = subs(equ, qnew, q);
equ = subs(equ, qnew, q);
equ = subs(equ, mc, zmc);
equ = subs(equ, md, zmd);
% equ = subs(equ, delta, [0.2,0.2]);

nzlc1 = sidx(solve(equ(5), zlc(1)), 2);     % bigger
nzlc2 = sidx(solve(equ(6), zlc(2)), 2);
nzmc1 = sidx(solve(equ(7), zmc(1)), 1);
nzmc2 = sidx(solve(equ(8), zmc(2)), 2);
nzld1 = sidx(solve(equ(9), zld(1)), 1);     % bigger
nzld2 = sidx(solve(equ(10), zld(2)), 1);
nzmd1 = sidx(solve(equ(11), zmd(1)), 2);
nzmd2 = sidx(solve(equ(12), zmd(2)), 2);

equ(5) = nzlc1 - zlc(1);
equ(6) = nzlc2 - zlc(2);
equ(7) = nzmc1 - zmc(1);
equ(8) = nzmc2 - zmc(2);
equ(9) = nzld1 - zld(1);
equ(10) = nzld2 - zld(2);
equ(11) = nzmd1 - zmd(1);
equ(12) = nzmd2 - zmd(2);

equ_now = subs(equ, irate, [0.5,0.3]);
wel11 = equ_now(17:20);
equ_now = equ_now(1:16);

[s3delta1, s3delta2, s3zmc1, s3zmc2, s3zlc1, s3zlc2, s3zmd1, s3zmd2, s3zld1, s3zld2, ...
    s3yc1, s3yc2, s3yd1, s3yd2, s3q1, s3q2] =...
    vpasolve(equ_now, [delta, zmc, zlc, zmd, zld, yc, yd, q], [ 65.477241797540354974388167328505, 80.217174036613074319664157255025, 76.835717286781466954304966552849, 76.835717286781466954304966552849, 443.91944687184207252598086515369, 443.91944687184207252598086515369, 83.105191145674486957941450791232, 83.175059966588357538713957956551, -601.45440547684333597546477670042, -601.96006474511285253996055168404, 0.13751686071987961921511259675287, 0.13751686071987961921511259675287, 0.13751686071987961921511259675287, 0.13751686071987961921511259675287, 1.2029749137273452313812528986736, 0.98507969004952223951961898952309]);

s3 = [s3delta1, s3delta2, s3zmc1, s3zmc2, s3zlc1, s3zlc2, s3zmd1, s3zmd2, s3zld1, s3zld2, ...
    s3yc1, s3yc2, s3yd1, s3yd2, s3q1, s3q2];

wel11 = double(subs(wel11, [delta, zmc, zlc, zmd, zld, yc, yd, q], s3));
y = sum(wel11)/4;

%%
[gg1, gg2] = meshgrid(0.01:0.05:0.3);
resu = zeros([length(gg1) length(gg1)]);
delta1res = zeros([length(gg1) length(gg1)]);
delta2res = zeros([length(gg1) length(gg1)]);

for ii = 1:length(gg1)
    for cc = 1:length(gg1)
        [resu(ii, cc), tmp] = pay11([gg1(ii,cc), gg2(ii,cc)]);
        delta1res(ii, cc) = tmp(1);
        delta2res(ii, cc) = tmp(2);
        disp([ii,cc]);
        disp([resu(ii, cc), tmp]);
    end
end
surf(resu);

%%
vpa(wel00, 6)
vpa(wel01, 6)
vpa(wel10, 6)
vpa(wel11, 6)

vpa(sum(wel00), 10)
vpa(sum(wel01), 10)
vpa(sum(wel10), 10)
vpa(sum(wel11), 10)
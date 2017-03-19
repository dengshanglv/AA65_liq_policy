clear;
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
theta = [0.8 1.25];

pr_mat = [0.7 0.3; 0.3 0.7];
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

xnewd = omega + kron((zmd+zld.*(1.+irate))./q, qnew./(1+gma_new));
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

[zconc1, zconc2, zcond1, zcond2] = solve(noweq([5,6,9,10]), [conc, cond]);

[nq1, nq2] = solve(noweq(1:2),q);
noweq = subs(noweq,q,[nq1,nq2]);

nconc1 = solve(noweq(3), conc(1));
nconc1 = nconc1(1);
nconc2 = solve(noweq(4), conc(2));
nconc2 = nconc2(1);
noweq = subs(noweq,conc,[nconc1, nconc2]);

ncond1 = solve(noweq(7), cond(1));
ncond1 = ncond1(1);
ncond2 = solve(noweq(8), cond(2));
ncond2 = ncond2(1);
noweq = subs(noweq,cond,[ncond1, ncond2]);

[nyc1, nyc2, nyd1, nyd2] = solve(noweq(11:14), [yc, yd]);
noweq = subs(noweq,[yc, yd],[nyc1, nyc2, nyd1, nyd2]);

eeq = [zconc1, zconc2, zcond1, zcond2] - [nconc1, nconc2, ncond1, ncond2];
zq1 = subs(nq1, [conc(1), cond(1)], [nconc1, ncond1]);
zq2 = subs(nq2, [conc(2), cond(2)], [nconc2, ncond2]);
eeq = subs(eeq, [yc, yd, q], [nyc1, nyc2, nyd1, nyd2, zq1, zq2]);

[smc1, smc2, smd1, smd2] = vpasolve(eeq, [mc, md],[ 1.9861012399000304028513282058772, 2.6862280239614762957398191825816, 1.8389866780949575711179509931154, 2.5522448490444860789513604239019]);

sq1 = subs(zq1, [mc, md], [smc1, smc2, smd1, smd2]);
sq2 = subs(zq2, [mc, md], [smc1, smc2, smd1, smd2]);
sconc1 = subs(nconc1, [mc, md], [smc1, smc2, smd1, smd2]);
sconc2 = subs(nconc2, [mc, md], [smc1, smc2, smd1, smd2]);
scond1 = subs(ncond1, [mc, md], [smc1, smc2, smd1, smd2]);
scond2 = subs(ncond2, [mc, md], [smc1, smc2, smd1, smd2]);
syc1 = subs(nyc1, [mc, md], [smc1, smc2, smd1, smd2]);
syc2 = subs(nyc2, [mc, md], [smc1, smc2, smd1, smd2]);
syd1 = subs(nyd1, [mc, md], [smc1, smc2, smd1, smd2]);
syd2 = subs(nyd2, [mc, md], [smc1, smc2, smd1, smd2]);
s00 = [sconc1,sconc2,scond1,scond2,sq1,sq2,syc1,syc2,syd1,syd2,smc1,smc2,smd1,smd2];

mindelta1 = (pi*smc1+(1-pi)*smd1)/sq1-1;
mindelta2 = (pi*smc2+(1-pi)*smd2)/sq2-1;
mindelta = [mindelta1, mindelta2];

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
eqa = simplify(eqa);
welr2 = eqa(16:19);
eqa = eqa(1:15);

nconc1 = sidx(solve(eqa(4), conc(1)), 1);  % smaller
nzlc2 = sidx(solve(eqa(5), zlc(2)), 2);    % bigger
nzmc2 = sidx(solve(eqa(6), zmc(2)), 1);    % bigger, so positive
ncond1 = sidx(solve(eqa(7), cond(1)), 1);  % smaller
nzld2 = sidx(solve(eqa(8), zld(2)), 1);    % smaller absolute value, bigger
nzmd2 = sidx(solve(eqa(9), zmd(2)), 1);    % bigger, so positive
[nmc1, nmd1] = solve(eqa(14:15), [mc(1), md(1)]);

eqa(4) = nconc1 - conc(1);
eqa(5) = nzlc2 - zlc(2);
eqa(6) = nzmc2 - zmc(2);
eqa(7) = ncond1 - cond(1);
eqa(8) = nzld2 - zld(2);
eqa(9) = nzmd2 - zmd(2);
eqa(14:15) = [nmc1, nmd1] - [mc(1), md(1)];

neqa = subs(eqa, delta(2), mindelta2*0.99);
[s1irate2, s1conc1, s1cond1, s1zmc2, s1zlc2, s1zmd2, s1zld2, s1q1, s1q2, s1yc1, s1yc2, ...
    s1yd1, s1yd2, s1mc1, s1md1] = ...
vpasolve(neqa, [irate(2), conc(1), cond(1), zmc(2), zlc(2), ...
    zmd(2), zld(2), q, yc, yd, mc(1), md(1)], [ 0.0027282396129125750667980485072276, 14.615367029293011283748301888421, 13.532777041965388195721546662308, 2.675647145383197273889626559612, 12.913443687029857926566586819173, 2.5408124600930483990693288252001, -16.225508380275794071450584432879, 0.92592796437080026026507572463579, 0.95219745611515476403747888555312, 0.09519147371668585614467830329925, 0.19363645190338261856974977731699, 0.09519147371668585614467830329925, 0.19363645190338261856974977731699, 1.9861012399000304028513282058772, 1.8389866780949575711179509931154]);

s1 = [s1irate2, s1conc1, s1cond1, s1zmc2, s1zlc2, s1zmd2, s1zld2, s1q1, s1q2, s1yc1, s1yc2, ...
    s1yd1, s1yd2, s1mc1, s1md1];

wel01 = subs(welr2, [irate(2), conc(1), cond(1), zmc(2), zlc(2), ...
    zmd(2), zld(2), q, yc, yd, mc(1), md(1)],s1);
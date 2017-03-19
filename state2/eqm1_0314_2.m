clear;
beta = 1;
gma = [0 0];
gma_new = [0 0];

pi = 0.4; % fix creditor's fraction.
omega = 1000;

tta = [0.8 0.9; 0.8 1; 0.8 1.1; 0.8 1.2; 0.9 0.9; 0.9 1; 0.9 1.1; ...
    0.9 1.2; 1 1; 1 1.1; 1 1.2; 1.1 1.1];

syms x1 x2 x3 y;
v = sqrt(x2);
u = log(x1);
udif = diff(u);
vdif = diff(v);
use = @(symfun, para)(cell2sym(arrayfun(@subs, repmat(symfun, [1 length(para)]), para, ...
    'UniformOutput',false)));

theta = sym('theta',[1 2],'positive');
% theta = [0.9 1.3];
pr_mat = [0.7 0.3; 0.3 0.7];
phi = [1,1];

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

xnewd = kron(omega.*theta + (zmd+zld.*(1.+irate))./q, qnew./(1+gma_new));
xnewc = kron((zmc+zlc.*(1.+irate))./q, qnew./(1+gma_new));

tmp = (use(udif, xnewc)+use(udif,subs(xnewc, zmc, zmc-mc))).*kron((1+irate)./q, qnew./(1+gma_new));
tmp = reshape(tmp, [2 2])';
right = beta*[dot(tmp(1,:), pr_mat(1,:)), dot(tmp(2,:), pr_mat(2,:))];

eq3 = 2*use(udif, omega-zmc-zlc) - right;

tmp2 = (use(udif, xnewc)+kron(phi.*use(vdif, yc), ones(1,...
    2)).*use(udif,subs(xnewc, zmc, zmc+mc))).*kron(1./q, qnew./(1+gma_new));
tmp2 = reshape(tmp2, [2 2])';
right2 = beta*[dot(tmp2(1,:), pr_mat(1,:)), dot(tmp2(2,:), pr_mat(2,:))];

eq4 = 2*use(udif, omega-zmc-zlc) - right2;

tmp = (use(udif, xnewd)+use(udif,subs(xnewd, zmd, zmd-md))).*kron((1+irate)./q, qnew./(1+gma_new));
tmp = reshape(tmp, [2 2])';
right = beta*[dot(tmp(1,:), pr_mat(1,:)), dot(tmp(2,:), pr_mat(2,:))];

eq5 = 2*use(udif, -zmd-zld) - right;

tmp2 = (use(udif, xnewd)+kron(phi.*use(vdif, yd), ones(1,...
    2)).*use(udif,subs(xnewd, zmd, zmd+md))).*kron(1./q, qnew./(1+gma_new));
tmp2 = reshape(tmp2, [2 2])';
right2 = beta*[dot(tmp2(1,:), pr_mat(1,:)), dot(tmp2(2,:), pr_mat(2,:))];

eq6 = 2*use(udif, -zmd-zld) - right2;

tmp3 = (use(u, subs(xnewc, zmc, zmc+mc))-use(u,xnewc)).*kron(phi, ones(1,2));
tmp3 = reshape(tmp3, [2 2])';
right3 = beta*[dot(tmp3(1,:), pr_mat(1,:)), dot(tmp3(2,:), pr_mat(2,:))];
eq7 = yc - right3;

tmp3 = (use(u, subs(xnewd, zmd, zmd+md))-use(u,xnewd)).*kron(phi, ones(1,2));
tmp3 = reshape(tmp3, [2 2])';
right3 = beta*[dot(tmp3(1,:), pr_mat(1,:)), dot(tmp3(2,:), pr_mat(2,:))];
eq8 = yd - right3;


tmp2 = (kron(phi.*use(vdif, yc), ones(1,2)).*use(udif,subs(xnewc, zmc, zmc+mc))...
    - use(udif,subs(xnewc, zmc, zmc-mc))).*kron([1 1], qnew./(1+gma_new));
tmp2 = reshape(tmp2, [2 2])';
right2 = beta*[dot(tmp2(1,:), pr_mat(1,:)), dot(tmp2(2,:), pr_mat(2,:))];

eqn4 = right2;

tmp2 = (kron(phi.*use(vdif, yd), ones(1,2)).*use(udif,subs(xnewd, zmd, zmd+md))...
    - use(udif,subs(xnewd, zmd, zmd-md))).*kron([1 1], qnew./(1+gma_new));
tmp2 = reshape(tmp2, [2 2])';
right2 = beta*[dot(tmp2(1,:), pr_mat(1,:)), dot(tmp2(2,:), pr_mat(2,:))];

eqn6 = right2;

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
eq = [eqn1, eq3, eqn4, eq5, eqn6, eq7, eq8];
eq = subs(eq,irate,[0,0]);
eq = subs(eq, qnew, q);
eq = subs(eq, [zmc, zmd], [omega-zlc-conc, -zld-cond]);
eq = simplify(eq);

x0 = repmat(0.1, [1 14]);
lb = zeros([1 14]);

noweq = subs(eq, theta, [0.9, 1.1]);

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

[smc1, smc2, smd1, smd2] = vpasolve(eeq, [mc, md], [...
    76.95816214879973223169501490536, 76.95816214879973223169501490536,...
    55.539029269916441424211715011069, 55.577927605383623236611918035265]);

%%
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

welcr = simplify(subs(welc, [zmc, irate, qnew, theta], ...
    [omega-zlc-conc, 0, 0, q, 0.9, 1.1]));
welcr = double(subs(welcr,[conc, cond, q, yc, yd, mc, md],s00));

weldr = simplify(subs(weld, [zmd, irate, qnew, theta], ...
    [-zld-cond, 0, 0, q, 0.9, 1.1]));
weldr = double(subs(weldr,[conc, cond, q, yc, yd, mc, md],s00));
wel00 = [welcr, weldr];

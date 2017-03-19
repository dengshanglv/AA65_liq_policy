clear;
beta = 1;
gma = [0 0];
gma_new = [0 0];

pi = 0.5; % fix creditor's fraction.
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
zmc = sym('zmc', [1 2], 'real');
zmd = sym('zmd', [1 2], 'real');
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
% conc = sym('conc', [1 2], 'positive');
% cond = sym('cond', [1 2], 'positive');
% 
% eqn1 = eq1+eq2;
% eq = [eqn1, eq3, eqn4, eq5, eqn6, eq7, eq8];
% eq = subs(eq,irate,[0,0]);
% eq = subs(eq, qnew, q);
% eq = subs(eq, [zmc, zmd], [omega-zlc-conc, -zld-cond]);
% eq = simplify(eq);
% 
% x0 = repmat(0.1, [1 14]);
% lb = zeros([1 14]);
% 
% noweq = subs(eq, theta, [0.9, 1.1]);
% objf = @(x)(sum(double(subs(noweq,[conc, cond, q, yc, yd, mc, md],x)).^2));
% pm = createOptimProblem('fmincon', 'objective', objf,'x0',x0, 'lb', lb);
% gs = GlobalSearch('StartPointsToRun', 'bounds-ineqs', 'Display', 'iter');
% [x00,f00,eflag00] = run(gs, pm);
% 
% welcr = simplify(subs(welc, [zmc, irate, qnew, theta], ...
%     [omega-zlc-conc, 0, 0, q, 0.9, 1.1]));
% welcr = double(subs(welcr,[conc, cond, q, yc, yd, mc, md],x00));
% 
% weldr = simplify(subs(weld, [zmd, irate, qnew, theta], ...
%     [-zld-cond, 0, 0, q, 0.9, 1.1]));
% weldr = double(subs(weldr,[conc, cond, q, yc, yd, mc, md],x00));
% wel00 = [welcr, weldr];
% mindelta = (pi*x00(11:12)+(1-pi)*x00(13:14))./x00(5:6)-1;

%% one i is + and the other is 0 (i1 0, i2 +)

conc = sym('conc', [1 2], 'positive');
cond = sym('cond', [1 2], 'positive');

how = [1, 2]; % first state set to 0
eqa(1) = eq1(how(1)) + eq2(how(1));
eqa([2 3]) = [eq1(how(2)), eq2(how(2))];
eqa([4 5]) = eq3;
eqa(6) = eq4(how(2));
eqa([7 8]) = eq5;
eqa(9) = eq6(how(2));
eqa([10 11 12 13]) = [eq7 eq8];
eqa([14 15]) = [eqn4(how(1)) eqn6(how(1))];

eqa(16:19) = [welc weld];
eqa(20:23) = lbond;

eqa = subs(eqa,irate(how(1)),0);
eqa = subs(eqa, qnew, q);
eqa = subs(eqa, [zmc(how(1)), zmd(how(1))], ...
    [omega-zlc(how(1))-conc(how(1)), -zld(how(1))-cond(how(1))]);
eqa = simplify(eqa);
eqa = subs(eqa, mc(how(2)), zmc(how(2)));
eqa = subs(eqa, md(how(2)), zmd(how(2)));
eqa = subs(eqa, delta(how(2)), 96);
eqa = simplify(eqa);
welr2 = eqa(16:19);
lowb = eqa(20:23);
eqa = eqa(1:15);

A = zeros([2 15]);
A(1, [6 7]) = 1;
A(2, [4 5]) = 1;

lb = repmat(-3000, [1 15]);
lb([1 2 3 4 6 8 9 10 11 12 13 14 15]) = 0;
b = [0; omega];
x0 = ones([1 15]);
x0(7) = -1;
x0([1 6]) = 0.05;

noweq = subs(eqa, theta, [0.9 1.1]);
welr2 = subs(welr2, theta, [0.9 1.1]);
lowb = subs(lowb, theta, [0.9, 1.1]);

conf = @(x)(-double(subs(lowb, [irate(how(2)), conc(how(1)), cond(how(1)), zmc(how(2)), zlc(how(2)), ...
    zmd(how(2)), zld(how(2)), q, yc, yd, mc(how(1)), md(how(1))],x)));
nonlinfcn = @(x)deal(conf(x),[]);
objf = @(x)(sum(double(subs(noweq,[irate(how(2)), conc(how(1)), cond(how(1)), zmc(how(2)), zlc(how(2)), ...
    zmd(how(2)), zld(how(2)), q, yc, yd, mc(how(1)), md(how(1))],x)).^2));
pm = createOptimProblem('fmincon', 'objective', objf,'x0',x0,...
    'Aineq',A,'bineq',b,'lb',lb);
gs = GlobalSearch('StartPointsToRun', 'bounds-ineqs', 'Display', 'iter');
[x01,f01, eflag01] = run(gs, pm);
wel01 = double(subs(welr2, [irate(how(2)), conc(how(1)), cond(how(1)), zmc(how(2)), zlc(how(2)), ...
    zmd(how(2)), zld(how(2)), q, yc, yd, mc(how(1)), md(how(1))],x01));


%% one i is + and the other is 0 (i1 +, i2 0)

conc = sym('conc', [1 2], 'positive');
cond = sym('cond', [1 2], 'positive');

how = [2, 1]; % first state set to 0
eqa(1) = eq1(how(1)) + eq2(how(1));
eqa([2 3]) = [eq1(how(2)), eq2(how(2))];
eqa([4 5]) = eq3;
eqa(6) = eq4(how(2));
eqa([7 8]) = eq5;
eqa(9) = eq6(how(2));
eqa([10 11 12 13]) = [eq7 eq8];
eqa([14 15]) = [eqn4(how(1)) eqn6(how(1))];

eqa(16:19) = [welc weld];
eqa(20:23) = lbond;

eqa = subs(eqa,irate(how(1)),0);
eqa = subs(eqa, qnew, q);
eqa = subs(eqa, [zmc(how(1)), zmd(how(1))], ...
    [omega-zlc(how(1))-conc(how(1)), -zld(how(1))-cond(how(1))]);
eqa = simplify(eqa);
eqa = subs(eqa, mc(how(2)), zmc(how(2)));
eqa = subs(eqa, md(how(2)), zmd(how(2)));
eqa = subs(eqa, irate(how(2)), 0.1);
eqa = simplify(eqa);
welr2 = eqa(16:19);
lowb = eqa(20:23);
eqa = eqa(1:15);

A = zeros([2 15]);
A(1, [6 7]) = 1;
A(2, [4 5]) = 1;

lb = repmat(-3000, [1 15]);
lb([1 2 3 4 6 8 9 10 11 12 13 14 15]) = 0;
b = [0; omega];
x0 = ones([1 15]);
x0(7) = -1;
x0([1 6]) = 0.05;

noweq = subs(eqa, theta, [0.9 1.1]);
welr2 = subs(welr2, theta, [0.9 1.1]);
lowb = subs(lowb, theta, [0.9, 1.1]);

conf = @(x)(-double(subs(lowb, [delta(how(2)), conc(how(1)), cond(how(1)), zmc(how(2)), zlc(how(2)), ...
    zmd(how(2)), zld(how(2)), q, yc, yd, mc(how(1)), md(how(1))],x)));
objf = @(x)(sum(double(subs(noweq,[delta(how(2)), conc(how(1)), cond(how(1)), zmc(how(2)), zlc(how(2)), ...
    zmd(how(2)), zld(how(2)), q, yc, yd, mc(how(1)), md(how(1))],x)).^2));
nonlinfcn = @(x)deal(conf(x),[]);
pm = createOptimProblem('fmincon', 'objective', objf,'x0',x0, ...
    'Aineq',A,'bineq',b,'lb',lb,'nonlcon',nonlinfcn);
gs = GlobalSearch('StartPointsToRun', 'bounds-ineqs', 'Display', 'iter');
ms = MultiStart('StartPointsToRun', 'bounds-ineqs', 'Display', 'iter','UseParallel', true);
[x10,f10, eflag10] = run(ms, pm, 20);

wel10 = double(subs(welr2, [delta(how(2)), conc(how(1)), cond(how(1)), zmc(how(2)), zlc(how(2)), ...
    zmd(how(2)), zld(how(2)), q, yc, yd, mc(how(1)), md(how(1))],x10));
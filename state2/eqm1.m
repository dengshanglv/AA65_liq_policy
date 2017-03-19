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

%% if i >> 0
equ = [eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8];
equ = subs(equ, qnew, q);
equ = subs(equ, mc, zmc);
equ = subs(equ, md, zmd);
equ = subs(equ, delta, [0.2,0.2]);

A = zeros([2 16]);
A(1, [11 15]) = 1;
A(2, [12 16]) = 1;

lb = repmat(-inf, [1 16]);
lb([1 2 13 14 15 16]) = 0;
b = zeros([2 1]);
x0 = ones([1 16]);
x0([11 12]) = -1;
x0([1 2 15 16]) = 0;

parpool;
result1 = repmat(-3, [length(tta) 19]);
parfor ii = 1:length(tta)
    noweq = subs(equ, theta, tta(ii,:));
    objf = @(x)(sum(double(subs(noweq,[irate, q, yc, yd, zlc, zld, zmc, zmd],x)).^2));
    pm = createOptimProblem('fmincon', 'objective', objf,'x0',x0, 'Aineq',A,'bineq',b,'lb',lb);
    gs = GlobalSearch;
    [x,f] = run(gs, pm);
    fprintf('the %.f loop at i>>0', ii);
    disp([f, x]);
    result1(ii, :) = [tta(ii,:), f, x];
end

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

result2 = repmat(-3, [length(tta) 17]);
parfor ii = 1:length(tta)
    noweq = subs(eq, theta, tta(ii,:));
    objf = @(x)(sum(double(subs(noweq,[conc, cond, q, yc, yd, mc, md],x)).^2));
    pm = createOptimProblem('fmincon', 'objective', objf,'x0',x0, 'lb', lb);
    gs = GlobalSearch;
    [x,f] = run(gs, pm);
    fprintf('the %.f loop at i==0', ii);
    disp([f, x]);
    result2(ii, :) = [tta(ii,:), f, x];
end

% res = result(:,4:end);
% conc_res = res(:, 1:2);
% cond_res = res(:, 3:4);
% q_res = res(:, 5:6);
% yc_res = res(:, 7:8);
% yd_res = res(:, 9:10);
% mc_res = res(:, 11:12);
% md_res = res(:, 13:14);
% 
% % test = subs(eq,[conc, cond, q, yc, yd, mc, md], ones([1 14]));
% % vpa(test',3)
% % % irrelevent to delta.
% 

% pi*mc_res./q_res
% (1-pi)*md_res./q_res

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

eqa = subs(eqa,irate(how(1)),0);
eqa = subs(eqa, qnew, q);
eqa = subs(eqa, [zmc(how(1)), zmd(how(1))], ...
    [omega-zlc(how(1))-conc(how(1)), -zld(how(1))-cond(how(1))]);
eqa = simplify(eqa);
eqa = subs(eqa, mc(how(2)), zmc(how(2)));
eqa = subs(eqa, md(how(2)), zmd(how(2)));
eqa = subs(eqa, delta(how(2)), 0.2);
eqa = simplify(eqa);

A = zeros([1 15]);
A(1, [6 7]) = 1;

lb = repmat(-inf, [1 15]);
lb([1 2 3 4 6 8 9 10 11 12 13 14 15]) = 0;
b = 1;
x0 = ones([1 15]);
x0(7) = -1;
x0([1 6]) = 0.05;

result3 = repmat(-3, [length(tta) 18]);
parfor ii = 1:length(tta)
    noweq = subs(eqa, theta, tta(ii,:));
    objf = @(x)(sum(double(subs(noweq,[irate(how(2)), conc(how(1)), cond(how(1)), zmc(how(2)), zlc(how(2)), ...
        zmd(how(2)), zld(how(2)), q, yc, yd, mc(how(1)), md(how(1))],x)).^2));
    pm = createOptimProblem('fmincon', 'objective', objf,'x0',x0, 'Aineq',A,'bineq',b,'lb',lb);
    gs = GlobalSearch;
    [x,f] = run(gs, pm);
    fprintf('the %.f loop at i1=0,i2>0', ii);
    disp([f, x]);
    result3(ii,:) = [tta(ii,:), f, x];
end

% test = subs(eqa,[irate(2), conc(1), cond(1), zmc(2), zlc(2), zmd(2), zld(2), ...
%     q, yc, yd, mc(1), md(1)], ones([1 15]));
% vpa(test',3)


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

eqa = subs(eqa,irate(how(1)),0);
eqa = subs(eqa, qnew, q);
eqa = subs(eqa, [zmc(how(1)), zmd(how(1))], ...
    [omega-zlc(how(1))-conc(how(1)), -zld(how(1))-cond(how(1))]);
eqa = simplify(eqa);
eqa = subs(eqa, mc(how(2)), zmc(how(2)));
eqa = subs(eqa, md(how(2)), zmd(how(2)));
eqa = subs(eqa, delta(how(2)), 0.2);
eqa = simplify(eqa);

A = zeros([1 15]);
A(1, [6 7]) = 1;

lb = repmat(-inf, [1 15]);
lb([1 2 3 4 6 8 9 10 11 12 13 14 15]) = 0;
b = 1;
x0 = ones([1 15]);
x0(7) = -1;
x0([1 6]) = 0.05;

result4 = repmat(-3, [length(tta) 18]);
parfor ii = 1:length(tta)
    noweq = subs(eqa, theta, tta(ii,:));
    objf = @(x)(sum(double(subs(noweq,[irate(how(2)), conc(how(1)), cond(how(1)), zmc(how(2)), zlc(how(2)), ...
        zmd(how(2)), zld(how(2)), q, yc, yd, mc(how(1)), md(how(1))],x)).^2));
    pm = createOptimProblem('fmincon', 'objective', objf,'x0',x0, 'Aineq',A,'bineq',b,'lb',lb);
    gs = GlobalSearch;
    [x,f] = run(gs, pm);
    fprintf('the %.f loop at i1>0,i2=0', ii);
    disp([f, x]);
    result4(ii,:) = [tta(ii,:), f, x];
end

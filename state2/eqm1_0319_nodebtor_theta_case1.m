clear;
beta = 1;
gma = [0 0];
gma_new = [0 0];

pi = 0.5; % fix creditor's fraction.
omega = 4.5;
% theta = sym('theta',[1 2],'positive');
theta = [0.9 1.1];

guess00 = ...
[ 1.916901868358666386305553056607, 1.916901868358666386305553056607, 1.0070403871175859304035000385095, 0.64196166253037009488759204091913, 0.53802887226187384164547345244174, 0.72056823455548175940342745123692, 0.13782646603130568239102037822635, 0.13782646603130568239102037822635, 0.30783264859519892892678005962144, 0.30783264859519892892678005962144]
guess1 = ...
[ 0.036955781650582146239888465388006, 1.916901868358666386305553056607, 1.241269661176348367893623288821, 0.29046397836173983728361028981942, 1.7898220693700208036069814159139, -0.65367436295367755823028233162959, 0.42091423523249262290041182728597, 0.71330584238904154133015468705186, 0.13782646603130568239102037822635, 0.12691825394933546344120254680522, 0.30783264859519892892678005962144]
guess2 = ...
[ 0.01534113627232469665351668065507, 1.916901868358666386305553056607, 0.64033100140951072171443121704198, 0.30044203110448676860106147665915, 1.7814759061872701566542002948244, -0.99365348651735221538725966581245, 0.54413222538720235493400105283554, 0.7213835651159114459900078631755, 0.13316664728075587617948427297706, 0.13782646603130568239102037822635, 0.30783264859519892892678005962144]

syms x1 x2 x3 y;
v = sqrt(x2);
u = log(x1);
udif = diff(u);
vdif = diff(v);
use = @(symfun, para)(cell2sym(arrayfun(@subs, repmat(symfun, [1 length(para)]), para, ...
    'UniformOutput',false)));

pr_mat = [0.7 0.3; 0.3 0.7];
phi = [1 1];

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

tmp = use(u,xnewd).*kron((1+irate)./q, qnew./(1+gma_new));
tmp = reshape(tmp, [2 2])';
right = beta*[dot(tmp(1,:), pr_mat(1,:)), dot(tmp(2,:), pr_mat(2,:))];

eq5 = use(udif, -zmd-zld) - right;       % Zld

% tmp2 = (use(udif, xnewd)+kron(phi.*use(vdif, yd), ones(1,...
%     2)).*use(udif,subs(xnewd, zmd, zmd+md))).*kron(1./q, qnew./(1+gma_new));
% tmp2 = reshape(tmp2, [2 2])';
% right2 = beta*[dot(tmp2(1,:), pr_mat(1,:)), dot(tmp2(2,:), pr_mat(2,:))];
% 
% eq6 = 2*use(udif, -zmd-zld) - right2;      % Zmd

tmp3 = (use(u, subs(xnewc, zmc, zmc+mc))-use(u,xnewc)).*kron(phi, ones(1,2));
tmp3 = reshape(tmp3, [2 2])';
right3 = beta*[dot(tmp3(1,:), pr_mat(1,:)), dot(tmp3(2,:), pr_mat(2,:))];
eq7 = yc - right3;                         % yc

tmp2 = (kron(phi.*use(vdif, yc), ones(1,2)).*use(udif,subs(xnewc, zmc, zmc+mc))...
    - use(udif,subs(xnewc, zmc, zmc-mc))).*kron([1 1], qnew./(1+gma_new));
tmp2 = reshape(tmp2, [2 2])';
right2 = beta*[dot(tmp2(1,:), pr_mat(1,:)), dot(tmp2(2,:), pr_mat(2,:))];

eqn4 = right2;                             % FOC for mc, when i==0

tmp4 = use(u, xnewc)+use(u, subs(xnewc, zmc, zmc-mc));
tmp4 = reshape(tmp4, [2 2])';
right4 = beta*[dot(tmp4(1,:), pr_mat(1,:)), dot(tmp4(2,:), pr_mat(2,:))];

welc = use(u, omega - zmc - zlc) + 0.5*(right4 + use(v,yc));

tmp5 = reshape(use(u, xnewd), [2 2])';
right4 = beta*[dot(tmp5(1,:), pr_mat(1,:)), dot(tmp5(2,:), pr_mat(2,:))];

weld = use(u,  - zmd - zld) + right4;

lbondd = omega.*theta.*q + (zmd+zld.*(1.+irate));
lbondc = (zmc+zlc.*(1.+irate));
lbond = [subs(lbondd, zmd, zmd-md), subs(lbondc, zmc, zmc-mc)];

%% if i == 0
conc = sym('conc', [1 2], 'positive');
cond = sym('cond', [1 2], 'positive');

eqn1 = eq1+eq2;
noweq = [eqn1, eq3, eqn4, eq5, eq7];
noweq = subs(noweq,irate,[0,0]);
noweq = subs(noweq, qnew, q);
noweq = subs(noweq, [zmc, zmd], [omega-zlc-conc, -zld-cond]);
noweq = simplify(noweq);

% [zmc1, zmc2] = solve(noweq([5,6]), mc);
% 
% nconc1 = sidx(solve(noweq(3), conc(1)),1);
% nconc2 = sidx(solve(noweq(4), conc(2)),1);
% ncond1 = sidx(solve(noweq(7), cond(1)),2);
% ncond2 = sidx(solve(noweq(8), cond(2)),2);
% 
% noweq(3) = conc(1) - nconc1;
% noweq(4) = conc(2) - nconc2;
% noweq(5) = mc(1) - zmc1;
% noweq(6) = mc(2) - zmc2;
% noweq(7) = cond(1) - ncond1;
% noweq(8) = cond(2) - ncond2;

[sconc1,sconc2,scond1,scond2,sq1,sq2,syc1,syc2,smc1,smc2] = ...
    vpasolve(noweq,[conc, cond, q, yc, mc], guess00);
% objf = @(x)(sum(double(subs(noweq, [conc, cond, q, yc, yd, mc, md], x)).^2));
% x0 = [ 14.5, 14.5, 13.2, 13.9, 1.14, 0.768, 0.137, 0.137, 0.137, 0.137, 2.32, 2.32, 2.11, 2.23];
% lb = zeros([14 1]);
% 
% pm = createOptimProblem('fmincon', 'objective', objf,'x0',x0, 'lb', lb);
% gs = GlobalSearch('Display', 'iter');
% [x,f] = run(gs, pm);

s00 = [sconc1,sconc2,scond1,scond2,sq1,sq2,syc1,syc2,smc1,smc2];

mindelta1 = pi*smc1/sq1-1;
mindelta2 = pi*smc2/sq2-1;
mindelta = [mindelta1, mindelta2];
welcr = simplify(subs(welc, [zmc, irate, qnew], ...
    [omega-zlc-conc, 0, 0, q]));
welcr = double(subs(welcr,[conc, cond, q, yc, mc],s00));

weldr = simplify(subs(weld, [zmd, irate, qnew], ...
    [-zld-cond, 0, 0, q]));
weldr = double(subs(weldr,[conc, cond, q, yc, mc],s00));
wel00 = [welcr, weldr];

%% i1 = 0, i2 > 0
conc = sym('conc', [1 2], 'positive');
cond = sym('cond', [1 2], 'positive');

eqa(1) = eq1(1) + eq2(1);
eqa([2 3]) = [eq1(2), eq2(2)];
eqa([4 5]) = eq3;
eqa(6) = eq4(2);
eqa([7 8]) = eq5;
eqa([9 10]) = eq7;
eqa(11) = eqn4(1);
eqa(12:15) = [welc weld];

eqa = subs(eqa,irate(1),0);
eqa = subs(eqa, qnew, q);
eqa = subs(eqa, [zmc(1), zmd(1)], ...
    [omega-zlc(1)-conc(1), -zld(1)-cond(1)]);
eqa = simplify(eqa);
eqa = subs(eqa, mc(2), zmc(2));
eqa = simplify(eqa);
welr2 = eqa(12:15);
eqa = eqa(1:11);

% % part1, conc1 yc1 mc1 when i=0. deterministic.
% part1 = subs(eqa([4 9 11]), q, [6, 6]); % q has nothing to do with solutions.
% guess1 = [ 14.6998, 1.68996, 0.102673];
% [z1conc1, z1mc1, z1yc1] = vpasolve(part1, [conc(1), mc(1), yc(1)], guess1);
% % part2, zlc2 zmc2 yc2 when i>0. i2 is a characteristic var.
% for ii = 1:50
%     part2 = subs(eqa([5 6 10]), [irate(2), q], [ii/100 6 6]); % only i2 can influence solutions.
%                                                               % not q1 q2.
%                                                               % i2 +, zlc2 +, zmc2 -, yc2 -
%     guess2 = [ 12.7621, 2.92534, 0.162424];
%     [z1zlc2, z1zmc2, z1yc2] = vpasolve(part2, [zlc(2), zmc(2), yc(2)], guess2);
%     z1res1(ii,:) = double([z1zlc2, z1zmc2, z1yc2]);
% end
% % part3, cond1, affected by q2/q1. negtively!
% for kk = 5:15
%     part3 = subs(eqa(7), q, [1 kk/10]);
%     guess3 = 0.2985;
%     z1cond1 = vpasolve(part3, cond(1),guess3);
%     z1res2(kk-4) = z1cond1;
% end
% % part4, zld2, affected by q1/q2 and i2, positively! probably a more general form of part3.
% % notably, zld2 is negtive, so larger means smaller in absolute value.
% for kk = 5:15
%     for ii = 1:50
%         % use q2/q1 as a the loop var.
%         part4 = subs(eqa(8), [q irate(2) zmd(2)], [1 kk/10 ii/100 0]);
%         guess4 = -0.29091627923804363171656938394774;
%         z1zld2 = vpasolve(part4, zld(2),guess4);
%         z1res3(kk-4,ii) = z1zld2;
%         z1wel01(kk-4,ii,:) = double(subs(welr2, [irate(2) q conc(1) yc(1) mc(1)...
%             zlc(2) zmc(2) yc(2) cond(1) zld(2) zmd(2)],...
%             [ii/100 1 kk/10 z1conc1 z1mc1 z1yc1 z1res1(ii,:) z1res2(kk-4) z1zld2 0]));
%     end
% end

% [xx,yy] = meshgrid(5:15,1:50);
% mesh(xx,yy,z1res3');
% z1wel01(:,:,5) = (z1wel01(:,:,1)+z1wel01(:,:,2))/2;
% z1wel01(:,:,6) = (z1wel01(:,:,3)+z1wel01(:,:,4))/2;
% z1wel01(:,:,7) = (z1wel01(:,:,5)+z1wel01(:,:,6))/2;
% mesh(xx,yy,z1wel01(:,:,5)'-mean2(z1wel01(:,:,5))); hold on;
% mesh(xx,yy,z1wel01(:,:,6)'-mean2(z1wel01(:,:,6))); hold on;
% mesh(xx,yy,z1wel01(:,:,7)'-mean2(z1wel01(:,:,7))); hold off;
% legend('c','d','all');

% nconc1 = sidx(solve(eqa(4), conc(1)), 1);  % smaller
% nzlc2 = sidx(solve(eqa(5), zlc(2)), 2);    % bigger
% nzmc2 = sidx(solve(eqa(6), zmc(2)), 2);    % bigger, so positive
% ncond1 = sidx(solve(eqa(7), cond(1)), 2);  % smaller
% nzld2 = sidx(solve(eqa(8), zld(2)), 1);    % bigger
% nmc1 = solve(eqa(11), mc(1));
% 
% eqa(4) = nconc1 - conc(1);
% eqa(5) = nzlc2 - zlc(2);
% eqa(6) = nzmc2 - zmc(2);
% eqa(7) = ncond1 - cond(1);
% eqa(8) = nzld2 - zld(2);
% eqa(11) = nmc1 - mc(1);

neqa = subs(eqa, delta(2), mindelta2 - 0.01);
neqa = subs(neqa, zmd(2), 0);
[s1irate2, s1conc1, s1cond1, s1zmc2, s1zlc2, s1zld2, s1q1, s1q2, s1yc1, s1yc2,s1mc1] = ...
vpasolve(neqa, [irate(2), conc(1), cond(1), zmc(2), zlc(2), ...
    zld(2), q, yc, mc(1)],guess1);
s1 = [s1irate2, s1conc1, s1cond1, s1zmc2, s1zlc2, s1zld2, s1q1, s1q2, s1yc1, s1yc2,s1mc1];

wel01 = subs(welr2, [irate(2), conc(1), cond(1), zmc(2), zlc(2),zld(2), q, yc, mc(1)],s1);
wel01 = subs(wel01, zmd(2), 0);

%%  i1>0 i2=0
conc = sym('conc', [1 2], 'positive');
cond = sym('cond', [1 2], 'positive');

eqa(1) = eq1(2) + eq2(2);
eqa([2 3]) = [eq1(1), eq2(1)];
eqa([4 5]) = eq3;
eqa(6) = eq4(1);
eqa([7 8]) = eq5;
eqa([9 10]) = eq7;
eqa(11) = eqn4(2);
eqa(12:15) = [welc weld];

eqa = subs(eqa,irate(2),0);
eqa = subs(eqa, qnew, q);
eqa = subs(eqa, [zmc(2), zmd(2)], [omega-zlc(2)-conc(2), -zld(2)-cond(2)]);
eqa = simplify(eqa);
eqa = subs(eqa, mc(1), zmc(1));
eqa = simplify(eqa);
welr2 = eqa(12:15);
eqa = eqa(1:11);

% nconc2 = sidx(solve(eqa(5), conc(2)), 1);  % smaller
% nzlc1 = sidx(solve(eqa(4), zlc(1)), 2);    % bigger
% nzmc1 = sidx(solve(eqa(6), zmc(1)), 2);    % bigger, so positive
% ncond2 = sidx(solve(eqa(8), cond(2)), 2);  % smaller
% nzld1 = sidx(solve(eqa(7), zld(1)), 1);    % bigger
% nmc2 = solve(eqa(11), mc(2));
% 
% eqa(5) = nconc2 - conc(2);
% eqa(4) = nzlc1 - zlc(1);
% eqa(6) = nzmc1 - zmc(1);
% eqa(8) = ncond2 - cond(2);
% eqa(7) = nzld1 - zld(1);
% eqa(11) = nmc2 - mc(2);

% for ii = 85:99
%     candelta = -ii/100;
%     neqa = subs(eqa, delta(1), candelta);
%     neqa = subs(neqa, zmd(1), 0);
%     guess = [ 0.047747017509509750990419202156917, 14.518816696398656412906897040188, 0.054558360147140907066899369290412, 2.152155851560617062578108275069, 13.301748597962019055776712763267, -0.081362652661085671368333359271424, 7.6862708984307752234932438395321, 7.7133124717271013400131017952609, 0.12303348491426073484524021304658, 0.13673226529299057933035815291063, 2.3179190370869255367205764488979];
%     [s2irate1, s2conc2, s2cond2, s2zmc1, s2zlc1, s2zld1, s2q1, s2q2, s2yc1, s2yc2,s2mc2] = ...
%     vpasolve(neqa, [irate(1), conc(2), cond(2), zmc(1), zlc(1), ...
%         zld(1), q, yc, mc(2)],guess);
%     s2 = [s2irate1, s2conc2, s2cond2, s2zmc1, s2zlc1, s2zld1, s2q1, s2q2, s2yc1, s2yc2,s2mc2];
% 
%     wel10 = subs(welr2, [irate(1), conc(2), cond(2), zmc(1), zlc(1), ...
%         zld(1), q, yc, mc(2)],s2);
%     wel10 = subs(wel10, zmd(1), 0);
%     res(ii-84,:) = [wel10,sum(wel10)/4,s2];
% end

neqa = subs(eqa, delta(1), mindelta1 - 0.01);
neqa = subs(neqa, zmd(1), 0);
[s2irate1, s2conc2, s2cond2, s2zmc1, s2zlc1, s2zld1, s2q1, s2q2, s2yc1, s2yc2,s2mc2] = ...
vpasolve(neqa, [irate(1), conc(2), cond(2), zmc(1), zlc(1), ...
    zld(1), q, yc, mc(2)],guess2);
s2 = [s2irate1, s2conc2, s2cond2, s2zmc1, s2zlc1, s2zld1, s2q1, s2q2, s2yc1, s2yc2,s2mc2];

wel10 = subs(welr2, [irate(1), conc(2), cond(2), zmc(1), zlc(1), ...
    zld(1), q, yc, mc(2)],s2);
wel10 = subs(wel10, zmd(1), 0);

%%
vpa(wel00, 6)
vpa(wel01, 6)
vpa(wel10, 6)

vpa(sum(wel00)/4, 10)
vpa(sum(wel01)/4, 10)
vpa(sum(wel10)/4, 10)

[sq1 sq2 sq2/sq1]
[s1q1 s1q2 s1q2/s1q1]
[s2q1 s2q2 s2q2/s2q1]

s00
s1
s2

mindelta
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

tmp2 = (use(udif, xnewc)+kron(theta.*use(vdif, yc), ones(1,...
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

tmp3 = use(u, subs(xnewc, zmc, zmc+mc))-use(u,xnewc);
tmp3 = reshape(tmp3, [2 2])';
right3 = beta*[dot(tmp3(1,:), pr_mat(1,:)), dot(tmp3(2,:), pr_mat(2,:))];
eq7 = yc - right3;                         % yc

tmp2 = (kron(theta.*use(vdif, yc), ones(1,2)).*use(udif,subs(xnewc, zmc, zmc+mc))...
    - use(udif,subs(xnewc, zmc, zmc-mc))).*kron([1 1], qnew./(1+gma_new));
tmp2 = reshape(tmp2, [2 2])';
right2 = [dot(tmp2(1,:), pr_mat(1,:)), dot(tmp2(2,:), pr_mat(2,:))];

eqn4 = right2;                             % FOC for mc, when i==0

tmp4 = use(u, xnewc)+use(u, subs(xnewc, zmc, zmc-mc));
tmp4 = reshape(tmp4, [2 2])';
right4 = beta*[dot(tmp4(1,:), pr_mat(1,:)), dot(tmp4(2,:), pr_mat(2,:))];

welc = use(u, omega - zmc - zlc) + 0.5*(right4 + use(v,yc));

tmp5 = reshape(use(u, xnewd), [2 2])';
right4 = beta*[dot(tmp5(1,:), pr_mat(1,:)), dot(tmp5(2,:), pr_mat(2,:))];

weld = use(u,  - zmd - zld) + right4;

%% if i == 0
conc = sym('conc', [1 2], 'positive');
cond = sym('cond', [1 2], 'positive');

eqn1 = eq1+eq2;
noweq = [eqn1, eq3, eqn4, eq5, eq7];
noweq = subs(noweq,irate,[0,0]);
noweq = subs(noweq, qnew, q);
noweq = subs(noweq, [zmc, zmd], [omega-zlc-conc, -zld-cond]);
noweq = simplify(noweq);

[sconc1,sconc2,scond1,scond2,sq1,sq2,syc1,syc2,smc1,smc2] = ...
    vpasolve(noweq,[conc, cond, q, yc, mc], [ 14.699811801620312593846229734541, 14.285702214684012081426303771732, 0.29846374303771044687517737433941, 0.30332288030702887647760055438811, 7.5008622276709884796392964455597, 7.70548745250447952104804783694, 0.10267308471444358142815654463766, 0.1759957789415155521994036192162, 1.6899587143919742564174116184134, 3.0913739883541854757101443968136]);

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

% eqa 1 2 3 mkt clear, 4 conc1 5 zlc2 6 zmc2 7 cond1 8 zld2 9 yc1 10 yc2
% 11 mc1
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
%         z1wel10(kk-4,ii,:) = double(subs(welr2, [irate(2) q conc(1) yc(1) mc(1)...
%             zlc(2) zmc(2) yc(2) cond(1) zld(2) zmd(2)],...
%             [ii/100 1 kk/10 z1conc1 z1mc1 z1yc1 z1res1(ii,:) z1res2(kk-4) z1zld2 0]));
%     end
% end
% 
% [xx,yy] = meshgrid(5:15,1:50);
% mesh(xx,yy,z1res3');
% z1wel10(:,:,5) = (z1wel10(:,:,1)+z1wel10(:,:,2))/2;
% z1wel10(:,:,6) = (z1wel10(:,:,3)+z1wel10(:,:,4))/2;
% z1wel10(:,:,7) = (z1wel10(:,:,5)+z1wel10(:,:,6))/2;
% mesh(xx,yy,z1wel10(:,:,5)'-mean2(z1wel10(:,:,5))); hold on;
% mesh(xx,yy,z1wel10(:,:,6)'-mean2(z1wel10(:,:,6))); hold on;
% mesh(xx,yy,z1wel10(:,:,7)'-mean2(z1wel10(:,:,7))); hold off;
% legend('c','d','all');

neqa = subs(eqa, delta(2), -0.81);
neqa = subs(neqa, zmd(2), 0);
[s1irate2, s1conc1, s1cond1, s1zmc2, s1zlc2, s1zld2, s1q1, s1q2, s1yc1, s1yc2,s1mc1] = ...
vpasolve(neqa, [irate(2), conc(1), cond(1), zmc(2), zlc(2), ...
    zld(2), q, yc, mc(1)],[ 0.042360174584666576213368838519326, 14.699811801620312593846229734541, 0.29854884273380714574507721543639, 2.9253379355997417384198871488503, 12.762093794163258411296088281678, -0.29091627923804363171656938394774, 7.5008196778229401302043465250112, 7.6982577252624782589997030232903, 0.10267308471444358142815654463766, 0.16242374827016772866982989081549, 1.6899587143919742564174116184134]);
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

neqa = subs(eqa, delta(1), -0.9);
neqa = subs(neqa, zmd(1), 0);
guess = [ 0.067025890072689376671256378005307, 14.285702214684012081426303771732, 0.30340218141405773482137269910707, 1.49883291319133177359946376062, 13.769135830743982921233299039136, -0.27963961202199695883812519355532, 7.4941645659566588679973188031002, 7.7054478019509650918761617645806, 0.086764571624980574985287030489505, 0.1759957789415155521994036192162, 3.0913739883541854757101443968136];
[s2irate1, s2conc2, s2cond2, s2zmc1, s2zlc1, s2zld1, s2q1, s2q2, s2yc1, s2yc2,s2mc2] = ...
vpasolve(neqa, [irate(1), conc(2), cond(2), zmc(1), zlc(1), ...
    zld(1), q, yc, mc(2)],guess);
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

% vpa(sum(wel00(1:2))/2, 10)
% vpa(sum(wel01(1:2))/2, 10)
% vpa(sum(wel00(3:4))/2, 10)
% vpa(sum(wel01(3:4))/2, 10)

[sq1 sq2 sq2/sq1]
[s1q1 s1q2 s1q2/s1q1]
[s2q1 s2q2 s2q2/s2q1]

[sq1 sq2 sq2/sq1]-[s1q1 s1q2 s1q2/s1q1]
[sq1 sq2 sq2/sq1]-[s2q1 s2q2 s2q2/s2q1]
clear;
% beta = 0.98;
% theta = [0.8 1.25];
% omega = 30;
syms beta omega positive;
theta = sym('theta',[1 2],'positive');

gma = [0 0];
gma_new = [0 0];

pi = 0.5; % fix creditor's fraction.

syms x1 x2 x3 y;
v = sqrt(x2);
u = log(x1);
udif = diff(u);
vdif = diff(v);
use = @(symfun, para)(cell2sym(arrayfun(@subs, repmat(symfun, [1 length(para)]), para, ...
    'UniformOutput',false)));

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

xnewd = kron(omega.*theta, [1 1]) + kron((zmd+zld.*(1.+irate))./q, qnew./(1+gma_new));
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

% tmp = (use(udif, xnewd)+use(udif,subs(xnewd, zmd, zmd-md))).*kron((1+irate)./q, qnew./(1+gma_new));
% tmp = reshape(tmp, [2 2])';
% right = beta*[dot(tmp(1,:), pr_mat(1,:)), dot(tmp(2,:), pr_mat(2,:))];
% 
% eq5 = 2*use(udif, -zmd-zld) - right;       % Zld
% 
% tmp2 = (use(udif, xnewd)+kron(phi.*use(vdif, yd), ones(1,...
%     2)).*use(udif,subs(xnewd, zmd, zmd+md))).*kron(1./q, qnew./(1+gma_new));
% tmp2 = reshape(tmp2, [2 2])';
% right2 = beta*[dot(tmp2(1,:), pr_mat(1,:)), dot(tmp2(2,:), pr_mat(2,:))];
% 
% eq6 = 2*use(udif, -zmd-zld) - right2;      % Zmd

tmp3 = (use(u, subs(xnewc, zmc, zmc+mc))-use(u,xnewc));
tmp3 = reshape(tmp3, [2 2])';
right3 = beta*[dot(tmp3(1,:), pr_mat(1,:)), dot(tmp3(2,:), pr_mat(2,:))];
eq7 = yc - right3;                         % yc

tmp2 = (kron(theta.*use(vdif, yc), ones(1,2)).*use(udif,subs(xnewc, zmc, zmc+mc))...
    - use(udif,subs(xnewc, zmc, zmc-mc))).*kron([1 1], qnew./(1+gma_new));
tmp2 = reshape(tmp2, [2 2])';
right2 = beta*[dot(tmp2(1,:), pr_mat(1,:)), dot(tmp2(2,:), pr_mat(2,:))];

eqn4 = right2;                             % FOC for mc, when i==0

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

%%
% subeq = subs([eq3(1) eq4(1) eq7(1) welc(1)], ...
%     [theta, qnew, mc, beta, omega], [0.8, 1.25, q, zmc, 0.98, 30]);
% 
% for ii = 1:50
%     for jj = 5:15
%         seq = subs(subeq, [irate(1), q], [ii/100, 1, jj/10]);
%         [nzmc1, nzlc1, nyc1] = vpasolve(seq(1:3),[zmc(1), zlc(1),yc(1)],[ 1.6591399010712415076836705651054, 13.636035142986937007390544345253, 0.10007744505829266329636751635906]);
%         subwel = subs(seq(4), [zmc(1), zlc(1),yc(1)], [nzmc1, nzlc1, nyc1]);
%         res(ii,jj-4,:) = double([subwel, nzmc1, nzlc1, nyc1]);
%     end
% end
% seq = subs(subeq, [irate(1), q], [0.01, 1, 0.8]);
% % seq(1) = zlc(1) - sidx(solve(seq(1),zlc(1)),2);
% % seq(2) = zmc(1) - sidx(solve(seq(2),zmc(1)),1);
% [nzmc1, nzlc1, nyc1] = vpasolve(seq(1:3),[zmc(1), zlc(1),yc(1)],[ 1.6591399010712415076836705651054, 13.636035142986937007390544345253, 0.10007744505829266329636751635906]);
% subwel = subs(seq(4), [zmc(1), zlc(1),yc(1)], [nzmc1, nzlc1, nyc1]);
% double([nzmc1, nzlc1, nyc1,subwel])

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

[sconc1,sconc2,scond1,scond2,sq1,sq2,syc1,syc2,syd1,syd2,smc1,smc2,smd1,smd2] = ...
    vpasolve(noweq,[conc, cond, q, yc, yd, mc, md], [ 14.615367029293011283748301888421, 14.409228195056342615288411571833, 14.421285183426302126527148737666, 14.585583891050439719061541210364, 0.48167389364034329486227468695656, 0.50259395694660883282502360890104, 0.09519147371668585614467830329925, 0.19473119541616853582582146315674, 0.09513840795958135409431963623444, 0.19461287699237335641715832611075, 1.9861012399000304028513282058772, 2.6862280239614762957398191825816, 1.9589287369624004692870222627505, 2.7179095134452263425656883915935]);

% objf = @(x)(sum(double(subs(noweq, [conc, cond, q, yc, yd, mc, md], x)).^2));
% x0 = [ 14.5, 14.5, 13.2, 13.9, 1.14, 0.768, 0.137, 0.137, 0.137, 0.137, 2.32, 2.32, 2.11, 2.23];
% lb = zeros([14 1]);
% 
% pm = createOptimProblem('fmincon', 'objective', objf,'x0',x0, 'lb', lb);
% gs = GlobalSearch('Display', 'iter');
% [x,f] = run(gs, pm);

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
nzmc2 = sidx(solve(eqa(6), zmc(2)), 2);    % bigger, so positive
nmc1 = solve(eqa(14), mc(1));

eqa(4) = nconc1 - conc(1);
eqa(5) = nzlc2 - zlc(2);
eqa(6) = nzmc2 - zmc(2);
eqa(14) = nmc1 - mc(1);

% neqa = subs(eqa, delta(2), 1.7);
% [s1irate2, s1conc1, s1cond1, s1zmc2, s1zlc2, s1zmd2, s1zld2, s1q1, s1q2, s1yc1, s1yc2, ...
%     s1yd1, s1yd2, s1mc1, s1md1] = ...
% vpasolve(neqa, [irate(2), conc(1), cond(1), zmc(2), zlc(2), ...
%     zmd(2), zld(2), q, yc, yd, mc(1), md(1)], [ 0.00084081792675536975040532940339449, 14.615367029293011283748301888421, 14.415924731023672473196058005587, 2.6829607906708327276127523313331, 12.907292597772275845351494707303, 2.7130957819505741441081649417686, -17.291427333940530773926953082146, 0.48435411984165812152782005299598, 0.50596091822657597157272944912948, 0.09519147371668585614467830329925, 0.19439304937304555656865421536061, 0.095135485914521411700612083609967, 0.19426892970948519898698146894221, 1.9861012399000304028513282058772, 1.9581566461984012222391673030071]);
% 
% s1 = [s1irate2, s1conc1, s1cond1, s1zmc2, s1zlc2, s1zmd2, s1zld2, s1q1, s1q2, s1yc1, s1yc2, ...
%     s1yd1, s1yd2, s1mc1, s1md1];
% 
% wel01 = subs(welr2, [irate(2), conc(1), cond(1), zmc(2), zlc(2), ...
%     zmd(2), zld(2), q, yc, yd, mc(1), md(1)],s1);

%%
report = repmat(-2, [151 19]);
parfor ii = 20:171
    neqa = subs(eqa, delta(2), ii/100);
    [s1irate2, s1conc1, s1cond1, s1zmc2, s1zlc2, s1zmd2, s1zld2, s1q1, s1q2, s1yc1, s1yc2, ...
        s1yd1, s1yd2, s1mc1, s1md1] = ...
    vpasolve(neqa, [irate(2), conc(1), cond(1), zmc(2), zlc(2), ...
        zmd(2), zld(2), q, yc, yd, mc(1), md(1)], [ 0.00084081792675536975040532940339449, 14.615367029293011283748301888421, 14.415924731023672473196058005587, 2.6829607906708327276127523313331, 12.907292597772275845351494707303, 2.7130957819505741441081649417686, -17.291427333940530773926953082146, 0.48435411984165812152782005299598, 0.50596091822657597157272944912948, 0.09519147371668585614467830329925, 0.19439304937304555656865421536061, 0.095135485914521411700612083609967, 0.19426892970948519898698146894221, 1.9861012399000304028513282058772, 1.9581566461984012222391673030071]);

    s1 = [s1irate2, s1conc1, s1cond1, s1zmc2, s1zlc2, s1zmd2, s1zld2, s1q1, s1q2, s1yc1, s1yc2, ...
        s1yd1, s1yd2, s1mc1, s1md1];

    wel01 = subs(welr2, [irate(2), conc(1), cond(1), zmc(2), zlc(2), ...
        zmd(2), zld(2), q, yc, yd, mc(1), md(1)],s1);
    report(ii, :) = double([s1irate2, s1conc1, s1cond1, s1zmc2, s1zlc2, s1zmd2, s1zld2, s1q1, s1q2, s1yc1, s1yc2, ...
        s1yd1, s1yd2, s1mc1, s1md1, wel01]);
end
saverep = report;
%%
clear rep;
report = saverep;
report = report(~ismember(report, repmat(-2, [1 19]), 'rows'),:);
report = report(~ismember(report, zeros([1 19]), 'rows'),:);
rep(:,1) = (report(:,14) + report(:,15))/2;        % real balance 1
rep(:,2) = (report(:,4) + report(:,6))/2;           % r. b. 2
rep(:,3:4) = (report(:,10:11) + report(:,12:13))/2; % st-2 output
rep(:,5) = rep(:,3) + (theta(1)*omega+omega)/2;    %gdp1
rep(:,6) = rep(:,4) + (theta(2)*omega+omega)/2;    %gdp2
rep(:,7) = (report(:,16)+report(:,17))/2;
rep(:,8) = (report(:,18)+report(:,19))/2;
rep(:,9) = (rep(:,7)+rep(:,8))/2;
rep2(:,1) = (0.7*report(:,8)/(1+gma(1)) + 0.3*report(:,9)/(1+gma(2)))./report(:,8);
rep2(:,2) = (0.3*report(:,8)/(1+gma(1)) + 0.7*report(:,9)/(1+gma(2)))./report(:,9);
rep2 = rep2 - 1;
rep = [rep2, rep];

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
nzmc1 = sidx(solve(eqa(6), zmc(1)), 2);    % bigger, so positive

nmc2 = solve(eqa(14), mc(2));

eqa(5) = nconc2 - conc(2);
eqa(4) = nzlc1 - zlc(1);
eqa(6) = nzmc1 - zmc(1);
eqa(14) = nmc2 - mc(2);

% eqa_now = subs(eqa, delta(1), 0.57);
% [s2irate1, s2conc2, s2cond2, s2zmc1, s2zlc1, s2zmd1, s2zld1, s2q1, s2q2, s2yc1, s2yc2, ...
%     s2yd1, s2yd2, s2mc2, s2md2] = ...
% vpasolve(eqa_now, [irate(1), conc(2), cond(2), zmc(1), zlc(1), ...
%     zmd(1), zld(1), q, yc, yd, mc(2), md(2)],[ 0.00040591635374116399079628081950222, 14.518816696398656412906897040188, 13.872725335519556168614844775128, 2.3164320737315599356018322774018, 13.164513831250045756194515237321, 1.9706154316542557251408031415962, -14.720957830020692206273051663327, 1.3653017533075846053320494964962, 0.80422898404089370923912909234216, 0.13660833294436349965044767307515, 0.13673226529299057933035815291063, 0.13078873324027073807007166530579, 0.12186069885165503830868361871485, 2.3179190370869255367205764488979, 2.0373037743609396405487375039048]);
% s2 = [s2irate1, s2conc2, s2cond2, s2zmc1, s2zlc1, s2zmd1, s2zld1, s2q1, s2q2, s2yc1, s2yc2, ...
%     s2yd1, s2yd2, s2mc2, s2md2];
% 
% wel10 = subs(welr2, [irate(1), conc(2), cond(2), zmc(1), zlc(1), ...
%     zmd(1), zld(1), q, yc, yd, mc(2), md(2)],s2);

report2 = repmat(-2, [1 19]);
parfor ii = 1:57
    eqa_now = subs(eqa, delta(1), ii/100);
    [s2irate1, s2conc2, s2cond2, s2zmc1, s2zlc1, s2zmd1, s2zld1, s2q1, s2q2, s2yc1, s2yc2, ...
        s2yd1, s2yd2, s2mc2, s2md2] = ...
    vpasolve(eqa_now, [irate(1), conc(2), cond(2), zmc(1), zlc(1), ...
        zmd(1), zld(1), q, yc, yd, mc(2), md(2)],[ 0.00040591635374116399079628081950222, 14.518816696398656412906897040188, 13.872725335519556168614844775128, 2.3164320737315599356018322774018, 13.164513831250045756194515237321, 1.9706154316542557251408031415962, -14.720957830020692206273051663327, 1.3653017533075846053320494964962, 0.80422898404089370923912909234216, 0.13660833294436349965044767307515, 0.13673226529299057933035815291063, 0.13078873324027073807007166530579, 0.12186069885165503830868361871485, 2.3179190370869255367205764488979, 2.0373037743609396405487375039048]);
    s2 = [s2irate1, s2conc2, s2cond2, s2zmc1, s2zlc1, s2zmd1, s2zld1, s2q1, s2q2, s2yc1, s2yc2, ...
        s2yd1, s2yd2, s2mc2, s2md2];

    wel10 = subs(welr2, [irate(1), conc(2), cond(2), zmc(1), zlc(1), ...
        zmd(1), zld(1), q, yc, yd, mc(2), md(2)],s2);
    report2(ii, :) = double([s2, wel10]);
end
saverep2 = report2;
%%
clear rep3;
report2 = saverep2;
report2 = report2(~ismember(report2, repmat(-2, [1 19]), 'rows'),:);
report2 = report2(~ismember(report2, zeros([1 19]), 'rows'),:);
rep3(:,1) = (report2(:,4) + report2(:,6))/2;           % r. b. 1
rep3(:,2) = (report2(:,14) + report2(:,15))/2;        % real balance 2
rep3(:,3:4) = (report2(:,10:11) + report2(:,12:13))/2; % st-2 output
rep3(:,5) = rep3(:,3) + (theta(1)*omega+omega)/2;    %gdp1
rep3(:,6) = rep3(:,4) + (theta(2)*omega+omega)/2;    %gdp2
rep3(:,7) = (report2(:,16)+report2(:,17))/2;
rep3(:,8) = (report2(:,18)+report2(:,19))/2;
rep3(:,9) = (rep3(:,7)+rep3(:,8))/2;
rep4(:,1) = (0.7*report2(:,8)/(1+gma(1)) + 0.3*report2(:,9)/(1+gma(2)))./report2(:,8);
rep4(:,2) = (0.3*report2(:,8)/(1+gma(1)) + 0.7*report2(:,9)/(1+gma(2)))./report2(:,9);
rep4 = rep4 - 1;
rep3 = [rep4, rep3];

% %%
% vpa(wel00, 6)
% vpa(wel01, 6)
% vpa(wel10, 6)
% 
% vpa(sum(wel00)/4, 10)
% vpa(sum(wel01)/4, 10)
% vpa(sum(wel10)/4, 10)

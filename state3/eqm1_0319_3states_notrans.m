clear;
beta = 1;
gma = [0 0 0];
gma_new = [0 0 0];

pi = 0.5; % fix creditor's fraction.
omega = 5;
% theta = sym('theta',[1 2],'positive');
theta = [1 1.05 1.1];
phi = [1 1 1];
NS = 3; % number of states.
pr_mat = [0.6 0.2 0.2 0.2 0.6 0.2 0.2 0.2 0.6];

guess00 = ...
[ 0.054847813137390982333957465884081, 0.050702534851471585715860319482897, 0.047316140597481193815296465869318, 2.4941770382768850524501437474731, 2.5024675948487238456863380402754, 2.5092403833567046294874657475026, 0.39919472148997269616758606170949, 0.4008775581691319047153760205178, 0.40120624737391557983341318353737, 0.13719293033107384658428659053138, 0.13735984641789253643216193474488, 0.13700419852068345335561559055283, 2.3961273354483329828819413207588, 2.3961273354483329828819413207588, 2.3961273354483329828819413207588, 0.3847908107439986611584750745268, 0.3847908107439986611584750745268, 0.3847908107439986611584750745268, 0.13782646603130568239102037822635, 0.13782646603130568239102037822635, 0.13782646603130568239102037822635];
guess1 = ...
[ 0.000049094183534151009639894043200152, 2.4114832380195587015312630023722, 2.5059183084592670632322814596045, 0.39339100915833805829113086946176, 2.1927856484896820024370394243769, 0.40783575579828698370339844989306, -2.9230539810415883929292144580553, 0.041299226760587117618227769011617, 0.035479216202359325751177142838184, 0.13271329319097929883023759690886, 0.1324602822143816286259547007431, 0.13177658886755127124205908773396, 0.13116424375499268040374394538242, 0.39354200434662717751343038465057, 0.40763106191078459303873390339754]
guess2 = ...
[ 0.000034667686349544878810265898139829, 2.413819495021815904644190803745, 2.5152188714384434853568085704155, 0.39352147883546564921451665902298, 2.194989554928233797956294203685, 0.40760322329391267696505615956595, -2.913502611873836809720879334148, 0.041305822591887657207493844062934, 0.035480816769870304999500312919766, 0.13270327759628620710780509959688, 0.13247445176200721315797672037043, 0.1317656048495903530939913531403, 0.1311762755129774552690772801198, 0.39342082063305715591961755095892, 0.40786292657402357527119153002588]

syms x1 x2 x3 y;
v = sqrt(x2);
u = log(x1);
% eta = 1.05;
% u = x1^(1-eta)/(1-eta);

udif = diff(u);
vdif = diff(v);
use = @(symfun, para)(cell2sym(arrayfun(@subs, repmat(symfun, [1 length(para)]), para, ...
    'UniformOutput',false)));

q = sym('q', [1 NS], 'real');
qnew = sym('qnew', [1 NS], 'real');
% qnew = q;

delta = sym('delta', [1 NS], 'real');
irate = sym('irate', [1 NS], 'real');
zmc = sym('zmc', [1 NS], 'positive');
zmd = sym('zmd', [1 NS], 'positive');
zlc = sym('zlc', [1 NS], 'real');
zld = sym('zld', [1 NS], 'real');
mc = sym('mc', [1 NS], 'positive');
md = sym('md', [1 NS], 'positive');
yc = sym('yc', [1 NS], 'positive');
yd = sym('yd', [1 NS], 'positive');

eq1 = pi.*zmc + (1-pi).*zmd - (1+delta).*q;
eq2 = delta.*q + pi.*zlc + (1-pi).*zld;

xnewd = kron(omega.*theta, ones([1 NS]))...
    + kron((zmd+zld.*(1.+irate))./q, qnew./(1+gma_new));
xnewc = kron((zmc+zlc.*(1.+irate))./q, qnew./(1+gma_new));

foczlc = (use(udif, xnewc)+use(udif,subs(xnewc, zmc, zmc-mc))).*kron((1+irate)./q, qnew./(1+gma_new));
eq3 = 2*use(udif, omega-zmc-zlc) - beta*sum(reshape(foczlc.*pr_mat,NS,NS));  % Zlc 

foczmc = (use(udif, xnewc)+kron(phi.*use(vdif, yc),...
    ones([1 NS])).*use(udif,subs(xnewc, zmc, zmc+mc))).*kron(1./q, qnew./(1+gma_new));
eq4 = 2*use(udif, omega-zmc-zlc) - beta*sum(reshape(foczmc.*pr_mat,NS,NS)); % Zmc

foczld = (use(udif, xnewd)+use(udif,subs(xnewd, zmd, zmd-md))).*kron((1+irate)./q, qnew./(1+gma_new));
eq5 = 2*use(udif, -zmd-zld) - beta*sum(reshape(foczld.*pr_mat,NS,NS));       % Zld

foczmd = (use(udif, xnewd)+kron(phi.*use(vdif, yd), ...
    ones([1 NS])).*use(udif,subs(xnewd, zmd, zmd+md))).*kron(1./q, qnew./(1+gma_new));
eq6 = 2*use(udif, -zmd-zld) - beta*sum(reshape(foczmd.*pr_mat,NS,NS));      % Zmd

ycr = (use(u, subs(xnewc, zmc, zmc+mc))-use(u,xnewc)).*kron(phi, ones([1 NS]));
eq7 = yc - beta*sum(reshape(ycr.*pr_mat,NS,NS));                         % yc

ydr = (use(u, subs(xnewd, zmd, zmd+md))-use(u,xnewd)).*kron(phi, ones([1 NS]));
eq8 = yd - beta*sum(reshape(ydr.*pr_mat,NS,NS));                         % yd


focmc = (kron(phi.*use(vdif, yc), ones([1 NS])).*use(udif,subs(xnewc, zmc, zmc+mc))...
    - use(udif,subs(xnewc, zmc, zmc-mc))).*kron(ones([1 NS]), qnew./(1+gma_new));
eqn4 = sum(reshape(focmc.*pr_mat,NS,NS));                 % FOC for mc, when i==0

focmd = (kron(phi.*use(vdif, yd), ones([1 NS])).*use(udif,subs(xnewd, zmd, zmd+md))...
    - use(udif,subs(xnewd, zmd, zmd-md))).*kron(ones([1 NS]), qnew./(1+gma_new));
eqn6 = sum(reshape(focmd.*pr_mat,NS,NS));                             % FOC for md, when i==0

welt2c = use(u, xnewc)+use(u, subs(xnewc, zmc, zmc-mc));
welc = use(u, omega - zmc - zlc) + 0.5*(beta*sum(reshape(welt2c.*pr_mat,NS,NS)) + use(v,yc));

welt2d = use(u, xnewd)+use(u, subs(xnewd, zmd, zmd-md));
weld = use(u,  - zmd - zld) + 0.5*(beta*sum(reshape(welt2d.*pr_mat,NS,NS)) + use(v,yd));

%% if i == 0
conc = sym('conc', [1 NS], 'positive');
cond = sym('cond', [1 NS], 'positive');

eqn1 = eq1+eq2;
noweq = [eqn1, eq3, eqn4, eq5, eqn6, eq7, eq8];
noweq = subs(noweq,irate,zeros([1 NS]));
noweq = subs(noweq, qnew, q);
noweq = subs(noweq, [zmc, zmd], [omega-zlc-conc, -zld-cond]);
noweq = simplify(noweq);

%% 
% tpri = reshape(noweq, 3, []);
% part1 = tpri(1, [2 3 6]);
% zguess1 = [ 2.3961273354483329828819413207588, 0.3847908107439986611584750745268, 0.13782646603130568239102037822635];
% [z1conc1, z1mc1, z1yc1] = vpasolve(subs(part1, q, [1 1 1]), [conc(1), mc(1), yc(1)], zguess1);
% part2 = tpri(2, [2 3 6]);
% zguess2 = [ 2.3961273354483329828819413207588, 0.3847908107439986611584750745268, 0.13782646603130568239102037822635];
% [z1conc2, z1mc2, z1yc2] = vpasolve(subs(part2, q, [1 1 1]), [conc(2), mc(2), yc(2)], zguess2);
% part3 = tpri(3, [2 3 6]);
% zguess3 = [ 2.3961273354483329828819413207588, 0.3847908107439986611584750745268, 0.13782646603130568239102037822635];
% [z1conc3, z1mc3, z1yc3] = vpasolve(subs(part3, q, [1 1 1]), [conc(3), mc(3), yc(3)], zguess3);
% 
% part4 = tpri(1, [4 5 7]);
% % z1guess4 = [ 2.3961273354483329828819413207588, 0.3847908107439986611584750745268, 0.13782646603130568239102037822635];
% [z1cond1, z1md1, z1yd1] = vpasolve(subs(part4, q, [0.4840 0.4401 0.4428]), [cond(1), md(1), yd(1)], zguess1);
% 
% noweq = subs(noweq([1:3,10:15,19:21]), [conc, mc, yc], [z1conc1, z1conc2, z1conc3, z1mc1, z1mc2, z1mc3, z1yc1, z1yc2, z1yc3]);
% %%
% objf = @(x)(double(subs(noweq,[q,cond,md,yd],x)));
% lb = zeros([1 12]);
% s = lsqnonlin(objf, double([1,1,1,z1conc1, z1conc2, z1conc3, z1mc1, z1mc2, z1mc3, z1yc1, z1yc2, z1yc3]), lb);
% 
% z1guess = [ 0.05484781313739098233395746588408, 0.050702534851471585715860319482897, 0.047316140597481193815296465869318, 2.4941770382768850524501437474731, 2.5024675948487238456863380402754, 2.5092403833567046294874657475026, 0.39919472148997269616758606170949, 0.4008775581691319047153760205178, 0.40120624737391557983341318353737, 0.13719293033107384658428659053138, 0.13735984641789253643216193474488, 0.13700419852068345335561559055283];
% [z1q1,z1q2,z1q3,z1cond1,z1cond2,z1cond3,z1md1,z1md2,z1md3,z1yd1,z1yd2,z1yd3] = vpasolve(noweq,[q,cond,md,yd],z1guess);

%%
[sq1,sq2,sq3,scond1,scond2,scond3,smd1,smd2,smd3,syd1,syd2,syd3,sconc1,sconc2,sconc3,...
    smc1,smc2,smc3,syc1,syc2,syc3] = vpasolve(noweq,[q cond md yd conc mc yc], guess00);

s00 = [sq1,sq2,sq3,scond1,scond2,scond3,smd1,smd2,smd3,syd1,syd2,syd3,sconc1,sconc2,sconc3,...
    smc1,smc2,smc3,syc1,syc2,syc3];

mindelta = (pi*[smc1 smc2 smc3] + (1-pi)*[smd1 smd2 smd3])./[sq1 sq2 sq3] - 1;

welcr = simplify(subs(welc, [zmc, irate, qnew], ...
    [omega-zlc-conc, 0, 0, 0, q]));
welcr = double(subs(welcr,[q cond md yd conc mc yc],s00));

weldr = simplify(subs(weld, [zmd, irate, qnew], ...
    [-zld-cond, 0, 0, 0, q]));
weldr = double(subs(weldr,[q cond md yd conc mc yc],s00));
wel00 = [welcr, weldr];

%% i1 = 0, i2 = 0, i3>0
conc = sym('conc', [1 NS], 'positive');
cond = sym('cond', [1 NS], 'positive');

eqa = [eq1(1:2) + eq2(1:2) eq1(3) eq2(3) eq3 eq4(3) eq5 eq6(3) eq7 eq8 eqn4(1:2) eqn6(1:2)];
eqa = [eqa welc weld];

eqa = subs(eqa,irate(1:2),[0 0]);
eqa = subs(eqa, qnew, q);
eqa = subs(eqa, [zmc(1:2), zmd(1:2)], ...
    [omega-zlc(1:2)-conc(1:2), -zld(1:2)-cond(1:2)]);
eqa = subs(eqa, mc(3), zmc(3));
eqa = subs(eqa, md(3), zmd(3));
eqa = simplify(eqa);
welr2 = eqa(23:28);
eqa = eqa(1:22);

neqa = subs(eqa, delta(3), mindelta(3) - 0.01);

[s1irate2, s1conc1, s1cond1, s1zmc2, s1zlc2, s1zmd2, s1zld2, s1q1, s1q2, s1yc1, s1yc2, ...
    s1yd1, s1yd2, s1mc1, s1md1] = ...
    vpasolve(neqa, [irate(3) q cond(1:2) zld(3) md(1:2) zmd(3) yd conc(1:2) zlc(3) mc(1:2) zmc(3) yc], guess1);
s1 = [s1irate2, s1conc1, s1cond1, s1zmc2, s1zlc2, s1zmd2, s1zld2, s1q1, s1q2, s1yc1, s1yc2, ...
    s1yd1, s1yd2, s1mc1, s1md1];

wel01 = subs(welr2, [irate(2), conc(1), cond(1), zmc(2), zlc(2), ...
    zmd(2), zld(2), q, yc, yd, mc(1), md(1)],s1);

%%
% report = repmat(-2, [151 19]);
% parfor ii = 20:171
%     neqa = subs(eqa, delta(2), ii/100);
%     [s1irate2, s1conc1, s1cond1, s1zmc2, s1zlc2, s1zmd2, s1zld2, s1q1, s1q2, s1yc1, s1yc2, ...
%         s1yd1, s1yd2, s1mc1, s1md1] = ...
%     vpasolve(neqa, [irate(2), conc(1), cond(1), zmc(2), zlc(2), ...
%         zmd(2), zld(2), q, yc, yd, mc(1), md(1)], [ 0.00084081792675536975040532940339449, 14.615367029293011283748301888421, 14.415924731023672473196058005587, 2.6829607906708327276127523313331, 12.907292597772275845351494707303, 2.7130957819505741441081649417686, -17.291427333940530773926953082146, 0.48435411984165812152782005299598, 0.50596091822657597157272944912948, 0.09519147371668585614467830329925, 0.19439304937304555656865421536061, 0.095135485914521411700612083609967, 0.19426892970948519898698146894221, 1.9861012399000304028513282058772, 1.9581566461984012222391673030071]);
% 
%     s1 = [s1irate2, s1conc1, s1cond1, s1zmc2, s1zlc2, s1zmd2, s1zld2, s1q1, s1q2, s1yc1, s1yc2, ...
%         s1yd1, s1yd2, s1mc1, s1md1];
% 
%     wel01 = subs(welr2, [irate(2), conc(1), cond(1), zmc(2), zlc(2), ...
%         zmd(2), zld(2), q, yc, yd, mc(1), md(1)],s1);
%     report(ii, :) = double([s1irate2, s1conc1, s1cond1, s1zmc2, s1zlc2, s1zmd2, s1zld2, s1q1, s1q2, s1yc1, s1yc2, ...
%         s1yd1, s1yd2, s1mc1, s1md1, wel01]);
% end
% saverep = report;
% %%
% clear rep;
% report = saverep;
% report = report(~ismember(report, repmat(-2, [1 19]), 'rows'),:);
% report = report(~ismember(report, zeros([1 19]), 'rows'),:);
% rep(:,1) = (report(:,14) + report(:,15))/2;        % real balance 1
% rep(:,2) = (report(:,4) + report(:,6))/2;           % r. b. 2
% rep(:,3:4) = (report(:,10:11) + report(:,12:13))/2; % st-2 output
% rep(:,5) = rep(:,3) + (theta(1)*omega+omega)/2;    %gdp1
% rep(:,6) = rep(:,4) + (theta(2)*omega+omega)/2;    %gdp2
% rep(:,7) = (report(:,16)+report(:,17))/2;
% rep(:,8) = (report(:,18)+report(:,19))/2;
% rep(:,9) = (rep(:,7)+rep(:,8))/2;
% rep2(:,1) = (0.7*report(:,8)/(1+gma(1)) + 0.3*report(:,9)/(1+gma(2)))./report(:,8);
% rep2(:,2) = (0.3*report(:,8)/(1+gma(1)) + 0.7*report(:,9)/(1+gma(2)))./report(:,9);
% rep2 = rep2 - 1;
% rep = [rep2, rep];

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

% lastsolve = zeros([1 12]);
% for ii = 1:50
%     for kk = 7:15
%         meceq = subs(eqa(4:15), [q, irate(1)], [1 kk/10 ii/100]);
%         mguess2 = [ 14.518816696398656412906897040188, 15.521266917582283687598302035455, 2.2816862308201021959132540359903, 13.193675465651312968663143389622, 1.910792588523776650437635702663, -14.296365089389625405385585744402, 0.13371755848114330370726366701934, 0.13673226529299057933035815291063, 0.13043866817086366093963733295564, 0.12963069652458473441171958686002, 2.3179190370869255367205764488979, 2.3838684847482660651634896033397];
%         [z2conc2, z2cond2, z2zmc1, z2zlc1, z2zmd1, z2zld1, z2yc1, z2yc2, z2yd1, z2yd2, z2mc2, z2md2] = ...
%         vpasolve(meceq, [conc(2), cond(2), zmc(1), zlc(1), zmd(1), zld(1), yc, yd, mc(2), md(2)],mguess2);
%         z2 = [z2conc2, z2cond2, z2zmc1, z2zlc1, z2zmd1, z2zld1, z2yc1, z2yc2, z2yd1, z2yd2, z2mc2, z2md2];
%         z2wel = subs(welr2, [irate(1), q, [conc(2), cond(2), zmc(1), zlc(1),...
%             zmd(1), zld(1), yc, yd, mc(2), md(2)]],[ii/100 1 kk/10 z2]);
%         if length(z2)==0 || ~isreal(z2wel)
%             [z2conc2, z2cond2, z2zmc1, z2zlc1, z2zmd1, z2zld1, z2yc1, z2yc2, z2yd1, z2yd2, z2mc2, z2md2] = ...
%             vpasolve(meceq, [conc(2), cond(2), zmc(1), zlc(1), zmd(1), zld(1), yc, yd, mc(2), md(2)],lastsolve);
%             z2 = [z2conc2, z2cond2, z2zmc1, z2zlc1, z2zmd1, z2zld1, z2yc1, z2yc2, z2yd1, z2yd2, z2mc2, z2md2];            
%         end
%         if length(z2)==0 || ~isreal(z2wel)
%             disp([ii,kk]);          
%         else
%             lastsolve = z2;
%             z2alloc10(ii,kk-6,:) = z2;
%             z2wel10(ii,kk-6,:) = z2wel;
%         end
%     end
% end
% [xx,yy] = meshgrid(7:15,1:50);
% z2wel10(find(z2wel10==0)) = nan;
% z2wel10(:,:,5) = (z2wel10(:,:,1)+z2wel10(:,:,2))/2;
% z2wel10(:,:,6) = (z2wel10(:,:,3)+z2wel10(:,:,4))/2;
% z2wel10(:,:,7) = (z2wel10(:,:,5)+z2wel10(:,:,6))/2;
% mesh(xx,yy,z2wel10(:,:,5)-mean2(z2wel10(:,:,5))); hold on;
% mesh(xx,yy,z2wel10(:,:,6)-mean2(z2wel10(:,:,6))); hold on;
% mesh(xx,yy,z2wel10(:,:,7)-mean2(z2wel10(:,:,7))); hold off;
% legend('c','d','all');

% nconc2 = sidx(solve(eqa(5), conc(2)), 1);  % smaller
% nzlc1 = sidx(solve(eqa(4), zlc(1)), 2);    % bigger
% nzmc1 = sidx(solve(eqa(6), zmc(1)), 2);    % bigger, so positive
% 
% nmc2 = solve(eqa(14), mc(2));
% 
% eqa(5) = nconc2 - conc(2);
% eqa(4) = nzlc1 - zlc(1);
% eqa(6) = nzmc1 - zmc(1);
% eqa(14) = nmc2 - mc(2);

eqa_now = subs(eqa, delta(1), mindelta1 - 0.005);
[s2irate1, s2conc2, s2cond2, s2zmc1, s2zlc1, s2zmd1, s2zld1, s2q1, s2q2, s2yc1, s2yc2, ...
    s2yd1, s2yd2, s2mc2, s2md2] = ...
vpasolve(eqa_now, [irate(1), conc(2), cond(2), zmc(1), zlc(1), ...
    zmd(1), zld(1), q, yc, yd, mc(2), md(2)],guess2);
s2 = [s2irate1, s2conc2, s2cond2, s2zmc1, s2zlc1, s2zmd1, s2zld1, s2q1, s2q2, s2yc1, s2yc2, ...
    s2yd1, s2yd2, s2mc2, s2md2];

wel10 = subs(welr2, [irate(1), conc(2), cond(2), zmc(1), zlc(1), ...
    zmd(1), zld(1), q, yc, yd, mc(2), md(2)],s2);
% 
% report2 = repmat(-2, [1 19]);
% parfor ii = 1:57
%     eqa_now = subs(eqa, delta(1), ii/100);
%     [s2irate1, s2conc2, s2cond2, s2zmc1, s2zlc1, s2zmd1, s2zld1, s2q1, s2q2, s2yc1, s2yc2, ...
%         s2yd1, s2yd2, s2mc2, s2md2] = ...
%     vpasolve(eqa_now, [irate(1), conc(2), cond(2), zmc(1), zlc(1), ...
%         zmd(1), zld(1), q, yc, yd, mc(2), md(2)],[ 0.00040591635374116399079628081950222, 14.518816696398656412906897040188, 13.872725335519556168614844775128, 2.3164320737315599356018322774018, 13.164513831250045756194515237321, 1.9706154316542557251408031415962, -14.720957830020692206273051663327, 1.3653017533075846053320494964962, 0.80422898404089370923912909234216, 0.13660833294436349965044767307515, 0.13673226529299057933035815291063, 0.13078873324027073807007166530579, 0.12186069885165503830868361871485, 2.3179190370869255367205764488979, 2.0373037743609396405487375039048]);
%     s2 = [s2irate1, s2conc2, s2cond2, s2zmc1, s2zlc1, s2zmd1, s2zld1, s2q1, s2q2, s2yc1, s2yc2, ...
%         s2yd1, s2yd2, s2mc2, s2md2];
% 
%     wel10 = subs(welr2, [irate(1), conc(2), cond(2), zmc(1), zlc(1), ...
%         zmd(1), zld(1), q, yc, yd, mc(2), md(2)],s2);
%     report2(ii, :) = double([s2, wel10]);
% end
% saverep2 = report2;
% %%
% clear rep3;
% report2 = saverep2;
% report2 = report2(~ismember(report2, repmat(-2, [1 19]), 'rows'),:);
% report2 = report2(~ismember(report2, zeros([1 19]), 'rows'),:);
% rep3(:,1) = (report2(:,4) + report2(:,6))/2;           % r. b. 1
% rep3(:,2) = (report2(:,14) + report2(:,15))/2;        % real balance 2
% rep3(:,3:4) = (report2(:,10:11) + report2(:,12:13))/2; % st-2 output
% rep3(:,5) = rep3(:,3) + (theta(1)*omega+omega)/2;    %gdp1
% rep3(:,6) = rep3(:,4) + (theta(2)*omega+omega)/2;    %gdp2
% rep3(:,7) = (report2(:,16)+report2(:,17))/2;
% rep3(:,8) = (report2(:,18)+report2(:,19))/2;
% rep3(:,9) = (rep3(:,7)+rep3(:,8))/2;
% rep4(:,1) = (0.7*report2(:,8)/(1+gma(1)) + 0.3*report2(:,9)/(1+gma(2)))./report2(:,8);
% rep4(:,2) = (0.3*report2(:,8)/(1+gma(1)) + 0.7*report2(:,9)/(1+gma(2)))./report2(:,9);
% rep4 = rep4 - 1;
% rep3 = [rep4, rep3];

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



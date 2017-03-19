clear;
beta = 1;
gma = [0 0];
gma_new = [0 0];

pi = 0.5; % fix creditor's fraction.
omega = 5;
% theta = sym('theta',[1 2],'positive');
theta = [1 1.1];

guess00 = ...
[ 2.4149281933810972699369717206347, 2.4123698075423505085169704129497, 2.5066082890310675711845706798981, 2.5198438446549838201942167607317, 0.039231758793917579439228799733619, 0.033893173901332835644406413159319, 0.13782646603130568239102037822635, 0.13782646603130568239102037822635, 0.1370683654284613110711383757952, 0.13676638583843601348038320036835, 0.38781001479864357005512049564093, 0.3873991671168349007409172774835, 0.40091860933248548272172119799857, 0.40239303267965088168172860765563]
guess1 = ...
[ 0.000042886531202840368879256530145945, 2.4149337224715125356117793335052, 2.5065796848296869229397475686883, 0.38737352351942472871549321413946, 2.2002473801101968205501811475799, 0.40236175759298409983111340848453, -2.9221720983893608433325056403922, 0.03924329634940027072423654890327, 0.033905281416622402882141064905858, 0.13782646603130568239102037822635, 0.13781317116024876659977817379549, 0.13706895657594706738786682976246, 0.13675424435801729888934994574063, 0.38781090270770920778439454296072, 0.40091529458138076084236380173326]
guess2 = ...
[ 0.000031104542582347677420808051923517, 2.4123762974764241629116237948576, 2.5198102696241473926290338088787, 0.38779203229499093297344552630858, 2.197269068705858009886055082591, 0.4008903419979219973480159827991, -2.9074534811100180109174798836005, 0.039248980944376464645018354049097, 0.03390671644971422222967119813182, 0.1378168234378101667547931029553, 0.13782646603130568239102037822635, 0.13705845461950390554710951901956, 0.13676578820597293887630261597066, 0.38740020932646914345510359671693, 0.40238639327069721439046486360745]

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

xnewd = kron([1 1], qnew) + kron(omega.*theta, [1 1])...
    + kron((zmd+zld.*(1.+irate))./q, qnew./(1+gma_new));
xnewc = kron([1 1], qnew) + kron((zmc+zlc.*(1.+irate))./q, qnew./(1+gma_new));

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
% 
% [zconc1, zconc2] = solve(noweq([5,6]), conc);
% 
% nconc1 = sidx(solve(noweq(3), conc(1)),1);
% nconc2 = sidx(solve(noweq(4), conc(2)),1);
% 
% noweq(3) = conc(1) - nconc1;
% noweq(4) = conc(2) - nconc2;
% noweq(5) = conc(1) - zconc1;
% noweq(6) = conc(2) - zconc2;

[sconc1,sconc2,scond1,scond2,sq1,sq2,syc1,syc2,syd1,syd2,smc1,smc2,smd1,smd2] = ...
    vpasolve(noweq,[conc, cond, q, yc, yd, mc, md], guess00);
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


% for ii = 1:50
%     for kk = 7:15
%         meceq = subs(eqa(4:15), [q, irate(2)], [1 kk/10 ii/100]);
%         mguess1 = [ 14.518816696398656412906897040188, 12.753164276949615230189953130185, 2.3126163266185756635206439650412, 13.167719791035983305479487807889, 2.0336501247729826621059453774297, -15.904257927097334843836969986482, 0.13673226529299057933035815291063, 0.13629038776488363002450503477062, 0.13092702578064287122563262629679, 0.12158276851846112621526361091364, 2.3179190370869255367205764488979, 1.9725625998425480857622955043783];
%         [z1conc1, z1cond1, z1zmc2, z1zlc2, z1zmd2, z1zld2, z1yc1, z1yc2, z1yd1, z1yd2, z1mc1, z1md1] = ...
%         vpasolve(meceq, [conc(1), cond(1), zmc(2), zlc(2), zmd(2), zld(2), yc, yd, mc(1), md(1)],mguess1);
%         z1 = [z1conc1, z1cond1, z1zmc2, z1zlc2, z1zmd2, z1zld2, z1yc1, z1yc2, z1yd1, z1yd2, z1mc1, z1md1];
%         z1wel = subs(welr2, [irate(2), q, conc(1), cond(1), zmc(2), zlc(2), ...
%             zmd(2), zld(2), yc, yd, mc(1), md(1)],[ii/100 1 kk/10 z1]);
%         z1alloc01(ii,kk-6,:) = z1;
%         z1wel01(ii,kk-6,:) = z1wel;
%     end
% end
% [xx,yy] = meshgrid(7:15,1:50);
% z1wel01(:,:,5) = (z1wel01(:,:,1)+z1wel01(:,:,2))/2;
% z1wel01(:,:,6) = (z1wel01(:,:,3)+z1wel01(:,:,4))/2;
% z1wel01(:,:,7) = (z1wel01(:,:,5)+z1wel01(:,:,6))/2;
% mesh(xx,yy,z1wel01(:,:,5)-mean2(z1wel01(:,:,5))); hold on;
% mesh(xx,yy,z1wel01(:,:,6)-mean2(z1wel01(:,:,6))); hold on;
% mesh(xx,yy,z1wel01(:,:,7)-mean2(z1wel01(:,:,7))); hold off;
% legend('c','d','all');

% nconc1 = sidx(solve(eqa(4), conc(1)), 1);  % smaller
% nzlc2 = sidx(solve(eqa(5), zlc(2)), 2);    % bigger
% nzmc2 = sidx(solve(eqa(6), zmc(2)), 2);    % bigger, so positive
% nmc1 = solve(eqa(14), mc(1));
% 
% eqa(4) = nconc1 - conc(1);
% eqa(5) = nzlc2 - zlc(2);
% eqa(6) = nzmc2 - zmc(2);
% eqa(14) = nmc1 - mc(1);

neqa = subs(eqa, delta(2), mindelta2 - 0.005);
[s1irate2, s1conc1, s1cond1, s1zmc2, s1zlc2, s1zmd2, s1zld2, s1q1, s1q2, s1yc1, s1yc2, ...
    s1yd1, s1yd2, s1mc1, s1md1] = ...
vpasolve(neqa, [irate(2), conc(1), cond(1), zmc(2), zlc(2), ...
    zmd(2), zld(2), q, yc, yd, mc(1), md(1)], guess1);
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



clear;
load('report.mat');
report = report(~ismember(report, repmat(-2,[1 19]), 'rows'),:);

% [irate(2), conc(1), cond(1), zmc(2), zlc(2), zmd(2), zld(2), q, yc, yd, mc(1), md(1)]
omega = 30;
theta = [0.8,1.25];
gma = [0,0];


rep(:,1) = (report(:,14) + report(:,15))/2;        % real balance 1
rep(:,2) = (report(:,4) + report(:,6))/2;           % r. b. 2
rep(:,3:4) = (report(:,10:11) + report(:,12:13))/2; % st-2 output
rep(:,5) = rep(:,3) + (theta(1)*omega+omega)/2;    %gdp1
rep(:,6) = rep(:,4) + (theta(2)*omega+omega)/2;    %gdp2
rep(:,7:8) = (report(:,16:17) + report(:,18:19))/2;
rep(:,9) = (rep(:,7)+rep(:,8))/2;
rep2(:,1) = (0.7*report(:,8)/(1+gma(1)) + 0.3*report(:,9)/(1+gma(2)))./report(:,8);
rep2(:,2) = (0.3*report(:,8)/(1+gma(1)) + 0.7*report(:,9)/(1+gma(2)))./report(:,9);
rep2 = rep2 - 1;

rep = [rep2, rep];
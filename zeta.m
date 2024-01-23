clc;
close all;
clearvars;

%%
RE=input('RE: ');
MI=input('MI: ');
s=RE+1i*MI;

t=input('t: ');

n=1e4;
addend=NaN(n,1);
for k=1:n
    addend(k,1)=1/(k^s);
end

RE_addend=real(addend);
signum_RE_addend=sign(RE_addend);
cumulative_signum_RE_addend=cumsum(signum_RE_addend);

MI_addend=imag(addend);
signum_MI_addend=sign(MI_addend);
cumulative_signum_MI_addend=cumsum(signum_MI_addend);

signum_addend=sign(addend);
RE_signum_addend=real(signum_addend);
MI_signum_addend=imag(signum_addend);

cumulative_signum_addend=cumsum(signum_addend);
RE_cumulative_signum_addend=real(cumulative_signum_addend);
MI_cumulative_signum_addend=imag(cumulative_signum_addend);

%%
s_name=num2str(s,'%.6f');
folder=fullfile(pwd,['zeta @ ',num2str(t,'%.6f'),' @ ',s_name]);
warning off;
mkdir(folder);
warning on;

figure;
hold on;
plot(addend,"k");
plot(real(addend(1,1)),imag(addend(1,1)),"rh");
plot(real(addend(end,1)),imag(addend(end,1)),"gh");
hold off;
xlabel('Real axis');
ylabel('Imaginary axis');
title(['Riemann \zeta(s) function addends for s=',s_name]);
savefig(fullfile(folder,'Addend.fig'));

figure;
hold on;
plot(signum_RE_addend);
plot(signum_MI_addend);
plot(RE_signum_addend);
plot(MI_signum_addend);
hold off;
ylim("padded");
legend(...
    'signum@RE@addend',...
    'signum@MI@addend',...
    'RE@signum@addend',...
    'MI@signum@addend');
xlabel('n (addend index)');
ylabel('Signum');
title(['Addend signum for s=',s_name]);
savefig(fullfile(folder,'Addend signum.fig'));

figure;
hold on;
plot(cumulative_signum_RE_addend);
plot(cumulative_signum_MI_addend);
plot(RE_cumulative_signum_addend);
plot(MI_cumulative_signum_addend);
hold off;
ylim("padded");
legend(...
    'cumulative@signum@RE@addend',...
    'cumulative@signum@MI@addend',...
    'RE@cumulative@signum@addend',...
    'MI@cumulative@signum@addend');
xlabel('n (addend index)');
ylabel('Cumulative signum');
title(['Cumulative addend signum for s=',s_name]);
savefig(fullfile(folder,'Cumulative addend signum.fig'));

figure;
hold on;
plot(signum_addend,"k");
plot(real(signum_addend(1,1)),imag(signum_addend(1,1)),"rh");
plot(real(signum_addend(end,1)),imag(signum_addend(end,1)),"gh");
hold off;
axis square;
xlabel('Real axis');
ylabel('Imaginary axis');
title(['Addend signum for s=',s_name]);
savefig(fullfile(folder,'Addend signum_2.fig'));

figure;
hold on;
plot(cumulative_signum_addend,"k");
plot(real(cumulative_signum_addend(1,1)),imag(cumulative_signum_addend(1,1)),"rh");
plot(real(cumulative_signum_addend(end,1)),imag(cumulative_signum_addend(end,1)),"gh");
hold off;
axis square;
xlabel('Real axis');
ylabel('Imaginary axis');
title(['Cumulative addend signum for s=',s_name]);
savefig(fullfile(folder,'Cumulative addend signum_2.fig'));

figure;
hold on;
plot(RE_addend(1:100,1));
plot(MI_addend(1:100,1));
hold off;
ylim("padded");
legend(...
    'RE@addend',...
    'MI@addend');
xlabel('n (addend index)');
ylabel('Part');
title(['Addend part for s=',s_name]);
savefig(fullfile(folder,'Addend part.fig'));

%%
n=14;
x0=(-MI*(1/t)*log(1:n))';
gamma=(RE*(1/abs(t))*log(1:n))';

s_functional=1-s;
RE_functional=real(s_functional);
MI_functional=imag(s_functional);
x0_functional=(-MI_functional*(1/t)*log(1:n))';
gamma_functional=(RE_functional*(1/abs(t))*log(1:n))';
s_name_functional=num2str(s_functional,'%.6f');

granularity=-1000:0.1:1000;

Cauchy_cumulative_distribution_function=NaN(n,size(granularity,2));
Cauchy_cumulative_distribution_function(1,granularity<0)=0;
Cauchy_cumulative_distribution_function(1,granularity>=0)=1;
Cauchy_cumulative_distribution_function_functional=NaN(n,size(granularity,2));
Cauchy_cumulative_distribution_function_functional(1,granularity<0)=0;
Cauchy_cumulative_distribution_function_functional(1,granularity>=0)=1;
for k=2:n
    Cauchy_cumulative_distribution_function(k,:)=(1/pi)*atan((granularity-x0(k))/gamma(k))+0.5;
    Cauchy_cumulative_distribution_function_functional(k,:)=(1/pi)*atan((granularity-x0_functional(k))/gamma_functional(k))+0.5;
end

alternating_Cauchy_cumulative_distribution_function=NaN(n,size(granularity,2));
alternating_Cauchy_cumulative_distribution_function_functional=NaN(n,size(granularity,2));
for k=1:n
    alternating_Cauchy_cumulative_distribution_function(k,:)=Cauchy_cumulative_distribution_function(k,:)*((-1)^(k+1));
    alternating_Cauchy_cumulative_distribution_function_functional(k,:)=Cauchy_cumulative_distribution_function_functional(k,:)*((-1)^(k+1));
end

figure;
hold on;
xline(0,'k');
yline(0,'k');
yline(0.5,'m:');
yline(-0.5,'m:');
yline(1,'m:');
yline(-1,'m:');
plot(granularity(granularity<0),alternating_Cauchy_cumulative_distribution_function(1,granularity<0),'k');
plot(granularity(granularity>=0),alternating_Cauchy_cumulative_distribution_function(1,granularity>=0),'k');
for k=2:n
    plot(granularity,alternating_Cauchy_cumulative_distribution_function(k,:),'k');
    plot(granularity,alternating_Cauchy_cumulative_distribution_function_functional(k,:),'r');
end
hold off;
ylim("padded");
xlabel('x');
ylabel('y');
legend(...
    '','','','','','','','',...
    'Alternating Cauchy CDF for s',...
    'Alternating Cauchy CDF for 1-s');
title(['Alternating Cauchy CDF for s=',s_name]);
savefig(fullfile(folder,'Alternating Cauchy CDF.fig'));

Cauchy_probability_density_function=NaN(n,size(granularity,2));
Cauchy_probability_density_function(1,granularity~=0)=0;
Cauchy_probability_density_function_functional=NaN(n,size(granularity,2));
Cauchy_probability_density_function_functional(1,granularity~=0)=0;
for k=2:n
    Cauchy_probability_density_function(k,:)=1./(pi*gamma(k)*(1+((granularity-x0(k))/gamma(k)).^2));
    Cauchy_probability_density_function_functional(k,:)=1./(pi*gamma_functional(k)*(1+((granularity-x0_functional(k))/gamma_functional(k)).^2));
end

alternating_Cauchy_probability_density_function=NaN(n,size(granularity,2));
alternating_Cauchy_probability_density_function_functional=NaN(n,size(granularity,2));
for k=1:n
    alternating_Cauchy_probability_density_function(k,:)=Cauchy_probability_density_function(k,:)*((-1)^(k+1));
    alternating_Cauchy_probability_density_function_functional(k,:)=Cauchy_probability_density_function_functional(k,:)*((-1)^(k+1));
end

Cauchy_amplitude=1./(pi*gamma);
Hippias_radius=1./(2*gamma);
Cauchy_amplitude_functional=1./(pi*gamma_functional);
Hippias_radius_functional=1./(2*gamma_functional);

figure;
hold on;
xline(0,'k');
for k=1:n
    plot(granularity,alternating_Cauchy_probability_density_function(k,:),'k');
end
plot(x0,Cauchy_amplitude,'k--');
plot(x0,-Cauchy_amplitude,'k--');
plot(x0,Hippias_radius,'k.');
for k=1:n
    plot(granularity,alternating_Cauchy_probability_density_function_functional(k,:),'r');
end
plot(x0_functional,Cauchy_amplitude_functional,'r--');
plot(x0_functional,-Cauchy_amplitude_functional,'r--');
plot(x0_functional,Hippias_radius_functional,'r.');
hold off;
ylim("padded");
xlabel('x');
ylabel('y');
legend(...
    '',...
    '','','','','','','','','','','','','','','Alternating Cauchy amplitude for s','','Hippias radius for s',...
    '','','','','','','','','','','','','','','Alternating Cauchy amplitude for 1-s','','Hippias radius for 1-s');
title(['Alternating Cauchy PDF for s=',s_name]);
savefig(fullfile(folder,'Alternating Cauchy PDF.fig'));

angular_granularity=(0:0.1:100)*(pi/200);
angular_granularity=angular_granularity(1,1:end-1);

versiera_x=NaN(n,size(angular_granularity,2));
versiera_y=NaN(n,size(angular_granularity,2));
versiera_x_functional=NaN(n,size(angular_granularity,2));
versiera_y_functional=NaN(n,size(angular_granularity,2));
for k=1:n
    versiera_x(k,:)=gamma(k)*tan(angular_granularity);
    versiera_y(k,:)=gamma(k)*((cos(angular_granularity).^2));
    versiera_x_functional(k,:)=gamma_functional(k)*tan(angular_granularity);
    versiera_y_functional(k,:)=gamma_functional(k)*((cos(angular_granularity).^2));
end
inflection_x=gamma/sqrt(3);
inflection_y=gamma*0.75;
inflection_x_functional=gamma_functional/sqrt(3);
inflection_y_functional=gamma_functional*0.75;

inflection_x_reference=(inflection_x(end)/RE)*0.5;
inflection_y_reference=(inflection_y(end)/RE)*0.5;
inflection_x_reference_functional=(inflection_x_functional(end)/RE_functional)*0.5;
inflection_y_reference_functional=(inflection_y_functional(end)/RE_functional)*0.5;

granularity=0:0.1:1000;

versiera_secant=NaN(n,size(granularity,2));
versiera_secant_functional=NaN(n,size(granularity,2));
for k=1:n
    versiera_secant(k,:)=-0.25*sqrt(3)*granularity+gamma(k);
    versiera_secant_functional(k,:)=-0.25*sqrt(3)*granularity+gamma_functional(k);
end

figure;
tiledlayout(1,2);

nexttile;
hold on;
for k=1:n
    plot(versiera_x(k,:),versiera_y(k,:),'k');
    plot(granularity,versiera_secant(k,:),'r:');
end
plot(inflection_x,inflection_y,'r');
plot(inflection_x(end),inflection_y(end),'r*');
plot(inflection_x_reference,inflection_y_reference,'co');
plot(gamma,gamma,'m');
plot(gamma,gamma/2,'g');
plot(granularity,tan(pi/3)*granularity,'b');
xline(gamma(end),'k');
yline(gamma(end),'k');
hold off;
axis square;
if gamma(end)>=gamma_functional(end)
    xlim([0 1.5*gamma(end)]);
    ylim([0 1.5*gamma(end)]);
else
    xlim([0 1.5*gamma_functional(end)]);
    ylim([0 1.5*gamma_functional(end)]);
end
xlabel('x');
ylabel('y');
legend(...
    '','','','','','','','','','','','','','',...
    '','','','','','','','','','','','',...
    'Agnesi versiera for s',...
    'versiera secant',...
    'versiera inflection',...
    ['versiera inflection at n=',num2str(n)],...
    ['versiera inflection at n=',num2str(n),', RE=1/2'],...
    'y=x',...
    'y=x/2',...
    '\theta=\pi/3');
title(['Agnesi versiera for s=',s_name]);

nexttile;
hold on;
for k=1:n
    plot(versiera_x_functional(k,:),versiera_y_functional(k,:),'k');
    plot(granularity,versiera_secant_functional(k,:),'r:');
end
plot(inflection_x_functional,inflection_y_functional,'r');
plot(inflection_x_functional(end),inflection_y_functional(end),'r*');
plot(inflection_x_reference_functional,inflection_y_reference_functional,'co');
plot(gamma_functional,gamma_functional,'m');
plot(gamma_functional,gamma_functional/2,'g');
plot(granularity,tan(pi/3)*granularity,'b');
xline(gamma_functional(end),'k');
yline(gamma_functional(end),'k');
hold off;
axis square;
if gamma(end)>=gamma_functional(end)
    xlim([0 1.5*gamma(end)]);
    ylim([0 1.5*gamma(end)]);
else
    xlim([0 1.5*gamma_functional(end)]);
    ylim([0 1.5*gamma_functional(end)]);
end
xlabel('x');
ylabel('y');
legend(...
    '','','','','','','','','','','','','','',...
    '','','','','','','','','','','','',...
    'Agnesi versiera for 1-s',...
    'versiera secant',...
    'versiera inflection',...
    ['versiera inflection at n=',num2str(n)],...
    ['versiera inflection at n=',num2str(n),', RE=1/2'],...
    'y=x',...
    'y=x/2',...
    '\theta=\pi/3');
title(['Agnesi versiera for 1-s=',s_name_functional]);

savefig(fullfile(folder,'Agnesi versiera.fig'));

%%
n=100;

s_name_reference='1/2+14.134725i';

gamma=(RE*(1/abs(t))*log(1:n))';
gamma_functional=(RE_functional*(1/abs(t))*log(1:n))';
gamma_reference=(0.5*(1/abs(t))*log(1:n))';

x0=(-MI*(1/t)*log(1:n))';
x0_functional=(-MI_functional*(1/t)*log(1:n))';
x0_reference=(-14.134725*(1/t)*log(1:n))';

x0_abs=abs(x0);
x0_abs_functional=abs(x0_functional);
x0_abs_reference=abs(x0_reference);

variance_lower_bound=2*(gamma.^2);
variance_lower_bound_functional=2*(gamma_functional.^2);
variance_lower_bound_reference=2*(gamma_reference.^2);

STD=sqrt(variance_lower_bound);
STD_functional=sqrt(variance_lower_bound_functional);
STD_reference=sqrt(variance_lower_bound_reference);

figure;
tiledlayout(2,2);

nexttile;
hold on;
plot(variance_lower_bound,'ko');
plot(variance_lower_bound_reference,'c*');
plot(variance_lower_bound_functional,'r.');
hold off;
axis square;
xlabel('n (addend index)');
ylabel('Variance lower bound');
legend(...
    'for s',...
    ['for ',s_name_reference],...
    'for 1-s',...
    'Location','eastoutside');
title(['Variance lower bound for s=',s_name]);

nexttile;
hold on;
plot(x0_abs,variance_lower_bound,'ko');
plot(x0_abs_reference,variance_lower_bound_reference,'c*');
plot(x0_abs_functional,variance_lower_bound_functional,'r.');
hold off;
axis square;
xlim([0 n]);
xlabel('|x_0|');
ylabel('Variance lower bound');
legend(...
    'for s',...
    ['for ',s_name_reference],...
    'for 1-s',...
    'Location','eastoutside');
title(['Variance lower bound for s=',s_name]);

nexttile;
hold on;
plot(STD,'ko');
plot(STD_reference,'c*');
plot(STD_functional,'r.');
hold off;
axis square;
xlabel('n (addend index)');
ylabel('STD lower bound');
legend(...
    'for s',...
    ['for ',s_name_reference],...
    'for 1-s',...
    'Location','eastoutside');
title(['STD lower bound for s=',s_name]);

nexttile;
hold on;
plot(x0_abs,STD,'ko');
plot(x0_abs_reference,STD_reference,'c*');
plot(x0_abs_functional,STD_functional,'r.');
hold off;
axis square;
xlim([0 n]);
xlabel('|x_0|');
ylabel('STD lower bound');
legend(...
    'for s',...
    ['for ',s_name_reference],...
    'for 1-s',...
    'Location','eastoutside');
title(['STD lower bound for s=',s_name]);

savefig(fullfile(folder,'Variance-STD.fig'));

%%
numbers=(1:n);

STD_score_n=false(n,n);
alternating_STD_score_n=zeros(n,n);
STD_score_x0_abs=false(n,n);
alternating_STD_score_x0_abs=zeros(n,n);
for k=1:n
    STD_score_n(k,:)=(numbers>=(k-STD(k)))&(numbers<=(k+STD(k)));
    alternating_STD_score_n(k,STD_score_n(k,:))=((-1)^(k+1));
    STD_score_x0_abs(k,:)=(numbers>=(x0_abs(k)-STD(k)))&(numbers<=(x0_abs(k)+STD(k)));
    alternating_STD_score_x0_abs(k,STD_score_x0_abs(k,:))=((-1)^(k+1));
end

STD_score_n=double(STD_score_n);
STD_count_n=sum(STD_score_n,1);
alternating_STD_count_n=sum(alternating_STD_score_n,1);

STD_score_x0_abs=double(STD_score_x0_abs);
STD_count_x0_abs=sum(STD_score_x0_abs,1);
alternating_STD_count_x0_abs=sum(alternating_STD_score_x0_abs,1);

%%
cd(folder);
clear alternating_* angular_granularity Cauchy_* folder granularity inflection_* k n numbers s_name s_name_functional versiera_*;
save;
cd ..;
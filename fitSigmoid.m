function [fsig,fitx,y,itd_thresh,beta,MSE] = fitSigmoid(x,y,DPT,pflg)
%This function fits a sigmoid through a set of anchor points [x,y], and 
%returns a threshold corresponding to DPT

%x: should be in units of log(ITD), ex. x = log([10 20 40 80 160 320])
%y: the corresponding d-prime estimate for x
%DPT: target d-prime
%pflg: if 1 plot, if 0 do not plot

%fsig: interpolated y with a resolution of 1 us
%fitx: interpolated x-axis
%y: original y input with cutoff applied
%itd_thresh: reulting ITD threshold corresponding to DPT
%beta: parameter estimates for the sigmoid
%MSE: error of fit

if nargin < 4 && nargout == 5
    disp('Too many output arguments')
end

%Define sigmoid following eq.3 from Moncada-Torres et al. (2018)
sfun = @(b,x) b(1) + (b(2)-b(1))./(1+10.^((b(3)-x).*b(4)));

%Define starting values for nonlinear regression
llim = 0; %lower limit
ulim = 4.65; %upper limit
p50 = mean(x); %50% point
slp = (max(y) - min(y))/(max(x) - min(x)); %slope
beta0 = [llim ulim p50 slp]; 

%Apply limits to y-values
for i = 1:length(y)
    if y(i) < llim
        y(i) = llim;
    elseif y(i) > ulim
        y(i) = ulim;
    end
end

%Correct non-monotonicity at large delays
[mx_val,mx_ind] = max(y);
for i = mx_ind:length(y)
    y(i) = mx_val;
end

[beta,~,~,~,MSE] = nlinfit(x,y,sfun,beta0);

fitx = log([exp(x(1)):1:exp(x(end))]);
for i = 1:length(fitx)
    fsig(i) = beta(1) + (beta(2)-beta(1))./(1+10.^(beta(3)-fitx(i)).*beta(4)); 
end


[~,ind] = min(abs(fsig - DPT));
itd_thresh = exp(fitx(ind));


if pflg == 1
    figure
    plot(fitx,fsig,'b','linewidth',1.5)
    hold on
    plot(x,y,'ko','markerfacecolor','k','markersize',10)
    xticks(x)
    xlbls = exp(x);
    xticklabels({num2str(xlbls(1)),num2str(xlbls(2)),num2str(xlbls(3)),...
        num2str(xlbls(4)),num2str(xlbls(5)),num2str(xlbls(6))})
    xlabel('ITD \mus')
    ylabel('d-prime')
    set(gca,'fontsize',14)
end

end

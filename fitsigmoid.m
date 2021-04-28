function [beta,MSE,fsig,fsx,y,itd_thresh] = fitsigmoid(x,y,pflg,DPT)
%x should be in units of log(ITD), ex. x = log([10 20 40 80 160 320])
%y is the corresponding d-prime estimate
%beta are the parameter estimates for the sigmoid
%fsig is the fitted sigmoid with a resolution of 1 us
%DPT is the target d-prime

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

%beta0(4) = slope

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

% %reimpose limits if necessary
% if beta(1) < beta0(1)
%     beta(1) = beta0(1);
% end
% 
% if beta(2) > beta0(2)
%     beta(2) = beta0(2);
% end

fsx = log([exp(x(1)):1:exp(x(end))]);
for i = 1:length(fsx)
    fsig(i) = beta(1) + (beta(2)-beta(1))./(1+10.^(beta(3)-fsx(i)).*beta(4)); 
end

if nargin == 4
    [~,ind] = min(abs(fsig - DPT));
    itd_thresh = exp(fsx(ind));
end

if pflg == 1
    figure
    plot(fsx,fsig,'b','linewidth',1.5)
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

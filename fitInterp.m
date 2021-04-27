function [fsig,fsx,y,itd_thresh] = fitInterp(x,y,pflg,DPT)
%x should be in units of log(ITD), ex. x = log([10 20 40 80 160 320])
%y is the corresponding d-prime estimate
%beta are the parameter estimates for the sigmoid
%fsig is the fitted sigmoid with a resolution of 1 us
%DPT is the target d-prime

if nargin < 4 && nargout == 5
    disp('Too many output arguments')
end

%Apply limits to y-values
llim = 0; %lower limit
ulim = 4.65; %upper limit
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

fsx = log([exp(x(1)):1:exp(x(end))]);

fsig = interp1(x,y,fsx,'pchip');

if nargin == 4
    [val,ind] = min(abs(fsig - DPT));
    if val > 0.1
        ind = length(fsig);
    end
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

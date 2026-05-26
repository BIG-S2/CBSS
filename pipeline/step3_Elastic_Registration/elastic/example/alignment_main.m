% Align two functions
load growth_female_vel.mat
lambda = 0;
option.parallel = 0;
option.closepool = 0;
option.showplot = 0;
option.smooth = 0;
option.sparam = 5;
if option.parallel == 1
    if matlabpool('size') == 0
        matlabpool open
    end
end
% Parameters

%fprintf('\n lambda = %5.1f \n', lambda);

a = size(t,1);
if (a ~=1)
    t = t';
end
binsize = mean(diff(t));
[M, N] = size(f);
f0 = f;

if option.smooth == 1
    for r = 1:option.sparam
        for i = 1:N
            f(2:(M-1),i) = (f(1:(M-2),i)+2*f(2:(M-1),i) + f(3:M,i))/4;
        end
    end
end

if option.showplot == 1
    figure(1); clf;
    plot(t, f, 'linewidth', 1);
    title('Original data', 'fontsize', 16);
    pause(0.1);
end

% Compute the q-function of the plot
[~,fy] = gradient(f,binsize,binsize);
q = fy./sqrt(abs(fy)+eps);
q1 = q(:,1)';
f_1 = f(:,1);
q2 = q(:,end-1)';
f_2 = f(:,end-1);

% align two functions;
[G,T] = DynamicProgrammingQ2(q1/norm(q1),t,q2/norm(q2),t,t,t);
gam0 = interp1(T,G,t);
gam_f = (gam0-gam0(1))/(gam0(end)-gam0(1));
f_2_aligned = interp1(t, f_2, (t(end)-t(1)).*gam_f + t(1))';


figure;
plot(gam_f,'linewidth',2);
set(gca,'FontSize',14);
title('\gamma','fontsize',16);

figure;
plot(t,f_1,'linewidth',1);
hold on, plot(t,f_2,'g','linewidth',1);
hold on, plot(t,f_2_aligned,'--r','linewidth',2);
legend('f1','f2','aligned f2');
set(gca,'FontSize',14);
title('Align f2 to f1','fontsize',16);



%% Problem 8
% find the center of the orbit [f] with respect to a given set of other
% functions
mnq = mean(q,2);
dqq = sqrt(sum((q - mnq*ones(1,N)).^2,1));
[~, min_ind] = min(dqq);
mq = q(:,min_ind);
mf = f(:,min_ind);

gam = zeros(N,size(q,1));
parfor k = 1:N
    q_c = q(:,k,1)'; mq_c = mq';
    [G,T] = DynamicProgrammingQ2(mq_c/norm(mq_c),t,q_c/norm(q_c),t,t,t);
    gam0 = interp1(T,G,t);
    gam(k,:) = (gam0-gam0(1))/(gam0(end)-gam0(1));  % slight change on scale
end
gamI = SqrtMeanInverse(gam);
gamI_dev = gradient(gamI, 1/(M-1));
mf_aligned = interp1(t, mf, (t(end)-t(1)).*gamI + t(1))';
mq_aligned = gradient(mf, binsize)./sqrt(abs(gradient(mf, binsize))+eps);

% gamma to find the center of the orbit;
figure;
plot(gamI,'linewidth',2);
set(gca,'FontSize',14);
title('f(\gamma)','fontsize',16);

% Plot of the original f and the center of [f]
plot(t,mf,'linewidth',2);
hold on, plot(t,mf_aligned,'--r','linewidth',2);
legend('f','center [f]');
set(gca,'FontSize',14);
title('Center [f]','fontsize',16);


%% Problem 9 
% compute the mean
mf = mf_aligned;
mq = mq_aligned;


fprintf('Computing Karcher mean of %d functions in SRVF space...\n',N);
ds = inf;
MaxItr = 30;
qun = zeros(1,MaxItr);
for r = 1:MaxItr
   % fprintf('updating step: r=%d\n', r);
%     if r == MaxItr
%         fprintf('maximal number of iterations is reached. \n');
%     end
    
    % Matching Step
    clear gam gam_dev;
    % use DP to find the optimal warping for each function w.r.t. the mean
    gam = zeros(N,size(q,1));
    gam_dev = zeros(N,size(q,1));
    parfor k = 1:N
        q_c = q(:,k,1)'; mq_c = mq(:,r)';
        [G,T] = DynamicProgrammingQ2(mq_c/norm(mq_c),t,q_c/norm(q_c),t,t,t);
        gam0 = interp1(T,G,t);
        gam(k,:) = (gam0-gam0(1))/(gam0(end)-gam0(1));  % slight change on scale
        gam_dev(k,:) = gradient(gam(k,:), 1/(M-1));
        f_temp(:,k) = interp1(t, f(:,k,1), (t(end)-t(1)).*gam(k,:) + t(1))';
		q_temp(:,k) = gradient(f_temp(:,k), binsize)./sqrt(abs(gradient(f_temp(:,k), binsize))+eps);
    end
    q(:,:,r+1) = q_temp;
    f(:,:,r+1) = f_temp;
    
    ds(r+1) = sum(simps(t, (mq(:,r)*ones(1,N)-q(:,:,r+1)).^2)) + ...
        lambda*sum(simps(t, (1-sqrt(gam_dev')).^2));
    
    % Minimization Step
    % compute the mean of the matched function
    mq(:,r+1) = mean(q(:,:,r+1),2);
    
    qun(r) = norm(mq(:,r+1)-mq(:,r))/norm(mq(:,r));
    if qun(r) < 1e-2 || r >= MaxItr
        break;
    end
end

if lambda == 0
    fprintf('additional run when lambda = 0\n');
    r = r+1;
    parfor k = 1:N
        q_c = q(:,k,1)'; mq_c = mq(:,r)';
        [G,T] = DynamicProgrammingQ2(mq_c/norm(mq_c),t,q_c/norm(q_c),t,t,t);
        gam0 = interp1(T,G,t);
        gam(k,:) = (gam0-gam0(1))/(gam0(end)-gam0(1));  % slight change on scale
        gam_dev(k,:) = gradient(gam(k,:), 1/(M-1));
    end
    gamI = SqrtMeanInverse(gam);
    gamI_dev = gradient(gamI, 1/(M-1));
    mq(:,r+1) = interp1(t, mq(:,r), (t(end)-t(1)).*gamI + t(1))'.*sqrt(gamI_dev');
	for k = 1:N
		q(:,k,r+1) = interp1(t, q(:,k,r), (t(end)-t(1)).*gamI + t(1))'.*sqrt(gamI_dev');
	    f(:,k,r+1) = interp1(t, f(:,k,r), (t(end)-t(1)).*gamI + t(1))';  
	    gam(k,:) = interp1(t, gam(k,:), (t(end)-t(1)).*gamI + t(1));
	end   
end

%% Aligned data & stats
fn = f(:,:,r+1);
qn = q(:,:,r+1);
q0 = q(:,:,1);
mean_f0 = mean(f0, 2);
std_f0 = std(f0, 0, 2);
mean_fn = mean(fn, 2);
std_fn = std(fn, 0, 2);
mqn = mq(:,r+1);
fmean = mean(f0(1,:))+cumtrapz(t,mqn.*abs(mqn));

gam = gam.';
[~,fy] = gradient(gam,binsize,binsize);
psi = sqrt(fy+eps);

option.showplot = 1;

if option.showplot == 1
    figure; clf;
    plot(t, f(:,:,1), 'linewidth', 1);
    title('Original data', 'fontsize', 16);
    pause(0.1);
    
    figure; clf;
    plot((0:M-1)/(M-1), gam, 'linewidth', 1);
    axis square;
    title('Warping functions', 'fontsize', 16);
    
    figure; clf;
    plot(t, fn, 'LineWidth',1);
    title(['Warped data, \lambda = ' num2str(lambda)], 'fontsize', 16);
    
    figure; clf;
    plot(t, mean_f0, 'b-', 'linewidth', 1); hold on;
    plot(t, mean_f0+std_f0, 'r-', 'linewidth', 1);
    plot(t, mean_f0-std_f0, 'g-', 'linewidth', 1);
    title('Original data: Mean \pm STD', 'fontsize', 16);
    
    figure; clf;
    plot(t, mean_fn, 'b-', 'linewidth', 1); hold on;
    plot(t, mean_fn+std_fn, 'r-', 'linewidth', 1);
    plot(t, mean_fn-std_fn, 'g-', 'linewidth', 1);
    title(['Warped data, \lambda = ' num2str(lambda) ': Mean \pm STD'], 'fontsize', 16);
    
    figure; clf;
    plot(t, fmean, 'g','LineWidth',1);
    title(['f_{mean}, \lambda = ' num2str(lambda)], 'fontsize', 16);
end
if option.parallel == 1 && option.closepool == 1
    matlabpool close
end
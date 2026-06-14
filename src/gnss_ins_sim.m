function gnss_ins_sim()
% GNSS_INS_SIM  Loosely-coupled GNSS/INS navigation with an 8-state EKF.
%
%   Personal project (Ali Murtaza). A planar (North-East) strapdown
%   inertial navigator is mechanised from noisy, biased IMU data and
%   fused with 1 Hz GPS position fixes by a loosely-coupled Extended
%   Kalman Filter. The filter state augments position, velocity and
%   heading with the gyro and accelerometer biases, so the IMU errors
%   are *observed* through GPS rather than left to accumulate.
%
%   A 25 s GPS dropout is injected to show inertial coasting (dead
%   reckoning) once the biases have been estimated, versus an open-loop
%   INS that shares the identical mechanisation but receives no GPS.
%
%   State x = [pN pE vN vE psi bg bax bay].   No toolboxes required.

clc; close all;
here = fileparts(mfilename('fullpath'));
outdir = fullfile(here, '..', 'assets');
if ~exist(outdir,'dir'); mkdir(outdir); end
rng(7);

%% ---------------------------------------------------- True trajectory
dt = 0.01; T = 120; t = (0:dt:T).'; N = numel(t);
spd = 12;                                   % ground speed [m/s]
psi = deg2rad(40*sin(2*pi*t/60));           % heading sweep [rad]
psidot = gradient(psi, dt);
vN_t = spd*cos(psi);  vE_t = spd*sin(psi);
pN_t = cumtrapz(t, vN_t);  pE_t = cumtrapz(t, vE_t);
aN_t = gradient(vN_t, dt); aE_t = gradient(vE_t, dt);
% True specific force resolved in body axes: a_body = R(-psi)*[aN;aE]
ax_b =  cos(psi).*aN_t + sin(psi).*aE_t;
ay_b = -sin(psi).*aN_t + cos(psi).*aE_t;

%% ----------------------------------------------------------- IMU model
bg_true = deg2rad(0.5);                      % gyro bias [rad/s]
ba_true = [0.05; 0.03];                      % accel bias [m/s^2]
ng = deg2rad(0.10); na = 0.02;              % gyro/accel white noise (1-sigma)
gyro = psidot + bg_true + ng*randn(N,1);
acc  = [ax_b, ay_b] + ba_true.' + na*randn(N,2);

%% ----------------------------------------------------------- GPS model
gps_dt = 1.0; gps_sigma = 2.5;              % 1 Hz, 2.5 m (1-sigma)
gpsMeas = nan(N,2);
gps_idx = find(abs(mod(t,gps_dt)) < dt/2);
drop_lo = 60; drop_hi = 85;                 % GPS dropout window [s]
for j = 1:numel(gps_idx)
    k = gps_idx(j);
    if t(k) >= drop_lo && t(k) <= drop_hi; continue; end
    gpsMeas(k,:) = [pN_t(k), pE_t(k)] + gps_sigma*randn(1,2);
end

%% ------------------------------------------------- Open-loop INS (ref)
psi_ins = psi(1); vins = [vN_t(1); vE_t(1)]; pins = [0;0];
PINS = zeros(N,2);

%% --------------------------------------------------------- EKF set-up
x = [0; 0; vN_t(1); vE_t(1); psi(1); 0; 0; 0];          % biases start at 0
P = diag([gps_sigma, gps_sigma, 0.5, 0.5, deg2rad(5), ...
          deg2rad(1), 0.1, 0.1].^2);
Qd = diag([1e-7, 1e-7, (na*dt)^2, (na*dt)^2, (ng*dt)^2, ...
           (deg2rad(0.003))^2*dt, (0.003)^2*dt, (0.003)^2*dt]);
H  = [1 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0];
Rk = gps_sigma^2*eye(2);
XEKF = zeros(N,8);

for k = 1:N
    % --- Open-loop INS mechanisation (no GPS), for comparison
    psi_ins = psi_ins + gyro(k)*dt;
    Rni = [cos(psi_ins) -sin(psi_ins); sin(psi_ins) cos(psi_ins)];
    a_ins = Rni*acc(k,:).';
    pins = pins + vins*dt + 0.5*a_ins*dt^2;
    vins = vins + a_ins*dt;
    PINS(k,:) = pins.';

    % --- EKF predict (nonlinear strapdown with estimated biases)
    cp = cos(x(5)); sp = sin(x(5));
    axc = acc(k,1) - x(7);  ayc = acc(k,2) - x(8);   % bias-corrected accel
    aN = cp*axc - sp*ayc;   aE = sp*axc + cp*ayc;
    xdot = [x(3); x(4); aN; aE; gyro(k)-x(6); 0; 0; 0];
    x = x + xdot*dt;

    % Jacobian of continuous dynamics
    A = zeros(8);
    A(1,3)=1; A(2,4)=1;
    A(3,5)= -sp*axc - cp*ayc;  A(3,7)= -cp;  A(3,8)=  sp;
    A(4,5)=  cp*axc - sp*ayc;  A(4,7)= -sp;  A(4,8)= -cp;
    A(5,6)= -1;
    Fk = eye(8) + A*dt;
    P = Fk*P*Fk.' + Qd;

    % --- EKF update on GPS availability
    if ~isnan(gpsMeas(k,1))
        y = gpsMeas(k,:).' - H*x;
        S = H*P*H.' + Rk;
        K = P*H.'/S;
        x = x + K*y;
        P = (eye(8) - K*H)*P;
        P = 0.5*(P+P.');               % keep symmetric
    end
    XEKF(k,:) = x.';
end

%% --------------------------------------------------------- Error stats
errINS = hypot(PINS(:,1)-pN_t, PINS(:,2)-pE_t);
errEKF = hypot(XEKF(:,1)-pN_t, XEKF(:,2)-pE_t);
drop = t>=drop_lo & t<=drop_hi;
settled = t>20 & ~drop;                         % after bias convergence
fprintf('\n=== GNSS/INS Loosely-Coupled Navigation (8-state EKF) ===\n');
fprintf('  EKF horizontal RMS (GPS available, settled): %.2f m\n', rms(errEKF(settled)));
fprintf('  EKF max error during %.0f s dropout:         %.2f m\n', drop_hi-drop_lo, max(errEKF(drop)));
fprintf('  EKF error after re-acquisition:              %.2f m\n', errEKF(end));
fprintf('  Estimated gyro bias  %.3f deg/s (true %.3f)\n', rad2deg(XEKF(end,6)), rad2deg(bg_true));
fprintf('  Estimated accel bias [%.3f %.3f] m/s^2 (true [%.3f %.3f])\n', ...
        XEKF(end,7), XEKF(end,8), ba_true(1), ba_true(2));
fprintf('  Open-loop INS final drift:                   %.1f m\n', errINS(end));

%% --------------------------------------------------------------- Plots
cT=[.2 .2 .2]; cI=[.85 .33 .10]; cE=[.10 .45 .85]; cG=[.10 .65 .30];

f1=figure('Color','w','Position',[80 80 760 700]); hold on; grid on; box on; axis equal;
plot(pE_t,pN_t,'-','Color',cT,'LineWidth',2.2);
plot(PINS(:,2),PINS(:,1),'-','Color',cI,'LineWidth',1.6);
plot(XEKF(:,2),XEKF(:,1),'--','Color',cE,'LineWidth',1.8);
gp=~isnan(gpsMeas(:,1));
scatter(gpsMeas(gp,2),gpsMeas(gp,1),10,cG,'filled','MarkerFaceAlpha',.5);
xlabel('East [m]'); ylabel('North [m]');
title('Trajectory: truth vs inertial-only vs GNSS/INS EKF','FontWeight','bold');
legend('true path','INS only (drift)','EKF estimate','GPS fixes','Location','best');
exportgraphics(f1, fullfile(outdir,'01_trajectory.png'),'Resolution',200);

f2=figure('Color','w','Position',[80 80 920 420]); hold on; grid on; box on;
area([drop_lo drop_hi],[1 1]*max(errINS)*1.05,'FaceColor',[1 .8 .4], ...
     'FaceAlpha',.3,'EdgeColor','none');
plot(t,errINS,'-','Color',cI,'LineWidth',1.6);
plot(t,errEKF,'-','Color',cE,'LineWidth',1.8);
xlabel('time [s]'); ylabel('horizontal position error [m]');
title('Position error - inertial-only vs EKF (GPS dropout shaded)','FontWeight','bold');
legend('GPS dropout','INS only','EKF','Location','northwest'); xlim([0 T]);
exportgraphics(f2, fullfile(outdir,'02_position_error.png'),'Resolution',200);

f3=figure('Color','w','Position',[80 80 920 700]);
subplot(2,1,1); hold on; grid on; box on;
area([drop_lo drop_hi],[1 1]*max(errEKF(drop))*1.1,'BaseValue',0,...
     'FaceColor',[1 .8 .4],'FaceAlpha',.3,'EdgeColor','none');
plot(t,errEKF,'-','Color',cE,'LineWidth',1.8);
ylabel('EKF error [m]'); xlim([40 110]);
title('EKF coasting through GPS dropout and re-acquisition','FontWeight','bold');
legend('GPS dropout','EKF error','Location','northwest');
subplot(2,1,2); hold on; grid on; box on;
plot(t,rad2deg(XEKF(:,6)),'LineWidth',1.6);
yline(rad2deg(bg_true),'--k');
plot(t,XEKF(:,7),'LineWidth',1.6); plot(t,XEKF(:,8),'LineWidth',1.6);
yline(ba_true(1),':k'); yline(ba_true(2),':k');
xlabel('time [s]'); ylabel('estimated IMU biases');
title('IMU bias estimates converging to truth','FontWeight','bold');
legend('gyro bias [deg/s]','gyro truth','accel bias_x [m/s^2]','accel bias_y [m/s^2]','Location','east');
exportgraphics(f3, fullfile(outdir,'03_dropout_and_bias.png'),'Resolution',200);

save(fullfile(outdir,'gnss_ins_results.mat'),'t','pN_t','pE_t','XEKF','PINS','errEKF','errINS');
fprintf('\nFigures and results written to %s\n', outdir);
end

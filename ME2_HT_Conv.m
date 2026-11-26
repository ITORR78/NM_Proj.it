clc; clear;

% Hohmann Transfer - Modified Euler (ME2) - 11/19
% Explicit Trapezoidal Method & Convergence

% Intial orbital parameters

% Initial/final orbital radii [km]
r1 = 5000;
r2 = 20000;

a = (r1+r2)/2; % Semi-major axis

% Standard gravitational parameter of Earth [km^3/s^2]
mu = 398600.4418;

% Initial circular orbital velocity
vc1 = sqrt(mu/r1);

% Transfer orbit velocity at periapsis
vt1 = sqrt(mu*(2/r1-1/a));

% First velocity change
dv1 = vt1-vc1;

% Final circular orbital velocity
vc2 = sqrt(mu/r2);

% Analytical transfer velocity at apoapsis
vt2 = sqrt(mu*(2/r2-1/a));

% Analytical second velocity change
dv2_exact = vc2-vt2;

% Transfer time (half elliptical orbit)
T = pi*sqrt(a^3/mu);


% main function

f = @(u) [u(3);
    u(4);
    -mu*u(1)/(u(1)^2+u(2)^2)^(3/2);
    -mu*u(2)/(u(1)^2+u(2)^2)^(3/2)];


% Original Modified Euler

% Number of steps
n = 2000;

% Step size
h = T/n;

% Initial state after first velocity burn
u_initial = [r1; 0; 0; vt1];

u = u_initial;

% Store results
U = zeros(4,n+1);
U(:,1) = u;

% Modified Euler Method (Explicit Trapezoidal)
for i = 1:n

    % First slope
    k1 = f(u);

    % Second slope
    k2 = f(u+h*k1);

    % Next approximation
    u = u+(h/2)*(k1+k2);

    % Store results
    U(:,i+1) = u;

end

% Final position and velocity
xf = u(1);
yf = u(2);

vxf = u(3);
vyf = u(4);

% Final orbital radius
rf = sqrt(xf^2+yf^2);

% Desired final circular velocity
vDesired = [0; -vc2];

% Second velocity change
dv2 = norm(vDesired-u(3:4));


% Original results

fprintf('\nORIGINAL MODIFIED EULER SIMULATION\n');
fprintf('----------------------------------\n');

fprintf('Transfer time: %.3f seconds\n',T);
fprintf('First velocity change: %.4f km/s\n',dv1);
fprintf('Final radius: %.4f km\n',rf);
fprintf('Second velocity change: %.4f km/s\n',dv2);


% Convergence study

% Number of steps to test
steps = 50*2.^(0:8);

numCases = length(steps);

% Preallocate arrays
h_arr = zeros(numCases,1);
rf_arr = zeros(numCases,1);
dv2_arr = zeros(numCases,1);

radErr = zeros(numCases,1);
velErr = zeros(numCases,1);

% Reset desired final velocity
vDesired = [0; -vc2];

% Loop through different step sizes
for k = 1:numCases

    n = steps(k);

    % Step size
    h = T/n;

    % Reset initial state
    u = u_initial;

    % Modified Euler Method (Explicit Trapezoidal)
    for i = 1:n

        % First slope
        k1 = f(u);

        % Second slope
        k2 = f(u+h*k1);

        % Next approximation
        u = u+(h/2)*(k1+k2);

    end

    % Final radius
    rf = sqrt(u(1)^2+u(2)^2);

    % Second velocity change
    dv2 = norm(vDesired-u(3:4));

    % Store results
    h_arr(k) = h;
    rf_arr(k) = rf;
    dv2_arr(k) = dv2;

    % Calculate absolute errors
    radErr(k) = abs(rf-r2);
    velErr(k) = abs(dv2-dv2_exact);

end


% Convergence order

% Preallocate convergence order arrays
radOrder = NaN(numCases,1);
velOrder = NaN(numCases,1);

for k = 2:numCases

    % Radius error convergence order
    radOrder(k) = log(radErr(k-1)/radErr(k)) / ...
                  log(h_arr(k-1)/h_arr(k));

    % Velocity error convergence order
    velOrder(k) = log(velErr(k-1)/velErr(k)) / ...
                  log(h_arr(k-1)/h_arr(k));

end


% Table of values

ResultsTable = table(steps',h_arr,rf_arr,dv2_arr,...
    radErr,velErr,radOrder,velOrder,...
    'VariableNames',{'Steps_n','StepSize_s',...
    'FinalRadius_km','dv2_km/s',...
    'RadiusError_km','dv2Error_km/s',...
    'RadiusOrder','VelocityOrder'});

fprintf('\nMODIFIED EULER CONVERGENCE STUDY\n');
fprintf('----------------------------------\n');

disp(ResultsTable);


% Convergence graph

figure;

loglog(steps,radErr,'-o',...
    'LineWidth',1.5,'MarkerSize',6);

hold on;

loglog(steps,velErr,'-s',...
    'LineWidth',1.5,'MarkerSize',6);

% Second-order reference line
reference = radErr(1)*(steps/steps(1)).^(-2);

loglog(steps,reference,'--k','LineWidth',1.5);

grid on;

xlabel('Number of Steps (n)');
ylabel('Absolute Error');

title('Modified Euler (ME2) Convergence Study');

legend('Final Radius Error [km]',...
    '\Delta v_2 Error [km/s]',...
    'Second-Order Reference',...
    'Location','southwest');
# MAE157A
Code containing the numerical calculations of the flight parameters of our capstone rocket by solving the Tsiolkovsky Rocket Eqn. 

%% TRAJECTORY MASTER CODE

%% AEROTECH F-50 Thrust Profile
clc; clear; close all;
%M_TO = m_prop + m_dry + m_rocket; %mass at TO

Thrust = [ 
0.012 51.377
0.023 61.197
0.026 66.117
0.044 66.564
0.082 69.685
0.152 73.264
0.208 75.053
0.237 77.279
0.254 76.832
0.272 77.726
0.307 77.726
0.330 76.832
0.336 78.621
0.342 76.832
0.354 79.590
0.363 76.385
0.371 77.756
0.395 76.385
0.447 75.937
0.523 73.711
0.652 68.344
0.810 60.302
0.828 62.539
0.836 58.076
0.901 53.603
1.079 37.074
1.158 29.480
1.196 25.464
1.246 16.976
1.301 9.380
1.430 0.000]; %Digitized thrust profile


%Add 4.3 seconds of zero thrust to the end of the thrust profile 
delay = 4.3;
t_free = Thrust(end,1)+0.01:0.01:Thrust(end,1) + delay;
test = zeros(length(t_free),1);
Thrust2 = [transpose(t_free) test];
Thrust = vertcat(Thrust, Thrust2);  %Now can be used for launch->ejection

%Isp
Isp = (76.83/1.43)/(0.0393/1.43 * 9.81); %Used for m_dot calculations later

%% Variable initializations

%Values from spec sheet
m_prop = 0.0393; %kg, will be updated post-script
m_dry = 0.0550; %kg, mass of motor sans propellant
m_frame = 0.302805; %kg, nosecone and frame and fins and internal structure
m_pl = 0.062369 + 0.01048932 + 0.0042; %kg, mass of the egg, altimeter, parachute
M_TO = 0.4333; %mass at takeoff

%Calculating mass profile across burnout
g = 9.81; %m/s^2
m_dot = Thrust(:,2)/(Isp*g); %m = T/Isp*g
its = length(Thrust); %Total number of iterations required
m_out = zeros(its,1); %m_out is calculating how much mass is exiting
m_bo = zeros(its,1); %mass profile across burnout
m_bo(1) = M_TO; %starting mass is MTOW
for i = 2:its
    dt = Thrust(i,1) - Thrust(i-1,1); %Calculate time step
    m_out(i) = m_dot(i)*dt; %Total ejected mass over time step
    m_bo(i) = m_bo(i-1) - m_out(i); %Mass profile across burn
end
M_BO = m_bo(end); %total mass end of burning phase

%Miscellaneous variables
theta_deg = 10;
theta_rad = degtorad(theta_deg); %Launch angle (rad)
CD = 0.80; %CD of aircraft, subject to change upon wind tunnel testing
A = 0.004561281; %m^2 assuming a fairing radius of 1.5in

%% Numerical calculation of Tsiolkovsky Rocket Eqn

%LAUNCH TO EJECTION
du_bo = zeros(its,1); %differential velocity
t_bo = Thrust(:,1); %time profile
T_bo = Thrust(:,2); %thrust profile
dt = zeros(its,1);
D_bo = zeros(its,1);
T = zeros(its,1);
grv = zeros(its,1);
V_bo = zeros(its,1);
rho_bo = zeros(its,1);
dist_bo = zeros(its,1);

%Initial conditions
V_bo(1) = 0; %Starting at rest
du_bo(1) = 0; %Starting at rest
rho_bo(1) = 1.225; %Sea-level altitude
dist_bo(1) = 0; %Starting on the ground

for i = 2:its
    dt(i) = t_bo(i) - t_bo(i-1); %Calculate your time step
    D_bo(i) = ((CD * 0.5 * rho_bo(i-1) .* V_bo(i-1).^2 * A) ./ m_bo(i)) .* dt(i); %Calculate drag using only i-1 values %%ORIGINALLY m_bo(i-1)
    T(i) = Thrust(i-1,2)./m_bo(i) .* dt(i); %Calculate Thrust
    grv(i) = g * cos(theta_rad) * dt(i); %Gravity term should be fairly constant with exception of time step
    du_bo(i) = T(i) - D_bo(i) - grv(i); %Velocity increment
%     if (du_bo(i) < 0) && (T(i) > 0) %Used to account for low thrust early on
%         du_bo(i) = T(i) - grv(i); 
%         if du_bo(i) < 0
%             du_bo(i) = 0;
%         end
%     end
    V_bo(i) = V_bo(i-1) + du_bo(i); %Vf = V0 + dv
    dist_bo(i) = dist_bo(i-1) + (V_bo(i)) .* dt(i); %X = x0 + v*t
    if dist_bo(i) < dist_bo(i-1) %After apogee gravity acts downward
        theta_rad = 0;
    end
    rho_bo(i) = 1.225*exp(-2.9e-5 * cos(theta_rad) * dist_bo(i).^ 1.15);
end

% %Altitude (y)
alt_bo = dist_bo * cos(theta_rad); %Altitude = Dist * cos(theta)

%EJECTION TO LANDING
t_free = Thrust(end,1) + 0.01:0.01:Thrust(end,1) + 150; %50 seconds to fall
it = length(t_free); %Total number of time steps for ejection->landing
CD_free = 1.5; %Parachute
A_free = 0.580644; %Assuming a square parachute 2.5 ft x 2.5 ft;
du_free = zeros(it,1);
V_free = zeros(it,1);
rho_free = zeros(it,1);
M_free = M_BO; %Constant falling mass, need to adjust for altimeter
dist_free = zeros(it,1);
D_free = zeros(it,1);

%Initial conditions
V_free(1) = V_bo(end); %Start with velocity at ejection
rho_free(1) = rho_bo(end); %Start with density at ejection
dt_free = 0.01; %Constant time steps
dist_free(1) = dist_bo(end); %Start with distanec at ejection
V_end = 0; %Used to find model's terminal velocity
TOF = it;

%Calculate ejection->landing
for j = 2:it
    D_free(j) = (CD_free * 0.5 * rho_free(j-1) .* V_free(j-1).^2 * A_free) / M_free;
    if V_free(j-1) > 0
        D_free(j) = -D_free(j);
    end
    du_free(j) = (-g + D_free(j))* dt_free;
    V_free(j) = V_free(j-1) + du_free(j);
    if (du_free(j) < 0.001) && (V_free(j) < V_end)
        V_end = V_free(j); %Find closest empirical V_terminal
    end
    dist_free(j) = dist_free(j-1) + (V_free(j)) .* dt_free; %X = X0 + V*t;
    if dist_free(j) * cos(theta_rad) <= 0 %If altitude < 0
        if j < TOF
            TOF = j; %index that it hits the ground
        end
        dist_free(j) = 0; %It hit the ground
        V_free(j) = 0; %It hit the ground
    end
    rho_free(j) = 1.225*exp(-2.9e-5 * cos(theta_rad) * dist_free(j).^ 1.15);
end

%True distance
alt_free = dist_free * cos(theta_rad);

%Check that payload converges to terminal velocity
V_term = -sqrt(g*M_free/(0.5*rho_free(end)*A_free*CD_free)); %Theoretical terminal velocity
ratio = (V_end / V_term)*100; %Comparing theoretical to model's
X = ['Your terminal velocity is ', num2str(ratio), '% the theoretical value of ', num2str(-V_term), ' m/s.'];
disp(X); %Sanity check lol

%% FINAL PLOTS
ms2mph = @(x) x*2.23694;
%Concatenation of launch->ejection and ejection->landing
time_final = vertcat(Thrust(:,1),transpose(t_free));
V_final = vertcat(V_bo, V_free);
V_final_mph = ms2mph(V_final);
dist_final = vertcat(dist_bo, dist_free);
alt_final = vertcat(alt_bo, alt_free);


%Velocity profile
% figure
subplot(2,2,3)
plot(time_final, V_final);
xlabel('Time (s)');
ylabel('Velocity (m/s)');
title('Velocity profile');
xlim([0 time_final(end)]);
hold on;
plot(Thrust(end,1), V_bo(end), 'rx', 'MarkerSize', 20);
hold off
X = ['Maximum velocity: ', num2str(max(V_final) * 2.23694), ' mph'];
disp(X);

%Altitude profile
% figure
subplot(2,2,[1 2])
plot(time_final, alt_final * 3.28084);
xlabel('Time (s)');
ylabel('Altitude (ft)');
title('Altitude profile');
xlim([0 time_final(end)]);
max_alt = max(alt_final);
X = ['Apogee: ', num2str(max(alt_final) * 3.28084), ' ft'];
hold on
plot(Thrust(end,1), alt_bo(end) * 3.28084, 'rx', 'MarkerSize', 20);
hold off;
disp(X);
% legend('Launch day profile', 'Numerically calculated profile');

%Acceleration
% figure
subplot(2,2,4)
plot(time_final(2:end), diff(V_final));
xlabel('Time (s)');
ylabel('Acceleration (m/s^2)');
title('Acceleration profile');
xlim([0 time_final(end)]);
ylim([-1 2.5])
%For altitude you can add an X at ejection and maybe dashed line to
%landing?

%% Get launch velocity

index_apogee = find(alt_final == max(alt_final)); %Make sure minimum is before apogee
alt_final_2 = alt_final - 4/3.2808; %Scale such that 4 ft = 0 position
index = find(alt_final_2(1:index_apogee,1) == min(abs(alt_final_2(1:index_apogee,1)))); 
Rail_V = V_final(index);
X = ['Rail velocity: ', num2str(Rail_V), ' m/s'];
disp(X);
%You get 1.2479 at alt_final = 3133 where 4/3.2808 = 1.2192

%% TOF

index_1 = length(Thrust(:,1));
index_TOF = index_1 + TOF;
Duration = time_final(index_TOF);
X = ['TOF: ', num2str(Duration), ' s'];
disp(X);





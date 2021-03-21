%%
% ELEC 4700 Assignment 3 
% Agberien Stephen 

%% Part 1
% This section involves modifying the Monte-Carlo simulator from assignment
% 1, and removing the bottle-neck, by adding a voltage across the x-axis of
% the semiconductor crystal. This voltage results in an electric field
% forming within the semiconductor. The applied voltage was set to 0.1 V,
% as given in the assignment instructions. 

nElectrons = 1000; 
nTime = 1000;
time = zeros(1, nTime);
m0 = 9.10938215e-31; %rest mass of an electron (kg)
mn = 0.26*m0; %effective mass of an electron (kg)
tau = 0.2e-12; %mean time between collisions (s)
Kb = 1.38064852e-23; %boltzmann constant (J/K)
T = 300; %temperature (K)
vth = sqrt(2*Kb*T/mn); %thermal velocity
dt = (100e-9)/vth/100;
l = vth*tau; %mean free path
temperature = zeros(1, nTime);
Pscat = 1 - exp(-dt/tau); 
x = linspace(0, 200, 400)*10^(-9); %x axis
y = linspace(0, 100, 400)*10^(-9); %y axis

voltagex = 0.1; 
voltagey = 0;
Ex = voltagex/(200e-9);
Ey = voltagey/(100e-9);
%%
% The electric field seen by the electrons is equal to:
E = sqrt(Ex^2 + Ey^2)

%%
% The force on each electron is equal to:
F = E * (1.60217653e-19)

%%
% The acceleration on each electron is equal to:
a = F/mn

%%
voltagex = 0.4; 
voltagey = 0;
Ex = voltagex/(200e-9);
Ey = voltagey/(100e-9);

Px = zeros(nElectrons, nTime);
Py = zeros(nElectrons, nTime); 

for n = 1 : nElectrons
    Px(n, 1) = x(randi(400));
    Py(n, 1) = y(randi(400));
end

Vx = zeros(nElectrons, nTime); 
Vy = zeros(nElectrons, nTime);
accelx = zeros(nElectrons, nTime); 
accely = zeros(nElectrons, nTime);

MaxwellBoltzmannVdist = makedist('Normal', 'mu', 0, 'sigma', sqrt(Kb*T/mn));

for k = 1 : nElectrons
    Vx(k, :) = random(MaxwellBoltzmannVdist);
    Vy(k, :) = random(MaxwellBoltzmannVdist);
    accelx(k, :) = Ex * (-1.60217653e-19/mn); %F = ma
    accely(k, :) = Ey * (-1.60217653e-19/mn);
end

avgV = sqrt(sum(Vx(:, 1).^2)/nElectrons + sum(Vy(:, 1).^2/nElectrons));

for j = 1 : nElectrons
    for w = 2 : nTime
        Vx(j, w) = Vx(j, w-1) + accelx(j, w-1) * dt;
        Vy(j, w) = Vy(j, w-1) + accely(j, w-1) * dt;
        if isnan(Px(j, w-1))
            if left == 1
                Px(j, w) = 0 + Vx(j, w)*dt;
            end
            if right == 1
                Px(j, w) = 200e-9 + Vx(j, w)*dt;
            end
        else
            Px(j, w) = Px(j, w-1) + Vx(j, w)*dt;
        end
        
        if Px(j, w) > 200e-9
            left = 1; 
            right = 0;
            Px(j, w) = NaN;
        end
        if Px(j, w) < 0
            left = 0;
            right = 1;
            Px(j, w) = NaN;
        end
        
        Py(j, w) = Py(j, w-1) + Vy(j, w)*dt; 
        if Py(j, w) > 100e-9
            Py(j, w) = 100e-9; 
            Vy(j, w:end) = -Vy(j, w);
        end
        if Py(j, w) < 0
            Py(j, w) = 0;
            Vy(j, w:end) = -Vy(j, w);
        end
        if Pscat > rand()
            Vx(j, w:end) = random(MaxwellBoltzmannVdist);
            Vy(j, w:end) = random(MaxwellBoltzmannVdist);
        end
    end
end

for i = 1:nTime
        temperature(i) = (sum(Vx(:, i).^2) + sum(Vy(:, i).^2))*mn/Kb/2/nElectrons;
        if i > 1
            time(i) = time(i-1) + dt;
        end
end

n = 100e15;
for i = 1:nTime
        Jd(i) = sqrt(sum(Vx(:, i).^2)/nElectrons + sum(Vy(:, i).^2/nElectrons))*(1.60217653e-19);
        if i > 1
            time(i) = time(i-1) + dt;
        end
end

for g = 1:10
    figure(1)
    plot(Px(g, :), Py(g, :))
    xlabel('x (nm)')
    ylabel('y (nm)')
    title('10 Electron Trajectories in N-type Si Semiconductor Crystal')
    hold on
end
%
% The electron drift current density is calculated by multiplying the
% charge of an electron, the number of electrons, the electron mobility in
% the material, and the electric field. The carrier velocity can be
% calculated by multiplying the electron mobility and the electric field.
% Therefore, the electron drift current density can be simplified as the
% product of the carrier charge, number of carriers, and the average
% carrier velocity. A plot of the drift current over time is shown in the
% following figure and a noticeable property is that the current
% drastically increases at first and then levels out. This effect can be
% due to the initialization of the velocities. Considering the electric 
% field begins to accelerate the electrons as time goes on, it is
% reasonable that the current follows an increasing trend. 

figure(2)
plot(time, Jd)
title('Drift Current vs. Time')
xlabel('time')
ylabel('current (A)')

densityM = [Px(:, 1000), Py(:, 1000)];
figure(3)
hist3(densityM, [200 100])
title('Electron Density Map (Final Positions)')
xlabel('x')
ylabel('y')
zlabel('# of electrons')

tempx = zeros(ceil(200), ceil(100));
tempy = zeros(ceil(200), ceil(100));
tempn = zeros(ceil(200), ceil(100)); 

for z = 1:nElectrons
    x = floor(Px(z, 1000)/1e-9);
    y = floor(Py(z, 1000)/1e-9);
    if (x == 0 || isnan(x)) 
        x = 1; 
    end 
    if (y == 0 || isnan(y))
        y = 1; 
    end
    tempy(x, y) = tempy(x, y) + Vy(z, 1000)^2;
    tempx(x, y) = tempx(x, y) + Vx(z, 1000)^2;
    tempn(x, y) = tempn(x, y) + 1;
end

temp2 = (tempx + tempy) .* mn ./ Kb ./ 2 ./ tempn;
temp2(isnan(temp2)) = 0; 
temp2 = temp2';

figure(4)
xtemp = linspace(1, 200, 200);
ytemp = linspace(1, 100, 100);
pcolor(xtemp, ytemp, temp2)
shading interp
colormap(jet)
title('Temperature Density Map (Final Velocities)')
xlabel('x')
ylabel('y')
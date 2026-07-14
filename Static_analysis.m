clear all; clc;
%% QZS Figure
close all;

% 设置下采样因子
downsample_factor = 10; % 例如，减少到原来的 1/10

% Data_0429_1 = csvread("Static_Testing/0429_QZS_1_50mm.csv",1);
Data_0429_1 = csvread("Static_Testing/Stat_QZS_1.csv",1);
deform_1    = Data_0429_1(:,1);
force_1     = Data_0429_1(:,2);
% 对变形数据进行下采样
% deform_1_downsampled = downsample(deform_1, downsample_factor);
% % 对力数据进行下采样
% force_1_downsampled = downsample(force_1, downsample_factor);

Data_0429_2 = csvread("Static_Testing/Stat_QZS_2.csv",1);
deform_2    = Data_0429_2(:,1);
% deform_2_downsampled = downsample(deform_2, downsample_factor);
force_2     = Data_0429_2(:,2);
% force_2_downsampled = downsample(force_2, downsample_factor);


Data_0429_3 = csvread("Static_Testing/Stat_QZS_3.csv",1);
deform_3    = Data_0429_3(:,1);
% deform_3_downsampled = downsample(deform_3, downsample_factor);
force_3     = Data_0429_3(:,2);
% force_3_downsampled = downsample(force_3, downsample_factor);


Data_0429_4 = csvread("Static_Testing/Stat_QZS_4.csv",1);
deform_4    = Data_0429_4(:,1);
% deform_4_downsampled = downsample(deform_4, downsample_factor);
force_4     = Data_0429_4(:,2);
% force_4_downsampled = downsample(force_4, downsample_factor);

%------- Spring Structure Maths Derivation, static 1 -----
a_1 = 50; %mm
b_1 = 80; %mm
x0_1 = 20;

k_1 = (250/1E3*9.8)/(80); % N/mm

x1_1 = [0:1:60]; % Displacement

EF_0_1 = sqrt(3)*sqrt(b_1^2 - (b_1-x0_1/2)^2 )*ones(size(x1_1));
EF_1_1 = sqrt(3)*sqrt(b_1^2 - ((2*b_1-x0_1-x1_1)/2).^2  );

delta_1 = EF_1_1-EF_0_1;

F_EF_1 = k_1.*delta_1;
% F_EB_1 = k_1*delta_1.*EF_1_1/(2*b_1);
F_EB_1 = F_EF_1.*cos(pi/6);

psi_1 = asin((2*b_1-x0_1-x1_1)/(2*b_1));
F_EB1_1 = F_EB_1./(2*cos(psi_1));

F_B1_1 = F_EB1_1.*sin(psi_1);

F_y_1 = 6*F_B1_1; % 3springs, 6 supports


%------- Spring Structure Maths Derivation, static 2 -----

a_2 = 50; %mm
b_2 = 80; %mm
x0_2 = 30;

k_2 = (400/1E3*9.8)/(80); % N/mm

x1_2 = [0:1:60]; % Displacement

EF_0_2 = sqrt(3)*sqrt(b_2^2 - (b_2-x0_2/2)^2 )*ones(size(x1_2));
EF_1_2 = sqrt(3)*sqrt(b_2^2 - ((2*b_2-x0_2-x1_2)/2).^2  );

delta_2 = EF_1_2-EF_0_2;

F_EF_2 = k_2.*delta_2;
% F_EB_2 = k_2*delta_2.*EF_1_2/(2*b_2);
F_EB_2 = F_EF_2.*cos(pi/6);

psi_2 = asin((2*b_2-x0_2-x1_2)/(2*b_2));
F_EB1_2 = F_EB_2./(2*cos(psi_2));

F_B1_2 = F_EB1_2.*sin(psi_2);

F_y_2 = 6*F_B1_2; % 3springs, 6 supports


%------- Spring Structure Maths Derivation, static 3 -----
a_3 = 50; %mm
b_3 = 80; %mm
x0_3 = 30;

k_3 = (800/1E3*9.8)/(80); % N/mm

x1_3 = [0:1:80]; % Displacement

EF_0_3 = sqrt(3)*sqrt(b_3^2 - (b_3-x0_3/2)^2 )*ones(size(x1_3));
EF_1_3 = sqrt(3)*sqrt(b_3^2 - ((2*b_3-x0_3-x1_3)/2).^2  );

delta_3 = EF_1_3-EF_0_3;

F_EF_3 = k_3.*delta_3;
% F_EB_3 = k_3*delta_3.*EF_1_3/(2*b_3);
F_EB_3 = F_EF_3.*cos(pi/6);

psi_3 = asin((2*b_3-x0_3-x1_3)/(2*b_3));
F_EB1_3 = F_EB_3./(2*cos(psi_3));

F_B1_3 = F_EB1_3.*sin(psi_3);

F_y_3 = 6*F_B1_3; % 3springs, 6 supports

%------ Plotting
color_full = ["#2e974e" ,  "#e25508" , "#2e7ebb" ,"#5f5f5f" , "#7262ac" , "#d92523"];
colors_2 = ["#b8e3b2","#fdc38d" ,"#b7d4ea", "#cecece", "#cfcfe5",   "#fcab8f", "#adadae", "#c2e6f7", "#e7c6db" , "#b3c7c8"]; 
rgb_array = cellfun(@(hex) sscanf(hex(2:end), '%2x%2x%2x', [1 3]), colors_2, 'UniformOutput', false);
rgb_array = cell2mat(rgb_array);
% disp(rgb_array)

figure()
hold on
plot(deform_4,force_4,'^','Color',colors_2(3)) % 04x05x50
plot(x1_1,F_y_1/1.3,'Color',color_full(3),'LineWidth',1.5)

plot(deform_1,force_1,'o','Color',colors_2(1)) % 05x04x35
plot(x1_2,F_y_2*2.25,'Color',color_full(1),'LineWidth',1.5)

plot(deform_3,force_3,'s','Color',colors_2(2)) % 05x05x35
plot(x1_3*6/7,F_y_3*2.3,'Color',color_full(2),'LineWidth',1.5)




% plot(deform_2_downsampled,force_2_downsampled,'.-','Color',colors_2(4))% 05x04x25



% Highlight
% x1 = [25,45,45,25];
% y1 = [1.5,1.5,2.5,2.5];
% plot1 = fill(x1,y1,'b');
% plot1.EdgeAlpha=0;
% plot1.FaceAlpha=0.3;
% 
% x2 = [32,45,45,32];
% y2 = [6,6,7,7];
% plot2 = fill(x2,y2,'g');
% plot2.EdgeAlpha=0;
% plot2.FaceAlpha=0.3;
% 
% x3 = [35,42,42,35];
% y3 = [12.7,12.7,13.7,13.7];
% plot3 = fill(x3,y3);
% plot3.EdgeAlpha=0;
% plot3.FaceAlpha=0.3;

% legend('05x04x35','05x05x35','04x05x50','Maths',fontsize=12);
legend('Spring 1 - Test','Spring 1 - Model','Spring 2 - Test','Spring 2 - Model','Spring 3 - Test','Spring 3 - Model')
axis([0 60 0 18])
grid on
box on;

% title('Force-Deflection Curve for Horizontal Spring QZS Structure',fontsize=12)
xlabel('Displacement [mm]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
ylabel('Loading [N]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')



%% 0816 - 3 Units
% Data_0816       = csvread("Static_Testing/0816_3Unit_1.csv",1);
% ----- Switch to 0816_3Unit_2: with constraints
Data_0816       = csvread("Static_Testing/0816_3Unit_2.csv",1);

deform_0816     = Data_0816(:,5);
force_0816      = 2*Data_0816(:,4);

len = size(force_0816,1);
stage_1 = round(49.4/max(deform_0816)*len);
stage_2 = round(79/max(deform_0816)*len);

figure()
plot(deform_0816,force_0816,'.-')
hold on
plot(deform_0816(1:stage_1),force_0816(1:stage_1),'.-','Color',color_full(3))
plot(deform_0816(stage_1:stage_2),force_0816(stage_1:stage_2),'.-','Color',color_full(1))
plot(deform_0816(stage_2:end),force_0816(stage_2:end),'.-','Color',color_full(2))

grid on
box on
% axis([0 200 0 12])
% title('Force-Deflection Curve for Three Units with Limits','FontSize',15)
xlabel('Displacement [mm]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
ylabel('Loading [N]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')











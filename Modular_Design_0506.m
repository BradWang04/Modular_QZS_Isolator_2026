clear all; clc; close all;
color_full = ["#2e974e" ,  "#e25508" , "#2e7ebb" ,"#5f5f5f" , "#7262ac" , "#d92523"];
colors_2 = ["#b8e3b2","#fdc38d" ,"#b7d4ea", "#cecece", "#cfcfe5",   "#fcab8f", "#adadae", "#c2e6f7", "#e7c6db" , "#b3c7c8"]; 

%% Experimental Data - Modular Design Schemes

%----- Case 1 -----
Data_Case1_Single_Config_1 = csvread("5-Modular_Design_Experiment/Case1-Single-Config_1.csv");
Data_Case1_Single_Config_2 = csvread("5-Modular_Design_Experiment/Case1-Single-Config_2.csv");
Data_Case1_Multiple = csvread("5-Modular_Design_Experiment/Case1_Multiple.csv");

% ----- Combined series -------
Data_050635_50_50_1     = csvread("5-Modular_Design_Experiment/050635_50-50_1.csv",1);
Data_050635_50_60_1     = csvread("5-Modular_Design_Experiment/050635_50-60_1.csv",1);
Data_050635_60_60_1     = csvread("5-Modular_Design_Experiment/050635_60-60_1.csv",1);
Data_050635_60_70_1     = csvread("5-Modular_Design_Experiment/050635_60-70_1.csv",1);
Data_050635_60_70_1     = csvread("5-Modular_Design_Experiment/050635_60-70_1.csv",1);
Data_050635_60t_70b_2   = csvread("5-Modular_Design_Experiment/050635_60t-70b_2.csv",1);
Data_050635_70t_60b_3   = csvread("5-Modular_Design_Experiment/050635_70t-60b_3.csv",1);




%% Experimental Data - ALL
% 设置下采样因子
downsample_factor = 10; % 例如，减少到原来的 1/10

% % Data_0429_1 = csvread("Static_Testing/0429_QZS_1_50mm.csv",1);
% Data_0429_1 = csvread("Static_Testing/Stat_QZS_1.csv",1);
% deform_1    = Data_0429_1(:,1);
% force_1     = Data_0429_1(:,2);
% % 对变形数据进行下采样
% % deform_1_downsampled = downsample(deform_1, downsample_factor);
% % % 对力数据进行下采样
% % force_1_downsampled = downsample(force_1, downsample_factor);
% 
% Data_0429_2 = csvread("Static_Testing/Stat_QZS_2.csv",1);
% deform_2    = Data_0429_2(:,1);
% % deform_2_downsampled = downsample(deform_2, downsample_factor);
% force_2     = Data_0429_2(:,2);
% % force_2_downsampled = downsample(force_2, downsample_factor);
% 
% 
% Data_0429_3 = csvread("Static_Testing/Stat_QZS_3.csv",1);
% deform_3    = Data_0429_3(:,1);
% % deform_3_downsampled = downsample(deform_3, downsample_factor);
% force_3     = Data_0429_3(:,2);
% % force_3_downsampled = downsample(force_3, downsample_factor);
% 
% 
% Data_0429_4 = csvread("Static_Testing/Stat_QZS_4.csv",1);
% deform_4    = Data_0429_4(:,1);
% % deform_4_downsampled = downsample(deform_4, downsample_factor);
% force_4     = Data_0429_4(:,2);
% % force_4_downsampled = downsample(force_4, downsample_factor);

Data_0607_3Unit_1   = csvread("Static_Testing/0607_3Unit_1.csv",1);
deform_5    = Data_0607_3Unit_1(:,1);
force_5     = Data_0607_3Unit_1(:,2);

% Data_0607_3Unit_2   = csvread("Static_Testing/0607_3Unit_2.csv",1);
% deform_5    = Data_0607_3Unit_2(:,1);
% force_5     = Data_0607_3Unit_2(:,2);
% 
% % Data_0816_3Unit_1   = csvread("Static_Testing/0816_3Unit_1.csv",1);
% % deform_6    = Data_0816_3Unit_1(:,1);
% % force_6     = Data_0816_3Unit_1(:,2);
% % 
% % Data_0816_3Unit_2   = csvread("Static_Testing/0816_3Unit_2.csv",1);
% % deform_7    = Data_0816_3Unit_2(:,1);
% % force_7     = Data_0816_3Unit_2(:,2);
% 
% Data_0605_QZS_singleUnit_Limits = csvread("Static_Testing/0605_QZS_singleUnit_Limits.csv",1);
% deform_8    = Data_0605_QZS_singleUnit_Limits(:,1);
% force_8     = Data_0605_QZS_singleUnit_Limits(:,2);
% 
% Data_0605_QZS_singleUnitWithLimits = csvread("Static_Testing/0605_QZS_singleUnitWithLimits.csv",1);
% deform_9    = Data_0605_QZS_singleUnitWithLimits(:,1);
% force_9    = Data_0605_QZS_singleUnitWithLimits(:,2);
% 
% Data_Stat_QZS_1 = csvread("Static_Testing/Stat_QZS_1.csv",1);
% deform_10   = Data_Stat_QZS_1(:,1);
% force_10   = Data_Stat_QZS_1(:,2);
% 
% Data_Stat_QZS_2 = csvread("Static_Testing/Stat_QZS_2.csv",1);
% deform_11   = Data_Stat_QZS_2(:,1);
% force_11    = Data_Stat_QZS_2(:,2);
% 
% Data_Stat_QZS_3 = csvread("Static_Testing/Stat_QZS_3.csv",1);
% deform_12   = Data_Stat_QZS_3(:,1);
% force_12    = Data_Stat_QZS_3(:,2);
% 
% Data_Stat_QZS_4 = csvread("Static_Testing/Stat_QZS_4.csv",1);
% deform_13   = Data_Stat_QZS_4(:,1);
% force_13    = Data_Stat_QZS_4(:,2);
% 

%-------- Data for [Combined Series]  ------------
deform_cs_1     = Data_050635_50_50_1(:,1);
force_cs_1      = Data_050635_50_50_1(:,2);

deform_cs_2     = Data_050635_50_60_1(:,1);
force_cs_2      = Data_050635_50_60_1(:,2);

deform_cs_3     = Data_050635_60_60_1(:, 1);
force_cs_3      = Data_050635_60_60_1(:, 2);

deform_cs_4     = Data_050635_60_70_1(:, 1);
force_cs_4      = Data_050635_60_70_1(:, 2);

deform_cs_5     = Data_050635_60t_70b_2(:, 1);
force_cs_5      = Data_050635_60t_70b_2(:, 2);

deform_cs_6     = Data_050635_70t_60b_3(:, 1);
force_cs_6      = Data_050635_70t_60b_3(:, 2);


%========= DATA visualization - [Combined Series] ========



%% ========= DATA visualization - Old ========
figure()
hold on
% plot(deform_5,force_5);
plot(deform_cs_1, force_cs_1, 'DisplayName', '050635-50-50.csv');
plot(deform_cs_2, force_cs_2, 'DisplayName', '050635-50-60.csv');
plot(deform_cs_3, force_cs_3, 'DisplayName', '050635-60-60.csv');
plot(deform_cs_4, force_cs_4, 'DisplayName', '050635-60-70.csv');
plot(deform_cs_5, force_cs_5, 'DisplayName', '050635-60t-70b.csv');
plot(deform_cs_6, force_cs_6, 'DisplayName', '050635-70t-60b.csv');

% 添加图例
legend('show', 'Location', 'best');

% 添加标题和坐标轴标签
title('Force vs Deformation');
xlabel('Deformation');
ylabel('Force');
grid on
% 保持图形窗口打开
hold off




%% Modular design: Case 1
% Same stiffness, different length 
% Spring selection:
% 04-50-5 & 05-70-6

a_1     = 50; %mm
b_1     = 80; %mm
l_0_1   = 100;
k_1     = (400/1E3*9.8)/(80); % N/mm
z_t_1 = [0:1:50]; % Displacement
l_1 = ( l_0_1^2 - 9/4*z_t_1.^2 + 3*z_t_1.*(3*b_1^2 - l_0_1^2)^0.5 ).^0.5; % Eq. 3-6
F_ver_1 = k_1*l_1.*(l_1 - l_0_1).*sqrt(3*b_1^2 - l_1.^2)/b_1^2;

a_2     = 50; %mm
b_2     = 80; %mm
l_0_2   = 20;
k_2     = (460/1E3*9.8)/(80); % N/mm
z_t_2 = [0:1:60]; % Displacement
l_2 = ( l_0_2^2 - 9/4*z_t_2.^2 + 3*z_t_2.*(3*b_2^2 - l_0_2^2)^0.5 ).^0.5; % Eq. 3-6
F_ver_2 = k_2*l_2.*(l_2 - l_0_2).*sqrt(3*b_2^2 - l_2.^2)/b_2^2;


rgb_array = cellfun(@(hex) sscanf(hex(2:end), '%2x%2x%2x', [1 3]), colors_2, 'UniformOutput', false);
rgb_array = cell2mat(rgb_array);

% Single
figure()
hold on
% deform_cs_1, force_cs_1
plot(Data_Case1_Single_Config_1(:,1),Data_Case1_Single_Config_1(:,2),'o','Color',colors_2(1),'LineWidth',1.5)

plot(z_t_1,F_ver_1,'Color',color_full(1),'LineWidth',1.5)

plot(Data_Case1_Single_Config_2(:,1)*0.85,Data_Case1_Single_Config_2(:,2)*1.05,'o','Color',colors_2(2),'LineWidth',1.5)
plot(z_t_2,F_ver_2,'Color',color_full(2),'LineWidth',1.5)

legend('Config 1, Exp','Config 1, Num','Config 2, Exp','Config 2, Num', 'FontName', 'Times New Roman', 'FontSize', 16)
box on;
grid on
axis([0 50 0 10])
xlabel('Displacement [mm]', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold')
ylabel('Loading [N]', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold')
%% Multiple
figure()
hold on
% deform_cs_1, force_cs_1

plot(Data_Case1_Multiple(:,1),Data_Case1_Multiple(:,2),'o','Color',colors_2(1),'LineWidth',1.5)
% plot(deform_cs_1,force_cs_1,'o','Color',colors_2(1),'LineWidth',1.5)
plot(z_t_1,F_ver_1,'Color',color_full(1),'LineWidth',1.5)


% plot(deform_2*1.2+30,force_2/4*1.05+1.5,'o','Color',colors_2(2),'LineWidth',1.5)
plot(z_t_2+max(z_t_1),F_ver_2+F_ver_1(end),'Color',color_full(2),'LineWidth',1.5)
legend('Case 1, Experimental','Case 1, Numerical', 'FontName', 'Times New Roman', 'FontSize', 16)

box on;
grid on
xlabel('Displacement [mm]', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold')
ylabel('Loading [N]', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold')






%%  Case 2
a_3     = 50; %mm
b_3     = 80; %mm
l_0_3   = 30;
k_3     = (640/1E3*9.8)/(80); % N/mm
z_t_3 = [0:1:60]; % Displacement
l_3 = ( l_0_3^2 - 9/4*z_t_3.^2 + 3*z_t_3.*(3*b_3^2 - l_0_3^2)^0.5 ).^0.5; % Eq. 3-6
F_ver_3 = k_3*l_3.*(l_3 - l_0_3).*sqrt(3*b_3^2 - l_3.^2)/b_3^2;

a_4     = 50; %mm
b_4     = 80; %mm
l_0_4   = 30;
k_4     = (980/1E3*9.8)/(80); % N/mm
z_t_4 = [0:1:60]; % Displacement
l_4 = ( l_0_4^2 - 9/4*z_t_4.^2 + 3*z_t_4.*(3*b_4^2 - l_0_4^2)^0.5 ).^0.5; % Eq. 3-6
F_ver_4 = k_4*l_4.*(l_4 - l_0_4).*sqrt(3*b_4^2 - l_4.^2)/b_4^2;

rgb_array = cellfun(@(hex) sscanf(hex(2:end), '%2x%2x%2x', [1 3]), colors_2, 'UniformOutput', false);
rgb_array = cell2mat(rgb_array);

% Temp - Data
def_1 = [0.171184023	0.684736091	1.283880171	2.225392297	2.39657632	2.995720399	3.423680456	3.851640514	4.279600571	4.793152639	5.306704708	5.649072753	6.419400856	7.104136947	7.617689016	8.559201141	9.415121255	10.01426534	10.95577746	11.64051355	12.15406562	13.35235378	14.29386591	15.14978602	15.92011412	16.5192582	17.46077033	18.05991441	19.08701854	20.37089872	20.88445078	21.39800285	22.33951498	22.93865906	23.9657632	25.33523538	26.36233951	27.81740371	28.58773181	29.78601997	30.98430813	31.75463623	32.61055635	33.29529244	33.98002853	34.92154066	35.69186876	36.20542083	37.31811698	38.34522111	39.20114123	40.4850214	41.51212553	42.62482168	43.90870185	44.67902996	45.53495007	46.73323823	47.93152639	49.38659058	50.49928673	51.7831669	52.89586305	53.92296719	54.7788873	55.89158345	57.17546362	58.28815977	58.88730385	59.74322397];
for_1 = [0.056179775	0.252808989	0.47752809	0.758426966	0.95505618	1.179775281	1.376404494	1.43258427	1.685393258	1.966292135	2.247191011	2.584269663	2.865168539	3.202247191	3.539325843	3.93258427	4.157303371	4.382022472	4.606741573	4.831460674	5	5.168539326	5.337078652	5.533707865	5.730337079	5.95505618	6.123595506	6.460674157	6.685393258	6.853932584	7.02247191	7.134831461	7.303370787	7.443820225	7.584269663	7.752808989	7.865168539	8.061797753	8.174157303	8.286516854	8.426966292	8.483146067	8.623595506	8.707865169	8.679775281	8.764044944	8.820224719	8.820224719	8.820224719	8.820224719	8.820224719	8.820224719	8.792134831	8.848314607	8.735955056	8.707865169	8.595505618	8.56741573	8.539325843	8.511235955	8.483146067	8.483146067	8.45505618	8.370786517	8.202247191	8.08988764	7.97752809	7.921348315	7.808988764	7.696629213];

def_2 = [0.085911982	0.429079953	0.601223886	1.11621582	1.716479795	1.974695695	2.489047689	2.918127641	3.261455597	3.776127561	4.119615503	4.548695455	4.891863426	5.235191382	5.750343301	6.093511272	6.694735158	7.124295066	7.640086925	8.326422867	8.670710734	9.186502593	9.787406509	10.13137441	10.9039423	11.24871012	11.84993401	12.1950218	12.45259776	12.10830989	10.47518232	11.0764062	11.59235805	8.843014652	7.468102977	12.7960857	13.05366166	13.31187756	14.16939752	14.51304545	15.02899729	15.71549322	16.40134921	16.74547709	17.77466103	18.63250097	19.40410895	20.00453291	20.60479689	21.80516485	22.5770928	23.34918074	23.9496047	24.46427667	25.06454064	26.00749263	26.35050062	27.20770061	28.23560468	28.92114069	29.77802072	30.46323676	31.1486128	32.00501287	32.86157292	33.88915701	34.74555708	35.94384524	36.6289013	38.34090152	37.74159745	38.8542936	39.79516578	41.07904596	42.44803818	44.07460637	45.01531857	45.95603077	47.58179903	48.00847921	49.29171944	50.57495967	51.51551189	53.22655219	53.56876025	55.19292866	54.42340048	56.30434493	57.84388124	58.78475342	59.04056955	59.72482568];
for_2 = [0.056074766	0.196261682	0.364485981	0.61682243	0.813084112	1.065420561	1.205607477	1.401869159	1.570093458	1.76635514	1.962616822	2.158878505	2.299065421	2.46728972	2.747663551	2.887850467	3.252336449	3.53271028	3.925233645	4.205607477	4.542056075	4.934579439	5.242990654	5.523364486	5.91588785	6.336448598	6.700934579	7.177570093	7.317757009	6.981308411	5.775700935	6.140186916	6.560747664	4.738317757	3.785046729	7.514018692	7.654205607	7.906542056	8.186915888	8.411214953	8.831775701	9.140186916	9.336448598	9.644859813	10.00934579	10.34579439	10.57009346	10.79439252	10.99065421	11.35514019	11.63551402	11.94392523	12.1682243	12.36448598	12.56074766	12.81308411	12.92523364	13.14953271	13.28971963	13.42990654	13.59813084	13.68224299	13.79439252	13.87850467	13.99065421	14.07476636	14.1588785	14.1588785	14.21495327	14.24299065	14.21495327	14.21495327	14.10280374	14.10280374	14.01869159	14.07476636	13.93457944	13.79439252	13.71028037	13.48598131	13.37383178	13.26168224	13.09345794	12.95327103	12.92523364	12.56074766	12.70093458	12.3364486	12.14018692	12.02803738	11.85981308	11.77570093];

figure()% single
hold on
plot(def_1,for_1*1.05,'o','Color',colors_2(3),'LineWidth',1.5)
plot(z_t_3,F_ver_3,'Color',color_full(3),'LineWidth',1.5)

plot(def_2,for_2*1.05,'o','Color',colors_2(4),'LineWidth',1.5)
plot(z_t_4,F_ver_4,'Color',color_full(4),'LineWidth',1.5)

legend('Config 3, Exp','Config 3, Num','Config 4, Exp','Config 4, Num', 'FontName', 'Times New Roman', 'FontSize', 16)
box on;
grid on
xlabel('Displacement [mm]', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold')
ylabel('Loading [N]', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold')


%%
figure()
hold on
plot(def_1,for_1*1.1,'o','Color',colors_2(3),'LineWidth',1.5)
plot(z_t_3,F_ver_3,'Color',color_full(3),'LineWidth',1.5)

plot(def_2+max(z_t_3),for_2*1.05+F_ver_3(end)*1.15,'o','Color',colors_2(4),'LineWidth',1.5)
plot(z_t_4+max(z_t_3),F_ver_4+F_ver_3(end),'Color',color_full(4),'LineWidth',1.5)
% legend('Config 3, Exp','Config 3, Num','Config 4, Exp','Config 4, Num', 'FontName', 'Calibri', 'FontSize', 16)
legend('Case 2, Experimental','Case 2, Numerical', 'FontName', 'Times New Roman', 'FontSize', 16)

box on;
grid on
xlabel('Displacement [mm]', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold')
ylabel('Loading [N]', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold')


%%  Case 3
a_5     = 50; %mm
b_5     = 80; %mm
l_0_5   = 60;
k_5     = (640/1E3*9.8)/(80); % N/mm
z_t_5 = [0:1:60]; % Displacement
l_5 = ( l_0_3^2 - 9/4*z_t_3.^2 + 3*z_t_3.*(3*b_3^2 - l_0_3^2)^0.5 ).^0.5; % Eq. 3-6
F_ver_5 = k_3*l_3.*(l_3 - l_0_3).*sqrt(3*b_3^2 - l_3.^2)/b_3^2;

a_6     = 50; %mm
b_6     = 80; %mm
l_0_6   = 60;
k_6     = (640/1E3*9.8)/(80); % N/mm
z_t_6 = [0:1:60]; % Displacement
l_6 = ( l_0_3^2 - 9/4*z_t_3.^2 + 3*z_t_3.*(3*b_3^2 - l_0_3^2)^0.5 ).^0.5; % Eq. 3-6
F_ver_6 = k_3*l_3.*(l_3 - l_0_3).*sqrt(3*b_3^2 - l_3.^2)/b_3^2;

Case_3 = csvread("5-Modular_Design_Experiment/Case_3_single.csv");
case3_deform = Case_3(:,1);
case3_force = Case_3(:,2);


figure()
hold on
% plot(def_1,for_1*1.1,'o','Color',colors_2(5),'LineWidth',1.5)
plot(case3_deform,case3_force,'o','Color',colors_2(5),'LineWidth',1.5)

plot(z_t_5,F_ver_5,'Color',color_full(5),'LineWidth',1.5)
legend('Config 5, Exp' ,'Config 5, Num' ,'FontName', 'Times New Roman', 'FontSize', 20)
box on;
grid on;
xlabel('Displacement [mm]', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold')
ylabel('Loading [N]', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold')

%% Case 3 - Multiple config

Case_3_double = csvread("5-Modular_Design_Experiment/Config3_double_2.csv");
case3_d_deform    = Case_3_double(:,1);
case3_d_force     = Case_3_double(:,2);


figure()
hold on
plot(case3_d_deform*0.9,case3_d_force,'o','Color',colors_2(6),'LineWidth',1.5)

plot(z_t_5*2,F_ver_5,'Color',color_full(6),'LineWidth',1.5)
% legend('Config 5 parallel, Exp' ,'Config 5 parallel, Num' ,'FontName', 'Calibri', 'FontSize', 20)
legend('Case 3, Experimental','Case 3, Numerical', 'FontName', 'Times New Roman', 'FontSize', 16)

box on;
grid on;
axis([0 100 0 10])
xlabel('Displacement [mm]', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold')
ylabel('Loading [N]', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold')






%%
rgb_array = cellfun(@(hex) sscanf(hex(2:end), '%2x%2x%2x', [1 3]), colors_2, 'UniformOutput', false);
rgb_array = cell2mat(rgb_array);


figure()
subplot(3,2,1)
hold on
plot(z_t_1,F_ver_1,'Color',color_full(1),'LineWidth',1.5)
plot(z_t_2,F_ver_2,'Color',color_full(2),'LineWidth',1.5)
% plot(z_t_2+max(z_t_1),F_ver_2+F_ver_1(end),'Color',color_full(2),'LineWidth',1.5)
legend('Config 1','Config 2')
box on;
xlabel('Displacement [mm]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
ylabel('Loading [N]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')

subplot(3,2,2)
hold on
plot(z_t_1,F_ver_1,'Color',color_full(1),'LineWidth',1.5)
% plot(z_t_2,F_ver_2,'Color',color_full(2),'LineWidth',1.5)
plot(z_t_2+max(z_t_1),F_ver_2+F_ver_1(end),'Color',color_full(2),'LineWidth',1.5)
legend('Config 1','Config 2')
box on;
xlabel('Displacement [mm]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
ylabel('Loading [N]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')


subplot(3,2,3)
hold on
plot(z_t_3,F_ver_3,'Color',color_full(3),'LineWidth',1.5)
plot(z_t_4,F_ver_4,'Color',color_full(4),'LineWidth',1.5)
% plot(z_t_2+max(z_t_1),F_ver_2+F_ver_1(end),'Color',color_full(2),'LineWidth',1.5)
legend('Config 3','Config 4')
box on;
xlabel('Displacement [mm]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
ylabel('Loading [N]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')

subplot(3,2,4)
hold on
plot(z_t_3,F_ver_3,'Color',color_full(3),'LineWidth',1.5)
% plot(z_t_4,F_ver_4,'Color',color_full(4),'LineWidth',1.5)
plot(z_t_4+max(z_t_3),F_ver_4+F_ver_3(end),'Color',color_full(2),'LineWidth',1.5)
legend('Config 3','Config 4')
box on;
xlabel('Displacement [mm]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
ylabel('Loading [N]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')


subplot(3,2,5)
hold on
plot(z_t_5,F_ver_5,'Color',color_full(3),'LineWidth',1.5)
plot(z_t_6,F_ver_6,'Color',color_full(4),'LineWidth',1.5)
% plot(z_t_6+max(z_t_5),F_ver_6+F_ver_5(end),'Color',color_full(2),'LineWidth',1.5)
legend('Config 5','Config 6')
box on;
xlabel('Displacement [mm]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
ylabel('Loading [N]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')

subplot(3,2,6)
hold on
plot(z_t_5,F_ver_5,'Color',color_full(3),'LineWidth',1.5)
% plot(z_t_6,F_ver_6,'Color',color_full(4),'LineWidth',1.5)
plot(z_t_6+max(z_t_5),F_ver_6+F_ver_5(end),'Color',color_full(2),'LineWidth',1.5)
legend('Config 5','Config 6')
box on;
xlabel('Displacement [mm]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
ylabel('Loading [N]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')















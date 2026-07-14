%% 多个静态QZS实验数据
clear all; clc;

color_full = ["#5f5f5f",  "#7262ac","#2e7ebb" ,"#2e974e" ,"#e25508","#d92523"   ];
colors_2 = ["#cecece",  "#cfcfe5","#b7d4ea","#b8e3b2","#fdc38d" ,"#fcab8f"]; 
%% Read measurement data for various configurations
% Geometry_50_60

filename_1 = 'Multiple_QZS/Geometry_50_60/b50_s04-5-20.xlsx';
data_1 = readtable(filename_1, 'Sheet', 'b50_s04-5-20');
load_b50_1 = data_1.LoadValue;
defl_b50_1 = data_1.PositionValue;

filename_2 = 'Multiple_QZS/Geometry_50_60/b50_s04-5-25.xlsx';
data_2 = readtable(filename_2, 'Sheet', 'b50_s04-5-25');
load_b50_2 = data_2.LoadValue;
defl_b50_2 = data_2.PositionValue;

filename_3 = 'Multiple_QZS/Geometry_50_60/b50_s05-4-25.xlsx';
data_3 = readtable(filename_3, 'Sheet',  'b50_s05-4-25');
load_b50_3 = data_3.LoadValue;
defl_b50_3 = data_3.PositionValue;

filename_4 = 'Multiple_QZS/Geometry_50_60/b50_s05-5-25.xlsx';
data_4 = readtable(filename_4, 'Sheet',  'b50_s05-5-25');
load_b50_4 = data_4.LoadValue;
defl_b50_4 = data_4.PositionValue;

filename_5 = 'Multiple_QZS/Geometry_50_60/b50_s05-6-35.xlsx';
data_5 = readtable(filename_5, 'Sheet',  'b50_s05-6-35');
load_b50_5 = data_5.LoadValue;
defl_b50_5 = data_5.PositionValue;

figure()
grid on
hold on
plot(defl_b50_1,load_b50_1, 'LineWidth', 1.2);
plot(defl_b50_2,load_b50_2, 'LineWidth', 1.2);
plot(defl_b50_3,load_b50_3, 'LineWidth', 1.2);
plot(defl_b50_4,load_b50_4, 'LineWidth', 1.2);
plot(defl_b50_5,load_b50_5, 'LineWidth', 1.2);

xlabel('Displacement [mm]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
ylabel('Loading [N]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
legend('0.4-5-20','0.4-5-25','0.5-4-25','0.5-5-25','0.5-6-35','FontName', 'Calibri', 'FontSize', 12, 'FontWeight', 'bold')
title('QZS Configuration: b=50mm', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
hold off
box on


%%
filename_b60_1 =    'Multiple_QZS/Geometry_50_60/b60_s04-5-20.xlsx';
data_b60_1 = readtable(filename_b60_1, 'Sheet', 'b60_s04-5-20');
load_b60_1 = data_b60_1.LoadValue;
defl_b60_1 = data_b60_1.PositionValue;

filename_b60_2 =    'Multiple_QZS/Geometry_50_60/b60_s04-5-25.xlsx';
data_b60_2 = readtable(filename_b60_2, 'Sheet', 'b60_s04-5-25');
load_b60_2 = data_b60_2.LoadValue;
defl_b60_2 = data_b60_2.PositionValue;

filename_b60_3 =    'Multiple_QZS/Geometry_50_60/b60_s05-4-25.xlsx';
data_b60_3 = readtable(filename_b60_3, 'Sheet', 'b60_s05-4-25');
load_b60_3 = data_b60_3.LoadValue;
defl_b60_3 = data_b60_3.PositionValue;

filename_b60_4 =    'Multiple_QZS/Geometry_50_60/b60_s05-5-25.xlsx';
data_b60_4 = readtable(filename_b60_4, 'Sheet', 'b60_s05-5-25');
load_b60_4 = data_b60_4.LoadValue;
defl_b60_4 = data_b60_4.PositionValue;

filename_b60_5 =    'Multiple_QZS/Geometry_50_60/b60_s05-6-35.xlsx';
data_b60_5 = readtable(filename_b60_5, 'Sheet', 'b60_s05-6-35');
load_b60_5 = data_b60_5.LoadValue;
defl_b60_5 = data_b60_5.PositionValue;

figure()
grid on
hold on
plot(defl_b60_1,load_b60_1, 'LineWidth', 1.2);
plot(defl_b60_2,load_b60_2, 'LineWidth', 1.2);
plot(defl_b60_3,load_b60_3, 'LineWidth', 1.2);
plot(defl_b60_4,load_b60_4, 'LineWidth', 1.2);
plot(defl_b60_5,load_b60_5, 'LineWidth', 1.2);

xlabel('Displacement [mm]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
ylabel('Loading [N]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
legend('0.4-5-20','0.4-5-25','0.5-4-25','0.5-5-25','0.5-6-35','FontName', 'Calibri', 'FontSize', 12, 'FontWeight', 'bold')
title('QZS Configuration: b=60mm', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
hold off
box on








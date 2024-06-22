% =========================================================================
% 作者：wy
% 日期：2023年10月2日
% 程序作用：仿真PAM的错误比特概率
% =========================================================================
clear;clc;close all;
%% 参数
M_SNR_dB = -5:8:35;                                % SNR

%% PAM 仿真
disp('开始仿真 PAM，请等待...');
PAM2 = SignalConstellation(2,'PAM');
PAM4 = SignalConstellation(4,'PAM');
PAM8 = SignalConstellation(8,'PAM');
PAM16 = SignalConstellation(16,'PAM');
PAM32 = SignalConstellation(32,'PAM');

disp('2PAM ...');
BEP_DoublyFlat_2PAM = BitErrorProbability(M_SNR_dB,PAM2.SymbolMapping/sqrt(2),PAM2.BitMapping);
disp('4PAM ...');
BEP_DoublyFlat_4PAM = BitErrorProbability(M_SNR_dB,PAM4.SymbolMapping/sqrt(2),PAM4.BitMapping);
disp('8PAM ...');
BEP_DoublyFlat_8PAM = BitErrorProbability(M_SNR_dB,PAM8.SymbolMapping/sqrt(2),PAM8.BitMapping);
disp('16PAM ...');
BEP_DoublyFlat_16PAM = BitErrorProbability(M_SNR_dB,PAM16.SymbolMapping/sqrt(2),PAM16.BitMapping);
%% 绘图
LineWidth = 1.4;
MarkerSize= 10;
figure();
semilogy(M_SNR_dB,BEP_DoublyFlat_2PAM, '-d','Color',0.85*[0,0,0],'LineWidth',LineWidth,'MarkerSize',MarkerSize);
hold on;grid on;
semilogy(M_SNR_dB,BEP_DoublyFlat_4PAM, '-*','Color',0.85*[1,0,0],'LineWidth',LineWidth,'MarkerSize',MarkerSize);
semilogy(M_SNR_dB,BEP_DoublyFlat_8PAM, '-o','Color',0.70*[0,1,0],'LineWidth',LineWidth,'MarkerSize',MarkerSize);
semilogy(M_SNR_dB,BEP_DoublyFlat_16PAM,'-+','Color',0.75*[0,0,1],'LineWidth',LineWidth,'MarkerSize',MarkerSize);
xlabel('SNR (dB)');
ylabel('BEP');
set(gca,'FontName','Times New Roman','FontSize',12,'LooseInset', [0,0,0,0]);



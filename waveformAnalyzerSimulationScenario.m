% Скрипт для запуска waveformAnalyxer
clc
clear
addpath quadriga_src/
addpath waveform/


waveformanalyzerObject1 = WaveformAnalyzer();
waveformanalyzerObject1.plotPowerSpectrumDensity();
waveformanalyzerObject1.plotPayloadConstellation();
waveformanalyzerObject1.calcWaveformParameters();
waveformanalyzerObject1.calcEvmPerformance();

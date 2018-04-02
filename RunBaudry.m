% Run results of Baudry for each experiment
pkg load statistics
graphics_toolkit gnuplot

% Plot results of 4.1
clear
cd Baudry;
filename = "4.1.mat";
load(strcat("Data/", filename));
MixCombi;
PlotResults;
PlotEntropy;

% Plot results of 4.2
clear
filename = "4.2.mat";
load(strcat("Data/", filename));
MixCombi;
PlotResults;
PlotEntropy;

clear
filename = "4.3.mat";
load(strcat("Data/", filename));
MixCombi;
PlotResults;
PlotEntropy;

clear
filename = "4.4.1.mat";
load(strcat("Data/", filename));
MixCombi;
PlotResults;
PlotEntropy;

clear
filename = "4.4.2.mat";
load(strcat("Data/", filename));
MixCombi;
PlotResults_3D;
PlotEntropy;

clear
filename = "GvHD+.mat";
load(strcat("Data/", filename));
MixCombi;
PlotResults_4D;
PlotEntropy;

clear
filename = "GvHD-.mat";
load(strcat("Data/", filename));
MixCombi;
PlotResults_4D;
PlotEntropy;

cd ..

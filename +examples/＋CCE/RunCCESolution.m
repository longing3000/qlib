%% Direct evolution
clear all;
clc;
cd /Users/ylp/Documents/AcademicLife/code/qlib
sol=model.phy.Solution.CCESolution.EnsembleCCESolution('EnsembleCCE_RoyCoord.xml');
% sol=model.phy.Solution.CCESolution.SingleSampleCCESolution('SingleSampleCCE_RoyCoord.xml');
% sol=model.phy.Solution.CCESolution.DECCESolution('DECCE_RoyCoord.xml');
sol.perform();
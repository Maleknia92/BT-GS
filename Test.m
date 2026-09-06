% This is a quick test to check the consistency of this BT-GS
% implementation with your MATLAB version. Simply run this script
% to check that BT-GS operates correctly on your system.
%
% It minimizes the Funny function, which is defined by
% f(x) = x^2, starting from x_0 = 100.
%
% Since the objective function is convex, we set Type = 1.
%
% If the run is successful, the Table of Results should appear in the Command Window.
%
% Check also the output structure in the Workspace.

clc

clear

x0=100;

Type=0;

output=BTGS(@Funny,x0,Type);












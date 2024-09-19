clear
clc
close all

filename = "FrontWheelData.m";
run(filename)
test = TestDisp(s);
hasPassed = test.solveAndCheckProblem();
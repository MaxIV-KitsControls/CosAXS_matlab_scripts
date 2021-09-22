function [dI,ampl_heat] = subtractHeat(q, dI, dI_heat, qRange)
% Oskar Berntsson, 2021
% Fits data representing the heating of the solvent in a specified range.
% This gives a scale factor (ampl_heat) which is then used to scale the
% pure solvent heating data and subtract the appropriate amount of solvent
% heating from the entire signal.
% Inputs
% q: nqx1
% dI: nqxnt, the data
% dI_heat; nqxn_heat, there can be more than one heat component
% qRange: 1x2, the q-range in which dI_heat is scaled to dI.
% Outputs
% dI: heat subtracted data
% ampl_heat: the heat scale factor
%%
ampl_heat = dI_heat( q>=qRange(1)&q<qRange(2),: ) \ dI( q>=qRange(1)&q<=qRange(2),: );
dI = dI - dI_heat*ampl_heat;
end


function mb = mass_balance(md)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

smb = [];
discharge = [];
icevolume = [];
time = [];

figure();

for i = 1:12:length(md.results.TransientSolution)
    if i ~= length(md.results.TransientSolution)
        year_smb = (md.results.TransientSolution(i:i+11).TotalSmb) ;
        year_smb = year_smb %./12; %.* (1e-9) .* md.constants.yts; % from kg s^-1 to Gt/yr not 1e-12???
        year_smb = sum(year_smb);

        year_discharge = md.results.TransientSolution(i:i+11).GroundinglineMassFlux;
        year_discharge = year_discharge%./12;
        year_discharge = sum(year_discharge);

        smb = [smb year_smb];
        discharge = [discharge year_discharge];

        time = [time md.results.TransientSolution(i).time];
        %if i ~= 1
            icevolume = [icevolume ((md.results.TransientSolution(i).IceVolume - md.results.TransientSolution(i+11).IceVolume) ./(1e9))];
        %end
    end

end

plot(time,smb,'DisplayName','smb');
hold on
plot(time,discharge,'DisplayName','discharge');
hold on
plot(time,smb-discharge,'DisplayName','smb-discharge');
hold on
plot(time,icevolume,'DisplayName','Icevolume change');
legend();
hold off
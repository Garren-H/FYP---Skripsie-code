% code which is used to test that the data in the 4D tensor corresponds to
% the data in the excel sheet.
Data_extraction_and_structuring; %calls the script which structure the data
fprintf('\nList of solvents:\n')
idx = 1:numel(Solvent.Name);
idx = idx';
for i = 1:idx(end)
    fprintf('%s - %s\n', num2str(i), Solvent.Name(i))
end
User_solvent = input('Please enter the number corresponding to the solvent: ', "s");
while sum(idx == str2num(User_solvent)) ~= 1
    
    fprintf('\nInvalid selection\n');
    User_solvent = input('Please enter the number corresponding to the solvent: ', "s");
    
end

User_solvent = Solvent.Name(str2num(User_solvent));

fprintf('\nList of solutes:\n')
idx = 1:numel(Solute.Name);
idx = idx';
for i = 1:idx(end)
    fprintf('%s - %s\n', num2str(i), Solute.Name(i))
end
User_solute = input('Please enter the number corresponding to the solute: ', "s");
while sum(idx == str2num(User_solute)) ~= 1
    
    fprintf('\nInvalid selection\n');
    User_solute = input('Please enter the number corresponding to the solute: ', "s");
    
end

User_solute = Solute.Name(str2num(User_solute));

fprintf('\nList of compositions:\n')
idx = 1:numel(mole_fraction_solvent);
idx = idx';
for i = 1:idx(end)
    fprintf('%s - %s\n', num2str(i), num2str(mole_fraction_solvent(i)))
    
end
User_mole_fract_solvent = input('Please enter the number corresponding to the composition: ', "s");
while sum(idx == str2num(User_mole_fract_solvent)) ~= 1
    
    fprintf('\nInvalid selection\n');
    User_mole_fract_solvent = input('Please enter the number corresponding to the composition: ', "s");
    
end

User_mole_fract_solvent = mole_fraction_solvent(str2num(User_mole_fract_solvent));

fprintf('\nList of temperature:\n')
unique_temps = unique(NewData.Temperature);
idx = 1:numel(unique_temps);
idx = idx';
for i = 1:idx(end)
    fprintf('%s - %s\n', num2str(i), num2str(unique_temps(i)))
    
end
User_temperature = input('Please enter the number corresponding to the temperature: ', "s");
while sum(idx == str2num(User_temperature)) ~= 1
    
    fprintf('\nInvalid selection\n');
    User_temperature = input('Please enter the number corresponding to the temperature: ', "s");
    
end

User_temperature = unique_temps(str2num(User_temperature));

idx1_sheet = ( (NewData.Name1 == User_solvent) & (NewData.Name2 == User_solute)  & ...
    (NewData.Temperature == User_temperature) & (NewData.Mole_fraction_Solvent == User_mole_fract_solvent) );

idx2_sheet = ( (NewData.Name1 == User_solute) & (NewData.Name2 == User_solvent)  & ...
    (NewData.Temperature == User_temperature) & (NewData.Mole_fraction_Solute == User_mole_fract_solvent) );



if sum(idx1_sheet) == 0 & sum(idx2_sheet) == 0
    
    fprintf('\nSheet Data:')
    fprintf('\nNo data for this spesific selection was found\n\n')
    
    fprintf('\nTensor data:')
    idx_solute = find(Solute.Name == User_solute);
    fprintf('\nSolute: %s', Solute.Name(idx_solute))
    idx_solvent = find(Solvent.Name == User_solvent);
    fprintf('\nSolvent: %s', Solvent.Name(idx_solvent))
    idx_comp = find(mole_fraction_solvent == User_mole_fract_solvent);
    fprintf('\nMole fraction of solvent: %f', mole_fraction_solvent(idx_comp))
    idx_temp = find(Temperature_values == User_temperature);
    fprintf('\nTemperature: %f', Temperature_values(idx_temp))
    fprintf('\nSpeed of sound: %f', tensor_4D(idx_solvent, idx_solute, idx_comp, idx_temp))
    
    fprintf('\nTensor data:')
    idx_solute1 = find(Solute.Name == User_solvent);
    fprintf('\nSolute: %s', Solute.Name(idx_solute1))
    idx_solvent1 = find(Solvent.Name == User_solute);
    fprintf('\nSolvent: %s', Solvent.Name(idx_solvent1))
    idx_comp1 = find(mole_fraction_solvent == User_mole_fract_solvent);
    fprintf('\nMole fraction of solvent: %f', mole_fraction_solvent(idx_comp1))
    idx_temp1 = find(Temperature_values == User_temperature);
    fprintf('\nTemperature: %f', Temperature_values(idx_temp1))
    fprintf('\nSpeed of sound: %f\n\n', tensor_4D(idx_solvent1, idx_solute1, idx_comp1, idx_temp1))
    
elseif sum(idx1_sheet) == 1
    
    fprintf('\nSheet Data:')
    fprintf('\nSolute: %s \nSolvent: %s \nMole fraction of solvent: %f \nTemperature: %f \nSpeed of Sound: %f\n\n', ...
        NewData.Name2(idx1_sheet), NewData.Name1(idx1_sheet), NewData.Mole_fraction_Solvent(idx1_sheet), ...
        NewData.Temperature(idx1_sheet), NewData.Speedofsound(idx1_sheet))
    fprintf('\nTensor data:')
    idx_solute = find(Solute.Name == User_solute);
    fprintf('\nSolute: %s', Solute.Name(idx_solute))
    idx_solvent = find(Solvent.Name == User_solvent);
    fprintf('\nSolvent: %s', Solvent.Name(idx_solvent))
    idx_comp = find(mole_fraction_solvent == User_mole_fract_solvent);
    fprintf('\nMole fraction of solvent: %f', mole_fraction_solvent(idx_comp))
    idx_temp = find(Temperature_values == User_temperature);
    fprintf('\nTemperature: %f', Temperature_values(idx_temp))
    fprintf('\nSpeed of sound: %f', tensor_4D(idx_solvent, idx_solute, idx_comp, idx_temp))
    
elseif sum(idx2_sheet) == 1
    
    fprintf('\nSheet Data:');
    fprintf('\nSolute = %s \nSolvent: %s \nMole fraction of solvent: %f \nTemperature: %f\nSpeed of Sound: %f\n\n', ...
        NewData.Name1(idx2_sheet), NewData.Name2(idx2_sheet), NewData.Mole_fraction_Solute(idx2_sheet), ...
        NewData.Temperature(idx2_sheet), NewData.Speedofsound(idx2_sheet))
    fprintf('\nTensor data:')
    idx_solute = find(Solute.Name == User_solvent);
    fprintf('\nSolute: %s', Solute.Name(idx_solute))
    idx_solvent = find(Solute.Name == User_solute);
    fprintf('\nSolvent: %s', Solvent.Name(idx_solvent))
    idx_comp = find(mole_fraction_solvent == User_mole_fract_solvent);
    fprintf('\nMole fraction of solvent: %f', mole_fraction_solvent(idx_comp))
    idx_temp = find(Temperature_values == User_temperature);
    fprintf('\nTemperature: %f', Temperature_values(idx_temp))
    fprintf('\nSpeed of sound: %f', tensor_4D(idx_solute, idx_solvent, idx_comp, idx_temp))
    
end

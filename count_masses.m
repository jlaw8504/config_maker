function mass_number = count_masses(filename)
% MASS_NUMBER counts the number of masses in a .cfg file or a .out file

%open the file for reading only
fid_in = fopen(filename);
%instantiate the mass_number variable
mass_number = 0;
%% loop through the file
% assign tline so that the lines can be looped through
tline = fgetl(fid_in);
while ischar(tline) == 1
    %remove leading/trailing whitespace from tline
    clean_line = strtrim(tline);
    if strncmpi(clean_line,'mass ',5) == 1
        mass_number = mass_number + 1;
    end
    % loop to the next line
    tline = fgetl(fid_in);
end
%close the file
fclose(fid_in);
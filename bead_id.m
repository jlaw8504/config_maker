function output = bead_id(infile)

% this code adds a number of cohesin molecules to the linear chain

% assign parameters
mass_mass = 3.38889e-020; % mass of a standard bead
mass_sep = 1e-008; % standard distances betweem masses
spring_rest = 1e-008; % spring distance at rest
spring_const = 0.226195; % standard spring constant
spring_damp = 10; % how much factor by which the gammma spring is weaker than the normal spring constant
time = 'test'; % this is used to find the last time step
colors = 'test'; % used to assign mass colors
condensin_mass_color = 2; % color of condensin beads
DNA_mass_color = [1 4]; % color of DNA beads
cohesin_mass_color = 5; % color of cohesin beads
mass_condensin = []; % this is used to look for condensin beads
mass_cohesin = []; % this is used to look for cohesin beads
springs_condensin = []; % this is used to look for springs bound to condensin
springs_cohesin = []; % this is used to look for springs bound to cohesin
access_range = 5; % size of binding box used to calculate accessiblity in terms of mass sep

% open up the file
fid_in = fopen(infile);

% assign tline so that the lines can be looped through
tline = fgetl(fid_in);

m = 1; % used to assign mass coords
n = 1; % used to assign springs

% loop through all the lines to find the masses and appropriate springs
while ischar(tline)
    % find the mass coordinates
    if size(strfind(tline,'mass '),1) ~= 0
        % split the string into pieces to parse coordinates
        b = strsplit(tline);
        % save the initial coordinates in a cell array
        mass_coords_in(m,1) = str2double(b{5});
        mass_coords_in(m,2) = str2double(b{6});
        mass_coords_in(m,3) = str2double(b{7});
        % record the color of each of the beads
        if size(b,2) == 7
            % rows without an eighth entry correpond to a color of 1 (red)
            mass_coords_in(m,4) = 1;
        elseif size(b,2) == 8
            % rows with an eighth entry have their color logged appropriately
            mass_coords_in(m,4) = str2double(b{8});
        else
        end
        mass_last(1,1) = str2double(b{3}); % assigns the last mass value in the file
        mass_coords_in(m,5) = str2double(b{3}); % puts the mass number in the table
        % increase the counter by 1
        m = m+1;
    else
    end
    if size(strfind(tline,'spring '),1) ~= 0
        % split the string into pieces to parse coordinates
        b = strsplit(tline);
        
        % save the desired springs into a cell array
        springs(n,1) = str2double(b{3});
        springs(n,2) = str2double(b{4});
        
        % increase the counter by 1
        n = n+1;
    else
    end
    if size(strfind(tline,'Time '),1) ~= 0
        % assign the time line to the variable time, this will always pick the last one
        if strcmp(time,'test')==1
            time_init = tline;
        else
        end
        
        time = tline;
    else
    end
    if size(strfind(tline,'MassColors'),1) ~= 0
        % assign the color line to the variable colors
        colors = tline;
    else
    end
    % loop to the next line
    tline = fgetl(fid_in);
end

% close the file
fclose('all');

% this next part logs the colors of all the beads in order
if size(strfind(colors,'MassColors'),1) ~= 0
    % open up the file (again)
    fid_in = fopen(infile);
    % assign tline so that the lines can be looped through (starts it over)
    tline = fgetl(fid_in);
    while ischar(tline)
        % looks for the beginning of the last time step
        if strcmp(colors,tline) == 1
            % loop to the next line
            tline = fgetl(fid_in);
            % this assigns the coordinates (which are listed backwards for some reason)
            for d = 0:mass_last(1,1)
                % split the line into components
                % log all of the coordinates into a cell matrix
                mass_colors(mass_last(1,1)+1-d,1) = str2double(tline);
                % loop to the next line
                tline = fgetl(fid_in);
            end
        else
            % loop to the next line if it isn't the right one
            tline = fgetl(fid_in);
        end
    end
else
end

% close the file
fclose('all');

% if this is not the original file, it looks for the last time step
if size(strfind(time,'Time '),1) ~= 0
    % open up the file (again)
    fid_in = fopen(infile);
    % assign tline so that the lines can be looped through (starts it over)
    tline = fgetl(fid_in);
    % loop through all the lines
    while ischar(tline)
        % looks for the beginning of the last time step
        if strcmp(time,tline) == 1
            % loop to the next line
            tline = fgetl(fid_in);
            % this assigns the coordinates (which are listed backwards for some reason)
            for d = 0:mass_last(1,1)
                % split the line into components
                e = strsplit(tline);
                % log all of the coordinates into a cell matrix
                mass_coords(mass_last(1,1)+1-d,1) = str2double(e{1});
                mass_coords(mass_last(1,1)+1-d,2) = str2double(e{2});
                mass_coords(mass_last(1,1)+1-d,3) = str2double(e{3});
                mass_coords(mass_last(1,1)+1-d,4) = mass_colors(mass_last(1,1)+1-d,1);
                mass_coords(mass_last(1,1)+1-d,5) = mass_last(1,1)-d;
                % loop to the next line
                tline = fgetl(fid_in);
            end
        else
            % loop to the next line if it isn't the right one
            tline = fgetl(fid_in);
        end
    end
else
    % if this is the orginial document, then the coords were parsed earlier
    mass_coords = mass_coords_in;
end

% close the file
fclose('all');



% assign counters for the following section
m = 1; % used to assign DNA
n = 1; % used to assign condensin
p = 1; % used to assign cohesin

% separate the condensin into a separate file
for z = 1:size(mass_coords)
    if max(mass_coords(z,4) == DNA_mass_color) == 1
        % assign everything with DNA colors to the DNA matrix
        mass_DNA(m,1:5) =  mass_coords(z,1:5);
        m = m+1;
    elseif max(mass_coords(z,4) == condensin_mass_color) == 1
        % assign everything with the condensin color to the condensin matrix
        mass_condensin(n,1:5) =  mass_coords(z,1:5);
        n = n+1;
    elseif max(mass_coords(z,4) == cohesin_mass_color) == 1
        % assign everything with the cohesin color to the condensin matrix
        mass_cohesin(p,1:5) =  mass_coords(z,1:5);
        p = p+1;
    else
    end
end

% time to find the accessiblity of all DNA beads
for z = 1:size(mass_DNA)
    accessibility(z,1) = calculate_access(mass_coords,mass_DNA(z,1:5),access_range*mass_sep);
end

% place the accessibility into the mass_DNA table
mass_DNA(:,6) = accessibility(:,1);

m = 1; % used to assign condensin springs

% create a list of springs that contain condensin
if size(mass_condensin,1)>0 % make sure we have condensin
    for z = 1:size(springs,1)
        if max(springs(z,2) == mass_condensin(:,5)) == 1
            % if the spring contains a condensin bead, log it in a table
            springs_condensin(m,:) = springs(z,:);
            % increase the counter by 1
            m = m + 1;
        else
        end
    end
else
end

m = 1; % used to assign cohesin springs

% create a list of springs that contain cohesin
if size(mass_cohesin,1)>0 % make sure we have condensin
    for z = 1:size(springs,1)
        if max(springs(z,2) == mass_cohesin(:,5)) == 1
            % if the spring contains a cohesin bead, log it in a table
            springs_cohesin(m,:) = springs(z,:);
            % increase the counter by 1
            m = m + 1;
        else
        end
    end
else
end

m = 1; % used to assign springs that do not involve cohesin or condensin

if size(mass_condensin,1) > 0 && size(mass_cohesin,1) > 0
    % if there is condensin and cohesin
    for z = 1:size(springs,1)
        if max(springs(z,1) == mass_condensin(:,5)) == 0 &&...
                max(springs(z,2) == mass_condensin(:,5)) == 0 &&...
                max(springs(z,1) == mass_cohesin(:,5)) == 0 &&...
                max(springs(z,2) == mass_cohesin(:,5)) == 0
            % if the spring is not bound to cohesin or condensin, log it in the table
            springs_DNA(m,:) = springs(z,:);
            % increase the counter by 1
            m = m + 1;
        else
        end
    end
elseif size(mass_condensin,1) > 0 && size(mass_cohesin,1) == 0
    % if there is condensin, but no cohesin
    for z = 1:size(springs,1)
        if max(springs(z,1) == mass_condensin(:,5)) == 0 &&...
                max(springs(z,2) == mass_condensin(:,5)) == 0
            % if the spring is not bound to condensin, log it in the table
            springs_DNA(m,:) = springs(z,:);
            % increase the counter by 1
            m = m + 1;
        else
        end
    end
elseif size(mass_condensin,1) == 0 && size(mass_cohesin,1) > 0
    % if there is cohesin, but no condensin
    for z = 1:size(springs,1)
        if max(springs(z,1) == mass_cohesin(:,5)) == 0 &&...
                max(springs(z,2) == mass_cohesin(:,5)) == 0
            % if the spring is not bound to cohesin, log it in the table
            springs_DNA(m,:) = springs(z,:);
            % increase the counter by 1
            m = m + 1;
        else
        end
    end
else % if there is no condensin or cohesin
    % all springs are DNA springs
    springs_DNA(:,:) = springs(:,:);
end

m = 1; % reset counter to assign unbound DNA

% make a table of all unbound DNA

if size(springs_condensin,1)>0
    % if there is a condensin spring, select the DNA beads not bound to condensin
    for z = 1:size(mass_DNA,1)
        % loop through all the DNA
        if max(mass_DNA(z,5) == springs_condensin(:,1)) == 0
            % if this mass is not involved in any springs with condensin, log it
            mass_DNA_unb(m,1:6) = mass_DNA(z,1:6);
            m = m+1; % increase the counter by 1
        else
        end
    end
else
    % if there is no condensin, then all are unbound
    mass_DNA_unb = mass_DNA;
end

for z = 1:size(mass_condensin,1)/11
    % lets find all the A, B1, and B2 beads
    condensin_A(z,1:5) = mass_condensin(1+(11*(z-1)),1:5);
    condensin_B1(z,1:5) = mass_condensin(9+(11*(z-1)),1:5);
    condensin_B2(z,1:5) = mass_condensin(10+(11*(z-1)),1:5);
end

if size(mass_condensin,1)>0
    % preassign this with 0.5 so that it can be checked later
    DNA_alpha = zeros(size(condensin_A,1),5)+0.5;
    
    % use the mass numbers of the A, B1, and B2 beads to find out what they are bound to
    for z = 1:size(condensin_A,1) % loop through the condensins
        for n = 1:size(springs_condensin,1) % loop through the springs
            % log all the DNA beads bound to A
            if condensin_A(z,5) == springs_condensin(n,2) && max(springs_condensin(n,1) == mass_condensin(:,5))==0
                DNA_alpha(z,1:5) = mass_coords(springs_condensin(n,1)+1,1:5);
            else
            end
            % log all the DNA beads bound to B1
            if condensin_B1(z,5) == springs_condensin(n,2) && max(springs_condensin(n,1) == mass_condensin(:,5))==0
                DNA_beta1(z,1:5) = mass_coords(springs_condensin(n,1)+1,1:5);
            else
            end
            % log all the DNA beads bound to B2
            if condensin_B2(z,5) == springs_condensin(n,2) && max(springs_condensin(n,1) == mass_condensin(:,5))==0
                DNA_beta2(z,1:5) = mass_coords(springs_condensin(n,1)+1,1:5);
            else
            end
        end
    end
    
    y = 1; % used to count loops
    
    for z = 1:size(condensin_A,1) % loop through the condensins
        
        end_check1 = 0; % used to look for the DNA end
        end_check2 = 0; % used to look for the DNA end
        
        bead_test1 = DNA_alpha(z,5); % assign the alpha bead as the start of the search
        bead_test2 = DNA_alpha(z,5); % assign the alpha bead as the start of the search
        
        s = 1; % used to log the coordinates
        
        % loop up through the DNA springs to find the beta1 bead or an end
        while end_check1 == 0 && bead_test1 ~= DNA_beta1(z,5)
            
            m = 0; % used as a check
            
            for n = 1:size(springs_DNA,1)
                if springs_DNA(n,1) == bead_test1
                    % assign the new bead to bind to
                    bead_new1 = springs_DNA(n,2);
                    
                    m = 1; % check changed to 1
                else
                end
            end
            
            if m == 1;
                bead_test1 = bead_new1; % reassigns the test bead
                loop_log1(s,1,z) = bead_test1; % log it in a table
                s = s+1; % increase the counter by 1
            else
                end_check1 = 1; % ends the while loop because an end was reached
            end
            
        end
        
        loop_size(z,1) = s-1; % shows how big the loop/end is
        
        s = 1; % used to log the coordinates
        
        % loop down through the DNA springs to find the beta1 bead or an end
        while end_check2 == 0 && bead_test2 ~= DNA_beta1(z,5)
            
            m = 0; % used as a check
            
            for n = 1:size(springs_DNA,1)
                if springs_DNA(n,2) == bead_test2
                    % assign the new bead to bind to
                    bead_new2 = springs_DNA(n,1);
                    
                    m = 1; % check changed to 1
                else
                end
            end
            
            if m == 1;
                bead_test2 = bead_new2; % reassigns the test bead
                loop_log2(s,1,z) = bead_test2; % log it in a table
                s = s+1; % increase the counter by 1
            else
                end_check2 = 1; % ends the while loop because an end was reached
            end
            
        end
        
        loop_size(z,2) = s-1; % shows how big the loop/end is
        
        if end_check1 == 0 && end_check2 == 0
            if loop_size(z,1) < loop_size(z,2)
                % if the first loop is smaller, log the first loop
                loops{y,1} = loop_log1(1:loop_size(z,1)-1,1,z);
                y = y + 1; %increase the counter by 1
            elseif loop_size(z,1) > loop_size(z,2)
                % if the second loop is smaller, log the second loop
                loops{y,1} = loop_log2(1:loop_size(z,2)-1,1,z);
                y = y + 1; %increase the counter by 1
            else
            end
        elseif end_check1 == 0 && end_check2 == 1
            % if the second test hit an end, log the first loop
            loops{y,1} = loop_log1(1:loop_size(z,1)-1,1,z);
            y = y + 1; %increase the counter by 1
        elseif end_check1 == 1 && end_check2 == 0
            % if the first test hit an end, log the second loop
            loops{y,1} = loop_log2(1:loop_size(z,2)-1,1,z);
            y = y + 1; %increase the counter by 1
        else
            % if both loops hit ends, then it is not a loop and nothing is logged
        end
        
    end
else
end

condensin_loop_beads = []; % assign an empty matrix for concatinating

for z = 1:size(loops,1) % for all the loops between condensin ends
    % concatinate all loops between condensin heads
    condensin_loop_beads = cat(1,condensin_loop_beads,loops{z,1});
end

m = 1; % counter to log beads that are not between condensin heads
s = 1; % counter to log beads that are between condensin heads

for z = 1:size(mass_DNA_unb) % loop through all unbound DNA
     if max(mass_DNA_unb(z,5) == condensin_loop_beads)==0
         % if it is not contained between condensin heads, log it
         mass_DNA_free(m,1:6) = mass_DNA_unb(1:6);
         m = m+1; % increase the counter by 1
     else
         mass_DNA_con_loops(s,1:6) = mass_DNA_unb(1:6);
         s = s+1; % increase the counter by 1
     end
end

output.mass_coords = mass_coords; % a list of all the masses
output.springs = springs; % a list of all the springs
output.mass_DNA = mass_DNA; % a list of all the DNA masses
output.mass_condensin = mass_condensin; % a list of all the condensin masses
output.mass_cohesin = mass_cohesin; % a list of all the cohesin masses
output.springs_DNA = springs_DNA; % a list of all the DNA springs
output.springs_condensin = springs_condensin; % a list of all the condensin springs
output.springs_cohesin = springs_cohesin; % a list of all the cohesin springs
output.condensin_A = condensin_A; % a list of all the A ends of condensin
output.condensin_B1 = condensin_B1; % a list of all the B1 ends of condensin
output.condensin_B2 = condensin_B2; % a list of all the B2 ends of condensin
output.DNA_alpha = DNA_alpha; % a list of all the DNA beads that the A ends of condensin are attached to
output.DNA_beta1 = DNA_beta1; % a list of all the DNA beads that the B1 ends of condensin are attached to
output.DNA_beta2 = DNA_beta2; % a list of all the DNA beads that the B2 ends of condensin are attached to
output.mass_DNA_unb = mass_DNA_unb; % a list of all the unbound DNA
output.mass_DNA_free = mass_DNA_free; % a list of DNA that is not between condensin heads
output.mass_DNA_con_loops = mass_DNA_con_loops; % a list of DNA that is between condensin heads
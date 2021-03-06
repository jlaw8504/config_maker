%Make a config file for chromoShake with two linear chains, linked with
%cohesin and several condensin mediate loops.  Wtih loops at either end to
%prevent slipping off of cohesin rings.



%% Open file for writing
config_file = 'test_cen.cfg';
fid = fopen(config_file, 'w');

%% Write out header based on default chromoShake values
%Note: this is currently hardcoded
string = 'meta collision_scheme 1\nmeta temperature_Celsius 25\nmeta viscosity_centiPoise 1\nmeta effective_damping_radius 8e-09\nmeta dna_modulus_gigaPascal 2\nmeta dna_radius_nanometers 0.6\nmeta damping_radius_factor 0.8\nstructure {\n  random_force 2.78554e-11\n  mass_damping 4.44973e+09\n  mass_radius 4.5e-09\n  time_step 2e-09\n  collision_spring_constant 0.0565487\n  spring_damping_factor 0\n  random_number_seed 42\n  color 1\n';
fprintf(fid, string);


%% Hard code variables
line_beads = 171;
final_loop_beads = 510;
cohesin_number = 10;
bead_per_cohesin = 16;
mass_sep = 1e-8; %10 nm separation
mass_mass = 3.38889e-20;
spring_rest =  1e-8; %10 nm
spring_const = 0.226195;
hinge_const = 4.0715e-12;

%% Dumbbell 1 data (pinned to spindle)
%use discrete_dumbbell.m function to generate positions
mass_number = line_beads + 2*final_loop_beads;
[x_pos,y_pos] = discrete_dumbbell(line_beads,final_loop_beads, mass_sep);
z_pos = zeros(size(x_pos));

%spring list for chain1
%want connection between all consectutive masses
spring_list = [0:mass_number-2;1:mass_number-1]';
%condensin list
% cond_list = [0 14;36 133; 156 170];
cond_list = [0, final_loop_beads-1; mass_number-final_loop_beads-1, mass_number-1];
%hinge list, want three consecutive numbers until end
hinge_list = [0:mass_number-3;1:mass_number-2;2:mass_number-1]';

%% Dumbbell 2 data
%want to duplicate dumbbell 1 but shift everything by 10 nm in z
x2_pos = x_pos;
y2_pos = y_pos;
z2_pos = z_pos  + 1e-08;

%shift all other data by mass_number
spring_list2 = spring_list + mass_number;
cond_list2 = cond_list + mass_number;
hinge_list2 = hinge_list + mass_number;

%% Cohesin data
%Generate cohesin template based on mass separation
[base_y, base_z] = discrete_circle(bead_per_cohesin,mass_sep);
%Generate desired number of cohesins around the loops
coh_index1 = (1:cohesin_number/2)*(final_loop_beads)/(cohesin_number);
coh_index2 = repmat(size(x_pos,2),size(coh_index1)) - coh_index1;
coh_indexes = [coh_index1, coh_index2];
%Use offsets to alter y positions of cohesin rings
coh_y = repmat(y_pos(coh_indexes),[16 1]);
%use repmat to repeat x and z positions of rings
coh_x = repmat(base_y',[1 cohesin_number]);
coh_z = repmat(base_z',[1 cohesin_number]);
% get positions in x of beads of interest and add coh_x to it to position
% cohesin beads in x dimension
offset_x = repmat(x_pos(coh_indexes),[16 1]);
coh_x = coh_x + offset_x;
%% Print out chain 1
%Spring and hinge idx
spr_idx = 1;
hinge_idx = 1;
for i = 1:mass_number
    fprintf(fid,'mass %d\t %.6g\t %.6g %.6g %.6g\n',i-1,...
        mass_mass, x_pos(i),y_pos(i), z_pos(i));
    if (i-1) == spring_list(spr_idx,2)
        fprintf(fid,'spring %d %d %.1g %.6g\n',spring_list(spr_idx,1),...
            spring_list(spr_idx,2), spring_rest, spring_const);
    %increment spr_idx
    spr_idx = spr_idx +1;
    end
    if (i-1) == hinge_list(hinge_idx,3)
        fprintf(fid,'hinge %d %d %d %.5g\n', hinge_list(hinge_idx,1),...
            hinge_list(hinge_idx,2), hinge_list(hinge_idx,3), hinge_const);
    end

end
%% Print out chain 2
spr_idx = 1;
hinge_idx = 1;
for i = 1:mass_number
    fprintf(fid,'mass %d\t %.6g\t %.6g %.6g %.6g %d\n',i-1+mass_number,...
        mass_mass, x2_pos(i),y2_pos(i), z2_pos(i),2);
    if (i-1) == spring_list(spr_idx,2)
        fprintf(fid,'spring %d %d %.1g %.6g\n',spring_list2(spr_idx,1),...
            spring_list2(spr_idx,2), spring_rest, spring_const);
    %increment spr_idx
    spr_idx = spr_idx +1;
    end
    if (i-1) == hinge_list(hinge_idx,3)
        fprintf(fid,'hinge %d %d %d %.5g\n', hinge_list2(hinge_idx,1),...
            hinge_list2(hinge_idx,2), hinge_list2(hinge_idx,3), hinge_const);
    end

end
%% Print out Cohesins
spr_idx = 1;
hinge_idx = 1;
mass_idx = 2*mass_number;
for col = 1:cohesin_number
    for row = 1:bead_per_cohesin
        if row == 1
            %print out mass
        fprintf(fid,'mass %d\t %.6g\t %.6g %.6g %.6g %d\n',...
            mass_idx,mass_mass,coh_x(row,col),coh_y(row,col),coh_z(row,col),5);
        %update mass idx
        mass_idx = mass_idx + 1;
        
        elseif row > 1 && row ~= bead_per_cohesin
            %print out mass
            fprintf(fid,'mass %d\t %.6g\t %.6g %.6g %.6g %d\n',...
            mass_idx,mass_mass,coh_x(row,col),coh_y(row,col),coh_z(row,col),5);
        %print out spring
            fprintf(fid,'spring %d %d %.1g %.6g\n',mass_idx-1,mass_idx,...
            spring_rest, spring_const);
        %update mass idx
        mass_idx = mass_idx + 1;
        
        elseif row > 2 && row ~= bead_per_cohesin
            %print out mass
            fprintf(fid,'mass %d\t %.6g\t %.6g %.6g %.6g %d\n',...
            mass_idx,mass_mass,coh_x(row,col),coh_y(row,col),coh_z(row,col),5);
        %print out spring
        fprintf(fid,'spring %d %d %.1g %.6g\n',mass_idx-1,mass_idx,...
            spring_rest, spring_const);
        %print out hinge
        fprintf(fid,'hinge %d %d %d %.5g\n',mass_idx-2,mass_idx-1,mass_idx,...
                hinge_const);
            %update mass idx
            mass_idx = mass_idx + 1;
            
        elseif row == bead_per_cohesin
            %print out mass
            fprintf(fid,'mass %d\t %.6g\t %.6g %.6g %.6g %d\n',...
            mass_idx,mass_mass,coh_x(row,col),coh_y(row,col),coh_z(row,col),5);
        %print out spring
        fprintf(fid,'spring %d %d %.1g %.6g\n',mass_idx-1,mass_idx,...
            spring_rest, spring_const);
        %print out hinge
        fprintf(fid,'hinge %d %d %d %.5g\n',mass_idx-2,mass_idx-1,mass_idx,...
                hinge_const);
            %print out final springs
            fprintf(fid,'spring %d %d %.1g %.6g\n',mass_idx,mass_idx-(bead_per_cohesin -1),...
            spring_rest, spring_const);
            %print out second-to-last hinge
             fprintf(fid,'hinge %d %d %d %.5g\n',mass_idx-1,mass_idx,mass_idx-(bead_per_cohesin -1),...
                hinge_const);
            %print out last hinge
            fprintf(fid,'hinge %d %d %d %.5g\n',mass_idx,mass_idx-(bead_per_cohesin -1)...
                ,mass_idx-(bead_per_cohesin -2),hinge_const);
            %update the mass idx
            mass_idx = mass_idx + 1;
        end
        
        
    end
  
end

%%Print out condensins
cond_total = [cond_list;cond_list2];
for k = 1:size(cond_total,1)
    fprintf(fid,'spring %d %d %.1g %.6g\n',cond_total(k,1),...
            cond_total(k,2), spring_rest, spring_const);
end
fprintf(fid,'}\n');
fclose(fid);

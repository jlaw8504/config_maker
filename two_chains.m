%Make a config file for chromoShake with two linear chains, linked with
%cohesin and several condensin mediate loops.  Wtih loops at either end to
%prevent slipping off of cohesin rings.



%% Open file for writing
config_file = 'test.cfg';
fid = fopen(config_file, 'w');

%% Write out header based on default chromoShake values
%Note: this is currently hardcoded
string = 'meta collision_scheme 1\nmeta temperature_Celsius 25\nmeta viscosity_centiPoise 1\nmeta effective_damping_radius 8e-09\nmeta dna_modulus_gigaPascal 2\nmeta dna_radius_nanometers 0.6\nmeta damping_radius_factor 0.8\nstructure {\n  random_force 2.78554e-11\n  mass_damping 4.44973e+09\n  mass_radius 4.5e-09\n  time_step 2e-09\n  collision_spring_constant 0.0565487\n  spring_damping_factor 0\n  random_number_seed 42\n  color 1\n';
fprintf(fid, string);


%% Hard code variables
mass_number = 171;
cohesin_number = 10;
bead_per_cohesin = 16;
mass_sep = 1e-8; %10 nm separation
mass_mass = 3.38889e-20;
spring_rest =  1e-8; %10 nm
spring_const = 0.226195;
hinge_const = 4.0715e-12;

%% Chain 1 data (chain1 pinned to spindle)
%increment distance only in xdir
%pre-allocate
x_pos = zeros([mass_number,1]);
y_pos = x_pos;
z_pos = x_pos;
for n = 1:mass_number
    x_pos(n) = (n-1) * mass_sep;
end

%spring list for chain1
%want connection between all consectutive masses
spring_list = [0:mass_number-2;1:mass_number-1]';
%condensin list
cond_list = [0 14;36 133; 156 170];
%hinge list, want three consecutive numbers until end
hinge_list = [0:mass_number-3;1:mass_number-2;2:mass_number-1]';

%% Chain2 data
%want to duplicate chain1 but shift everything by 10 nm in y
x2_pos = x_pos;
y2_pos = y_pos + 1e-08;
z2_pos = z_pos;

%shift all other data by mass_number
spring_list2 = spring_list + mass_number;
cond_list2 = cond_list + mass_number;
hinge_list2 = hinge_list + mass_number;

%% Cohesin data
%evenly space cohesin in a ring
[coh_x, coh_z] = discrete_circle(bead_per_cohesin,mass_sep);
coh_y = zeros([size(coh_x,1),1]);

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
fprintf(fid,'}\n');
fclose(fid);

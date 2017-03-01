function [x,y] = discrete_circle(num_sides, mass_sep)

%list of angles in radians
t = (1/(2*num_sides):1/num_sides:1)*2*pi;
%Calc interior angle, this is used to calculate radius of circle
int_angle = ((num_sides-2)*180)/num_sides;
half_angle = int_angle/2;
radius = (mass_sep/2)/cos(deg2rad(half_angle));
x = cos(t)*radius;
y = sin(t)*radius;
end
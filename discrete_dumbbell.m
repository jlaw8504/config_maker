function [x_pos,y_pos] = discrete_dumbbell(line_beads,final_loop_beads, mass_sep)
%% Generate central line
x_pos = zeros([1, line_beads]);
y_pos = 0:mass_sep:(line_beads-1)*mass_sep;

%% Generate final loops
[x_loop, y_loop] = discrete_circle(final_loop_beads,mass_sep);
%flip the x_loop and y loop for ordering purposes
x_loop_flip = fliplr(-x_loop);
%shift the circle such that the first position is equal to the final x and
%y positoins of the omega shape + mass_separation
x_diff_first = x_pos(1) - x_loop_flip(length(x_loop_flip));
x_loop_first = x_loop_flip + x_diff_first + mass_sep;
y_diff_first = y_pos(1) - y_loop(length(y_loop));
y_loop_first = -y_loop + y_diff_first;
%position final loop
x_diff_last = x_pos(length(x_pos)) - x_loop_flip(1);
x_loop_last = x_loop_flip + x_diff_last + mass_sep;
y_diff_last = y_pos(length(y_pos)) - y_loop(1);
y_loop_last = -y_loop + y_diff_last;

%% Concatenate final loops to x and y pos
x_pos = [x_loop_first, x_pos, x_loop_last];
y_pos = [y_loop_first, y_pos, y_loop_last];
end
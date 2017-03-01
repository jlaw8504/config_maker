function [x_pos,y_pos] = discrete_omega(circle_beads,line_beads,final_loop_beads, mass_sep)
%% Generate central circle
[x_pos,y_pos] = discrete_circle(circle_beads,mass_sep);

%% attach flanking beads in straight line in y
%get final x position
final_x = x_pos(length(x_pos))+10;
%repmat x_position for flanking lines
flank_line_x = repmat(final_x, [1 line_beads]);
%create y postions based on first circle y-position
first_flank_line_y = y_pos(1):mass_sep:y_pos(1)+((line_beads-1) * mass_sep);
%flip the array order left to right
first_flank_line_y = fliplr(first_flank_line_y);
%create y positions based on last circle y-position
last_flink_line_y = y_pos(length(y_pos)):-mass_sep:...
    y_pos(length(y_pos)) - ((line_beads-1)*mass_sep);
%concatenate the x and y positions of the line to the circle
x_pos = [flank_line_x, x_pos, flank_line_x];
y_pos = [first_flank_line_y, y_pos, last_flink_line_y];

%% Generate final loops
[x_loop, y_loop] = discrete_circle(final_loop_beads,mass_sep);
%flip the x_loop and y loop for ordering purposes
x_loop_flip = fliplr(-x_loop);
%shift the circle such that the first position is equal to the final x and
%y positoins of the omega shape + mass_separation
x_diff_first = x_pos(1) - x_loop_flip(length(x_loop_flip));
x_loop_first = x_loop_flip + x_diff_first + mass_sep;
y_diff_first = y_pos(1) - y_loop(length(y_loop));
y_loop_first = y_loop + y_diff_first;
%position final loop
x_diff_last = x_pos(length(x_pos)) - x_loop_flip(1);
x_loop_last = x_loop_flip + x_diff_last + mass_sep;
y_diff_last = y_pos(length(y_pos)) - y_loop(1);
y_loop_last = y_loop + y_diff_last;

%% Concatenate final loops to x and y pos
x_pos = [x_loop_first, x_pos, x_loop_last];
y_pos = [y_loop_first, y_pos, y_loop_last];
end
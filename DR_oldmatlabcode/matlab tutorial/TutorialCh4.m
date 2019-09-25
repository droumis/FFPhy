% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Chapter 4 - Matlab Tutorial
% by Suzanne Baker
%
% Included in Chapter 4:
%
% figure . plot and plot options . clf . subplot . axis, title, xlabel, ylabel . LineWidth 
% text and gtext . ginput . legend  . plot3 and rotate3d
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% 1. figure - this command creates a new figure window.  You can specify the figure number by
%       typing figure(num).  
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% 2. plot - there are a lot of options that go along with plot, so its good to take a look at:
clear all % clear all variables in the matlab workspace
close all % closes all figures in the matlab workspace  
help plot

%       Here I will go through options I find handy.

pi_group=[0:pi/100:2*pi]; %just make an array of numbers

figure %creates a new figure window
plot(pi_group) % creates plot with default plotting styles
% if you are plotting 1 group of numbers, the x axis will default to 1, 2, 3...

figure(2)
plot(pi_group,pi_group,'r:') % r makes the line red, : makes it dotted, there are many options

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% 3. clf - clears current figure
figure(2)
clf

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% 4. subplot - this allows many little plots within one figure

figure(2)
subplot(2,1,1)
plot(pi_group)
subplot(2,1,2)
plot(pi_group,pi_group)

% now its easier to compare the x axes of the 2 figures.  You also can control the axes limits with axis:

figure(2)
subplot(2,1,1)
axis([0 200 0 2*pi])
subplot(2,1,2)
axis([0 2*pi 0 2*pi])


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% 5. axis, title, xlabel, ylabel

% some interesting axis options to play with

figure(3)
plot(pi_group,sin(pi_group))

% a. execute then look at figure 3
axis equal

% b. execute then look at figure 3
axis tight 

% c. execute then look at figure 3
axis off

% title, xlabel, ylabel

title('Insert Snappy figure title here')
axis on %if axis off was the last axis command, then you won't see axes labeling
xlabel('0->2pi')
ylabel('sine(0->2pi)')

close all %this will close all current figures to help clean up the situation

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% 6. LineWidth
figure(1)
clf
gg=plot(pi_group,sin(pi_group),'c.-');
hold on% allows another plot to be added without overwriting the current plot
hh=plot(pi_group,cos(pi_group),'g--');
hold on
ii=plot(pi_group,tan(pi_group),'y-');
axis([0 2*pi -1 1])

% There are many options that can be set using the set command.  Assigning each
% specific plot to a specific variable allows you to assign such options.  To see
% a list of possibilities, use the get command

get(ii)

% I don't use most of these and don't know what a lot of them do, but I do use:

set(ii,'LineWidth',4)

% Line Width can also be changed from the figure window:
%   -click on the arrow below the menu bar (between the printer and the A character)
%   -click once on a line and it will select that line, if you click twice, a box will 
%       appear with many options... these options should be simple enough to understand
%   -double click on the axis to get axis options
%   -click on A to add text to the picture, the arrow and line drawing commands work
%       similar to Word

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% 7. text and gtext - text specifies coordinates for a piece of text to be plotted
text(0,0,'origin')

% gtext puts the string up where the user clicked (after executing these lines, go into
%   figure 1, crosshairs should appear, use your mouse to click anywhere.
figure(1)
gtext('you clicked here')

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% 8. ginput - this command returns coordinates of a point where the user clicks on the graph
%   Check out ginput options

help ginput

figure(1)
[x y]=ginput(1); % the 1 makes its so the user only clicks once, you can see the x and y
%                   coordinates chosen at the command prompt
disp(['x = ' num2str(x)])
disp(['y = ' num2str(y)])

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% 9. legend - creates a legend for such a hectic plot.  Check out the options:

help legend

legend('cyan - sin','green - cos','yellow - tan',3)
% Using the select arrow in the figure window, the legend can be moved around as well.

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 10. plot3 and rotate3d - plots 3d data

help plot3
close all
figure(1)
plot3(cos(pi_group),sin(pi_group),pi_group)
figure(1)
rotate3d % this is a cool command which lets you rotate your image around with your mouse



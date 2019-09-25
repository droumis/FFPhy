% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Chapter 3 - Matlab Tutorial
% by Suzanne Baker
%
% Included in Chapter 3:
%
% input . paths . functions . 
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% 1. paths - Matlab will run programs that are either in the current directory (found by typing pwd at the 
% command prompt) or already set in the path.  To see what directories Jarrod has already set
% in the default paths or add more directories, go to File-SetPath.  If you know a program exists
% but are having problems running it, check your paths.  As a quick exercise, because the default
% spm right now is spm99, lets try running spm.  At your command prompt just type spm, and notice
% it is most likely spm99.  Quit that.  Got to SetPath and choose "Add with subfolders.
% In the box that appears, type in (without the %):

%/usr/bic/matlab6p1/toolbox/spm2b

% 8 new folders should appear at the top of your path.  Click close (you can save it so everytime
% you open matlab this will be your personal path).  Now at the command prompt, when you type in spm
% spm2 should start up.

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% 2. input - this is a command which allows user input of a specific variable when the program 
% is running.  Following a string to prompt he user can input a string or numbers.

t = input('What are your initials: ','s');

% so now if you check t, you will see it has been saved as characters that you just typed in.

% Another option is:

s = input('How many trials will you perform \n type 1 for 1 \n type 2 for 2 \n type 3 for 3 \n :');

% (notice the use of \n to make a new line)

%       

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Functions are those handy things that your fellow matlabbers have already been writing for the lab.
% SPM code is also written as a bunch of functions.  

% Essentially a function is just a program that you can call from the command line.  However
% Matlab will only keep specific variables from a function in its working memory when a function
% is called (the list of variables can be seen by typing whos at the prompt).  From the function 
% example you'd find if you type in 'help function' at the command prompt (don't execute this):

           function [avg,stdev] = stat(x)
           %STAT Interesting statistics.
           n = length(x);
           avg = mean(x,n);
           stdev = sqrt(sum((x-mean(x,n)).^2)/n);
           
           
% If I called this function from the command line, I'd type something like:

%> [variable1, variable2]=stat(timeseries)

% In this situation, timeseries, variable1 and variable2 are held in working memory, but not n, mean, 
% stdev, or x.  This is an important point.  Often if I'm debugging a function I'm writing, I will
% either highlite the program and hit F9 until I know its running properly (this way I have access
% to all variables to help in the troubleshooting process) or you can put disp commands output
% variable values as the function runs.

% If you were to implement this into a program, the first line of the program should be:

           function [mean,stdev] = stat(x)

% The lines that follow that are commented (preceeded by %) will be displayed if you type
% 'help stat' at the command line. 

% Copy this short program into a new .m file (File-New) and save it as stat.m (functions 
% should be saved as the name of the function in your directory.  Execute the following
% line (highlite and hit F9), so that you have a timeseries to feed into the stat program.

timeseries=[1 2 4 6 3 9 5 6 12 4 23 1 2 4 6 2 2 1 6 3];

% Now run the stat program from the command prompt to see if it works.  

% I have put a simple program written as a function that calls spm2 functions called rot_1vol_spm.m
% in the directory with the tutorials.  This is a good fucntion to run and understand, especially
% because it has spm2 functions in it to read spm .img files and write .img files.  I suggest
% you run rot_1vol_spm.m and take a look at the code.

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Here is an exercise if you want to test if you understand functions.  Create a function called
% scale.m where you input a matrix of any size, the current min and max and new min and max.  

% So if I have the matrix [ 1 2 3 4 5 20] and the current range of the matrix could be 0 to 25, 
% but I want the matrix to be scaled between 0 and 1, then the output matrix answer should be
% [0.0400    0.0800    0.1200    0.1600    0.2000    0.8000]... that is just a simple example.  
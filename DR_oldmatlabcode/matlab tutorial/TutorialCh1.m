%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Chapter 1 - Matlab Tutorial
% by Suzanne Baker
%
% Included in Chapter 1:
%
% . ; . who and whos . clear . ls . pwd . chdir and cd . help . lookfor . save and load . how information is stored in variables .  matrix operations
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% 1. % - First thing to notice is a percent sign (%) allows you to comment lines.
%       I can also put a % after something to add a comment.  Anything that won't 
%       be executed as code appears as green on my screen.

a=5 % this makes a equal to the number 5


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 2. F9 - The next important thing to note is if you highlite a line with your
%       and while the line is highlited you hit the F9 key, that will run the lines 
%       highlited as if they were typed after the prompt in the command window.  
%       This is handy if you don't want to run the entire program or
%       if you just want to test a portion of the program.  If there is output
%       it will put shown in the command window.


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% 3. ; - Following any line with a ; will keep an output from being shown on 
%       the screen.  This is extremely important for manipulating HUGE matrices.
%       Just to show you really fast:
%           i. highlite b= line and hit F9 
            b=[1 2 3 4 5 6 7 8 9];
%           ii. Check your command window, it should show:
%                >>             b=[1 2 3 4 5 6 7 8 9];
%                >> 
%           iii. highlite the c= line and then hit F9
            c=[10 20 30 40 50 60]
%           iv. Check your commant window and it should show:
%                >>             c=[10 20 30 40 50 60]
%
%                c =
%
%                    10    20    30    40    50    60
%
%                >> 

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% 4. who and whos - on the command line, type both who <enter> and whos <enter>
%       Who and whos allow you to see the variables that are actively being stored.
%       Thus far I would get for who:
%>> who

%Your variables are:

%a  b  c  

%        And for whos:
        
%>> whos
%  Name      Size           Bytes  Class

%  a         1x1                8  double array
%  b         1x9               72  double array
%  c         1x6               48  double array

%Grand total is 16 elements using 128 bytes


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% 5. clear - Clear allows you to clear all or specific variables.  Execute this:

d=5;
e=6;
f=7;
clear b
whos

%   That demonstrates clearing 1 variable.  You could also list the variables to be cleared

clear d f
whos

%   You can list as many variables you want.  You can also choose to clear all the variables
%       (which is a popular command to put at the beginning of a program).

clear all
whos

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% 6. ls - this will list all files in the directory you are in.

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% 7. pwd - this will give the current path you are in

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% 8. cd or chdir - these commands help you navigate. 

mypwd=pwd;

cd .. % takes you back a directory

chdir(['/raid0/despo/slbaker/042203/']) % you have to give your path

cd exp5_ge_2shot.fid % this will dump you into that directory

chdir(mypwd); % this puts you back in whatever original directory you started in


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% 9. help - on the command line, you can type help followed by a command name to read more about the command:
%       I will now sneakily add in important commands such as lookfor and save that you should read the 
%       help blurb about.  Help blurbs are for functions and they are the set of commented lines until the program
%       starts (little thing to note for when you are writing your own functions).

help lookfor
help save
help load

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% 10. save and load - the simplest example of save is :

ex1='string';
ex2=[1 2 3 4; 2 3 4 5; 3 4 5 6];
save test_filename ex1 ex2
ls % you should see a filename called test_filename.mat listed
clear all % no more variables exist now
whos
load test_filename
whos % now you can see that you have loaded your parameters back in


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% 11. This brings up different types of ways data can be stored (number, string, convert from string to number and vice versa, and structure).  I have already shown simple
%       numerical versions.  Here are a few more:

% ~~~~~~~~~~~~~~~~~~~~~~~~~NUMBER~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
simple = 5;
a = [1 2 3 9;4 5 6 12]

%   Other easy ways to create matrices:

biga = [1:20]

incUPa = [1:2:20]

incDOWNa = [100:-3:24]

%   You can check the size:

size(a)

%   You can access different parts of the matrix, execute each line independently to see the output

a(1,2)
a(2,:)
a(1,2:4)

%   Use : to access all values in that row/column or something like 2:4 to only access 2 through 4 values

%   You can also create multi dimensional data, which we will get into later.


% ~~~~~~~~~~~~~~~~~~~~~~~~~STRING or CHAR~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%   A string is another data type possibility which will be listed as char in the whos command:

b = 'left right'
c = ['tire';'plug';'food']
d = '142';

%   You can do exactly the same in terms of accessing values of the matrix, even though the matrix contains string values
%       Execute each line independently:

size(b)
b(2:6)
size(c)
c(1,:)
c(3,3)

% ~~~~~~~~~~~~~~~~~~~~~~~~~CONVERT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


%   You can convert a number value to string 

g = num2str(simple)
whos

%   You can convert a string value into a number (notice h is then listed as a double value after the whos command:

h = str2num(d);
whos

% ~~~~~~~~~~~~~~~~~~~~~~~~~STRUCTURE or CELL ARRAY~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear all

%   Little squiggly brackets creates a structure (a.k.a. cell array):
i={23};
whos

%   This is handy for saving a lot of different information under one variable and not worrying about combining
%       different string and number information as well as matrices of different dimensions.

g.subjectname='Mark DEsposito';
g.subjectbday='01/01/68';
g.timeseries=[1 5 7 8 9 1 2 3 5 4 2 1 7 3 7 2 3 1 4];
g.slice=[1 1 1 1;2 2 2 2;3 3 3 3];
g.hemi{1}=[1 2 56 3 3 7 2 4 2;1 2 56 3 3 7 2 4 2];
g.hemi{2}=[0 1 0  1 0 1 1 1 0;1 1 0  1 0 1 0 0 0];


%  So now see what happens when you type g, g.slice, and g.timeseries at the prompt of the command line.

%  To pull out information highlite and F9 each line individually:

g
g.timeseries(3:6)
g.slice(:,2)
g.hemi
temp=g.hemi{2}
temp(2,1)


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% 12. Lastly, some matrix functions.  Matlab has a great demo of this.  Go to the command window -> Help -> Demos.  
%       Choose Matrices in the left window and Basic Matrix Operations in the right window.  On the bottom right, 
%       Click Run Basic Matrix... Feel free to look at all of the slides, but the ones of most interest are:
%       Slides of interest: 2, 3, 4, 7, 8, 9, 10, 11, 12, 19, 20

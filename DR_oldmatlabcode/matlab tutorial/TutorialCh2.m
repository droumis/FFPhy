% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Chapter 2 - Matlab Tutorial
% by Suzanne Baker
%
% Included in Chapter 2:
%
% disp . if . switch . for . min/max/mean/median/std/floor/ceil/round . find . reshape
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% 1. disp - this is an easy one.  Its good for putting into programs here and there
%           for debugging
%           You can read what Matlab says about this command by executing:
help disp

a=[4 10 7 9 23];

disp(['The 3rd value of a is ' num2str(a(3))]);


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% 2. if - If statement allows you to execute a command based on the value of a variable.
%       == for equals                   ~= for not equal
%       >  for greather than            >= greater than or equal to
%       <  for less than                <= less than or equal to
%       &  and                          |  or 

%           You can read what Matlab says about this command by executing:
help if


%  Change the value of a to see what is displayed

clear a
a = 5;

if a == 1
    disp('a = 1')
elseif a > 2 & a < 3.5
    
    disp(['a is greater than 2 and less than 3.5, specifically a = ' num2str(a)])
    % here I had to put square brackets around the disp because there was more than just
    % a phrase in quotes.
    
elseif a == 4 | a == 5
    
    disp(['a is 4 or 5, specifically, a = ' num2str(a)]);
    
else
    
    disp(['a is not 1, 4, 5 or between 2 and 3.5, a = ' num2str(1)]);
    % else means its none of the previous choices
end


% you can do this with any number (including 0) elseif statements.  If you use else, you only use 1 else statement.

if a > 2
    b=pwd;
    disp(pwd)
end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% 3. switch - Switch is a more constrained version of if statements where you only consider if a variable is equal
%           to a value.  I don't use it much, but its simple enough to go over
%           Change the value of h to see what is displayed

%           You can read what Matlab says about this command by executing:
help switch


h = 3;

switch h
case 1
    disp('case is 1');
case 2
    disp('case is 2');
case {3,4}
    disp('case is 3 or 4');
otherwise
    disp('case is not 1, 2, 3, 4');
end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% 4. for - For loops are so important and so popular.

%           You can read what Matlab says about this command by executing:
help for

for i=1:4
    disp(['Loop ' num2str(i)]);
    disp(['i = ' num2str(i)]);
    j=i*4;
    disp(['j = ' num2str(j)]);
end


% Loops can count up or down by whatever increment you want:
% Things to consider -      1. what happens if I made incnum negative,
%                           2. what happens if incnum is some sort of decimal (like 2.2)
%                           3. can you make a loop that counts backwards
startnum=2;
endnum=40;
incnum=4;
count=0;

for value=startnum:incnum:endnum
    count=count+1;
    disp(['Loop ' num2str(count)]);
    disp(['value = ' num2str(value)])
end
    

% Lastly, nested loops are important.  In this example, I'm making a 3 dimensional matrix

clear all
for i=1:4
    for j=1:2
        for k=1:7
            mat(i,j,k)=i*100+(j-1)*10+k/10;
        end
    end
end

[x y z]=size(mat)
% in this example, x=maximum of i, y=maximum of j, z=maximum of k
mat

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 5. min, max, mean, median, std, floor, ceil, round - I think these are pretty self explanatory, 
% especially once you play with the numbers below, but if you are curious, just use the help command

temp=[1.1 2.2 3.3 4.4 5.5 6.6 7.7];

% this will use the mat matrix that was made above
% you can always assign the output of the matrix to a value, such as

output=min(temp);
disp([output])

min(mat)
min(min(mat))
min(min(min(mat)))

% max works in the same way, when you have such a big matrix, all of these can be combined in different ways

max(temp)
max(max(max(mat)))

% mean

mean(temp)

% median

median(temp)

% std

std(temp)

% floor rounds all the numbers in the matrix down

floor(temp)

% ceil rounds them all up

ceil(temp)

% round just rounds them

round(temp)

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 6. find - this is actually one of my favorites.  It finds the values in a matrix that meet a 
%       specific criteria (good for thresholding)
%       read the help and then I'll give examples

help find

% just to create a new matrix really fast:

clear all
img=[1:32]'*[1:32];
% just to give you an idea of the profile of this image

disp(['min value of img - ' num2str(min(min(img)))]);
disp(['max value of img - ' num2str(max(max(img)))]);
disp(['mean value of img - ' num2str(mean(mean(img)))]);

% lets say we want to use all voxels above a threshold of 800

[i j] = find(img > 800);

[sizei NoNeed]=size(i)

for k=1:sizei
    disp(['i = ' num2str(i(k)) ' and j = ' num2str(j(k))]);
end

% i and j are the indeces to the matrix.  Then when you want to look at the values in img that are > 800

for k=1:sizei
    if k==1
        disp(['actual value (greater than 800, we hope) - ' num2str(img(i(k),j(k)))])
    else
        disp(['actual value - ' num2str(img(i(k),j(k)))])
    end
    
end

% the find statement can be as complicated as you want, for example

clear all
a=[1 2 3 4 5 6 7 4 2 1 6 3 1 1 9 3 5 2 1];

h = find(a > 2 & a < 7)

[NoNeed sizeh]=size(h);
    
for i=1:sizeh
    disp(['Index #' num2str(h(i)) ' is a matrix value ' num2str(a(h(i))) ])
end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 7. reshape - reshape lets you changes the organization of a matrix.  For instance, if you get a
%       2 by 3 matrix, you could make it 6 by 1, 1 by 6 or 3 by 2.  Be very careful with this command.

help reshape

clear all
a=[1 2 3; 100 200 300];

[sizex sizey=size(a);
% make a 1 by 6
tempa=reshape(a,1,6)

% put a back to its original dimensions (you really have to be able to do this!!)
reshape(tempa,sizex, sizey)

% other reshaping options
reshape(a,6,1)
reshape(a,3,2)



function raise_me (hObject, event, me)
% Pulls the clicked object to the top of the ui stack; useful
% for raising partially masked objects to the front of a plot.
% GUI-shortcut for 'uistack':  Left-clicking brings to the top,
% right-clicking sends to bottom.

%   Last Modified By: sbm on Sat Aug 21 00:57:32 2004

myfig = gcf;
oldptr = get(myfig, 'Pointer');
set(myfig, 'Pointer', 'watch');

switch(get(myfig, 'SelectionType')),
    case ('normal'),
        uistack(me, 'top');
    case ('alt'),
        uistack(me, 'bottom');
end

set(myfig, 'Pointer', oldptr);
drawnow;

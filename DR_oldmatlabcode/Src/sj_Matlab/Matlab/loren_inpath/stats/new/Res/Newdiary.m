% NEWDIARY: closes old diary and begins new one.

diary off;

if (exist('c:\matlab5\session.log')),
  delete c:\matlab5\session.log;
end;

diary c:\matlab5\session.log;


function coremex(mexstring)
%COREMEX           Macro to compile CORE_ functions.
%   COREMEX(MEXSTRING) is an macro that calls:
%       mex MEXSTRING CORE_library.c CORE_mextools.c msvc_libmwlapack.lib
%   where MEXSTRING is "[OPTS] SOURCEFILE".  This compiles the MEX code in
%   SOURCEFILE and links to both CORE_library functions and BLAS/LAPACK.
%
%   'msvc_libmwlapack' is the static link definition for Microsoft
%   Visual C++; this file can be edited to work with Borland C++ Builder
%   link definitions instead by replacing 'msvc_libmwlapack' with
%   'libmwlapack'.

eval(['mex ' mexstring ' CORE_library.c CORE_mextools.c msvc_libmwlapack.lib']);  % MSVC
%eval(['mex ' mexstring ' CORE_library.c CORE_mextools.c libmwlapack.lib']);      % Borland
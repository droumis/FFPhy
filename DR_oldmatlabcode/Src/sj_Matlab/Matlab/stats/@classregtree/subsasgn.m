function [varargout] = subsasgn(varargin)
%SUBSASGN Subscripted reference for a CLASSREGTREE object.
%   Subscript assignment is not allowed for a CLASSREGTREE object.

%   Copyright 2006 The MathWorks, Inc.
%   $Revision: 1.1.6.2.2.1 $  $Date: 2007/01/16 03:38:08 $

error('stats:classregtree:subsasgn:NotAllowed',...
      'Subscripted assignments are not allowed.')
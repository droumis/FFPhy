function [varargout] = subsasgn(varargin)
%SUBSASGN Subscripted reference for a PIECEWISEDISTRIBUTION object.
%   Subscript assignment is not allowed for a PIECEWISEDISTRIBUTION object.

%   Copyright 2006-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2007/01/16 03:38:20 $

error('stats:piecewisedistribution:subsasgn:NotAllowed',...
      'Subscripted assignments are not allowed.')
%NEWVALUES = VECTORFILL(VALUES, FILLVALUE)
%NEWVALUES = VECTORFILL(VALUES, FILLVALUE,DIRECTION)
%NEWVALUES = VECTORFILL(VALUES, FILLVALUE,DIRECTION,INDECES)
%
%Takes the input vector VALUES, and changes the numbers at the 
%specified indeces that equal FILLVALUE to the nearest number that
%does not equal FILLVALUE.  DIRECTION can be either -1,0, or 1, and
%specifies whether to fill the numbers with the last non-FILLVALUE (-1)
%or the next non-FILLVALUE (1), or the nearest non-FILLVALUE (0 - default).  
%If no non-FILLVALUE exist in that direction (1 or -1)
%then those values remain unchanged. The user can also specify which indeces of VALUES
%operate on (default is all indeces).
%
%
function selectPlaceCells

% parameters
runAll=        0; 
newPassTime=   0; 
newPlaceField= 1; 
cleanup=       1;

if(runAll | exist('select-all.mat','file')~=2)
    selectAllCells;
end
if(newPassTime)
    disp('**** getting passtimes ****')
    fmaux.selectid= 'all';
    fmaux.select=['select-' fmaux.selectid];
    selectPasstimes
    deleteFile('stats-all-passtimes.mat');
end
if(newPlaceField | newPassTime)
    disp('**** determining place fields ****')
    fmaux.selectid= 'all-passtimes';
    fmaux.select=['select-' fmaux.selectid];
    selectPlaceFields
    deleteFile('stats-all-passtimes-placefields.mat');
end

if(newPlaceField | newPassTime | cleanup)
    disp('**** clean up ****')
    fmaux.selectid= 'all-passtimes-placefields';
    fmaux.select=['select-' fmaux.selectid];
    cleanSelect('select-placecells');
%    deleteFile('stats-placecells.mat');
end

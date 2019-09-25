function animaldirectory = animaldef(animalname,daynumber)

load rootdir;

switch lower(animalname)
    case 'dinah'
       animroot = fullfile(rootdir,'dinah');
       animdir = {'din0525'};
    case 'phineas'
       animroot = fullfile(rootdir,'phineas-new');
       animdir = {'','','','','','','','phin1118','','phin1119b'};
    case 'jabez'
       animroot = fullfile(rootdir,'jabez');
       animdir = {'','','','','','','jab0715','','jab0716','jab0717'};
    case 'kohath'
       animroot = fullfile(rootdir,'kohath');
       animdir = {'','','ko0916','ko0917','ko0918','ko0919','ko0920','ko0921','ko0922','ko0923','ko0924','ko0925'};
    case 'mephibosheth'
       animroot = fullfile(rootdir,'mephibosheth');
       animdir = {'meph0425','meph0426','meph0427','meph0428','meph0429','meph0430','meph0501','meph0502','meph0504','meph0506','meph0506_2','meph0507'};
    case 'nebuchadnezzar'
       animroot = fullfile(rootdir,'nebuchadnezzar');
       animdir = {'','neb1219','neb1220','neb1221','neb1222'};
    case 'obed'
       animroot = fullfile(rootdir,'obed');
       animdir = {'ob0220','ob0221','ob0222','ob0223','ob0224'};
    otherwise
        error(['Animal ',animalname, ' not defined.']);
end

animaldirectory = fullfile(animroot,animdir{daynumber});

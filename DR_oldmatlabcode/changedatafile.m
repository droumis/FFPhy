function tmpVal = changedatafile(prefix,varName)

files = dir([prefix, varName, '*']);
for f = 1:length(files)
    load(files(f).name);
    eval(['tmpVal = ',varName,';']);
    for C1 = 1:length(tmpVal)
        if (~isempty(tmpVal{C1}))
            for C2 = 1:length(tmpVal{C1})
                if isstruct(tmpVal{C1}{C2})
                    tmpVal{C1}{C2} = changeStruct(tmpVal{C1}{C2});
                elseif ~isempty(tmpVal{C1}{C2})
                    for C3 = 1:length(tmpVal{C1}{C2})
                        if (~isempty(tmpVal{C1}{C2}{C3}))
                            for C4 = 1:length(tmpVal{C1}{C2}{C3})
                                if isstruct(tmpVal{C1}{C2}{C3}{C4})
                                    tmpVal{C1}{C2}{C3}{C4} = changeStruct(tmpVal{C1}{C2}{C3}{C4});
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    eval([varName ' = tmpVal;']);
    save(files(f).name, varName);
end


function outStruct = changeStruct(inputStruct)
    
outStruct = inputStruct;
fieldcount = 0;
noerror = 1;
while (~isempty(inputStruct.fields) && noerror)
    [tmpfield inputStruct.fields] = strtok(inputStruct.fields);
    fieldcount = fieldcount + 1;
    try
        outStruct = setfield(outStruct,tmpfield,inputStruct.data(:,fieldcount));
    catch
        noerror = 0;
        disp(['cant create field ', tmpfield]);
    end
end
outStruct = rmfield(outStruct, 'data');
outStruct = rmfield(outStruct, 'fields');


    
                            
    



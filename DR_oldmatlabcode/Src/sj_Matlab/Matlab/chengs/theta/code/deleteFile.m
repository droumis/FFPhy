function deleteFile(fname)

if exist(fname, 'file')== 2
    delete(fname);
end

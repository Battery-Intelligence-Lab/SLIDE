function tab = slide_to_table(path)
% @file slide_to_table.m
% @brief binary writing/reading functions, defines `.slide` format
% @author Volkan Kumtepeli
% @date 05 Aug 2024
% @details Uses version 0.

fileID = fopen(path);
readed = uint8(fread(fileID))';

% Read and process metadata
metadata_length = 4;  % Version 0
meta_data = readed(1:metadata_length);

fileVersion = meta_data(1);

if(meta_data(2))
    fileEndianness = 'L';
else
    fileEndianness = 'B';
end

fileDataType = char(meta_data(3));
fileDataSize = meta_data(4); % How many bytes it is.

if(fileDataType == 'f') % floating point type
    if(fileDataSize == 4)
        matlabDataType = 'single'; % single precision
    elseif(fileDataSize == 8)
        matlabDataType = 'double'; % double precision
    else
        error('Float with unknown data size!\n');
    end

elseif(fileDataType == 'u') % unsigned integer
    matlabDataType = ['uint', num2str(8*fileDataSize)];
elseif(fileDataType == 'i') % signed integer
    matlabDataType = ['int', num2str(8*fileDataSize)];
end

[~,~,computerEndianness] = computer; % Get computer specs.

if(fileEndianness ~= computerEndianness)
    error(['File endianness (', fileEndianness, ', ) is different than computer endianness (', computerEndianness, ').' ...
        ' Endianness correction feature is not implemented yet!\n']);
end

fprintf('Now reading %s, with version %d, %s-endian and %s data.\n',path, fileVersion, fileEndianness, matlabDataType);

i_one = find(readed(1+metadata_length:end)==0, 1);
i_one = i_one + metadata_length;
header = char(readed(1+metadata_length:i_one-1));

headers = split(header,',');

Ncol = length(headers);

arr = typecast(readed(i_one+1:end), matlabDataType);
arr = reshape(arr, Ncol, [])';

tab = array2table(arr, 'VariableNames', headers);
end
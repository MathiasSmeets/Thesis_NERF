function data = upload_recording_matrix(data_file, precision, nch)

if nargin <2
    precision = 'int16';
end
if nargin < 3
    nch = [];
end

fid = fopen(data_file);
dd = fread(fid, precision);
if ~isempty(nch)
    data = reshape(dd,[nch size(dd,1)/nch]);
else
    data = dd;
end

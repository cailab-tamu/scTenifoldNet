% Convert a community list to a Pajek partition file
%
% Input
%   - com : community list per vertex
%   - name: file name (without extension)
%   - dir : destination directory
%
% Author: Erwan Le Martelot
% Date: 01/06/11

function [] = com2pajek(com, name, dir)

    dst_file = [name,'.clu'];
    if nargin == 3
        dst_file = [dir, '/', dst_file];
    end

    fid = fopen(dst_file, 'w');
    
    fprintf(fid, '*Vertices %d\n', length(com));
    for i=1:length(com)
        fprintf(fid, ' %d\n', com(i));
    end

    fclose(fid);

end

% InfoMap method from (Roswall and Bergstrom, 2008)
%
% Input
%   - adj: adjacency matrix
%
% Output
%   - com: communities (listed for each node)
%
% Author: Erwan Le Martelot
% Date: 16/06/11

function [com] = infomap(adj)

    % Set the path and command line name
    dir_path = '~/Documents/Dev/infomap_undir/';
    command = 'infomap';

    command = [dir_path,command];
    command = [command,' ',num2str(randi(99999999)),' adj.net',' 1'];

    % Get edges list and output to file
    adj2pajek(adj,'adj','.');
 
    % Call community detection algorithm
    tic;
    system(command);
    toc;
    
    % Load data and create data structure
    fid = fopen('adj.clu','rt');
    com = textscan(fid,'%f','Headerlines',1);
    fclose(fid);
    com = com{1};

    % Delete files
    delete('adj.net');
    delete('adj_map.net');
    delete('adj_map.vec');
    delete('adj.clu');
    delete('adj.map');
    delete('adj.tree');

end

%addFolderPath
%
%  addFolderPath( folder, nMaxSub, nMinSub )
%
%  Adds a folder and possibly all subfolders to the matlab path.  The big
%  advantage of using this over Matlab's menu-based method for adding
%  subfolders is that this function excludes the following subfolders 
%        Those starting with '.' or '@' or "+"
%        Those named: 'Scripts'
%       (In particular if avoids adding .git and its subfolders).
%  folder:  folder to add, default is: '.'
%  nMaxSub: 1 (default) will add 1 level of subfolders only
%           0 no subfolders
%           n will add up to n levels of subfolders
%  nMinSub: 0 (default) will add selected folder up to nMaxSub
%           1 will add only direct subfolders up to nMaxSub
%           n will add only n'th level subfolders to nMaxSub
%
%  Examples:
%    addFolderPath:            to add current folder and direct subfolders
%    addFolderPath('.', 0):    to add just current folder and no subfolders
%    addFolderPath(folder, 1, 1):  to add just direct subfolders of folder
%
%  Note: you'll need to save the path if you want Matlab to remember it for
%  the next session.  
%
%  ==============================
%  The way I set up a new Matlab install is:
%  1. Create home folder, ex:  $home/Documents/MATLAB-<version> and set
%     MATLABPATH environment variable to this.  Restart Matlab.
%  2. Cd to the folder containing addFolderPath.m, and call: addFolderPath
%  3. Next call: addFolderPath(<repo-path>) for each repo to add it
%       -- this will add just 1 level of subfolders, but can add more if needed
%          using: addFolderPath(<repo-path>, n)
%  4. Open Matlab menu: Home > Set Path, and save the pathdef.m file in the
%     home folder specified by MATLABPATH.
%  ==============================
%
%  Daniel Morris, Nov 2016
%
%  See also: userpath
%

%  Copyright (c) 2018, Daniel Morris
%  Licensed under the Simplified BSD License
function addFolderPath( folder, nMaxSub, nMinSub )

if nargin<1, folder = '.'; end
if nargin<2, nMaxSub = 1; end
if nargin<3, nMinSub = 0; end

curdir = pwd;
cd(folder);

if nMaxSub > 0
    files = dir('*');
    for its=1:numel(files)
        if files(its).isdir
            if ~strcmp(files(its).name(1),'.') && ...
                    ~strcmp(files(its).name(1),'@') && ...
                    ~strcmp(files(its).name(1),'+') && ...
                    ~strcmpi(files(its).name, 'scripts')
                addFolderPath( files(its).name, nMaxSub-1, nMinSub-1 );
            end
        end
    end
end

if nMinSub <= 0
    addpath( pwd );
    fprintf('Added: %s\n',pwd);
end

cd(curdir);

end

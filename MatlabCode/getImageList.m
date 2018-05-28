%getImageList
%
%  ilist = getImageList( folderList, presuffix, suffix )
%
%  Get list of images in folder or folders
%  Returns files matching: folder/*<presuffix>.<suffix>
%  folderList: single folder name or else cell array of folder names
%  ilist: cell array of images that match search
%  
%  Daniel Morris, Jan 2018
%
%  See also: 
%

%  Copyright (c) 2018, Daniel Morris
%  Licensed under the Simplified BSD License
function ilist = getImageList( folderList, presuffix, suffix )

if nargin < 3, error('requires 3 arguments'); end

getImages = @(folder, pre, post) dir([folder filesep '*' pre '.' post]);
getImageList = @(folder, pre, post) arrayfun(@(name) [folder filesep name.name], getImages(folder, pre, post), 'uniformoutput',false ); 
getAllImageLists = @(folderlist, pre, post) cellfun( @(folder) getImageList( folder, pre, post), folderlist, 'uniformoutput', false );

if iscell(folderList)
    ilist = getAllImageLists( folderList, presuffix, suffix );
    ilist = vertcat(ilist{:});
elseif ischar( folderList )
    ilist = getImageList( folderList, presuffix, suffix );
else
    error('folderList must be a char or cell array');
end

end



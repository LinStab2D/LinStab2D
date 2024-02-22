function [axs] = plotFlow(X,Y,vars,nr,nc,usedInd,varargin)
    % plotFlow(X,Y,vars,nr,nc,usedInd,varargin)
    %creates a subplot structure to make contour plots of flow quantities.
    % Inputs 
    %   X       : gridpoints x coordiantes
    %   Y       : gridpoints y coordiantes
    %   vars    : a N by 2 cell with each line containing  {F,label}, where
    %   F is a matrix with the same size of X and Y, to be ploted, and
    %   label is the title of the plot
    %   nr,nc   : number of rowns and columns of subplots. nr*nc >=
    %                   size(vars,1)
    %   usedInd (optional) : if the fields in vars do not contain all grid 
    %                       points (e.g., when masks are used), provide 
    %                       usedInd (typically mesh.usedInd to indicate
    %                       the grid points used.
    %   varargin (optional): extra arguments to pcolor, used to create
    %                        plots.
    % Outputs
    %   axs     : list with references to each of the subplots axis.
    
    fprintf( '--- Plotting fields...  '); tic();

    warning('off','MATLAB:contour:ConstantData'); % Disables warning; There are often constant fields to be ploted
    
    nVars   = size(vars,1);
    if exist('nr');  if isempty(nr); nr = ceil(sqrt(nVars)); end; end
    if exist('nc');  if isempty(nc); nc = nr               ; end; end

    if ~( exist('usedInd') && ~isempty(usedInd)) 
        usedInd = 1:numel(X);
    end

    for i=1:nVars
        axs(i) = subplot(nr,nc,i);
        var = nan(size(X));

        var(usedInd) = vars{i,1};

        pcolor(X,Y,var); shading interp
        colorbar;
        title(vars{i,2});
    end
        
    warning('on','MATLAB:contour:ConstantData'); % Re-enables warning
    fprintf( ' Done in %.0f seconds.\n',toc);

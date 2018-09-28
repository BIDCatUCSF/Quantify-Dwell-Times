%
%  Sort Interacting Surfaces Function for Imaris
%
%  Copyright PBeemiller September 2010.
%
%  To make this function available to Imaris, copy this file into the
%  XTensions folder in the Imaris installation folder, e.g., C:\Program
%  Files\Bitplane\Imaris x64 7.1.1\XTensions. After copying this file to
%  the XTensions folder and restarting Imaris, you can find this function
%  in the Surfaces Function menu, as well as in the Image Processing menu
%  under the Surfaces Functions group.
%
%    <CustomTools>
%      <Menu>
%       <Submenu name="Surfaces Functions">
%        <Item name="Quantify Dwell Times" icon="Matlab"
%        tooltip="Quantify the duration of interactions between surfaces.">
%          <Command>MatlabXT::xtsurfaceinteractionduration(%i)</Command>
%        </Item>
%       </Submenu>
%      </Menu> <SurpassTab>
%        <SurpassComponent name="bpSurfaces">
%          <Item name="Quantify Dwell Times" icon="Matlab"
%          tooltip="Quantify the duration of interactions between surfaces.">
%            <Command>MatlabXT::xtsurfaceinteractionduration(%i)</Command>
%          </Item>
%        </SurpassComponent>
%      </SurpassTab>
%    </CustomTools>
% 
%  Description: Meaures the duration of contacts that surfaces make with
%  other surfaces objects. Both surfaces must be tracked.
% 

function xtsurfaceinteractionduration(xObjectID)
    %
    % SURFACEDWELLTIMESXT Measure how long tracked surfaces dwell in
    % interactions with target surfaces.
    %
    % SURFACEDWELLTIMESXT identifies and quantifies interactions that a surface
    % makes with other surface objects. Surface objects are marked as
    % interacting based on whether the vertices of the surfaces come within a
    % specified distance of the vertices in another target Surfaces object. The
    % duration (in frames/time indexes) of all interactions formed by the
    % tracked objects in the Imaris Surfaces object are recorded track-by-track
    % and reported in an Excel spreadsheet.
    %
    %% Connect to the Imaris XT interface.
    % The line below can be used to connect to Imaris when running this
    % function as a script for troubleshooting, etc. Assign xObjectID in the
    % Base workspace or replace with an integer. The connectimaris function
    % needs to be in the MATLAB path.
    % xImarisApp = connectimaris(xObjectID);

    % The following code should be used when running the extension from Imaris.
    if ischar(xObjectID)
        xObjectID = round(str2double(xObjectID));
    end

    % If it's not already added, append ImarisLib to the dynamic Java class path.
    if isempty(regexp(javaclasspath('-dynamic'), 'ImarisLib.jar', 'once'))
        javaaddpath '.\ImarisLib.jar'
    end

    % Create a new object from the ImarisLib class.
    xImarisLib = ImarisLib;

    % Connect to MATLAB to the calling Imaris window with the specified ID.
    xImarisApp = xImarisLib.GetApplication(xObjectID);
    
    %% Prompt the user for the surface to analyze for interactions.
    % Create a list of all surpass children and a list of surfaces. The
    % surfaces list is presented to the user for selection (so that they don't
    % accidentally click a non-surface object). Then we use the surface name to
    % get the child index of the selected surface in the Surpass scene data
    % container.
    ImarisChildren = struct('name', {});
    ImarisSurfaces = struct('name', {});

    for m = 1:xImarisApp.GetSurpassScene.GetNumberOfChildren
        % Get the next child of the Surpass container.
        MthChild = xImarisApp.GetSurpassScene.GetChild(m - 1);

        % If the child is a surface, add it to the children and surfaces list.
        % Otherwise, add it to the children list only.
        if xImarisApp.GetFactory.IsSurfaces(MthChild)
            ImarisSurfaces(length(ImarisSurfaces) + 1).name = char(MthChild.GetName);
            ImarisChildren(length(ImarisChildren) + 1).name = char(MthChild.GetName);
        else
            ImarisChildren(length(ImarisChildren) + 1).name = char(MthChild.GetName);
        end
    end

    % Prompt the user for the Surpass surface they want to sort.
    sortListIdx = listdlg('ListString', {ImarisSurfaces.name}, ...
        'Name', 'Surface Selection', 'SelectionMode', 'single', ...
        'PromptString', 'Select a tracked surface to sort:');

    % Get the child number of the surface.
    sortSurpassIdx = find(strcmp({ImarisChildren.name}, ...
        ImarisSurfaces(sortListIdx).name));

    % Get the surface object.
    SortSurface = xImarisApp.GetFactory.ToSurfaces(...
        xImarisApp.GetSurpassScene.GetChild(sortSurpassIdx - 1));
    
    %% Prompt the user for the interacting surface.
    % Prompt the user for the interacting Surpass surface. If the user does not
    % want to sort this surface, it does not need to be sorted.
    partnerListIdx = listdlg('ListString', {ImarisSurfaces.name}, ...
        'Name', 'Surface Selection', 'SelectionMode', 'single', ...
        'PromptString', 'Select a target surface:', 'InitialValue', 1);

    % Get the child number of the target surfaces.
    partnerSurpassIdx = find(strcmp({ImarisChildren.name}, ...
        ImarisSurfaces(partnerListIdx).name));

    % Get the surface object.
    PartnerSurface = xImarisApp.GetFactory.ToSurfaces(...
        xImarisApp.GetSurpassScene.GetChild(partnerSurpassIdx - 1));
    
    %% Get the data set geometry.
    xMin = xImarisApp.GetDataSet.GetExtendMinX;
    zMin = xImarisApp.GetDataSet.GetExtendMinZ;

    xMax = xImarisApp.GetDataSet.GetExtendMaxX;
    zMax = xImarisApp.GetDataSet.GetExtendMaxZ;

    xSize = xImarisApp.GetDataSet.GetSizeX;
    zSize = xImarisApp.GetDataSet.GetSizeZ;
    
    %% Calculate the lateral and axial sampling resolution.
    rUnit = (xMax - xMin)/xSize;
    zUnit = (zMax - zMin)/zSize;
    % Use the sampling resolution to calculate a default gap cutoff to feed
    % into inputdlg.
    gapDefault = num2str(min([rUnit, zUnit]));
    
    %% Ask for the contact distance cutoff.
    Options.Resize='on';
    Options.WindowStyle='normal';
    Options.Interpreter='tex';
    gapParameter = inputdlg({'Gap cutoff distance (\mum):'}, ...
        'Contact Gap', [1 32], {gapDefault}, Options);
    gapDistance = str2double(gapParameter{1});
    
    %% Start a waitbar to track the processing.
    progress = waitbar(0, 'Getting surfaces to analyze', ...
        'Name', 'Quantify Surface Interactions');
    
    %% Get the tracked surfaces we want to analyze.
    % Get the tracked surfaces object information.
    sortSurfaceCount = SortSurface.GetNumberOfSurfaces;
    sortEdges = SortSurface.GetTrackEdges;
    sortTrackIds = SortSurface.GetTrackIds;

    % Pre-allocate a structure to hold the tracked surface data.
    SortStruct(sortSurfaceCount) = struct('Time', [], ...
        'Vertices', [], 'TrackId', []);

    % Now get the track data.
    for s = 1:sortSurfaceCount
        % Update the waitbar.
        waitbar(s/sortSurfaceCount, progress)

        % Find the current surface's entries in the edge list.
        [edgesRowID, ~] = find(sortEdges == s - 1, 1, 'first');

        % Get the parent (track) IDs from the first edge entry (if they exist,
        % 2nd entries inthe edges array will be identical).
        SortStruct(s).TrackId = sortTrackIds(edgesRowID);

        % Get the track time indexes for easy searches later.
        SortStruct(s).Time = SortSurface.GetTimeIndex(s - 1);
        SortStruct(s).Vertices = SortSurface.GetVertices(s - 1);
    end
    
    %% Get the binding partner surfaces.
    % Update the waitbar.
    waitbar(0, progress, 'Getting interaction partner surfaces')

    % Get the binding partner surfaces. If we are searching the surfaces for
    % homotypic interactions, we can just alias the surface to save the memory
    % (we don't modify the surface data structure(s), just read from them).
    if sortSurpassIdx == partnerSurpassIdx
        PartnerStruct = SortStruct;
    else
        % Get the tracked surfaces object information.
        partnerSurfaceCount = PartnerSurface.GetNumberOfSurfaces;
        partnerTrackIds = PartnerSurface.GetTrackIds;
        partnerEdges = PartnerSurface.GetTrackEdges;

        % Pre-allocate a structure to hold the tracked surface data.
        PartnerStruct(partnerSurfaceCount) = struct('Time', [], ...
            'Vertices', [], 'TrackId', []);

        % Now get the track data.
        for s = 1:partnerSurfaceCount
            % Update the waitbar.
            waitbar(s/partnerSurfaceCount, progress)

            % Find the current surface's entries in the edge list.
            [edgesRowID, ~] = find(partnerEdges == s - 1, 1, 'first');

            % Get the parent (track) IDs from the first edge entry.
            PartnerStruct(s).TrackId = partnerTrackIds(edgesRowID);

            % Get the track indexes.
            PartnerStruct(s).Time = PartnerSurface.GetTimeIndex(s - 1);
            PartnerStruct(s).Vertices = PartnerSurface.GetVertices(s - 1);
        end
    end
    
    %% Create a cell for recording all the interaction data.
    % Count the number of tracks to analyze.
    sortTrackCount = length(unique(sortTrackIds));

    % Count the number of parnter tracks.
    if sortSurpassIdx == partnerSurpassIdx
        partnerTrackCount = sortTrackCount;
    else
        partnerTrackCount = length(unique(partnerTrackIds));
    end

    interactionRecord = zeros(sortTrackCount, partnerTrackCount, ...
        xImarisApp.GetDataSet.GetSizeT, 'uint16');
    
    %% For each track, search for encounters with the target surfaces.
    % For each surface object in each track (each time point in the track), we
    % want to keep track of whether it contacts a target surface. After scoring
    % for contacts with the target surface, we determine whether the track
    % meets the criteria to be an interacting track.
    for s = 1:length(SortStruct)
        % Update the progress bar.
        waitbar(s/length(SortStruct), progress, 'Searching for interactions')

        % Convert the current surface's Imaris parent (track ID) to a value
        % that we can use to index into the interactionRecord array.
        sTrack = SortStruct(s).TrackId - 10^9 + 1;

        % Find target surfaces present at the current time point.
        targetIdxs = find([PartnerStruct(:).Time] == SortStruct(s).Time);

        if sortSurpassIdx == partnerSurpassIdx
            selfSurfaceIdx = find([SortStruct(:).TrackId] == SortStruct(s).TrackId);
            targetSelfIdx = ismember(targetIdxs, selfSurfaceIdx);
            targetIdxs(targetSelfIdx) = [];
        end

        % Get the surface's Delasunay triangulation and determine the bounding box.
        SortDelaunay = DelaunayTri(...
            double(unique(SortStruct(s).Vertices, 'rows')));
        bBoxMin = min(SortDelaunay.X, [], 1) - (gapDistance + eps('single'));
        bBoxMax = max(SortDelaunay.X, [], 1) + (gapDistance + eps('single'));

        % At each time point, find objects that come within the cutoff
        % distance.
        for t = 1:length(targetIdxs)        
            %% Check the nearby targets for encounters with the surface.
            % Get the potential binding partner's vertices.
            partnerVertices = double(PartnerStruct(targetIdxs(t)).Vertices);

            % Find the vertices greater than the bounding box min.
            gtCoords = bsxfun(@ge, partnerVertices, bBoxMin);
            gtVertices = sum(gtCoords, 2) == 3;

            % Find the vertices less than the bounding box max.
            ltCoords = bsxfun(@le, partnerVertices, bBoxMax);
            ltVertices = sum(ltCoords, 2) == 3;

            % Find the vertices that fall within the bounding box.
            bBoxVertices = gtVertices & ltVertices;

            % Mask the vertices to include only those that are within
            % the bounding box.
            testVertices = partnerVertices(bBoxVertices, :);

            % If there are vertices within the bounding box, we search for
            % encounters.
            if ~isempty(testVertices)
                % Check for target vertices residing within the
                % surface.
                enclosingSimplex = pointLocation(SortDelaunay, testVertices);
                enclosedVertexIdxs = ~isnan(enclosingSimplex);
                enclosedVerticesCount = sum(enclosedVertexIdxs);

                % Calculate the distance to the nearest point in the
                % surface to sort.
                [~, nearDistances] = nearestNeighbor(SortDelaunay, ...
                    testVertices(~enclosedVertexIdxs, :));

                % Determine if any points get within the cutoff
                % distance.
                partnerFound = (sum(nearDistances <= gapDistance) + ...
                    enclosedVerticesCount) > 0;

                % If the current target is close enough that we consider it
                % interacting, record it's ID in the interaction record.
                if partnerFound
                    % Get the current target's track parent ID.
                    tTrack = PartnerStruct(targetIdxs(t)).TrackId - 10^9 + 1;

                    % Record the interaction in the current surfaces
                    interactionRecord(sTrack, tTrack, ...
                        PartnerStruct(targetIdxs(t)).Time + 1) = 1;
                end
            end
        end
    end
    
    %% Calculate the duration the interactions for each track.
    dwellCell = cell(sortTrackCount, 1);
    for r = 1:size(interactionRecord, 1)
        % Update the progress bar.
        waitbar(r/size(interactionRecord, 1), progress, 'Measuring dwell times')

        % Extract the track's interactions from the record.
        trackInteractions = squeeze(interactionRecord(r, :, :));

        % Figure out which targets the track interacts with.
        interactionMask = sum(trackInteractions, 2);% >= 1;
        partnerIds = transpose(find(interactionMask));

        % Measure the dwell times of the interactions.
        for d = partnerIds
            encounterMask = transpose(trackInteractions(d, :));

            % Pad a vector the size of the encounter mask.
            encounterRuns = [0; single(encounterMask); 0];

            % Find the 0 (non-interacting/free) indices.
            freeIdxs = find(~encounterRuns);

            % Find the length of the gaps between free periods.
            runBreaks = 1 - diff(freeIdxs);

            % Insert the free gap durations into the padded encounter mask.
            encounterRuns(freeIdxs) = [runBreaks; 0];

            % Find the indices representing the start of the interactions.
            dwellStarts = encounterRuns < 0;

            % Count how much longer the interactions last at each time point.
            encounterRuns = cumsum(-encounterRuns(1:end - 1));

            dwellTimes = transpose(encounterRuns(dwellStarts));

            % Record the dwell times.
            dwellCell{r} = [dwellCell{r}, dwellTimes];
        end
    end
    %% Reorganize the data into Excel friendly form and write it to file.
    trackDwellCounts = cellfun(@length, dwellCell);
    dwellArray = nan(sortTrackCount, max(trackDwellCounts));
    for r = 1:sortTrackCount
        dwellArray(r, 1:trackDwellCounts(r)) = dwellCell{r};
    end

    % Construct the Excel file name.
    imsName = char(xImarisApp.GetCurrentFileName);
    imsBase = strrep(imsName, '.ims', '');
    excelName = [imsBase ' - ' char(SortSurface.GetName) ' Dwell Times on ' ...
        char(PartnerSurface.GetName) '.xlsx'];
    
    % It’s the nominal duration (number of time points) 
	%that an interaction lasted. So if a row has values of 1, 3, and 1 
	%in the columns, the cell formed interactions that lasted 1 frame, 
	%3 frames and 1 frame.	
	col_header = {'Track ID', 'Separate Continuous Interactions (# of timepoints)'};
    % Write the dwell time data.
	xlswrite(excelName, col_header,1 , 'A2')
    xlswrite(excelName, unique(sortTrackIds), 1, 'B2')
    if isempty(dwellArray)
        msgbox('No interactions found.', 'Quantify Dwell Times')
    else
        xlswrite(excelName, dwellArray, 1, 'C2')
    end
    %%
    close(progress)
end    

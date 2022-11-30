%Wind tunnel analysis script
%Reads in video file, corresponding cropping function, and height data from
%wind tunnel csv.
%this script integrates the fitted function to find the length of the half
%height point, which should be more accurate.

%Loop through directory, get list of video files and cropping functions
videos = dir('*.mp4');
crops = dir('*_hsv.m');
videonames = {videos.name};
cropnames = {crops.name};
%nfiles = length(videonames);
nfiles = size(videonames,2);

%read in wind tunnel csv
pheno = csvread('windtunnel.csv');


%waitbar loop through files
h = waitbar(0,'Please wait...');

for n = 1:nfiles
    %get nth video, cropping function in directory
    currentvideo = cell2mat(videonames(n));
    currentcrop = cell2mat(cropnames(n));
    video = VideoReader(currentvideo);
    mask = erase(currentcrop,'.m');
    
    %get half height value from pheno
    halfheight = pheno(n,2) / 2;
    %get quarter height value from pheno
    quarterheight = pheno(n,2) / 4;
    
    % bend analysis portion
    %get number of frames 
    nframes = video.NumberOfFrames;
    
    %preallocate coefficient matrix
    coef = zeros(nframes,2);
    centroidOut = zeros(nframes,2);
    %preallocate angle, displacement array
    %fields: frame, angle, heightathalfheightx,
    %xdisphalfheight,xdispquarterheight,xhalfheight,yhalfheight,curvelength,xhalfheightcalc,yhalfheightcalc,projectedarea
    angles = zeros(nframes,11);
    %pre allocate the initial area value
    initialArea = [];

    %loop through frames
    for k = 1:nframes
        %read in image
        image = read(video,k);
        resolution = 7.3; %pixels per cm
        %mask image using video specific masking function
        maskedimage = feval(mask,image); 
        %change everything else in background to 0
        maskedimage(:,1:550) = 0; %lower columns
        maskedimage(:,1475:1920) = 0; %higher columns
        maskedimage(1:230,:) = 0; %lower rows
        maskedimage(890:1080,:) = 0; %higher rows

        %get the unique indices for nonzero row, col indices in the masked image
        [maskRowIndices, maskColIndices] = find(maskedimage ~= 0); %get non zero row and col indices
        maskRowIndices = unique(maskRowIndices); %fewer rows (vertical part of plant)
        maskColIndices = unique(maskColIndices); %more cols (windward edge of plant)

        %get the max nonzero row index for each nonzero col index. 
        nmaskColInds = size(maskColIndices,1);
        maxRowIndices = zeros(size(maskColIndices,1),1);
        for i = 1:nmaskColInds
            nonZeroRows = find(maskedimage(:,maskColIndices(i)) ~= 0 ); %get all non zero indices within each column
            maxRowIndices(i) = max(nonZeroRows); %get the maximum (windward) non zero row index
            %maxRowIndices(i) = max(find(maskedimage(:,maskColIndices(i)) ~= 0 ));
        end 

        %get points to represent plant in wind tunnel
        %maxRowIndices and maskColIndices have the same dimensions. 
        %scatter(maskColIndices,maxRowIndices)
        %combine the col and row index arrays. 
        combined = [maxRowIndices,maskColIndices]; 
        %scatter((combined(:,1).*-1),(combined(:,2).*-1))
        %rotate points in combined array 180 degrees
        rotatedCombined = combined.*-1;
        %divide points by resolution
        rotatedCombined = rotatedCombined./resolution; 
        %scatter(rotatedCombined(:,1),rotatedCombined(:,2))

        %model the bending.  Use a power law function f(x) = a*x^b. First find
        %the min x and y values as the origin, then subtract these from the
        %negative values in the rotated combined matrices.  add 0.1 to each
        %array to avoid problems fitting the model at 0.  
        minx = min(rotatedCombined(:,1));
        miny = min(rotatedCombined(:,2)); 
        X = rotatedCombined(:,1) - (minx - 0.1);
        Y = rotatedCombined(:,2) - (miny - 0.1); 
        %power law fit
        f = fit(X,Y,'power1');

        %calculating parameters of plant bending
        %write the frame number
        angles(k,1) = k;
        %find where the bending is at half height 
        x = (cosd(50)*halfheight);
        Ymod = f(x); %get the modeled value of Y
        %calculate angle with respect to plant base
        Angmod = atand(Ymod/x); 
        angles(k,2) = Angmod;

        %Write the y value  for each frame in
        %field 3 of angles array. This is the y value at the fixed x point
        %where half height as hypoteneuse forms a 50 degree angle
        angles(k,3) = f(x);

        %write the x displacement (for a consistent value of y at halfheight, quarterheight)
        %fit a new function to the opposite dataset, where plant height is on
        %the x axis
        f2 = fit(Y,X,'power1');
        xdisp = f2(halfheight);
        angles(k,4) = xdisp;
        xdisp2 = f2(quarterheight);
        angles(k,5) = xdisp2; 

        %Get the x and y values for the halfheight point at the given Angmod
        xhalfheight = (sind(180-90-Angmod)*halfheight);
        yhalfheight = f(xhalfheight);
        angles(k,6) = xhalfheight;
        angles(k,7) = yhalfheight;
        
        %Calculate Centroid on masked image(with respect to base of plant at reference
        %deformation, then rotate coordinates
        [yC, xC] = ndgrid(1:size(maskedimage, 1), 1:size(maskedimage, 2));
        centroid = mean([xC(logical(maskedimage)), yC(logical(maskedimage))]);
        centroid = centroid.*-1; %rotate centroid 180 degrees
        centroid = centroid./resolution;  %divide points by resolution
        y_G = centroid(1)-miny; %correct centroid Y coordinate 
        
        


        %calculus method for obtaning x and y values for the halfheight point
        %based on fit curve length
        %Get the x and y values for the halfheight point. Halfheight point found by integrating to find the length of the curve
        %1st derivative of power function: d/dx = b*a*x.^b-1 
        %obtain xval range of points within the function for integration
        minx2 = min(X);
        maxx = max(X);
        %get coeffvalues from the fit
        coeff = coeffvalues(f);
        %define the integrand as an anonymous function of the 1st
        %derivative (integration to find arc length matlab documentation).
        %Takes form integral(sqrt(1+(dy/dx).^2))
        flen = @(x) sqrt(1+(coeff(2)*coeff(1)*x.^(coeff(2)-1)).^2);
        %get total length along curve fit by integrating this function from
        %minx to maxx
        len = integral(flen,minx2,maxx);
        angles(k,8) = len;
        %Now create the zeroing function.  The function glen will be zero
        %when the length of curve (estimated by flen between minx and some upper limit x) is at halfheight
        glen = @(x) ((integral(flen,minx2,x)) - halfheight); 
        %find the value of x (xhalfheight point) for which glen function is 0.
        %There will be warnings here, since some of these initial values will
        %produce an imaginary answer
        warning('off','all')

        if k <= 1300 
            xhalfheightcalc = fzero(glen,10);
        else 
            xhalfheightcalc = fzero(glen,100);
        end 

        %get the y half height by inputting this value into the function f
        yhalfheightcalc = f(xhalfheightcalc); 
        %write half height coordinates to angles array
        angles(k,9) = xhalfheightcalc;
        angles(k,10) = yhalfheightcalc; 
                       
        %get image area. Store plant area for only the first frame, then
        %calculate projection for each subsequent frame
        if k == 1
            nonzeroPixels = nnz(maskedimage); %number of nonzero pixels
            initialArea = nonzeroPixels / (resolution^2); %in cm2
            angles(k,7) = initialArea;
        else
            %area in frame is a projection of the initial area dependent on
            %the stem angle
            projectedArea = sind(Angmod)*initialArea;
            angles(k,11) = projectedArea; 
        end 


        %write coef values to matrix
        coef(k,:) = coeffvalues(f);
        
        %write y_G centroid y, x S.2 coordinate of halfway point to matrix
        centroidOut(k,:) = [y_G x];
        
    end 
    
    %write coeff, angle arrays to files
    potnum = erase(currentvideo,'.mp4');
    coeffname = char(sprintf('fit%s.txt',potnum));% change file name from coeff to fit (somehow R having issues in reading files that start with same letters)
    anglename = char(sprintf('bend%s.txt',potnum));
    centroidname = char(sprintf('centroid%s.txt',potnum));

    %write coef array to text
    dlmwrite(coeffname,coef);
    %write angles array to text
    dlmwrite(anglename,angles);
    %write centroid and x coords to text
    dlmwrite(centroidname,centroidOut);
    
    %waitbar
    waitbar(n / nfiles)
    
   
    
end 

close(h) 
%}
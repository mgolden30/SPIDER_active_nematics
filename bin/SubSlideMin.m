function  [Out] = SubSlideMin(ImgIn,WindowSize,Plot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Image processing - subtract sliding min (normilization) %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function subrtacts the sliding min in a window of size WindowSize
% [Out] = SubSlideMin(ImgIn,WindowSize,Plot)
%
%  ImgIn      -- Should be the image (not the file name)
%
%  WindowSize -- Size of the window the local min is calculated over
%                If no input is given a default window of 3x3 is used
%                The WindowSize should be roughly the particle diameter
%
%  Plot       -- (optional) 0 to not show results, 1 to show results
%                If no input is given no plot is shown
%
%
% Dependencies:
%               none
% 
% This code was originaly writen by Chris Crowley on 10/29/0218
% 
% Code was updated with image padding so the boarders are corectly handled
% by Chris Crowley on 7/6/2020
% 

%% Taking care of a few things first
    if ~exist('WindowSize','var')
         % If no input was given for WindowSize use the default value
          WindowSize = 3;
    end
    if ~exist('Plot','var')
         % If no input was given for Plot, don't plot
          Plot = 0;
    end
    if mod(WindowSize,2)==0
        WindowSize = WindowSize - 1;
        disp('...............................................................................')
        disp('This function requires an odd window size..')
        disp('    By choosing a window size that is even you are not symetricaly normalizing.')
        disp('    The window size has ben reduced by 1 to make it odd for this calculation.  ')
        disp('...............................................................................')
    end
%% Actually doing the sliding average
% this slides the window to the right then down


% Pad the image with NaN around the image. This is done so you can subtract
% the min at the edge of the image correctly
ImgIn_pad = double(padarray(ImgIn+1,[(WindowSize-1)/2 (WindowSize-1)/2],0,'both')); % the +1 is so there are no zero pixels
ImgIn_pad(ImgIn_pad==0)=nan; % since the lina above has a +1, the only zeros are from the pad
ImgIn_pad = ImgIn_pad-1; % remove the +1 to make it the original image intensity.

% calculate the min
[m,n]=size(ImgIn);
for j = (WindowSize-1)/2+1:(WindowSize-1)/2+n
    for i = (WindowSize-1)/2+1:(WindowSize-1)/2+m
        Window = ImgIn_pad(i-(WindowSize-1)/2:i+(WindowSize-1)/2,j-(WindowSize-1)/2:j+(WindowSize-1)/2);
        B(i-(WindowSize-1)/2,j-(WindowSize-1)/2) = uint16(min(min(Window)));
    end
end

%% Ploting 
    if Plot == 1
        subplottight(1,2,1)
        imagesc(ImgIn)
        colormap('gray')
        axis off;
        title('Original image')
        subplottight(1,2,2)
        imagesc(ImgIn-B)
        axis off;
        title('Output image')
    end

    Out = ImgIn-B;
end

%% Ploting function 
% this is to plot the images without huge ass margins
function h = subplottight(n,m,i)
    [c,r] = ind2sub([m n], i);
    ax = subplot('Position', [(c-1)/m, 1-(r)/n, 1/m, 0.95/n]);
    if(nargout > 0)
      h = ax;
    end
end

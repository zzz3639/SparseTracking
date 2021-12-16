function [trackmerge] = trackunion( track1, track2 )
% Merge 2 track sets into one
%   track1 and track2 should have same number of frames
%   trackset.no is T*1 array, T is number of frames
%   trackset.track is a matrix, [x,y,I,T,id], 2 position values, intensities, frame, molecule id

T=length(track1.no);
trackmerge.no = zeros(T,1);
trackmerge.no = track1.no + track2.no;
idmax = max(track1.track(:,5));
idmintrack2 = min(track2.track(:,5));
if idmintrack2==0
    track2.track(:,5) = track2.track(:,5) + idmax + 1;
else
    track2.track(:,5) = track2.track(:,5) + idmax;
end

trackmerge.track = [track1.track;track2.track];

end

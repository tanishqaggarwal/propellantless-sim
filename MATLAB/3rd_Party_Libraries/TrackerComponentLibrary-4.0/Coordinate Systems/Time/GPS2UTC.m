function [Jul1,Jul2]=GPS2UTC(Jul1,Jul2)
%%UTC2TCG  Convert from the timescale used by the Global Positioning System
%          (GPS) given as a two-part Julian date to universal coordinated
%          time (UTC), represented as a two-part pseudo Julian date.
%
%INPUTS:    Jul1, Jul2  Two parts of a Julian date given in GPS time.
%                       The units of the date are days. The full date is
%                       the sum of both terms. The date is broken into two
%                       parts to provide more bits of precision. It does
%                       not matter how the date is split.
%
%OUTPUTS:   Jul1, Jul2  The time as a pseudo-Julian date in UTC.
%
%The UTC date is only pseudo-Julian, because there is not a fixed number
%of seconds in a Julian day. The convention used in the IAU standard is
%that the Julian day matches the UTC day regardless of whether the UTC day
%is 86399, 86400 or 86401 SI seconds (depending on the presence of leap
%seconds).
%
%UTC began at 1960 January 1.0 (JD 2436934.5) and this function should not
%be called with an earlier date.
%
%This function just calls GPS2TAI and TAI2UTC.
%
%Many temporal coordinate systems standards are compared in [1].
%
%REFERENCES:
%[1] D. F. Crouse, "An Overview of Major Terrestrial, Celestial, and
%    Temporal Coordinate Systems for Target Tracking," Formal Report, Naval
%    Research Laboratory, no. NRL/FR/5344--16-10,279, 10 Aug. 2016, 173
%    pages.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

[Jul1,Jul2]=GPS2TAI(Jul1,Jul2);
[Jul1,Jul2]=TAI2UTC(Jul1,Jul2);
end

%LICENSE:
%
%The source code is in the public domain and not licensed or under
%copyright. The information and software may be used freely by the public.
%As required by 17 U.S.C. 403, third parties producing copyrighted works
%consisting predominantly of the material produced by U.S. government
%agencies must provide notice with such work(s) identifying the U.S.
%Government material incorporated and stating that such material is not
%subject to copyright protection.
%
%Derived works shall not identify themselves in a manner that implies an
%endorsement by or an affiliation with the Naval Research Laboratory.
%
%RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
%SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY THE NAVAL
%RESEARCH LABORATORY FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS
%OF RECIPIENT IN THE USE OF THE SOFTWARE.

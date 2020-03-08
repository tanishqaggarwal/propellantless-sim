function [MSLPoints,coeffData]=ellips2MSLHelmert(points,useNGAApprox,modelType,coeffData)
%%ELLIPS2MSLHELMERT  Given points in WGS-84 ellipsoidal coordinates,
%                    convert the ellipsoidal height components of the
%                    points into heights above mean-sea level (MSL)
%                    using Helmert's projection method with the EGM2008
%                    gravitational (and terrain correction) models. This
%                    differs from the more precise but seldom-used
%                    Pizzetti's projection method in that the curvature of
%                    the plumb line is completely ignored.
%
%INPUTS: points One or more points given in WGS-84 geodetic latitude and
%               longitude, in radians, and height, in meters for which the
%               corresponding MSL heights are desired. To convert N points,
%               points is a 3XN matrix with each column having the format
%               [latitude;longitude; height].
%  useNGAApprox If one wishes the intermediate tide-free geoid computation
%               to match the results of the National Geospatial
%               Intelligence Agency's code to three digits after the
%               decimal point, then this should be true. Setting this to
%               false might produce results that are marginally more
%               accurate. The default is false if this parameter is
%               omitted.
%     modelType An optional parameter specifying coefficient model to
%               load if coeffData is not provided. Possible values are
%               0 (The default if omitted) Load the EGM2008 model.
%               1 Load the EGM96 model.
%     coeffData A set of pre-loaded coefficients that can speed up the
%               computation by eliminating the need to compute them on the
%               fly. coeffData is as defined for the coeffData return value
%               of getEGMGeoidHeight.
%
%OUTPUTS: MSLPoints A 3XN array of the converted points, where each vector
%                   contains [latitude;longitude;MSLHeight]; The latitude
%                   and longitude are unchanged from the input data, but
%                   the heights have been converted from WGS-84 ellipsoidal
%                   heights to MSL heights.
%         coeffData The coeffData coefficients that can be passed to
%                   another call of ellips2MSLHelmert or MSL2EllipseHelmert
%                   to make it faster.
%
%As described in Chapter 5.5 of [1], the MSL height using Helmert's
%projection is just the ellipsoidal height minus the geoid height. This
%function calls getEGMGeoidHeight to get the geoidal height and
%subtracts it off.
%
%Helmert's projection calculated the height using a line from the point to
%the reference ellipsoid, subtracting off the distance from the point to
%the geoid along the line. Pizzetti's projection follows the curved plumb
%line down to the geoid. Thus, the latitude and longitude of the projected
%point on the geoid are not the same as that of the point being projected.
%The MSL height obtained is also slightly different. However, the
%difference is small and Pizzetti's projection is difficult to use, so
%Helmert's projection is almost always used.
%
%The inverse of this function is MSL2EllipseHelmert.
%
%REFERENCES:
%[1] B. Hofmann-Wellenhof and H. Moritz, Physical Geodesy, 2nd ed. 
%    SpringerWienNewYork, 2006.
%
%January 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(useNGAApprox))
    useNGAApprox=false;
end

if(nargin<3||isempty(modelType))
    modelType=0; 
end

if(nargin>3)
    [geoidHeight,coeffData]=getEGMGeoidHeight(points(1:2,:),1,useNGAApprox,modelType,coeffData);
else
    [geoidHeight,coeffData]=getEGMGeoidHeight(points(1:2,:),1,useNGAApprox,modelType);
end

MSLHeight=points(3,:)-geoidHeight';
numPoints=size(points,2);

MSLPoints=zeros(3,numPoints);
MSLPoints(1:2,:)=points(1:2,:);
MSLPoints(3,:)=MSLHeight;
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

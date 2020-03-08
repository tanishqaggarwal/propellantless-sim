function sgn=exactSignOfSum(S)
%%EXACTSIGNOFSUM Compute the exact sign of the sum of a number of finite
%                numbers. The sign computation occurs exactly, regardless
%                of whether the actual sum value would have overflowed or
%                suffered an accumulation of finite precision errors
%                whereby just adding them would have given the wrong
%                result.
%
%INPUTS: S An NX1 or 1XN vector whereby one desires the sign of the sum of
%          the elements.
%
%OUTPUTS: sgn This is 1 if the exact sum of the elements of S is
%             positive, 0 if it is zero and -1 if it is negative.
%
%The algorithm is taken from [1], where code is provided in an appendix.
%The code with minor changes and corrections is also available from Jon
%Rokne's web site at
%http://pages.cpsc.ucalgary.ca/~rokne/
%The implementation here uses the corrections and uses the sort algorithm
%in C++ standard template library rather than the sort algorithm provided
%by Rokne.
%
%The algorithm can be compiled for use in Matlab  using the 
%CompileCLibraries function.
%
%The algorithm is run in Matlab using the command format
%sgn=exactSignOfSum(S);
%
%REFERENCES:
%[1] H. Ratschek and J. Rokne, "Exact computation of the sign of a finite
%    sum," Applied Mathematics and Computation, vol. 99, no. 2-3, pp. 99-
%    127, 15 Mar. 1999.
%
%December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

error('This function is only implemented as a mexed C or C++ function. Please run CompileCLibraries.m to compile the function for use.')

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

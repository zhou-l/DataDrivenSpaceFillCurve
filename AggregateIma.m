function ImaOUT = AggregateIma(ImaIN,n,Fct)
% AggregateIma resize an image to a smaller one using 6 possible function
% (Fct). 
% in:
% - ImaIN: image N dimension
% - n: factor for aggregation. 
% - Fct: type of aggregation (1=mean (default); 2=mode; 3=std; 4=min;
% 5=max; 6=median)
% out = ImaOUT: aggregated image.
% ex: if ImaIN is a 100 by 100 image and n=4, ImaOUT will be a 25*25 image
% 
% Notice that it works with N dimension image and an extra dimension is
% temporally added to compute the aggregation. there is no pixel by pixel
% (time consuming) computation.
% Author: Martin Claverie, Univesity of Maryland
% date: November 2012
d=[];
ClassT = class(ImaIN);
DimIma = ndims(ImaIN);
for i=1:n
  for j=1:n
      if DimIma == 2
          d=cat(DimIma+1,d,ImaIN(i:n:end-(n-i),j:n:end-(n-j),:));
      else
          for k=1:n
             d=cat(DimIma+1,d,ImaIN(i:n:end-(n-i),j:n:end-(n-j),k:n:end-(n-k)));
          end
      end
  end
end
if nargin<3
  Fct=1;
end
switch Fct
  case 1
    eval(['ImaOUT =' ClassT '(mean(double(d),DimIma+1)) ;']); 
  case 2
    eval(['ImaOUT =' ClassT '(mode(double(d),DimIma+1)) ;']); 
  case 3
    eval(['ImaOUT =' ClassT '(std(double(d),[],DimIma+1)) ;']);
  case 4
    eval(['ImaOUT =' ClassT '(min((d),[],DimIma+1)) ;']); 
  case 5
    eval(['ImaOUT =' ClassT '(max((d),[],DimIma+1)) ;']); 
  case 6
    eval(['ImaOUT =' ClassT '(median((d),DimIma+1)) ;']); 
end

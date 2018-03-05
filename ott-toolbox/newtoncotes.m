function output=newtoncotes(varargin);
%this is a modified version fo the newton cotes formula for simpson's 
%1/3 and 3/8 rules. It requires at leas three points in the vertical
%direction of the matrix to work.

lengthx=length(varargin{1})-1;

num38ths=fix(lengthx/3);
remlen=lengthx-num38ths*3;

if remlen==1
    num38ths=num38ths-1;
    remlen=remlen+3;
end

vector13rds=0;
if remlen==4
    vector13rds=[1,3+num38ths*3];
end
if remlen==2
    vector13rds=num38ths*3+1;
end

vector38ths=(2*length(vector13rds)-1)+[1:3:num38ths*3]-1;


if length(varargin)==1
    output=3*sum(varargin{1}(vector38ths,:)+3*varargin{1}(vector38ths+1,:)+3*varargin{1}(vector38ths+2,:)+varargin{1}(vector38ths+3,:))*.125;
    if vector13rds(1)~=0
        output=output+sum(varargin{1}(vector13rds,:)+4*varargin{1}(vector13rds+1,:)+varargin{1}(vector13rds+2,:))*3.333333333333333e-1;
    end
else
    x1=repmat(varargin{1}(vector38ths+3)-varargin{1}(vector38ths),[1,size(varargin{2},2)]);
    
    output=sum(x1.*(varargin{2}(vector38ths,:)+3*varargin{2}(vector38ths+1,:)+3*varargin{2}(vector38ths+2,:)+varargin{2}(vector38ths+3,:)))*.125;
    if vector13rds(1)~=0
        x2=repmat(varargin{1}(vector13rds+2)-varargin{1}(vector13rds),[1,size(varargin{2},2)]);
        output=output+sum(x2.*(varargin{2}(vector13rds,:)+4*varargin{2}(vector13rds+1,:)+varargin{2}(vector13rds+2,:)))*1.6666666666666666e-1;
    end
end

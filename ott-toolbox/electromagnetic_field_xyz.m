function [structureoutput]=electromagnetic_field_xyz(xyz,nm,ab,pq,cd,relindex,tol)
% electromagnetic_field_xyz.m : Calculates the fields of any of the
%                                incident, scattered or internal
%                                beam shape coefficients.
%
% Usage:
%
% [structureoutput]=electromagnetic_field_xyz(xyz,[n;m],[a;b],[p;q],[c;d],relindex)
% or
% [structureoutput]=electromagnetic_field_xyz(xyz,[n;m],[a;b],[],[])
% or
% [structureoutput]=electromagnetic_field_xyz(xyz,[n;m],[a;b],[p;q])
% or
% [structureoutput]=electromagnetic_field_xyz(xyz,[n;m],[],[p;q])
%
% the output structure contains three set of vectors of the following:
% Eincident, Hincident.
% Escattered, Hscattered.
% Einternal, Hinternal.
% depending on the context of the input. e.g. if only ab has non-zero length
% then only the Eincident, Hincident fields will appear in the output.
%
% ab, pq, cd are the full or sparse column vectors of the te/tm modes.
%
% n,m are the mode indices, these can be in truncated form. The calculation
% will be quicker if a truncated n and m can be used.
%
% NOTE: If internal fields are calculated only the theta and phi components
% of E are continuous at the boundary. Conversely, only the r component of
% D is continuous at the boundary.
%
% PACKAGE INFO

verbose=0;
if ~exist('relindex','var')
    relindex=1;
end

lengthnm=length(nm)/2;

n=nm(1:length(nm)/2);
m=nm(length(nm)/2+1:end);

ci=combined_index(n,m);

try
    lengthab=length(ab)/2;
    ab=full(ab([ci;ci+lengthab]));
catch
    ab=[];
end

try
    lengthpq=length(pq)/2;
    pq=full(pq([ci;ci+lengthpq]));
catch
    pq=[];
end

try
    lengthcd=length(cd)/2;
    cd=full(cd([ci;ci+lengthcd]));
catch
    cd=[];
end

if ~exist('tol','var')
    tol=1e-15;
end
%calculate the space
lengthab=length(ab)/2;
lengthcd=length(cd)/2;
lengthpq=length(pq)/2;

[rv,tv,pv]=xyz2rtp(xyz(:,1),xyz(:,2),xyz(:,3));
%ab length can be bigger than cd or pq as it is the beam we start with.

%look for biggest and smallest elements in the matrix, kill elements
%less than tol of the max.

if lengthab==0
    if lengthcd==0
        behaviour=4;
        if lengthpq==0
            error('No entries detected!')
            return
        end
        
    else
        behaviour=5;
        if lengthpq==0
            behaviour=6;
        end
    end
else
    if lengthcd==0
        behaviour=1;
        if lengthpq==0
            behaviour=0;
        end
    else
        behaviour=2;
        if lengthpq==0
            behaviour=3;
        end
    end
end

if verbose
    behaviour
end

E1 = zeros(size(xyz));
E2 = E1;
E3 = E1;
H1 = E1;
H2 = E1;
H3 = E1;

un=unique(n);

%NumberOfLoops=length(calcvecels)
if verbose
    tic
end

switch behaviour
    case 0
        
        for ii=1:length(un)
            
            if length(m(n==un(ii)))>=1
                indx=find(n==un(ii));
                abv=ab([indx;indx+lengthab]);
                
                [M3,N3]=vswf(un(ii),m(n==un(ii)),2*pi*rv,tv,pv,3);
                
                for jj=1:length(abv)/2
                    E1=E1+(abv(jj))*M3(:,[0:2]*length(indx)+jj)+(abv(length(abv)/2+jj)).*N3(:,[0:2]*length(indx)+jj);
                    H1=H1+(abv(jj))*N3(:,[0:2]*length(indx)+jj)+(abv(length(abv)/2+jj)).*M3(:,[0:2]*length(indx)+jj);
                end
                
            end
            
            %             if verbose
            %                 if ni==10
            %                     disp(['Estimated time to completion: ' num2str(2*toc*length(calcvecels)/11) ' seconds!']);
            %                 end
            %             end
            
        end
        
    case 1
        
        for ii=1:length(un)
            
            if length(m(n==un(ii)))>=1
                indx=find(n==un(ii));
                abv=ab([indx;indx+lengthab]);
                pqv=pq([indx;indx+lengthab]);
                
                [M1,N1,~,~,M3,N3]=vswf(un(ii),m(n==un(ii)),2*pi*rv,tv,pv);
                for jj=1:length(abv)/2
                    E1=E1+full(abv(jj))*M3(:,[0:2]*length(indx)+jj)+full(abv(length(abv)/2+jj)).*N3(:,[0:2]*length(indx)+jj);
                    H1=H1+full(abv(jj))*N3(:,[0:2]*length(indx)+jj)+full(abv(length(abv)/2+jj)).*M3(:,[0:2]*length(indx)+jj);
                    
                    E2=E2+full(pqv(jj))*M1(:,[0:2]*length(indx)+jj)+full(pqv(length(abv)/2+jj)).*N1(:,[0:2]*length(indx)+jj);
                    H2=H2+full(pqv(jj))*N1(:,[0:2]*length(indx)+jj)+full(pqv(length(abv)/2+jj)).*M1(:,[0:2]*length(indx)+jj);
                    
                end
                
            end
            
            %             if verbose
            %                 if ii==1
            %                     disp(['Estimated time to completion: ' num2str(2*toc*length(calcvecels)/11) ' seconds!']);
            %                 end
            %             end
            
        end
    case 2
        
        for ii=1:length(un)
            
            if length(m(n==un(ii)))>=1
                indx=find(n==un(ii));
                
                abv=ab([indx;indx+lengthab]);
                pqv=pq([indx;indx+lengthab]);
                cdv=cd([indx;indx+lengthab]);
                
                [M1,N1,~,~,M3,N3]=vswf(un(ii),m(n==un(ii)),2*pi*rv,tv,pv);
                
                [M3a,N3a] = vswf(un(ii),m(n==un(ii)),2*pi*rv*(relindex),tv,pv,3);
                
                for jj=1:length(abv)/2
                    E1=E1+full(abv(jj))*M3(:,[0:2]*length(indx)+jj)+full(abv(length(abv)/2+jj)).*N3(:,[0:2]*length(indx)+jj);
                    H1=H1+full(abv(jj))*N3(:,[0:2]*length(indx)+jj)+full(abv(length(abv)/2+jj)).*M3(:,[0:2]*length(indx)+jj);
                    
                    E2=E2+full(pqv(jj))*M1(:,[0:2]*length(indx)+jj)+full(pqv(length(abv)/2+jj)).*N1(:,[0:2]*length(indx)+jj);
                    H2=H2+full(pqv(jj))*N1(:,[0:2]*length(indx)+jj)+full(pqv(length(abv)/2+jj)).*M1(:,[0:2]*length(indx)+jj);
                    
                    E3=E3+full(cdv(jj))*M3a(:,[0:2]*length(indx)+jj)+full(cdv(length(abv)/2+jj)).*N3a(:,[0:2]*length(indx)+jj);
                    H3=H3+full(cdv(jj))*N3a(:,[0:2]*length(indx)+jj)+full(cdv(length(abv)/2+jj)).*M3a(:,[0:2]*length(indx)+jj);
                    
                end
                
            end
            
            %             if verbose
            %                 if ni==10
            %                     disp(['Estimated time to completion: ' num2str(2*toc*length(calcvecels)/11) ' seconds!']);
            %                 end
            %             end
            
        end
        
    case 4
        
        for ii=1:length(un)
            
            if length(m(n==un(ii)))>=1
                indx=find(n==un(ii));
                
                pqv=pq([indx;indx+lengthpq]);
                
                [M1,N1]=vswf(un(ii),m(n==un(ii)),2*pi*rv,tv,pv,1);
                
                for jj=1:length(pqv)/2
                    
                    E2=E2+full(pqv(jj))*M1(:,[0:2]*length(indx)+jj)+full(pqv(length(pqv)/2+jj)).*N1(:,[0:2]*length(indx)+jj);
                    H2=H2+full(pqv(jj))*N1(:,[0:2]*length(indx)+jj)+full(pqv(length(pqv)/2+jj)).*M1(:,[0:2]*length(indx)+jj);
                    
                end
                
            end
            
            %             if verbose
            %                 if ni==10
            %                     disp(['Estimated time to completion: ' num2str(2*toc*length(calcvecels)/11) ' seconds!']);
            %                 end
            %             end
            %
        end
    case default
        disp('This behaviour is not implimented')
        return
        
        
end
if verbose
    toc
end

H1=-1i*H1; %LOOK HERE TO FIX STUFF
H2=-1i*H2; %LOOK HERE TO FIX STUFF
H3=-1i*H3; %LOOK HERE TO FIX STUFF

%res flipped because it's a meshgrid
if behaviour==2||behaviour==3||behaviour==5||behaviour==6
    [Exv3,Eyv3,Ezv3] = rtpv2xyzv(E3(:,1),E3(:,2),E3(:,3),rv,tv,pv);
    Ex3=Exv3;
    index= isnan(Ex3);
    Ex3(index)=0;
    Ey3=Eyv3;
    index= isnan(Ey3);
    Ey3(index)=0;
    Ez3=Ezv3;
    index= isnan(Ez3);
    Ez3(index)=0;
    
    [Hxv3,Hyv3,Hzv3] = rtpv2xyzv(H3(:,1),H3(:,2),H3(:,3),rv,tv,pv);
    Hx3=Hxv3;
    index= isnan(Hx3);
    Hx3(index)=0;
    Hy3=Hyv3;
    index= isnan(Hy3);
    Hx3(index)=0;
    Hz3=Hzv3;
    index= isnan(Hz3);
    Hx3(index)=0;
    
    structureoutput.Einternal=[squeeze(Ex3),squeeze(Ey3),squeeze(Ez3)];
    structureoutput.Hinternal=[squeeze(Hx3),squeeze(Hy3),squeeze(Hz3)];
    
    warning('Must scale grid for internal fields by relative refractive index.')
end

if behaviour==1|behaviour==2|behaviour==4|behaviour==5
    
    [Exv2,Eyv2,Ezv2] = rtpv2xyzv(E2(:,1),E2(:,2),E2(:,3),rv,tv,pv);
    Ex2=Exv2;
    index= isnan(Ex2);
    Ex2(index)=0;
    Ey2=Eyv2;
    index= isnan(Ey2);
    Ey2(index)=0;
    Ez2=Ezv2;
    index= isnan(Ez2);
    Ez2(index)=0;
    
    [Hxv2,Hyv2,Hzv2] = rtpv2xyzv(H2(:,1),H2(:,2),H2(:,3),rv,tv,pv);
    Hx2=Hxv2;
    index= isnan(Hx2);
    Hx2(index)=0;
    Hy2=Hyv2;
    index= isnan(Hy2);
    Hx2(index)=0;
    Hz2=Hzv2;
    index= isnan(Hz2);
    Hx2(index)=0;
    
    structureoutput.Escattered=[squeeze(Ex2),squeeze(Ey2),squeeze(Ez2)];
    structureoutput.Hscattered=[squeeze(Hx2),squeeze(Hy2),squeeze(Hz2)];
    
end
if behaviour<4
    
    [Exv1,Eyv1,Ezv1] = rtpv2xyzv(E1(:,1),E1(:,2),E1(:,3),rv,tv,pv);
    Ex1=Exv1;
    index= isnan(Ex1);
    Ex1(index)=0;
    Ey1=Eyv1;
    index= isnan(Ey1);
    Ey1(index)=0;
    Ez1=Ezv1;
    index= isnan(Ez1);
    Ez1(index)=0;
    
    [Hxv1,Hyv1,Hzv1] = rtpv2xyzv(H1(:,1),H1(:,2),H1(:,3),rv,tv,pv);
    Hx1=Hxv1;
    index= isnan(Hx1);
    Hx1(index)=0;
    Hy1=Hyv1;
    index= isnan(Hy1);
    Hx1(index)=0;
    Hz1=Hzv1;
    index= isnan(Hz1);
    Hx1(index)=0;
    
    structureoutput.Eincident=[squeeze(Ex1),squeeze(Ey1),squeeze(Ez1)];
    structureoutput.Hincident=[squeeze(Hx1),squeeze(Hy1),squeeze(Hz1)];
end

return
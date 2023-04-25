function outStruct = calc3Dstrain(X,Y,Z,DX,DY,DZ)
% Function for calculating the 3D strain and rotation field from
% 3-component surface displacement fields. This version inverts for strains
% withing a gridded box that is defined by the step size in variable
% "step". This requires X and Y locations of observations in meters, Z elevations in meters, and
% two different displacement components, DX, DY, and DZ in meters. Each of these
% variables should be in a grid, rather than in a vector.

% Initiate variables
[n,m]       = size(X);
i           = 1;
Xs          = [];
Ys          = [];
shearXY     = [];
shearXZ     = [];
shearYZ     = [];
dilatation  = [];
S1          = [];
S2          = [];
S3          = [];
cosx        = [];
cosy        = [];
trd         = [];
plg         = [];
rotation    = [];
W           = 0.25*eye(12); %Regularization coefficient for mimimum-moment regularization

% Define grid spacing size
step        = 4; 
hstep       = step/2;


for k=1:5:n-5
    for j =1:5:m-5
        tX      = X(k:k+4,j:j+4);
        tY      = Y(k:k+4,j:j+4);
        tZ      = Z(k:k+4,j:j+4);
        tDX     = DX(k:k+4,j:j+4);
        tDY     = DY(k:k+4,j:j+4);
        tDZ     = DZ(k:k+4,j:j+4);
        id      = find(isfinite(tDX)&isfinite(tDY)&isfinite(tDZ));
        if(isempty(id)|length(id)<25)
            S1      = [S1; NaN];
            S2      = [S2;NaN];
            S3      = [S3;NaN];
            shearXY = [shearXY; NaN];
            shearXZ = [shearXZ; NaN];
            shearYZ = [shearYZ; NaN];
            dilatation = [dilatation; NaN];
            
            trd     = [trd;NaN];
            plg     = [plg;NaN];
            Xs      =[Xs;NaN];
            Ys      = [Ys;NaN];
        else
            
            tX=tX(id); tY=tY(id); tDX=tDX(id); tDY=tDY(id); tZ= tZ(id); tDZ=tDZ(id);
            
            Xs      = [Xs;mean(tX)]; % Center X of each grid cell
            Ys      = [Ys;mean(tY)]; % Center Y of each grid cell
            
            numVals = length(tX)*3;
            odds    = [1:3:numVals];
            evens   = [2:3:numVals];
            ups     = [3:3:numVals];
            
         % Build Green's function (G) and data vector (d)

            G       = zeros(numVals,12);
            d       = zeros(numVals,1);
            G(odds,1)   = 1;
            G(evens,2)  = 1;
            G(ups,3)    = 1;
            G(odds,4)   = tX-mean(tX);
            G(odds,5)   = tY-mean(tY);
            G(odds,6)   = tZ-mean(tZ);
            G(evens,7)  = tX-mean(tX);
            G(evens,8)  = tY-mean(tY);
            G(evens,9)  = tZ-mean(tZ);
            G(ups,10)   = tX-mean(tX);
            G(ups,11)   = tY-mean(tY);
            G(ups,12)   = tZ-mean(tZ);
            
            
            d(odds)     = tDX;
            d(evens)    = tDY;
            d(ups)      = tDZ;
            
            model = lsqlin([G; W],[d;zeros(12,1)]);
            if (isnan(sum(model)))
                S1      = [S1; NaN];
                S2      = [S2;NaN];
                S3      = [S3;NaN];
                shearXY = [shearXY; NaN];
                shearXZ = [shearXZ; NaN];
                shearYZ = [shearYZ; NaN];
                dilatation = [dilatation; NaN];
                
                trd     = [trd;NaN];
                plg     = [plg;NaN];
            else
                L       = [model(4) model(5) model(6);
                    model(7) model(8) model(9);
                    model(10) model(11) model(12)];
                E = 0.5*[L(1,1)*2 L(1,2)+L(2,1) L(1,3)+L(3,1);
                    L(1,2)+L(2,1) L(2,2)*2 L(2,3)+L(3,2);
                    L(3,1)+L(1,3) L(3,2)+L(2,3) 2*L(3,3)];
                rotation = [rotation;0.5*L(2,1)-L(1,2)];

                [a,b]           = eig(E);
                b               = diag(b);
                id1             = find(abs(b)==max(abs(b)));
                id2             = find(abs(b)==min(abs(b)));
                maxStrain       = a(:,id1);
                minStrain       = a(:,id2);
                cosx            = [cosx;maxStrain(1)];
                cosy            = [cosy; maxStrain(2)];
                S1              = [S1;b(id1(1))];
                S2              = [S2;b(2)];
                S3              = [S3;b(id2(1))];
                shearXY         = [shearXY;E(1,2)];
                shearXZ         = [shearXZ;E(1,3)];
                shearYZ         = [shearYZ;E(2,3)];
                dilatation      = [dilatation;sum(b)];
                
                trd             = [trd;my_Cart2Sph(maxStrain(2),maxStrain(1),0)];
                plg             = [plg;asind(a(id1(1),3))];
            end
        end
    end    
end

outStruct = struct('Xs',Xs,'Ys',Ys,'S1',S1,'S2',S2,'S3',S3,'shearXY',shearXY,'shearXZ',shearXZ,'shearYZ',shearYZ,'dilatation',dilatation,'trd',trd,'plg',plg,'rotation',rotation);

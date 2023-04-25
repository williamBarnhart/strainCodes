function outStruct = calc2Dstrain(X,Y,DX,DY)
% Function for calculating the 2D strain and rotation field from
% 2-component surface displacement fields. This version inverts for strains
% withing a gridded box that is defined by the step size in variable
% "step". This requires X and Y locations of observations in meters, and
% two different displacement components, DX and DY in meters. Each of these
% variables should be in a grid, rather than in a vector.

% Initiate variables
[n,m]       = size(X);
i           = 1;
Xs          = [];
Ys          = [];
shear       = [];
dilatation  = [];
S1          = [];
S2          = [];
cosx        = [];
cosy        = [];
trd         = [];
rotation    = [];
W           = 0.25*eye(6); %Regularization coefficient for mimimum-moment regularization

% Define grid spacing size
step        = 4; 
hstep       = step/2;

% Solve for strain field of each grid cell
for k=1:hstep:n-hstep
    for j =1:hstep:m-hstep
        tX      = X(k:k+hstep,j:j+hstep);
        tY      = Y(k:k+hstep,j:j+hstep);
        tDX     = DX(k:k+hstep,j:j+hstep);
        tDY     = DY(k:k+hstep,j:j+hstep);
        id      = find(isfinite(tDX)&isfinite(tDY));
        tX=tX(id); tY=tY(id); tDX=tDX(id); tDY=tDY(id);
        
        Xs          = [Xs;mean(tX)]; %Center X of each grid cell
        Ys          = [Ys;mean(tY)]; %Center Y of each grid cell
        
        numVals     = length(tX)*2;
        odds        = [1:2:numVals];
        evens       = [2:2:numVals];
        
        % Build Green's function (G) and data vector (d)
        G           = zeros(numVals,6);
        d           = zeros(numVals,1);
        G(odds,1)   = 1;
        G(evens,2)  = 1;
        G(odds,3)   = tX-mean(tX);
        G(odds,4)   = tY-mean(tY);
        G(evens,5)  = tX-mean(tX);
        G(evens,6)  = tY-mean(tY);
        d(odds)     = tDX;
        d(evens)    = tDY;
        
        % Solve for Displacement gradient tensor (L) and reference frame
        model = lsqlin([G; W],[d;zeros(6,1)]);
        if (isnan(sum(model)))
            S1          = [S1; NaN];
            S2          = [S2;NaN];
            shear       = [shear; NaN];
            dilatation  = [dilatation; NaN];
            cosx        = [cosx;NaN];
            cosy        = [cosy; NaN];
        else
            % Build displacement gradient tensor
            L           = [model(3) model(4);
                model(5) model(6)];
            
            % Decompose L into strain (E) and rotation tensors)
            E           = 0.5*[L(1,1)*2 L(1,2)+L(2,1);L(1,2)+L(2,1) L(2,2)*2];
            rotation    = [rotation;0.5*L(2,1)-L(1,2)];
            
            % Get eigenvectors and eigenvalues of E (direction and
            % magnitudes of principal strains)
            [a,b]       = eig(E);
            b           = diag(b);
            id1         = find(abs(b)==max(abs(b)));
            id2         = find(abs(b)==min(abs(b)));
            maxStrain   = a(:,id1);
            minStrain   = a(:,id2);
            cosx        = [cosx;maxStrain(1)];
            cosy        = [cosy; maxStrain(2)];
            S1          = [S1;b(id1(1))];
            S2          = [S2;b(id2(1))];
            shear       = [shear;E(2)];
            dilatation  = [dilatation;sum(diag(E))];
            trd         = [trd;my_Cart2Sph(maxStrain(2),maxStrain(1),0)]; %trend of S1
        end
    end
end

outStruct = struct('Xs',Xs,'Ys',Ys,'S1',S1,'S2',S2,'shear',shear,'dilatation',dilatation,'trd',trd,'rotation',rotation);

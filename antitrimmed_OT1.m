% check image
clear
P=imread('FLIR3099.jpg'); %MUJER
%P=imread('FLIR2257.jpg'); %HOMBRE

Pgris=0.299*P(:,:,1)+0.587*P(:,:,2)+0.114*P(:,:,3);
Pgris(:,[285:320])=[];
Pgris(:,[1:55])=[];


P=Pgris;
ima=P;
p=5;
k1=25;
k2=25;


fileID = fopen('salida.txt','a+');
fprintf(fileID,'%s %s %2.0f %6.0f %6.0f ',"antiTrimmed_OT","P1_denoised.mrc",5,128,p);
fclose(fileID);

 
ima=double(ima);
copy=ima;         % make a copy
ima=ima(:);       % vectorize ima
mi=min(ima);      % deal with negative 
ima=ima-mi+1;     % and zero values

s=length(ima);

% create image histogram

m=max(ima)+1;
h=zeros(1,m);
hc=zeros(1,m);

for i=1:s
  if(ima(i)>0) h(ima(i))=h(ima(i))+1;end
end

ind=find(h);hl=length(ind);


long=(1:length(h));


h(h==0)=[];


long=ind;
T=length(h);




% %k-centrum%
lambda=zeros(1,T);
for i=1:T-k2
lambda(i)=0;
end
for i=T-k2+1:T
lambda(i)=1;
end


 D=abs(long(:)-long(:)');
 Dist=D;
 for i=1:T
     D(i,:)=D(i,:)*h(i);
 end

C=D;
D=unique(sort(D));
for i=2:length(D)
    D2(i-1)=D(i)-D(i-1);    
end
D(1)=[];

gamma=max(lambda(:));

alpha1=length(find(lambda==0));%kcentrum
alpha2=0;%kcentrum



beta=length(find(lambda==max(lambda)));

DD=gamma'*D2; %cada fila es lamba(fila)*D
DD=DD';
DD=DD(:);
LD=length(DD);

%[u x y d t]
CC=C';
fu=ones(1,T^2).*CC(:)';
fxfy=zeros(1,T^2+T);
for i=1:T-1
    lambdadif(i)=lambda(T-i+1)-lambda(T-i);
end
lambdadif(T)=lambda(1);

di=ones(1,T);
t=(k2)*ones(1);


ff=[fu,fxfy,di,t];



%%%\sum_k x_{jk}C_{jk}-d_{i}-t\leq 0
frestricu1=sparse(T,T^2);
 frestricx1=sparse(T,T^2);
 frestricy1=sparse(T,T); 
 frestricd1=sparse(T,T);
 frestrictt1=sparse(T,1);
 ve=zeros(T,T^2);
 

 
 frestricx1=[];
 frestrictt1=[];
  frestricd1=[];
for j=1:T
    vector=ones(1,T).*C(j,:);
    vector2=sparse(padarray(vector,[0 (j-1)*T],0,'pre'));
    vector2=sparse(padarray(vector2,[0 T^2-(T*j)],0,'post'));
    %ve=repmat(vector2,T,1);
    frestricx1(j,:)=vector2;

end
for i=1:T
    frestrictt1(i,1)=-1;
    i
end
frestricd1=speye(T)*(-1);

frestrict1=[frestricu1, frestricx1, frestricy1, frestricd1, frestrictt1];

clear frestricy1 frestricz1 frestricu1

%%%x_{ij}\leq y_j
frestricu2=sparse(T^2,T^2);
frestricx2=sparse(T^2,T^2);
frestricy2=sparse(T^2,T); 
frestricd2=sparse(T^2,T);
frestrictt2=sparse(T^2,1);

frestricx2=speye(T^2);


frestricy2=repmat(speye(T)*(-1),T,1);

frestrict2=[frestricu2, frestricx2, frestricy2, frestricd2, frestrictt2];

 %%%u_{ij}\leq y_j
frestricu3=sparse(T^2,T^2);
frestricx3=sparse(T^2,T^2);
frestricy3=sparse(T^2,T); 
frestricd3=sparse(T^2,T);
frestrictt3=sparse(T^2,1);

frestricu3=speye(T^2);


frestricy3=repmat(speye(T)*(-1),T,1);

frestrict3=[frestricu3, frestricx3, frestricy3, frestricd3, frestrictt3];

% %     frestrict1(aux,:)=[];
clear frestricy1 frestricz1 frestricu1


%%%sum_j y_j=p
frestricu4=sparse(1,T^2);
frestricx4=sparse(1,T^2);
frestricy4=sparse(1,T);
frestricd4=sparse(1,T);
frestrictt4=sparse(1,1);

 for i=1:T
     %frestricy(i,i)=-1; % es (i,i)
     frestricy4(1,i)=1;
 end


frestrict4=[frestricu4,frestricx4,frestricy4,frestricd4,frestrictt4];
clear frestricy1 frestricz1 frestricu2 frestricv2

%%sum_j x_{ij} = 1   \forall i=1,...,n
frestricu5=sparse(T,T^2);
frestricx5=sparse(T,T^2);
frestricy5=sparse(T,T);
frestricd5=sparse(T,T);
frestrictt5=sparse(T,1);
vector=ones(1,T);

for i=1:T
    vector2=sparse(padarray(vector,[0 (i-1)*T],0,'pre'));
    vector2=sparse(padarray(vector2,[0 T^2-(T*i)],0,'post'));
    frestricx5(i,:)=vector2;
    i
end

frestrict5=[frestricu5,frestricx5,frestricy5,frestricd5,frestrictt5];
clear frestricy2 frestricz2 vector2 vector

%%sum_j u_{ij} \leq 1   \forall i=1,...,n
frestricu6=sparse(T,T^2);
frestricx6=sparse(T,T^2);
frestricy6=sparse(T,T);
frestricd6=sparse(T,T);
frestrictt6=sparse(T,1);
vector=ones(1,T);

for i=1:T
    vector2=sparse(padarray(vector,[0 (i-1)*T],0,'pre'));
    vector2=sparse(padarray(vector2,[0 T^2-(T*i)],0,'post'));
    frestricu6(i,:)=vector2;

    i
end

frestrict6=[frestricu6,frestricx6,frestricy6,frestricd6,frestrictt6];
clear frestricy2 frestricz2 vector2 vector


%%% \sum_ij u_{ij}= k1  
frestricu7=ones(1,T^2);
frestricx7=sparse(1,T^2);
frestricy7=sparse(1,T);
frestricd7=sparse(1,T)
frestrictt7=sparse(1,1);

frestrict7=[frestricu7,frestricx7,frestricy7,frestricd7,frestrictt7];

 %%% u_{ij} \leq x_{ij}  
frestricu8=sparse(T^2,T^2);
frestricx8=sparse(T^2,T^2);
frestricy8=sparse(T^2,T);
frestricd8=sparse(T^2,T);
frestrictt8=sparse(T^2,1);

frestricu8=speye(T^2);
frestricx8=speye(T^2)*-1;

frestrict8=[frestricu8,frestricx8,frestricy8,frestricd8,frestrictt8];

%%% \sum_a / d_{i,a}> d_{i,j}  x_{i,a}  + y_j \leq 1   \forall i,j
frestricu9=sparse(T^2,T^2);
frestricx9=sparse(T^2,T^2);
frestricy9=sparse(T^2,T);
frestricd9=sparse(T^2,T);
frestrictt9=sparse(T^2,1);
 
frestricx9=[];
frestricx9bloque=sparse(T,T^2);
 for i=1:T
        frestricx9bloque=sparse(T,T^2);
     for j=1:T
          frestricx9bloque(j,(i-1)*T+find(C(i,:)>C(i,j)))=1; 
     end 
     frestricx9=[frestricx9;frestricx9bloque];
 end
 
frestricy9=repmat(speye(T),T,1);

frestrict9=[frestricu9, frestricx9, frestricy9, frestricd9, frestrictt9];

clear frestricx5 frestricy5 frestricu5 frestricv5




Aineq = [frestrict1;frestrict2;frestrict8;frestrict9];%;frestrict10];
bineq = [sparse(1,T)';sparse(1,T^2)';sparse(1,T^2)';ones(1,T^2)'];%;ones(1,T^2)'];%
Aeq = [frestrict4;frestrict5;frestrict7];
beq = [p;ones(1,T)';k1];
bineq2=ones(size(Aineq,1),1)*-inf;

clear frestrict1 frestrict2 frestrict3 frestrict4 frestrict5 


vectorulb=zeros(1,T^2);
vectoruub=ones(1,T^2);

vectorxlb=zeros(1,T^2);
vectorxub=ones(1,T^2);


vectorylb=zeros(1,T);
vectoryub=ones(1,T);

vectordlb=zeros(1,T);
vectordub=ones(1,T)*inf;

  vectortlb=ones(1,1)*-inf;
vectortub=ones(1,1)*inf;

lb    = [vectorulb';vectorxlb';vectorylb';vectordlb';vectortlb'];
ub    = [vectoruub';vectorxub';vectoryub';vectordub';vectortub'];
cadena='';

cadena(1:T^2)='I';
cadena(T^2+1:T^2+T^2)='I';
cadena(T^2+T^2+1:T^2+T^2+T)='I';
cadena(T^2+T^2+T+1:T^2+T^2+T+T)='C';
cadena(T^2+T^2+T+1+1:T^2+T^2+T+T+1)='C';
ctype = cadena;

   

options = cplexoptimset('cplex');
options.Display = 'on';

     options.NodeFileInd=3;

ff=ff';

 options.timelimit=7200;
[x, fval, exitflag, output] = cplexlp (ff', Aineq, bineq, Aeq, beq, lb, ub, [ ], options);
       
   
RLval=fval;
RLtime=output.time;

fileID = fopen('salida.txt','a+');
fprintf(fileID,'%6.2f %6.2f ',RLval,RLtime);
fclose(fileID);
 
  
  %%%%%%%%%%%%%%Problema entero%%%%%%%%%%%%%%%%
   
   
   
cplex = Cplex('mipex1');
cplex.Model.sense = 'minimize';
cplex.Model.obj   = [ff'];
cplex.Model.lb    = [ lb];
cplex.Model.ub    = [ub];
cplex.Model.ctype = cadena;
cplex.Model.A =  [Aineq;Aeq];
cplex.Model.lhs = [bineq2; beq];
cplex.Model.rhs = [  bineq;  beq];

cplex.Param.timelimit.Cur=7200;
      



cplex.solve();

x=cplex.Solution.x;
gap = cplex.Solution.miprelgap;
%   
opt = cplexoptimset('exportmodel', 'myModel.lp');

tol=1e-5;
   for i=T^2+T^2+1:T^2+T^2+T
       xx(i-T^2-T^2)=x(i);
   end
   centroides=find(xx>=tol);
   
   for i=1:length(centroides)
       x((centroides(i)-1)*T+centroides(i))=1;
   end

for i=1:T^2
     imagen1(i)=x(i);
end
 for i=T^2+1:T^2+T^2
     imagen2(i-T^2)=x(i);
 end
 imagen=imagen1+imagen2;

   for i=T^2+T^2+T+1:T^2+T^2+T+T
       xd(i-(T^2+T^2+T))=x(i);
   end
    for i=T^2+T^2+T+T+1:T^2+T+T+1
       xt(i-(T^2+T^2+T+T^2))=x(i);
    end
    valores=find(imagen>=tol);
  
    
    imagen=reshape(imagen,T,T);
    imagen=imagen';
    %centro=find(centroides);
    centro=long(centroides(:));
    
    for i=1:length(centroides)
        imagen(centroides(i),centroides(i))=1;
    end
    
    
       sumaasignaciones=sum(imagen');
     vec=[];
     for i=1:T
         if sumaasignaciones(i)==0
             for j=1:length(centroides)
                 vec(j)=C(i,centroides(j));
             end
             aux=find(C(i,:)==min(vec));
             longitud=length(aux);
             for ii=1:longitud
                 if sum(centroides==aux(ii))==1
                  imagen(i,aux(ii))=1;
                 end
             end
         end
     end
    

seg=zeros(size(P,1)*size(P,2),1);
for i=1:length(centro)
    aux(i)=i;
end
     for i=1:length(centro)
%for i=1:1
        im=ind(find(imagen(:,centroides(i))));
        for j=1:length(im)
        mask=find(ima==im(j));
            %seg(mask)=centro(i); 
            seg(mask)=aux(i); 
        end
     end
     
   seg3D=reshape(seg,size(P,1),size(P,2),size(P,3));
segmentacion=seg3D/max(seg3D(:)); %Normalizamos la segmentacion
figure; imshow(segmentacion,[])

%Para calcular el n√∫mero de pixeles en un cluster
numero_pixeles_cluster_blanco=sum(segmentacion(:)==max(segmentacion(:)));
numero_pixeles_cluster_negro=sum(segmentacion(:)==min(segmentacion(:)));

   
   
nodos=cplex.Solution.nodecnt;
 bb=cplex.Solution.objval-gap*cplex.Solution.objval;
 gapRL=(cplex.Solution.objval-RLval)/cplex.Solution.objval;
  fileID = fopen('salida.txt','a+');
 fprintf(fileID,' %6.2f %6.2f %6.2f %6.4f %6.4f %6.0f\r\n',cplex.Solution.objval,bb,cplex.Solution.time,gap,gapRL,nodos);
fclose(fileID); 

%%%%PARA MOSTRAR EL BORDE DEL CLUSTER DE COLOR ROJO%%%%%
 seg2=seg3D/max(seg3D(:));
 seg2(seg2(:,:)<1)=0;
 ed=edge(seg2);
 ed=logical(ed);
 
% 
   seg_or = cat(3, seg3D/max(seg3D(:)), seg3D/max(seg3D(:)), seg3D/max(seg3D(:)));
   seg_or=uint8(255*seg_or);
% 
 imshow(seg_or,[])
 red = cat(3, ones(size(ed)), zeros(size(ed)), zeros(size(ed)));
 hold on
 h = imshow(red);
 hold off
 set(h, 'AlphaData', ed)
 
 
fprintf('\n--- n˙mero de pÌxeles por cl˙ster ---\n')
for i = 1:p
    num_pixeles = sum(seg == i);
    fprintf('Cl˙ster %d: %d pÌxeles\n', i, num_pixeles);
end 
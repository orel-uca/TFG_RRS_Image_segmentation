% check image
clear
%P=imread('FLIR3099.jpg'); %MUJER
P=imread('FLIR2257.jpg'); %HOMBRE


Pgris=0.299*P(:,:,1)+0.587*P(:,:,2)+0.114*P(:,:,3);
Pgris(:,[285:320])=[];
Pgris(:,[1:55])=[];


   P=Pgris;
ima=P;
p=6;
k1=20;
k2=20;


  fileID = fopen('salida.txt','a+');
 fprintf(fileID,'%s %s %2.0f %6.0f %6.0f ',"trimmed_OT","P1_denoised.mrc",5,128,p);
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

%trimmed mean
kk=k1;
lambda=ones(1,T)*1;
for i=1:k1
lambda(i)=0;
end
for i=1:k2
lambda(T-(i-1))=0;
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
%D2(1)=[];
D(1)=[];

gamma=max(lambda(:));

alpha1=length(find(lambda==0));%kcentrum
alpha2=0;%kcentrum

beta=length(find(lambda==max(lambda)));

DD=gamma'*D2; %cada fila es lamba(fila)*D
DD=DD';
DD=DD(:);
LD=length(DD);

%[x y d t]
fxfy=zeros(1,T^2+T);

for i=1:T-k2-1
     lambdadif(i)=lambda(T-k2-i+1)-lambda(T-k2-i);
end
lambdadif(T-k2)=0;

dij=ones(T-k2,T);
ti=ones(1,T-k2);
for i=1:T-k2
    dij(i,:)=lambdadif(i).*dij(i,:);
    ti(i)=lambdadif(i)*ti(i)*i;
end
dij=dij'; dij=dij(:);



ff=[fxfy,dij',ti];


%%%\sum_k x_{jk}C_{jk}-d_{ij}-t_i\leq 0
 frestricx1=sparse((T-k2)*T,T^2);
 frestricy1=sparse((T-k2)*T,T); 
 frestricd1=sparse((T-k2)*T,(T-k2)*T);
 frestrictt1=sparse((T-k2)*T,T-k2);
 ve=zeros(T,T^2);
 

  frestricx1=[];
 frestrictt1=[];
for j=1:T
    vector=ones(1,T).*C(j,:);
    vector2=sparse(padarray(vector,[0 (j-1)*T],0,'pre'));
    vector2=sparse(padarray(vector2,[0 T^2-(T*j)],0,'post'));
    %ve=repmat(vector2,T,1);
    frestricx1aux(j,:)=vector2;

end
for i=1:T-k2
    frestricx1=[frestricx1;frestricx1aux];

                vectort=sparse(T,T-k2);
                vectort(:,i)=ones(T,1)*(-1);
    frestrictt1=[frestrictt1;vectort];
i
 end
frestricd1=speye((T-k2)*T)*-1;

 frestrict1=[frestricx1, frestricy1, frestricd1, frestrictt1];

clear frestricy1 frestricz1 frestricu1

%%%x_{ij}\leq y_j
 frestricx2=sparse(T^2,T^2);
 frestricy2=sparse(T^2,T); 
 frestricd2=sparse(T^2,(T-k2)*T);
 frestrictt2=sparse(T^2,T-k2);

frestricx2=speye(T^2);


   frestricy2=repmat(speye(T)*(-1),T,1);
   
 frestrict2=[frestricx2, frestricy2, frestricd2, frestrictt2];
 
% %     frestrict1(aux,:)=[];
clear frestricy1 frestricz1 frestricu1


%%%sum_j y_j=p
frestricx3=sparse(1,T^2);
frestricy3=sparse(1,T);
frestricd3=sparse(1,(T-k2)*T);
frestrictt3=sparse(1,T-k2);

 for i=1:T
     %frestricy(i,i)=-1; % es (i,i)
     frestricy3(1,i)=1;
 end


frestrict3=[frestricx3,frestricy3,frestricd3,frestrictt3];
clear frestricy1 frestricz1 frestricu2 frestricv2

%%sum_j x_{ij} \leq 1   \forall i=1,...,n
 frestricx4=sparse(T,T^2);
 frestricy4=sparse(T,T);
 frestricd4=sparse(T,(T-k2)*T);
 frestrictt4=sparse(T,T-k2);
 vector=ones(1,T);

for i=1:T
    vector2=sparse(padarray(vector,[0 (i-1)*T],0,'pre'));
    vector2=sparse(padarray(vector2,[0 T^2-(T*i)],0,'post'));
    frestricx4(i,:)=vector2;
    i;
end

frestrict4=[frestricx4,frestricy4,frestricd4,frestrictt4];
clear frestricy2 frestricz2 vector2 vector

%%% \sum_ij x_{ij}= n-(k2)  
frestricx5=ones(1,T^2);
frestricy5=sparse(1,T);
frestricd5=sparse(1,(T-k2)*T);
frestrictt5=sparse(1,T-k2);

frestrict5=[frestricx5,frestricy5,frestricd5,frestrictt5];


%%% \sum_a / d_{i,a}> d_{i,j}  x_{i,a} + y_j \leq 1   \forall i,j
frestricx6=sparse(T^2,T^2);
frestricy6=sparse(T^2,T);
frestricd6=sparse(T^2,(T-k2)*T);
frestrictt6=sparse(T^2,T-k2);

frestricx6=[];
frestricx6bloque=sparse(T,T^2);
for i=1:T
    frestricx6bloque=sparse(T,T^2);
    for j=1:T
        frestricx6bloque(j,(i-1)*T+find(C(i,:)>C(i,j)))=1; 
    end 
    frestricx6=[frestricx6;frestricx6bloque];
end
 
frestricy6=repmat(speye(T),T,1);
 
frestrict6=[frestricx6, frestricy6, frestricd6, frestrictt6];

clear frestricx5 frestricy5 frestricu5 frestricv5



Aineq = [frestrict1;frestrict2;frestrict4];
bineq = [sparse(1,(T-k2)*T)';sparse(1,T^2)';ones(1,T)'];
Aeq = [frestrict3;frestrict5];
beq = [p;T-k2];
bineq2=ones(size(Aineq,1),1)*-inf;

clear frestrict1 frestrict2 frestrict3 frestrict4 frestrict5 


vectorxlb=zeros(1,T^2);
vectorxub=ones(1,T^2);


vectorylb=zeros(1,T);
vectoryub=ones(1,T);

vectorulb=zeros(1,(T-k2)*T);
vectoruub=ones(1,(T-k2)*T)*inf;

  vectorvlb=ones(1,T-k2)*-inf;
vectorvub=ones(1,T-k2)*inf;

lb    = [vectorxlb';vectorylb';vectorulb';vectorvlb'];
ub    = [vectorxub';vectoryub';vectoruub';vectorvub'];
cadena='';
cadena(1:T^2)='C';

cadena(T^2+1:T^2+T)='C';
cadena(T^2+T+1:T^2+T+(T-k2)*T)='C';
cadena(T^2+T+(T-k2)*T+1:T^2+T+(T-k2)*T+T-k2)='C';
ctype = cadena;


%options = cplexoptimset();
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
 
  
  %%%%%%%%%%%Problema entero
  
  
  cadena='';
  cadena(1:T^2)='I';

cadena(T^2+1:T^2+T)='I';
cadena(T^2+T+1:T^2+T+(T-k2)*T)='C';
cadena(T^2+T+(T-k2)*T+1:T^2+T+(T-k2)*T+T-k2)='C';
ctype = cadena;

   
options = cplexoptimset('cplex');
options.Display = 'on';

options.NodeFileInd=3;

ff=ff';
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
   for i=T^2+1:T^2+T
       xx(i-T^2)=x(i);
   end
   centroides=find(xx>=tol);
   
   for i=1:length(centroides)
       x((centroides(i)-1)*T+centroides(i))=1;
   end

for i=1:T^2
     imagen(i)=x(i);
 end


    valores=find(imagen>=tol);
  
    
    imagen=reshape(imagen,T,T);
    imagen=imagen';
    
        imagenfval=imagen;
    fval2=0;
    for i=1:T
        if length(find(imagenfval(i,:)==1))>0
        fval2=fval2+C(i,find(imagenfval(i,:)==1));
        fval22(i)=C(i,find(imagenfval(i,:)==1));
        end
    end
    fval22=sort(fval22); fval22(fval22==0)=[];  
    for i=1:k1-p
    fval22(1)=[]; 
    end
    fval22=sum(fval22);
    
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
                  imagen(i,aux(ii))=1;  ii=longitud+1;
                 end
             end
         end
     end
    

seg=zeros(size(P,1)*size(P,2),1);
for i=1:length(centro)
    aux(i)=i;
end
 for i=1:length(centro)
    im=ind(find(imagen(:,centroides(i))));
    for j=1:length(im)
        mask=find(ima==im(j));
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

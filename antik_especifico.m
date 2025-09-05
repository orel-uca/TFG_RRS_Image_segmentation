% check image
clear
tic
P=imread('FLIR3099.jpg'); %MUJER
%P=imread('FLIR2257.jpg'); %HOMBRE


Pgris=0.299*P(:,:,1)+0.587*P(:,:,2)+0.114*P(:,:,3);
Pgris(:,[285:320])=[];
Pgris(:,[1:55])=[];



P=Pgris

  
ima=P;
porcentaje=40;
p=4;
k=20;

  fileID = fopen('salida.txt','a+');
 fprintf(fileID,'%s %s %2.0f %6.0f %6.0f',"antik_especifico","P1_denoised.mrc",5,128,p);
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


%%%anti k-centrum
porcentaje=floor(porcentaje*length(ind)/100);
k=porcentaje;
 for i=1:k
lambda(i)=1;
end
for i=k+1:T
lambda(i)=0;
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


alpha1=0;%antikcentrum
alpha2=length(find(lambda==0));%antikcentrum

beta=length(find(lambda==max(lambda)));

DD=gamma'*D2; %cada fila es lamba(fila)*D
DD=DD';
DD=DD(:);
LD=length(DD);

%[y z u]
CC=C';
f1=ones(1,T^2).*CC(:)';
f2=zeros(1,T);

ff=[f1,f2];


%%%x_{ij}\leq y_j
frestricx1=sparse(T^2,T^2);
frestricy1=sparse(T^2,T); 
frestricx1=speye(T^2);
frestricy1=repmat(speye(T)*(-1),T,1);
frestrict1=[frestricx1, frestricy1];
 
% %     frestrict1(aux,:)=[];
clear frestricy1 frestricz1 frestricu1


%%%sum_j y_j=p
frestricx2=sparse(1,T^2);
frestricy2=sparse(1,T);

 for i=1:T
     %frestricy(i,i)=-1; % es (i,i)
     frestricy2(1,i)=1;
 end


frestrict2=[frestricx2,frestricy2];%[[opMatrix(frestric)],A];
clear frestricy1 frestricz1 frestricu2 frestricv2

%%sum_j x_{ij} <= 1   \forall i=1,...,n
 frestricx3=sparse(T,T^2);
 frestricy3=sparse(T,T);
 vector=ones(1,T);

for i=1:T
    vector2=sparse(padarray(vector,[0 (i-1)*T],0,'pre'));
    vector2=sparse(padarray(vector2,[0 T^2-(T*i)],0,'post'));
    frestricx3(i,:)=vector2;
    i
end

 frestrict3=[frestricx3,frestricy3];
clear frestricy2 frestricz2 vector2 vector

%%% \sum_ij x_{ij}= n/2  
frestricx4=ones(1,T^2);
 frestricy4=sparse(1,T);

 frestrict4=[frestricx4,frestricy4];

 
clear frestricx5 frestricy5 frestricu5 frestricv5

Aineq = [frestrict1;frestrict3];
bineq = [sparse(1,T^2)';ones(1,T)'];%
Aeq = [frestrict2;frestrict4];
%beq = [p;floor(porcentaje2*T/100)];
beq = [p;k];
bineq2=ones(size(Aineq,1),1)*-inf;

clear frestrict1 frestrict2 frestrict3 frestrict4 frestrict5 


vectorxlb=zeros(1,T^2);
vectorxub=ones(1,T^2);


vectorylb=zeros(1,T);
vectoryub=ones(1,T);


lb    = [vectorxlb';vectorylb'];
ub    = [vectorxub';vectoryub'];
cadena='';
cadena(1:T^2)='C';

cadena(T^2+1:T^2+T)='C';
ctype = cadena;

options = cplexoptimset('cplex');
options.Display = 'on';
    
options.NodeFileInd=3;
ff=ff';
[x, fval, exitflag, output] = cplexlp (ff', Aineq, bineq, Aeq, beq, lb, ub, [ ], options);
           
       
 RLval=fval;
 RLtime=output.time;
  
 fileID = fopen('salida.txt','a+');
 fprintf(fileID,'%6.2f %6.2f',RLval,RLtime);
 fclose(fileID);
      
      
      %%%%Problema entero%%%%%
      
cadena='';
cadena(1:T^2)='I';

cadena(T^2+1:T^2+T)='I';
ctype = cadena;

options = cplexoptimset('cplex');
options.Display = 'on';

options.NodeFileInd=3;
ff=ff';
cplex = Cplex('');
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
     
     
seg=zeros(size(P,1)*size(P,2)*size(P,3),1);
clear aux
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
segmentacion=seg3D/max(seg3D(:));



 
 %%%% BORDES ROJOS EN ESCALA DE GRISES ENTRE TODOS LOS CLÚSTERES %%%%
seg_normal = double(seg3D) / max(seg3D(:));  
ed = edge(seg_normal, 'sobel');             
ed = logical(ed);

figure;
imshow(seg_normal, []);                    
hold on
red = cat(3, ones(size(ed)), zeros(size(ed)), zeros(size(ed)));  
h = imshow(red);
set(h, 'AlphaData', ed * 0.8);              
hold off

%numero de pixeles por cluster
 
fprintf('\n--- número de píxeles por clúster ---\n')
for i = 1:p
    num_pixeles = sum(seg == i);
    fprintf('Clúster %d: %d píxeles\n', i, num_pixeles);
end







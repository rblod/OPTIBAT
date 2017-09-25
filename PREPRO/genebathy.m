 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %Largeur terrace
% genebathy.m
function z=genebathy

X=0:499;%vecteur position

% interval number bars
ib=[0 5];
bnb=round((ib(2)-ib(1)).*rand(1) + ib(1));%number of bars
%slope range
ib=[0.02 0.05];
sl=((ib(2)-ib(1)).*rand(1,1) + ib(1)); %linear beach slope

%bar amplitude range
ib=[0.25 6];
bamp=sort(((ib(2)-ib(1)).*rand(bnb,1) + ib(1)));

%bar position range
ib=[1*length(X)/5 3*length(X)/5];
bx=sort(round((ib(2)-ib(1)).*rand(bnb,1) + ib(1))); 

%bar width range
ib=[50 600];
bw=sort(round((ib(2)-ib(1)).*rand(bnb,1) + ib(1))); 

z=-X.*sl;    
for ib=1:bnb
    bw(ib)=3.*bx(ib);
    while z(bx(ib))+bamp(ib)*exp(  -((X(bx(ib))-bx(ib)).^2)/(bw(ib)) )>z(bx(ib))/3;bamp(ib)=bamp(ib)/2;end 
z=z+bamp(ib)*exp(  -((X-bx(ib)).^2)/(bw(ib)) );
end

% figure(56);plot(zf)

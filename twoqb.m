%Example given in Zurek 2002
discord=[];
z=0:0.01:1;
theta=-pi/4:0.1:pi/4;
for r=1:1:length(z)
%Define density matrix, random for now
%rho=randn(4);
%Normalize and symmetrize
%rho=rho*rho'/trace(rho*rho');
rho=[0.5 0 0 z(r)/2; 0 0 0 0; 0 0 0 0; z(r)/2 0 0 0.5];


%von Neumann entropy
e = eig(rho);
hab = -e'*log2(e+(e==0));

%Partial trace
rhoa=[rho(1,1)+rho(2,2) rho(1,3)+rho(2,4);
      rho(1,3)+rho(2,4) rho(3,3)+rho(4,4)];
rhob=[rho(1,1)+rho(3,3) rho(2,1)+rho(4,3);
      rho(2,1)+rho(4,3) rho(2,2)+rho(4,4)];

%entropy of a
ha=-eig(rhoa)'*log2(eig(rhoa)+(eig(rhoa)==0));
%entropy of b
hb=-eig(rhob)'*log2(eig(rhob)+(eig(rhob)==0));
%mutual info
minfo=ha+hb-hab;

for t=1:1:length(theta)
%Measurement basis from eigenvectors
mb=[];
phi=pi;
mb(:,:,1)=kron(eye(2),[(cos(theta(t)))^2, exp(i*phi)*sin(theta(t))*cos(theta(t)); exp(i*phi)*sin(theta(t))*cos(theta(t)) (sin(theta(t)))^2]);
mb(:,:,2)=kron(eye(2),[(sin(theta(t)))^2, -exp(-i*phi)*sin(theta(t))*cos(theta(t));-exp(-i*phi)*sin(theta(t))*cos(theta(t)) (cos(theta(t)))^2]);

  had=[];
  for j=1:size(mb,3)
      %Eq. 9 in zurek
      rhoad=mb(:,:,j)*rho*mb(:,:,j)/(trace(mb(:,:,j)*rho));
      %Conditonal entropy given measurement
      e = eig(rhoad);
      had(j) = (trace(mb(:,:,j)*rho))*(-e'*log2(e+(e==0)));
  end

%Discord
discord(r,t)=ha-hab+sum(had);

end

end

mesh(sqrt(real(discord).^2+imag(discord).^2))

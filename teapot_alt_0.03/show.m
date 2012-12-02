comm='';

load out/time.dat;
sizet=size(time);
for i=0:1:sizet-1,
  comm=sprintf('load out/psi%.4d.dat;',i);
  eval(comm);
  comm=sprintf('load out/theta%.4d.dat;',i);
  eval(comm);
  comm=sprintf('load out/w%.4d.dat;',i);
  eval(comm);
%newplot;
%hold on;
  fprintf(1,'%.4d\n',i);

  subplot(2,2,1);
%  comm=sprintf('pcolor(w%.4d);',i);
  comm=sprintf('contour(w%.4d,16);',i);
  eval(comm);
  shading interp;
%  axis([1 200 1 20]);
  s=sprintf('Omega. T=%f, nom=%d of %d',time(i+1),i+1,sizet(1));
  title(s);

  subplot(2,2,2);
%  comm=sprintf('pcolor(theta%.4d);',i);
  comm=sprintf('contour(theta%.4d,16);',i);
  eval(comm);
  shading interp;
 % axis([1 200 1 20]);
%  s=sprintf('T=%f, nom=%d of %d\n',time(i+1),i+1,sizet);
  title('theta');

  subplot(2,2,3);
%  comm=sprintf('pcolor(psi%.4d);',i);
  comm=sprintf('contour(psi%.4d,16);',i);
  eval(comm);
  shading interp;
%  axis([1 200 1 20]);
%  s=sprintf('T=%f, nom=%d of %d\n',time(i+1),i+1,sizet);
  title('psi');

%  comm=sprintf('sizepsi=size(psi%.4d);',i);
%  eval(comm);
%  vx=zeros(sizepsi(1),sizepsi(2));
%  vy=zeros(sizepsi(1),sizepsi(2));

%  for k=2:sizepsi(1)-1,
%    for l=2:sizepsi(2)-1,
%      comm=sprintf('vx(k,l)=(psi%.4d(k+1,l)-psi%.4d(k-1,l))/0.4;',i,i);
%      eval(comm);
%      comm=sprintf('vy(k,l)=-(psi%.4d(k,l+1)-psi%.4d(k,l+1))/0.4;',i,i);
%      eval(comm);
%    end;
%  end;

%  subplot(2,2,4);
  %pcolor(vx);
  %shading interp;
%  k=1:sizepsi(1);
%  l=1:sizepsi(2);
%  quiver(k,l,vx,vy,2);
%  title('v');
    
% Compare tool

%  comm=sprintf('load out1/rho%.4d.dat;',i);
%  eval(comm);
%  comm=sprintf('plot(rho%.4d(5,:),%co-r%c);',i,39,39);
%  eval(comm);
%hold off;

% Saving tool

%fil=sprintf('png/gr%.4d',i);
%saveas(gcf,fil,'png');

comm=sprintf('print -dpng -r128 png/gr%.4d',i);
eval(comm);

%  pause;
  comm=sprintf('clear psi%.4d;',i);
  eval(comm);
  comm=sprintf('clear theta%.4d;',i);
  eval(comm);
  comm=sprintf('clear w%.4d;',i);
  eval(comm);
end;

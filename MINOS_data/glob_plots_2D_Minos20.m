chi_grid=load('pp.MINOS.gauss.chi.dat');   % Datos de considerar la JERAQUÍA NORMAL
cl_beam_68=load('MINOS_cl_68_beam_nh.dat');
cl_beam_90=load('MINOS_cl_90_beam_nh.dat');

[a2,b]=size(chi_grid);
a=sqrt(a2)
k=1;
l=1;
for i=1:a
    
   theta(i)= chi_grid(l,1) ;
   sin2th2(i)=(sin(2.0*theta(i)))^2;

    m=1;
   for j=1:a
    mass(m)=chi_grid(k,2)   ;
    grid(i,j)= chi_grid(k,3)  ;
   k=k+1;
   m=m+1;
   
   end
   l=l+a;
end


%contour3(X,Y,Z,[2,6,9]) %  graficar contour plot   
interpts=1000
del_x=(theta(a)-theta(1))/(interpts-1);
del_y=(mass(a)-mass(1))/(interpts-1);


[Xq,Yq]=meshgrid(theta(1):del_x:theta(a),mass(1):del_y:mass(a));

Zq = interp2(theta,mass,grid,Xq,Yq,'spline');


for i=1:interpts
 sin2_2th_int(i)= (sin(2.0*Xq(1,i)))^2;  
theta_int(i)=Xq(1,i);
 mass_int(i)=  Yq(i,1) ;
sin2_th_int(i)= (sin(Xq(1,i)))^2;  

end

%  figure
%  surf(mass_int,theta_int,Zq)

%%%%%  surface plot

%   figure
%   surf(mass,theta,grid)
 

  minimum_int = min(min(Zq));
minimum=min(min(grid))

[imin,jmin]=find(grid==minimum);
masss=mass(jmin)
s2_th=theta(imin)
grid(imin,jmin);

[imin_int,jmin_int]=find(Zq==minimum_int);


sin2_2th2_min=sin(2.0*theta(imin))^2;

sin2_2th2_min_int=sin2_2th_int(imin_int);


mass_min=mass(jmin);
mass_min_int=mass_int(jmin_int);

minteta=theta_int(imin_int);
minsin2_theta=sin2_th_int(imin_int);

chi2_min=grid(imin,jmin)
chi2_min_int=Zq(imin_int,jmin_int)


for i=1:a

chi2_dm(i)=min(grid(:,i))-minimum;
chi2_sin2th2(i)=min(grid(i,:))-minimum;

end

for i=1:interpts
chi2_dm_int(i)=min(Zq(:,i))-minimum_int;
chi2_sin2_2th_int(i)=min(Zq(i,:))-minimum_int;
 
end

%%%%%  chi2 vs dm 
figure
plot(mass,chi2_dm)
%  hold on 
%  plot(mass_int,chi2_dm_int)


ylim([0.0 10.0])

figure
plot(theta,chi2_sin2th2)
%  hold on
%  plot(sin2_2th_int,chi2_sin2_2th_int)

ylim([0.0 10.0])


Zq_zero=Zq-minimum_int;
%  figure
%  contour3(sin2_2th_int,mass_int,Zq_zero,[2,4,9])
m=1;
l=1;
k=1;
for i=1:interpts
    for j=1:interpts
        chi=Zq_zero(i,j);
        
        if (  (chi<=2.3) )
        clmass_1(k)=mass_int(j)*1000 ;   
        clsin2_2th_1(k)= sin2_2th_int(i)   ;
        clsin2_th_1(k)=sin2_th_int(i);
        cltheta_1(k)=theta_int(i);
        k=k+1;
        end
        
        
        %if (  (chi<=4.61) )
        if (  (chi<=4.3) )
        clmass_2(l)=mass_int(j)*1000 ;   
        clsin2_2th_2(l)= sin2_2th_int(i)   ;
        clsin2_th_2(l)=sin2_th_int(i);
        cltheta_2(l)=theta_int(i);

        l=l+1;
        end        
        
        if (  (chi<=6.630) )         
        clmass_3(m)=mass_int(j)*1000 ;   
        clsin2_2th_3(m)= sin2_2th_int(i)   ;
        clsin2_th_3(m)=sin2_th_int(i);
        cltheta_3(m)=theta_int(i);
        

        m=m+1;
        end               
    
        
     end
    
end

%  
 figure 
%  plot (cltheta_3,clmass_3,'* k')
%   hold on
 
 plot (cltheta_2,clmass_2,'* k')
  hold on
 plot (cltheta_1,clmass_1,'* g')
 plot(cl_beam_68(:,1),cl_beam_68(:,2)*1000,'* b')
 plot(cl_beam_90(:,1),cl_beam_90(:,2)*1000,'* r')


%  xlim([0.065 0.1])
%  ylim([2.2  2.8])
%  
%  
%  
xlim([0.3 0.75])
ylim([2.0  3.0]) % ylim([1.5  3.1])


%  figure
%  
%  plot (clsin2_th_3,clmass_3,'* k')
%  hold on
%  plot (clsin2_th_2,clmass_2,'* r')
%  hold on
%  plot (clsin2_th_1,clmass_1,'* g')

%  xlim([0.065 0.1])
%  ylim([2.2  2.8])



% chibin=load('chi2_bin.txt');
% figure
% bar(chibin(:,2))

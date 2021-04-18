


       subroutine  minos_2020( Y,min_chi2)
       implicit none
       include 'MINOS_20.inc'       
       	include 'UMAT4.inc'
       	include 'GLOB.inc'
       	
       	
       integer i,j,k
       real*8 Y(13)
       real*8  dm223,theta13,theta23,norm
       real*8  w12,w13,w14,w23,w24,w34
       real*8  dm212,dm234	
       real*8 theta12,theta14,theta24,theta34
       real*8 min_tot_th,min_tot_nosc,min_tot_data       
       real*8 pi,epsilon1, Ener
       real*8 int_min_data,int_min_nosc, int_min_th
       real*8 P_aver,LoE,Pmumu
       	real*8 min_th(40),min_nosc(40),min_data(40)
        real*8 min_chi2_pull,min_chi2,Prob_k(40)
        real*8 NOS_Minos(40),err_th(40),Ener_cent
        real*8 min_osc(40),min_mc(40),min_bf(40),
     c sigma(40),Enu(40),bin_min(40),bin_max(40),
     c mc_tot, osc_tot
        real*8 sigma_k

c i -> MinData20(i,2,40)      
c i=1 ->MINOS-POT-numu-minos-mas.dat
c i=2 ->MINOS-POT-numu-numubar.dat
c i=3 ->MINOS-best-fit.dat
c i=4 ->MINOS-bines.dat
c i=5 ->MINOS-bottom_sigma_limits.dat
c i=6 ->MINOS-data-points.dat
c i=7 ->MINOS-no-oscillation.dat
c i=8 -> MINOS-upper-sigma-limits.dat


        do i=1,39
        min_bf(i)=MinData20(3,2,i) 
        bin_min(i)=MinData20(4,1,i) 
        bin_max(i)=MinData20(4,2,i) 
        
        min_osc(i)=MinData20(6,2,i) 
        min_mc(i)=MinData20(7,2,i)  
        sigma(i)= abs(MinData20(8,2,i)-MinData20(5,2,i))/2.0
        sigma(i)=sigma(i)/min_osc(i)
        
c        write(*,*) sigma(i),sigma(i)/min_osc(i),i
    
         enddo



        
c            min_bf(:)=MinData20(3,2,:) 
c            bin_min(:)=MinData20(4,1,:) 
c            bin_max(:)=MinData20(4,2,:) 
c            
c            min_osc(:)=MinData20(6,2,:) 
c            min_mc(:)=MinData20(7,2,:)  
c            sigma(:)= abs(MinData20(8,2,:)-MinData20(5,2,:))/2.0
         
        
        
         pi=3.14159265
    
         
       dm212=  7.54d-5
       dm223=  2.480d-3
       dm234=  0.0d0
       theta12= 0.5873
       theta13= 0.1454
       theta23= 0.6642
       theta14= 0.0d0
       theta24= 0.0d0
       theta34= 0.0d0
! 	  norm=Y(10)
 	    norm=1.0d0 

       call umat(theta12,theta13,theta23,
     c		theta14,theta24,theta34)        
     
c         open(10,file='fakeMinos_NOS.dat')
c          open(11,file='RealMinos_NOS.dat')

         
         epsilon1=0.005d0
 
          do k=1,39
         Ener=bin_min(k)
         
c              if  (  k.eq.12 ) then 
c              Ener= 0.999* Ener
c              endif 
c              if  (  k.eq.11) then 
c              Ener= 0.996* Ener
c              endif 
c              if  (  k.eq.10 ) then 
c              Ener= 0.995* Ener
c              endif 
c             if  (  k.eq.9) then 
c              Ener= 0.99* Ener
c              endif          
c              if  (  k.eq.8) then 
c              Ener= 0.98* Ener
c              endif                
c             if  (  k.eq.7) then 
c              Ener= 0.97* Ener
c              endif           
c              if  (  k.eq.6) then 
c              Ener= 0.93 * Ener
c              endif                    
c              if  (  k.eq.5) then 
c              Ener= 0.80 * Ener
c              endif                     
c              if  (  k.eq.4) then 
c              Ener= 0.78 * Ener
c              endif                              
c              if  (  k.eq.3) then 
c              Ener= 0.90* Ener
c              endif                            
         
         P_aver=0.0d0
 
 	   do  while (Ener.le.bin_max(k))         
 	   
 	    LoE=735.0d0/Ener     
  
           w12=1.27d0*dm212*LoE 
           w23=1.27d0*dm223*LoE 
             w34=1.27d0*dm234*LoE 
             w13=w12+w23      
             w14=w13+w34      
             w24=w23+w34      
                       
            Pmumu= 1.0d0  -4.d0*       
     c        (cmm(1)*(sin(w12)**2)      
     c       +cmm(2)*(sin(w13)**2)      
     c       +cmm(3)*(sin(w23)**2)      
     c       +cmm(4)*(sin(w14)**2)        
     c       +cmm(5)*(sin(w24)**2)             
     c       +cmm(6)*(sin(w34)**2) )    
         P_aver=P_aver+Pmumu*(epsilon1)  
         Ener=Ener+epsilon1
  	   
  	   enddo
     
        P_aver=P_aver/(bin_max(k)-bin_min(k))
        
        Prob_k(k)=P_aver
        
         NOS_Minos(k)=  min_bf(k)/Prob_k(k)    ! min_osc(k)/Prob_k(k)       
          err_th(k)=abs(NOS_Minos(k)-min_mc(k)) /min_mc(k)       
          
c      c          print*,NOS_Minos(k),min_mc(k),err_th(k)
          
c  c          write(10,*)  k,NOS_Minos(k)
c  c          write(11,*)  k, min_mc(k)
c                   NOS_Minos(k)=min_mc(k)
          
         enddo 
         
       
c          close(10)
c          close(11)
     
    
        
        
        
       dm212=Y(1)
       dm223=1.0d0*Y(2)
       dm234=Y(3)
       theta12=Y(4)
       theta13=Y(5)
       theta23=Y(6)
       theta14=Y(7)
       theta24=Y(8)
       theta34=Y(9)
! 	norm=Y(10)
 	    norm=1.0d0 
         
         
      
c         dm212=8.0d-5
c         dm223=2.540d-3
c         dm234=0.0d0
c         theta12=0.5
c         theta13=0.1d0
c         theta23=pi/4.0d0
c         theta14=0.0d0
c         theta24=0.0d0
c         theta34=0.0d0
        

       call umat(theta12,theta13,theta23,
     c		theta14,theta24,theta34)               
     
     
     
        epsilon1=0.005d0

       min_tot_th=0.0d0
       min_tot_nosc=0.0d0
       min_tot_data=0.0d0        
       
        

         do k=1,39
c         NOS_Minos(k)=min_mc(k)

         
         
c        min_mc(k) 
        Ener=bin_min(k)
        int_min_th=0.0d0
        int_min_nosc=0.0d0
        int_min_data=0.0d0
        P_aver=0.0d0
        
          
           if  (  k.eq.12 ) then 
           Ener= 0.999* Ener
           endif 
           if  (  k.eq.11) then 
           Ener= 0.996* Ener
           endif 
           if  (  k.eq.10 ) then 
           Ener= 0.995* Ener
           endif 
          if  (  k.eq.9) then 
           Ener= 0.99* Ener
           endif          
           if  (  k.eq.8) then 
           Ener= 0.98* Ener
           endif                
          if  (  k.eq.7) then 
           Ener= 0.97* Ener
           endif           
           if  (  k.eq.6) then 
           Ener= 0.93 * Ener
           endif                  
           if  (  k.eq.5) then 
           Ener= 0.80 * Ener
           endif                     
           if  (  k.eq.4) then 
           Ener= 0.78 * Ener
           endif                              
           if  (  k.eq.3) then 
           Ener= 0.90* Ener
           endif                              
      
        
        Ener_cent=bin_min(k) + ((bin_max(k)-bin_min(k))/2.0d0) 	   

	    LoE=735.0d0/Ener_cent     
	    
 
          w12=1.27d0*dm212*LoE 
          w23=1.27d0*dm223*LoE 
            w34=1.27d0*dm234*LoE 
            w13=w12+w23      
            w14=w13+w34      
            w24=w23+w34              
           
           
               Prob_k(k) = 1.0d0  -4.d0*       
     c        (cmm(1)*(sin(w12)**2)      
     c       +cmm(2)*(sin(w13)**2)      
     c       +cmm(3)*(sin(w23)**2)      
     c       +cmm(4)*(sin(w14)**2)        
     c       +cmm(5)*(sin(w24)**2)             
     c       +cmm(6)*(sin(w34)**2) )    
            
            
        
        
      
	   do  while (Ener.le.bin_max(k))      
      
	    LoE=735.0d0/Ener     
 
          w12=1.27d0*dm212*LoE 
          w23=1.27d0*dm223*LoE 
            w34=1.27d0*dm234*LoE 
            w13=w12+w23      
            w14=w13+w34      
            w24=w23+w34      
           
           
           Pmumu= 1.0d0  -4.d0*       
     c        (cmm(1)*(sin(w12)**2)      
     c       +cmm(2)*(sin(w13)**2)      
     c       +cmm(3)*(sin(w23)**2)      
     c       +cmm(4)*(sin(w14)**2)        
     c       +cmm(5)*(sin(w24)**2)             
     c       +cmm(6)*(sin(w34)**2) )    
         
           int_min_data=int_min_data + min_osc(k)*(epsilon1)  

c             int_min_th= int_min_th + min_mc(k)*Pmumu*(epsilon1)   
            int_min_th= int_min_th +  NOS_Minos(k)*Pmumu*(epsilon1)   
       	
c             int_min_nosc= int_min_nosc + min_mc(k)*(epsilon1)        
             int_min_nosc= int_min_nosc +  NOS_Minos(k)*(epsilon1)        
        
        P_aver=P_aver+Pmumu*(epsilon1)  
        Ener=Ener+epsilon1
 	   
 	   enddo

    
       P_aver=P_aver/(bin_max(k)-bin_min(k))
       

c              min_th(k)=min_mc(k)*Prob_k(k)
c             min_th(k)=NOS_Minos(k)*Prob_k(k)

c                min_nosc(k)=NOS_Minos(k)
c                min_th(k)=NOS_Minos(k)*Prob_k(k)
c               min_data(k)=min_osc(k)    
      
      
       
              min_th(k)=int_min_th          !  / (bin_max(k)-bin_min(k))
              min_nosc(k)=int_min_nosc  ! /(bin_max(k)-bin_min(k))
              min_data(k)=int_min_data   !/(bin_max(k)-bin_min(k))
              Prob_k(k)=P_aver
         
         
         min_tot_th=min_tot_th+min_th(k)
         min_tot_nosc=min_tot_nosc+min_nosc(k)
         min_tot_data=min_tot_data+min_data(k)
         
c          write(*,*) , min_th(k),min_data(k), P_aver        
         
        enddo 
        
        
c         write(*,*)  min_tot_th,min_tot_nosc,min_tot_data
        
       min_chi2_pull=0.0d0
       
               do k=1,39
                 sigma_k= min_data(k) ! * sigma(k)  !min_th(k)  !   !*   
             
                min_chi2_pull=min_chi2_pull+
     c    (( min_th(k) -  min_data(k))  )**2   / (1.0*sigma_k)    
     
c               min_chi2_pull=min_chi2_pull+
c      c    (( min_th(k)-  min_data(k)) / min_th(k)  )**2
 
         
 
 
c             min_chi2_pull=min_chi2_pull+
c      c    (( min_th(k)-  min_data(k))     )**2     /(min_data(k)) 
     
c               min_chi2_pull=min_chi2_pull+
c      c   (   ( (min_data(k)/ min_nosc(k))-Prob_k(k) ) )**2
        
               enddo
       
c         call chis2_min_2020_SN(min_th,min_nosc,min_data,min_chi2_pull)
 
c          call chis2_min_2020_test(min_th,min_nosc,
c      c      min_data,min_chi2_pull)
  
c        call chis2_min_2020(min_th,min_nosc,min_data,min_chi2_pull)
       
c          min_chi2_pull=min_chi2_pull+
c      c  ((min_tot_th-min_tot_data)/min_tot_data)**2      
        
        min_chi2=min_chi2_pull
            print*, min_chi2
        
         
       end 
        


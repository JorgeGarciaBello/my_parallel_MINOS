        program minos_grid
        
        implicit none
        integer  i,j,points,k
        parameter (points=200 )
        
        real*8 chi_dayabay, chi_reno, chi_doubleCHOOZ
        real*8 chi_min_global,db_chi_min
        real*8 db_chi_square_spectral_analysis2_period
        real*8 DC_FUNC    
        real*8 reno_chi_square_spectral_analysis
        real*8 Y(13)      !Y=( t12, t13,t14,t23,t24,t34,d13,d24,d34,dm21,dm31,dm41)        
        real*8 delta_dm23
        real*8  var_dm23max,var_the13max,var_the23max
        real*8 var_the23min,var_the13min,var_dm23min
        real*8 sumth13, sumth23, sumdm23,  delta_th13,delta_th23
        real*8 var_th13(points), var_th23(points),var_dm23(points)
        real*8 chi2_grid(points,points,points)
        real*8 min_chi2,pi,sin2_th
     
        CHARACTER(30)  names

        
        names='30.minos2'
        
        
         call Read_MinData20                
c          call Read_MinNuMuAtmCV
c          call db_read_data()      ! Lee datos de Dayabay        
c          call readRENOData()    ! Lee datos de RENO   
         
         
         
       var_the23min=0.5
       var_the23max=1.0 

       var_the13min=-0.45
       var_the13max=0.45
 
       var_dm23max=9.0d-3
       var_dm23min=0.90d-3
       
       delta_th13=abs( var_the13max - var_the13min) /float(points-1)     
       delta_th23=abs(var_the23max - var_the23min)/float(points-1)
       delta_dm23=abs( var_dm23max - var_dm23min) /float(points-1)     

       
       
       sumth13=var_the13min!initial parameters
       sumth23=var_the23min
       sumdm23=var_dm23min
       	
       do i=1,points
       var_th13(i)=sumth13
       var_th23(i)=sumth23 
       
        var_dm23(i)=sumdm23
       
!        	print*,var_th13(i),var_th23(i),var_dm23(i)
       
       sumth13=sumth13+delta_th13
       sumth23=sumth23+delta_th23
       sumdm23=sumdm23+delta_dm23
       
       
       enddo
       
c        Y(1)=0.0d0*7.650d-5 ! dm2_12   
c        !Y(2)=2.5d-3 !dm2_23
c        Y(3)=0.0d0 !dm2_34
c        
c        Y(4)=0.0d0*0.584d0 !theta_12
c        !Y(5)=0.1  !theta_13
c       ! Y(6)=0.78  !theta_23
c        
c        Y(7)=0.0d0 !theta_14
c        Y(8)=0.0d0 !theta_24
c        Y(9)=0.0d0 !theta_34
c        Y(10)=0.0     
       
          Y(1) =7.54d-5
          Y(2) =2.480d-3
          Y(3) =0.0d0
          Y(4) =0.5873
          Y(5) =0.1454
          Y(6) =0.6642
          Y(7) =0.0d0
          Y(8) =0.0d0
          Y(9) =0.0d0
          Y(10)=0.0
         
!=======================================================       
         pi= 3.14159265
c          Y(2)=2.5d-3
c          Y(5)=0.0
c          Y(6)=pi/4.0d0             
c          call minos_2020_print( Y,min_chi2)
c          print*, min_chi2
!=======================================================       

         min_chi2=0.0d0
         do i=1,points
         do k=1,points
          
         
         Y(2)=var_dm23(i)
c          Y(5)=0.0d0
         Y(6)=var_th23(k)        
!           call daya_bay_cov(Y,db_chi_min)
!           chi2_grid(i,j,k) =db_chi_min
c           call renoChi2(Y,chi_reno)
               call minos_2020( Y,min_chi2)
c               call minos_2020_old( Y,min_chi2)
c               call minos_2020_plus( Y,min_chi2)


            chi2_grid(i,1,k) =min_chi2
         
c           write(*, '(4F20.8)') Y(2),Y(5),Y(6),chi2_grid(i,1,k) 
                 
c            write(*, *) Y(2),Y(6),chi2_grid(i,1,k) 
  
        enddo
        enddo    
            

! open(10,file='g.th23th13dm.50'//names)
       open(30,file='g.th23dm.'//names)
           
            
            
        do k=1,points
        do i=1,points
        
        sin2_th=sin(var_th23(k))**2

       write(30, '(4F20.8)') sin2_th, var_dm23(i),chi2_grid(i,1,k)
c          write(30, *)  sin2_th, var_dm23(i),chi2_grid(i,1,k)
c          write(*, *)  sin2_th, var_dm23(i),chi2_grid(i,1,k)

       
        enddo
        enddo
            
        close(30)            
            
            
       
       

       end  

c        program Read_MinData20
       
       subroutine Read_MinData20
       
c                        Read_MinNuMuDomBeamNFid
       
       implicit none
       
       include 'MINOS_20.inc'
c        
       integer LNBLNK
       integer i,j,k
c        real*8 MinData20(8,2,40)
        character(150) file_name(10),dummyCh,paths

       
         CALL chdir("/Users/roa/workfolder/sources_20/minos/")
        
        paths= "MINOS2/MINOS-data/ "
         paths=   paths(:LNBLNK(paths)) 

c         nflx_max= home(:LNBLNK(home)) //'skatm-solmax-sumphi.d'       
       
       

!descripcion del arreglo, 
!todos son arreglos de 2x25

c i -> MinData20(i,2,2)      
c i=1 ->MINOS-POT-numu-minos-mas.dat
c i=2 ->MINOS-POT-numu-numubar.dat
c i=3 ->MINOS-best-fit.dat
c i=4 ->MINOS-bines.dat
c i=5 ->MINOS-bottom_sigma_limits.dat
c i=6 ->MINOS-data-points.dat
c i=7 ->MINOS-no-oscillation.dat
c i=8 -> MINOS-upper-sigma-limits.dat


 
c        file_name(1)=paths// 'MINOS-numu-atm'
       
       file_name(1)=
     c  'MINOS-POT-numu-minos-mas.dat'       
         file_name(1)= paths(:LNBLNK(paths))// file_name(1)

       file_name(2)=
     c   'MINOS-POT-numu-numubar.dat'
              file_name(2)=paths(:LNBLNK(paths))//file_name(2)
       file_name(3)=
     c  'MINOS-best-fit.dat'
              file_name(3)=paths(:LNBLNK(paths))//file_name(3)     
       file_name(4)=
     c  'MINOS-bines.dat'
              file_name(4)=paths(:LNBLNK(paths))//file_name(4)     
       file_name(5)=
     c  'MINOS-bottom_sigma_limits.dat'
              file_name(5)=paths(:LNBLNK(paths))//file_name(5)     
       file_name(6)=
     c  'MINOS-data-points.dat'
              file_name(6)=paths(:LNBLNK(paths))//file_name(6)     
       file_name(7)=
     c  'MINOS-no-oscillation.dat'
              file_name(7)=paths(:LNBLNK(paths))//file_name(7) 
              
       file_name(8)=              
     c  'MINOS-upper-sigma-limits.dat'              
              file_name(8)=paths(:LNBLNK(paths))//file_name(8)               
       
       do i=1,8
       
         open(i,file= file_name(i))
         
         do j=1,39
         read(i,*) MinData20(i,1,j),MinData20(i,2,j)
 
  
 
         
         enddo
        
c        write(*,*)  file_name(i)       
           close(i)
          enddo
          
          
c             i=8
c           do j=1,39
c           
c            write(*,*)   MinData20(i,1,j),MinData20(i,2,j)
c           
c           enddo


      CALL chdir("/Users/roa/workfolder/sources_20/minos/sources") 
       
       
       end 

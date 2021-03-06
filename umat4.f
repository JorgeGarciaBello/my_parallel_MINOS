	subroutine umat(th12,th13,th23,th14,th24,th34)
!c  returns coefficients cij for calculating oscillation ij
!c  returns u(i,j) in common unt
!c
	implicit none
	include 'UMAT4.inc'


	real*8  u(4,4),s13
	real*8  s12,c12,c13,s23,c23,s14,c14,s24,c24,s34,c34
	real*8  th12,th13,th23,th14,th24,th34
	

!c calculate sin's and cos's
    	s12=sin(th12)
      	c12=cos(th12)
      	s13=sin(th13)
      	c13=cos(th13)
      	s23=sin(th23)
      	c23=cos(th23)
      	s14=sin(th14)
      	c14=cos(th14)
      	s24=sin(th24)
      	c24=cos(th24)
      	s34=sin(th34)
      	c34=cos(th34)      
            
      
!   calculate u's


      	u(1,1)=C12* C13 *C14
      	u(1,2)=C13* C14 *S12
      	u(1,3)=C14* S13
      	u(1,4)=S14
      	 
      	u(2,1)=C24*(-C23*S12 - C12*S13*S23) - C12*C13*S14*S24
      	u(2,2)=C24* (C12*C23 - S12*S13*S23) -    	
     c             C13*S12*S14*S24
      	u(2,3)=C13*C24*S23 - S13*S14*S24
      	u(2,4)= C14*S24    

      	u(3,1)=C34 *(-C12 *C23* S13 + S12* S23) - 	
     c       (C12* C13* C24* S14* S34) - 			
     c        ((-C23* S12 - C12* S13* S23)* S24 *S34)
      	u(3,2)=C34* (-C23* S12* S13 - C12* S23) - 	
     c       (C13* C24* S12* S14* S34) -  			
     c        ((C12* C23 - S12* S13* S23)* S24* S34)
     
      	u(3,3)=(C13* C23* C34) - (C24* S13* S14* S34) - 
     c            (C13* S23* S24* S34)
      	u(3,4)=C14* C24* S34
    
      	u(4,1)=(-C12* C13* C24* C34* S14) -  		
     c       (C34* (-C23*S12 - C12* S13 *S23)* S24) -	
     c       ((-C12* C23* S13 + S12* S23)* S34)
      
      	u(4,2)=(-C13 *C24* C34 *S12 *S14) - 		
     c       (C34 *(C12* C23 - S12 *S13 *S23) *S24) - 	
     c       ((-C23 *S12 *S13 - C12* S23) *S34)
   
      	u(4,3)=(-C24 *C34 *S13 *S14) -			
     c      (C13 *C34 *S23* S24) - (C13* C23* S34)

	u(4,4)=C14* C24 *C34
   

!c construct coefficients
 

      
      
       	cee(1)=u(1,1)*u(1,2)*u(1,1)*u(1,2)          
       	cee(2)=u(1,1)*u(1,3)*u(1,1)*u(1,3)        
       	cee(3)=u(1,2)*u(1,3)*u(1,2)*u(1,3)          
       	cee(4)=u(1,1)*u(1,4)*u(1,1)*u(1,4)            
       	cee(5)=u(1,2)*u(1,4)*u(1,2)*u(1,4)          
       	cee(6)=u(1,3)*u(1,4)*u(1,3)*u(1,4)       
       
      
      

       	cem(1)=u(1,1)*u(1,2)*u(2,1)*u(2,2)          
       	cem(2)=u(1,1)*u(1,3)*u(2,1)*u(2,3)        
       	cem(3)=u(1,2)*u(1,3)*u(2,2)*u(2,3)          
       	cem(4)=u(1,1)*u(1,4)*u(2,1)*u(2,4)            
       	cem(5)=u(1,2)*u(1,4)*u(2,2)*u(2,4)          
       	cem(6)=u(1,3)*u(1,4)*u(2,3)*u(2,4)             
      
      



      
       	cet(1)=u(1,1)*u(1,2)*u(3,1)*u(3,2)          
       	cet(2)=u(1,1)*u(1,3)*u(3,1)*u(3,3)        
       	cet(3)=u(1,2)*u(1,3)*u(3,2)*u(3,3)          
       	cet(4)=u(1,1)*u(1,4)*u(3,1)*u(3,4)            
       	cet(5)=u(1,2)*u(1,4)*u(3,2)*u(3,4)          
       	cet(6)=u(1,3)*u(1,4)*u(3,3)*u(3,4)

      
       	ces(1)=u(1,1)*u(1,2)*u(4,1)*u(4,2)          
       	ces(2)=u(1,1)*u(1,3)*u(4,1)*u(4,3)        
       	ces(3)=u(1,2)*u(1,3)*u(4,2)*u(4,3)          
       	ces(4)=u(1,1)*u(1,4)*u(4,1)*u(4,4)            
       	ces(5)=u(1,2)*u(1,4)*u(4,2)*u(4,4)          
       	ces(6)=u(1,3)*u(1,4)*u(4,3)*u(4,4)
      


      
      
       	cmm(1)=u(2,1)*u(2,2)*u(2,1)*u(2,2)          
       	cmm(2)=u(2,1)*u(2,3)*u(2,1)*u(2,3)        
       	cmm(3)=u(2,2)*u(2,3)*u(2,2)*u(2,3)          
       	cmm(4)=u(2,1)*u(2,4)*u(2,1)*u(2,4)            
       	cmm(5)=u(2,2)*u(2,4)*u(2,2)*u(2,4)          
       	cmm(6)=u(2,3)*u(2,4)*u(2,3)*u(2,4)             
            
      

      
      

       	cmt(1)=u(2,1)*u(2,2)*u(3,1)*u(3,2)          
       	cmt(2)=u(2,1)*u(2,3)*u(3,1)*u(3,3)        
       	cmt(3)=u(2,2)*u(2,3)*u(3,2)*u(3,3)          
       	cmt(4)=u(2,1)*u(2,4)*u(3,1)*u(3,4)            
       	cmt(5)=u(2,2)*u(2,4)*u(3,2)*u(3,4)          
       	cmt(6)=u(2,3)*u(2,4)*u(3,3)*u(3,4)             


       	cms(1)=u(2,1)*u(2,2)*u(4,1)*u(4,2)          
       	cms(2)=u(2,1)*u(2,3)*u(4,1)*u(4,3)        
       	cms(3)=u(2,2)*u(2,3)*u(4,2)*u(4,3)          
       	cms(4)=u(2,1)*u(2,4)*u(4,1)*u(4,4)            
       	cms(5)=u(2,2)*u(2,4)*u(4,2)*u(4,4)          
       	cms(6)=u(2,3)*u(2,4)*u(4,3)*u(4,4)             
                  

      


      
       	ctt(1)=u(3,1)*u(3,2)*u(3,1)*u(3,2)          
       	ctt(2)=u(3,1)*u(3,3)*u(3,1)*u(3,3)        
       	ctt(3)=u(3,2)*u(3,3)*u(3,2)*u(3,3)          
       	ctt(4)=u(3,1)*u(3,4)*u(3,1)*u(3,4)            
       	ctt(5)=u(3,2)*u(3,4)*u(3,2)*u(3,4)          
       	ctt(6)=u(3,3)*u(3,4)*u(3,3)*u(3,4)      
      
       	cts(1)=u(3,1)*u(3,2)*u(4,1)*u(4,2)          
       	cts(2)=u(3,1)*u(3,3)*u(4,1)*u(4,3)        
       	cts(3)=u(3,2)*u(3,3)*u(4,2)*u(4,3)          
       	cts(4)=u(3,1)*u(3,4)*u(4,1)*u(4,4)            
       	cts(5)=u(3,2)*u(3,4)*u(4,2)*u(4,4)          
       	cts(6)=u(3,3)*u(3,4)*u(4,3)*u(4,4)      

       	css(1)=u(4,1)*u(4,2)*u(4,1)*u(4,2)          
       	css(2)=u(4,1)*u(4,3)*u(4,1)*u(4,3)        
       	css(3)=u(4,2)*u(4,3)*u(4,2)*u(4,3)          
       	css(4)=u(4,1)*u(4,4)*u(4,1)*u(4,4)            
       	css(5)=u(4,2)*u(4,4)*u(4,2)*u(4,4)          
       	css(6)=u(4,3)*u(4,4)*u(4,3)*u(4,4)      
      

      return
      end

      subroutine lasinv(n,ss,rho,ia,is,thr,ww,wwi,jerr)                     
      real ss(n,n),rho(n),ww(n,n),wwi(n,n)                                  
      real, dimension (:,:), allocatable :: vv,xs                            
      real, dimension (:), allocatable :: s,x,z,ws                           
      integer, dimension (:), allocatable :: mm                              
      nm1=n-1                                                                
      allocate(vv(1:nm1,1:nm1),stat=jerr)                                    
      ierr=0                                                                 
      if(ia.eq.0) allocate(xs(1:nm1,1:n),stat=ierr)                          
      jerr=jerr+ierr                                                         
      allocate(s(1:nm1),stat=ierr)                                          
      jerr=jerr+ierr                                                        
      allocate(x(1:nm1),stat=ierr)                                          
      jerr=jerr+ierr                                                        
      allocate(z(1:nm1),stat=ierr)                                          
      jerr=jerr+ierr                                                        
      allocate(mm(1:nm1),stat=ierr)                                         
      jerr=jerr+ierr                                                        
      if(ia .ne. 0)goto 10021                                               
      allocate(ws(1:n),stat=ierr)                                           
      jerr=jerr+ierr                                                        
10021 continue                                                              
      if(jerr.ne.0) return                                                  
      shr=0.0                                                               
10030 do 10031 j=1,n                                                        
10040 do 10041 k=1,n                                                        
      if(j.eq.k)goto 10041                                                  
      shr=shr+abs(ss(j,k))                                                  
10041 continue                                                              
10042 continue                                                              
10031 continue                                                              
10032 continue                                                              
      if(shr .ne. 0.0)goto 10061                                            
      ww=0.0                                                                
      wwi=0.0                                                               
10070 do 10071 j=1,n                                                        
      ww(j,j)=ss(j,j)+rho(j)                                                
      wwi(j,j)=1.0/ww(j,j)                                                  
10071 continue                                                              
10072 continue                                                              
      return                                                                
10061 continue                                                              
      shr=thr*shr/nm1                                                       
      if(ia .eq. 0)goto 10091                                               
      if(is.eq.0) wwi=0.0                                                   
10100 do 10101 m=1,n                                                        
      call setup(m,n,ss,rho,ss,vv,s,z)                                      
      l=0                                                                   
10110 do 10111 j=1,n                                                        
      if(j.eq.m)goto 10111                                                  
      l=l+1                                                                 
      x(l)=wwi(j,m)                                                         
10111 continue                                                              
10112 continue                                                              
      call lasso7(z,nm1,vv,s,shr/sum(abs(vv)),x)                            
      call cleanup(m,ia,n,wwi,vv,x,s,z,mm)                                  
10101 continue                                                              
10102 continue                                                              
      return                                                                
10091 continue                                                              
      if(is .ne. 0)goto 10131                                               
      ww=ss                                                                 
      xs=0.0                                                                
      goto 10141                                                            
10131 continue                                                              
10150 do 10151 j=1,n                                                        
      xjj=-2.0*wwi(j,j)                                                     
      l=0                                                                   
10160 do 10161 k=1,n                                                        
      if(k.eq.j)goto 10161                                                  
      l=l+1                                                                 
      xs(l,j)=wwi(k,j)/xjj                                                  
10161 continue                                                              
10162 continue                                                              
10151 continue                                                              
10152 continue                                                              
10141 continue                                                              
10121 continue                                                              
10170 do 10171 j=1,n                                                        
      ww(j,j)=ss(j,j)+rho(j)                                                
10171 continue                                                              
10172 continue                                                              
10180 continue                                                              
10181 continue                                                              
      dlx=0.0                                                               
10190 do 10191 m=1,n                                                        
      x=xs(:,m)                                                             
      ws=ww(:,m)                                                            
      call setup(m,n,ss,rho,ww,vv,s,z)                                      
      call lasso7(z,nm1,vv,s,shr/sum(abs(vv)),x)                            
      call cleanup(m,ia,n,ww,vv,x,s,z,mm)                                   
      dlx=max(dlx,sum(abs(ww(:,m)-ws)))                                     
      xs(:,m)=x                                                             
10191 continue                                                              
10192 continue                                                              
      if(dlx.lt.shr)goto 10182                                              
      goto 10181                                                            
10182 continue                                                              
      call inv(n,ww,xs,wwi)                                                 
      return                                                                
      end                                                                   
      subroutine setup(m,n,ss,rho,ww,vv,s,r)                                
      real ss(n,n),rho(n),ww(n,n),vv(n-1,n-1),s(n-1),r(n-1)                 
      l=0                                                                   
10200 do 10201 j=1,n                                                        
      if(j.eq.m)goto 10201                                                  
      l=l+1                                                                 
      r(l)=rho(j)                                                           
      s(l)=ss(j,m)                                                          
      i=0                                                                   
10210 do 10211 k=1,n                                                        
      if(k.eq.m)goto 10211                                                  
      i=i+1                                                                 
      vv(i,l)=2.0*ww(k,j)                                                   
10211 continue                                                              
10212 continue                                                              
10201 continue                                                              
10202 continue                                                              
      return                                                                
      end                                                                   
      subroutine lasso7(rho,n,vv,s,thr,x)                                   
      real rho(n),vv(n,n),s(n),x(n)                                         
10220 do 10221 j=1,n                                                        
      s(j)=s(j)-dot_product(vv(:,j),x)                                      
10221 continue                                                              
10222 continue                                                              
10230 continue                                                              
10231 continue                                                              
      dlx=0.0                                                               
10240 do 10241 j=1,n                                                        
      xj=x(j)                                                               
      x(j)=0.0                                                              
      t=s(j)+vv(j,j)*xj                                                     
      if (abs(t)-rho(j).gt.0.0) x(j)=sign(abs(t)-rho(j),t)/vv(j,j)          
      if(x(j).eq.xj)goto 10241                                              
      del=x(j)-xj                                                           
      dlx=max(dlx,abs(del))                                                 
      s=s-del*vv(:,j)                                                       
10241 continue                                                              
10242 continue                                                              
      if(dlx.lt.thr)goto 10232                                              
      goto 10231                                                            
10232 continue                                                              
      return                                                                
      end                                                                   
      subroutine cleanup(m,ia,n,ww,vv,x,s,z,mm)                             
      real ww(n,n),vv(n-1,n-1),x(n-1),s(n-1),z(n-1)                         
      integer mm(n-1)                                                       
      if(ia .ne. 0)goto 10261                                               
      call fatmul(n-1,vv,x,s,z,mm)                                          
      goto 10271                                                            
10261 continue                                                              
      s=x                                                                   
10271 continue                                                              
10251 continue                                                              
      l=0                                                                   
10280 do 10281 j=1,n                                                        
      if(j.eq.m)goto 10281                                                  
      l=l+1                                                                 
      ww(j,m)=s(l)                                                          
      ww(m,j)=ww(j,m)                                                       
10281 continue                                                              
10282 continue                                                              
      return                                                                
      end                                                                   
      subroutine fatmul(n,vv,x,s,z,m)                                       
      parameter(fac=0.2)                                                    
      real vv(n,n),x(n),s(n),z(n)                                           
      integer m(n)                                                          
      l=0                                                                   
10290 do 10291 j=1,n                                                        
      if(x(j).eq.0.0)goto 10291                                             
      l=l+1                                                                 
      m(l)=j                                                                
      z(l)=x(j)                                                             
10291 continue                                                              
10292 continue                                                              
      if(l .le. int(fac*n))goto 10311                                       
      s=matmul(vv,x)                                                        
      goto 10321                                                            
10311 continue                                                              
10330 do 10331 j=1,n                                                        
      s(j)=dot_product(vv(j,m(1:l)),z(1:l))                                 
10331 continue                                                              
10332 continue                                                              
10321 continue                                                              
10301 continue                                                              
      return                                                                
      end                                                                   
      subroutine inv(n,ww,xs,wwi)                                           
      real ww(n,n),xs(n-1,n),wwi(n,n)                                       
      nm1=n-1                                                               
      xs=-2.0*xs                                                            
      wwi(1,1)=1.0/(ww(1,1)+dot_product(xs(:,1),ww(2:n,1)))                 
      wwi(2:n,1)=wwi(1,1)*xs(:,1)                                           
      wwi(n,n)=1.0/(ww(n,n)+dot_product(xs(:,n),ww(1:nm1,n)))               
      wwi(1:nm1,n)=wwi(n,n)*xs(:,n)                                         
10340 do 10341 j=2,nm1                                                      
      jm1=j-1                                                               
      jp1=j+1                                                               
      wwi(j,j)=1.0/(ww(j,j)+dot_product(xs(1:jm1,j),ww(1:jm1,j))  +dot_product(xs(j:nm1,j),ww(jp1:n,j)))
      wwi(1:jm1,j)=wwi(j,j)*xs(1:jm1,j)                                     
      wwi(jp1:n,j)=wwi(j,j)*xs(j:nm1,j)                                     
10341 continue                                                              
10342 continue                                                              
      return                                                                
      end                                                                   

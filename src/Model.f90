SUBROUTINE Model
   use param_Y
   use param_rheo
   use prec
   implicit none

   integer, parameter :: nlay=4
   integer :: i,j,ncc, ilay, i0
   real(kind=dp) :: dr
   real(kind=dp) :: eta_i, eta_0
   integer, dimension(nlay) :: n_layer(nlay)   
   real(kind=dp), dimension(nlay) :: Rlim, rho_lay, K_lay, mu_lay, vp_lay, vs_lay, eta_lay

!--Define constants and allocate
   O=real(ord)
   pi=acos(-1._dp)
   GG=6.67408e-11_dp

!--Define interior layer parameters ------

Rlim(4)=Rtot
rho_lay(4)=920._dp
mu_lay(4)=3.3e9_dp
K_lay(4)=10e9_dp
eta_lay(4)=etaice

Rlim(3)=Rlim(4)-bice
rho_lay(3)=rhoocean
mu_lay(3)=0.e0_dp
K_lay(3)=2e9_dp
eta_lay(3)=0.e0_dp

Rlim(2)=Rmant
rho_lay(2)=rhomant
mu_lay(2)=70e9_dp
K_lay(2)=200.e9_dp
eta_lay(2)=1e19_dp

rho_lay(1)=rhocore
mu_lay(1)=0.e0_dp
K_lay(1)=200.e9_dp
eta_lay(1)=0._dp
Rlim(1)=Rcore

write(*,*) Rlim(1), Rlim(2), Rlim(3), Rlim(4)
write(*,*) rho_lay(1), rho_lay(2), rho_lay(3), rho_lay(4)


n_layer(1)=10
n_layer(2)=64
n_layer(3)=10
n_layer(4)=10

ncore=n_layer(1)+1
nmant=ncore+n_layer(2)+1
ncrust=nmant+n_layer(3)+1

vs_lay=sqrt(mu_lay/rho_lay)
vp_lay=sqrt((K_lay+4._dp/3._dp*mu_lay)/rho_lay)


n=0
do i=1,nlay 
n=n+n_layer(i)+1
enddo 

 r(1)=0._dp
 rho(1)=rho_lay(1)
 vp(1)=vp_lay(1)
 vs(1)=vs_lay(1)
 eta(1)=eta_lay(1)
 mu(1)=mu_lay(1)
 Ka(1)=K_lay(1)

 dr=Rlim(1)/n_layer(1)

 do i=2,n_layer(1)+1
  r(i)=r(i-1)+dr
  rho(i)=rho_lay(1)
  vp(i)=vp_lay(1)
  vs(i)=vs_lay(1)
  eta(i)=eta_lay(1)
  mu(i)=mu_lay(1)
  Ka(i)=K_lay(1)
 enddo


ncc=n_layer(1)+1

 do ilay=2,nlay
  ncc=ncc+1
  r(ncc)=R(ncc-1)
  rho(ncc)=rho_lay(ilay)
  vp(ncc)=vp_lay(ilay)
  vs(ncc)=vs_lay(ilay)
  eta(ncc)=eta_lay(ilay)
  mu(ncc)=mu_lay(ilay)
  Ka(ncc)=K_lay(ilay)
  dr=(Rlim(ilay)-Rlim(ilay-1))/n_layer(ilay)
  
 do i=2,n_layer(ilay)+1
   ncc=ncc+1  
   r(ncc)=r(ncc-1)+dr
   rho(ncc)=rho_lay(ilay)
   vp(ncc)=vp_lay(ilay)
   vs(ncc)=vs_lay(ilay)
   eta(ncc)=eta_lay(ilay)
   mu(ncc)=mu_lay(ilay)
   Ka(ncc)=K_lay(ilay)
 enddo
enddo


!phi_melt=0.d0 !in %
phi_melt=phi_melt/1000.
cmelt=67._dp/15._dp
alphamelt=30._dp


gamma=5.

open(601, file=file_eta)
open(701, file=file_rmelt)
do i=ncore+1,nmant
read(601,*) r(i), eta(i)
read(701,*) r(i), rmelt(i), delta_eta(i), zeta_eta(i), phi_crit(i), delta_mu(i), zeta_mu(i), &
            phi_crit_mu(i), delta_K(i), zeta_K(i), phi_crit_K(i)
if (melt==1) then
errorfunc=(1-zeta_eta(i))*erf(((sqrt(pi))/(2*(1-zeta_eta(i))))*((1-phi_melt)/(1-phi_crit(i))) &
            *(1+((1-phi_melt)/(1-phi_crit(i)))**gamma))
errorfunc_mu=(1-zeta_mu(i))*erf(((sqrt(pi))/(2*(1-zeta_mu(i))))*((1-phi_melt)/(1-phi_crit_mu(i))) &
                *(1+((1-phi_melt)/(1-phi_crit_mu(i)))**gamma))
errorfunc_K=(1-zeta_K(i))*erf(((sqrt(pi))/(2*(1-zeta_K(i))))*((1-phi_melt)/(1-phi_crit_K(i))) &
                *(1+((1-phi_melt)/(1-phi_crit_K(i)))**gamma))

etamant_ast=((1+((1-phi_melt)/(1-phi_crit(i)))**(delta_eta(i)))/(1-errorfunc)**(2.5*(1-phi_crit(i))))
mu_mant_ast=1e1*(((1+((1-phi_melt)/(1-phi_crit_mu(i)))**(delta_mu(i)))/(1-errorfunc_mu)**(2.5*(1-phi_crit_mu(i)))))
K_mant_ast=1e9*(((1+((1-phi_melt)/(1-phi_crit_K(i)))**(delta_K(i)))/(1-errorfunc_K)**(2.5*(1-phi_crit_K(i)))))
if (rmelt(i) .ne. 0) then
eta(i)=etamant_ast
mu(i)=mu_mant_ast
Ka(i)=K_mant_ast
vs(i)=sqrt(mu(i)/rho(i))
vp(i)=sqrt((Ka(i)+4._dp/3._dp*mu(i))/rho(i))
endif
endif
enddo

!Do i=1,n
!write(*,*) r(i), rho(i), vp(i), vs(i), eta(i), mu(i), Ka(i), rmelt(i)
!Enddo

END SUBROUTINE Model

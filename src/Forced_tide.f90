!**************************************************************
!Calcul des réponses au potentiel de marée
!d'apres Takeushi et Saito (1972)
!Mathieu Le Feuvre / Gabriel Tobie
!**************************************************************
PROGRAM Forced_tide
   use param_Y
   use param_rheo
   use prec
   use interf

   implicit none

   integer :: i, n_in, iecc, iecc_max, imelt_min, imelt_max, imelt_step
   real :: ecc_in, bast
   real(kind=dp):: w_in, k2_out, Dt_out, Q_out, h2_out, l2_out, Ts_in, rtop
   real(kind=dp),allocatable:: r_in(:), rho_in(:), vp_in(:), vs_in(:), eta_in(:)
   real(kind=dp), allocatable :: mu_in(:), Ka_in(:)
   logical :: ecc_loop, map

   namelist /input/ bice, etaice, ecc_in, ecc_loop, bast, rtop, imelt_min, imelt_max, imelt_step, map, hrad, &
                    file_eta, file_rmelt

 pi=acos(-1._dp)

 ecc_loop=.false.
 map=.false.
 hrad=0.
 open(10,file='input.nml')
 read(10,nml=input)
 close(10)

 open(21,file='yi-real.dat')
 open(22,file='yi-imag.dat')
 open(301,file='k2_Q_Ptide.txt')
 open(302,file='Hmu_htide.txt')
 open(402,file='Melt/Melt_ecc_Ptot_Pmant_bice_etaice_bast.dat')
 open(500,file='Melt/melt_ptide_pmant_pice_phaselag.dat')
 open(550,file='Melt/melt_y5_y1_y3_k2_h2_l2.dat')
 open(600,file='Melt/R_rho_eta_mu_y1_y5_htide_phi_21Gyr_phi.dat')

!- Mass and radius -
Rtot=1561.e3_dp

! iron core
rhocore=5150.d0
Rcore=600.d3

! silicate mantle properties
Rmant=1421.d3
rhomant=3500.d0
vsmant=sqrt(70.D9/rhomant)

! ocean properties
rhoocean=1000.d0

!melt model ; 0=no melt model ; 1 = melt model
melt=1

! viscoelastic model assumption: 0=purely elastic; 1=Maxwell; 2=Andrade
visco=3
alpha_gam=0.30_dp

! orbital eccentricity
ecc=ecc_in
iecc_max=1
if (ecc_loop) iecc_max=302

! Tidal period
Ts_in=3.55
Ts_in=Ts_in*(24._dp*3600._dp)
w_in=2._dp*pi/Ts_in

phi_melt=0.d0 !in pour mille
do imelt=imelt_min,imelt_max,imelt_step
phi_melt=phi_melt+imelt
write(*,*) 'phi_melt=', phi_melt
rewind(601)
rewind(701)

Call Model

allocate(r_in(n),rho_in(n), vp_in(n), vs_in(n), eta_in(n))

r_in=r
rho_in=rho
vp_in=vp
vs_in=vs
eta_in=eta
mu_in=mu
Ka_in=Ka
n_in=n

do iecc=1,iecc_max,5
if (ecc_loop) then
ecc=dble(iecc)*0.05/302
write(*,*) 'ecc=', ecc
endif

call k2_Dt(r_in,rho_in,vp_in,vs_in,eta_in,mu_in,Ka_in,n_in,w_in,k2_out,Dt_out,Q_out, h2_out, l2_out)

write(21,*) n, Ts_in/24.d0/3600.d0
write(22,*) n, Ts_in/24.d0/3600.d0

write(*,*) 'eta_mantle/(10^18 Pa.s)=',(etamant/1.D18), vsmant, ' >  k2=', k2_out, 'h2=', h2_out, 'Q_out=', Q_out, 'l2=', l2_out,&
'dt=', Dt_out/3600./24., 'Eglob=', 21._dp/2._dp*k2_out/Q_out*(w_in*R_in(n_in))**5*(ecc)**2/GG/1.D9
write(301,*) etamant, vsmant, k2_out, Q_out, 21._dp/2._dp*k2_out/Q_out*(w_in*R_in(n_in))**5*(ecc)**2/GG/1.D9
write(302,*) etamant, vsmant, k2_out, Q_out, 21._dp/2._dp*k2_out/Q_out*(w_in*R_in(n_in))**5*(ecc)**2/GG/1.D9

write(500,*) phi_melt, Ptot/1.d9, Pmant/1.d9, Pice/1.d9, (REAL(Yi(5,nc)))/(imag(Yi(5,nc))),&
(real(Yi(1,nc)))/(imag(Yi(1,nc))), (real(Yi(3,nc)))/(imag(Yi(3,nc))), bice, etaice, bast
write(550,*) phi_melt, REAL(Yi(5,nc)), imag(Yi(5,nc)), real(Yi(1,nc)), imag(Yi(1,nc)), real(Yi(3,nc)), imag(Yi(3,nc)), &
k2_out, h2_out, l2_out
write(402,*) phi_melt, ecc, Ptot/1.d9, Pmant/1.d9, bice, etaice, bast, (phi_melt*(4*pi*rtop**3 - 4*pi*(rtop-bast)**3)/3)/1e9

do i=1,n
write(600,*) R(i), rho(i), eta(i),mu(i), yi(1,i), yi(5,i), real(htide(i)), phi_melt
enddo

write(*,*) ' '

if (map) call Tidal_heating_map_output

enddo

deallocate(r_in, rho_in, vp_in, vs_in, eta_in)
enddo

close(21)
close(22)
close(301)

END PROGRAM Forced_tide

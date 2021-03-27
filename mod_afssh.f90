Module mod_afssh
!! Hammes-Schiffer, Tully, JCP 101, 4657 (1994)
implicit none
real*8, parameter :: clight=2.99792458D10,av=6.0221367D23,hbar=1.05457266D-34
real*8, parameter :: kb=1.3806503d-23
real*8, parameter :: pascal_to_atm=9.86923267d-6,kcal_to_J=6.95d-21
real*8, parameter :: amu2kg=1.66053892d-27
real*8, parameter :: au2kg=9.10938291d-31,au2J=4.35974417d-18,au2m=5.2917721092d-11,au2s=2.418884326505d-17
real*8, parameter :: q_el=2.307075d-28
complex*16,parameter :: iota = (0.d0,1.d0)
real*8 pi,wave_to_J

!! Potential
real*8 A_pot,B_pot,C_pot,D_pot
real*8,allocatable:: omg(:),x_min(:)
real*8 k_coup(3)
real*8,allocatable :: mass(:)
real*8 mass_qm
real*8 omg_B,gamma_B,temperature,beta

!! On the fly diagonalization
real*8,allocatable :: KE_DVR(:,:),r_grid(:)
real*8 r_min,r_max,r_del
real*8 r_exp,r_cutoff

!! Output/Input
real*8 k2_val
real*8 cnt_frust,cnt_collapse,cnt_init,cnt_term
real*8,allocatable :: pop(:,:),pop_surf(:,:),pop_amp(:,:)

!! Classical
integer nclass,idistribution
real*8,allocatable :: x(:),v(:),acc(:)
real*8,allocatable :: x_old(:),v_old(:),acc_old(:),x_hop(:)
real*8 tim_hop
integer iforward,flag_terminate,flag_frust,flag_reactant,flag_hop
complex*16,allocatable :: delr(:,:),delp(:,:)
complex*16,allocatable :: delr_old(:,:),delp_old(:,:)

!! Quantum
integer nquant,state,nbasis,state_tentative
integer state_old
real*8,allocatable :: si_adiab(:,:),V_k(:),d_ij(:,:,:),vdotd(:,:),V_k_old(:)
real*8,allocatable :: Hamil_diab(:,:),delH_dels(:,:,:),delH_dels_ad(:,:,:)
real*8,allocatable :: pot(:,:),force(:,:),force_old(:,:)
complex*16,allocatable :: ci(:),ci_old(:)
real*8,allocatable :: si_adiab_prev(:,:)
complex*16,allocatable :: mat(:,:),mat_adiab(:,:)
real*8,allocatable :: hop_prob(:),W_overlap(:,:),hop_prob_net(:)

!! Evolution
integer n_traj,nsteps,nsteps_eq,nstep_write,iwrite,nstep_avg
real*8 dtc,total_time,curr_time,traj_num,tim_eq
real*8 energy,pot_en,KE_en,energy_cutoff,energy_old
real*8 ensq_avg,en_avg
integer nst_av
integer ihop,icollapse,iterminate,ifriction,iaverage
real*8 prob_tun

!! Parallelization
integer iparallel,iflow,iproc,iwait,ifolder

!! Misc
integer nold,cnt_rate
real*8 tim_tot,tim_ev_cl,tim_diag,tim_cl,tim_rattle,tim_pbc,tim_LJ_tot,tim_solv_solv,tim_check,tim_check2
real*8 tim_T_jk
integer,allocatable:: seed(:)
real*8,allocatable:: work(:)
complex*16,allocatable:: cwork(:)
integer,allocatable:: iwork(:),isuppz(:)

contains
!---------------------------------------------------------- 
!---------------------------------------------------------- 

subroutine setup
  implicit none
  character st_ch
  integer i,size_seed,seed2(2)
  real*8 rnd,c_0,c_e,kt

  pi=dacos(-1.d0)
  wave_to_J=2*pi*clight*hbar
  nold=0

  open(10,file="AFSSH.inp")
  read(10,*) iflow
  read(10,*) iproc
  read(10,*) iparallel
  read(10,*) iwait
  read(10,*) k2_val
  read(10,*) N_traj
  read(10,*) dtc
  read(10,*) tim_eq
  read(10,*) total_time
  read(10,*) iwrite
  read(10,*) nstep_write
  read(10,*) nstep_avg
  read(10,*) idistribution
  read(10,*) flag_frust
  read(10,*) energy_cutoff
  read(10,*) nclass
  read(10,*) nquant
  read(10,*) nbasis
  read(10,*) gamma_B
  read(10,*) temperature
  read(10,*) iforward
  read(10,*) icollapse
  read(10,*) ifriction
  read(10,*) seed2
  read(10,*) st_ch
  close(10)
  !----------------------------------------------------------

  if(st_ch.ne.'x') then
    write(6,*) "problem in reading input file"
    stop
  endif
  !---------------------------------------------------------- 

  energy_cutoff=energy_cutoff*wave_to_J
  kt=kb*temperature

  nsteps=nint(total_time/dtc)+1
  nsteps_eq=nint(tim_eq/dtc)+1
  beta=1.d0/(kb*temperature)
  !gamma_D=omg_B**2/gamma_B
  !lambda_D=V_reorg/4.d0

  !-----------------------------------------------------------------  
  i=nsteps/nstep_avg+1
  allocate(pop(nquant,i),pop_surf(nquant,i),pop_amp(nquant,i))
  allocate(x(nclass),v(nclass),acc(nclass))
  allocate(x_old(nclass),v_old(nclass),acc_old(nclass),x_hop(nclass))
  allocate(mass(nclass),omg(nclass),x_min(nclass))
  allocate(KE_DVR(nbasis,nbasis),r_grid(nbasis))
  allocate(delr(nquant,nclass),delp(nquant,nclass))
  allocate(delr_old(nquant,nclass),delp_old(nquant,nclass))
  allocate(si_adiab(nbasis,nquant),ci(nquant),V_k(nquant),V_k_old(nquant))
  allocate(Hamil_diab(nbasis,nbasis),delH_dels(nbasis,nbasis,nclass),delH_dels_ad(nquant,nquant,nclass))
  allocate(pot(nquant,nquant),force(nquant,nclass),force_old(nquant,nclass))
  allocate(mat(nbasis,nbasis),mat_adiab(nquant,nquant))
  allocate(d_ij(nquant,nquant,nclass),vdotd(nquant,nquant),hop_prob(nquant),W_overlap(nquant,nquant))
  allocate(hop_prob_net(nquant))
  allocate(ci_old(nquant),si_adiab_prev(nbasis,nquant))
  call random_seed(size=size_seed)
  allocate(seed(size_seed))
  do i=1,size_seed/2
    seed(i)=seed2(1)*(2*i+i*i-7)
  enddo
  do i=size_seed/2+1,size_seed
    seed(i)=seed2(2)*(i/2+34-i**3)
  enddo
  call random_seed(put=seed)
  call system_clock(count_rate=cnt_rate)
  !-----------------------------------------------------------------  

  if(iflow==2) then
    open(10,file="ifolder.inp")
    read(10,*) ifolder
    close(10)
    N_traj=N_traj/iparallel
    if(ifolder>1) then
      do i=1,(ifolder-1)*N_traj
        seed=seed+1
        call init_cond
      enddo
      call random_seed(put=seed)
    endif
  else
    ifolder=1
  endif

  tim_tot=0.d0

end subroutine setup
!---------------------------------------------------------- 

subroutine main
  implicit none
  integer i,j,k,n
  real*8 t1,t2
  integer if_reactant

  call files(0)

  call cpu_time(t1)

  call setup_parameters
  call initialize_averages
  !call draw_pes

  do i=1,N_traj
    traj_num=i
    call init_cond
    call equilibrate(if_reactant)
    if(if_reactant==1) then
      cnt_init=cnt_init+1
      iaverage=1;en_avg=0.d0;ensq_avg=0.d0
      call evolve(nsteps)
    endif
  enddo
  call write_average

  call cpu_time(t2);tim_tot=tim_tot+t2-t1
  call files(1)

end subroutine main
!---------------------------------------------------------- 

subroutine files(flag)
  implicit none
  integer,intent(in)::flag

  if(flag==0) then
    open(10,file="output")
    open(11,file="output_cl")
    open(12,file="output_qm")
    open(13,file="output_hop")
    open(14,file="output_overlap")
    open(15,file="output_dec")

    open(100,file="pop.out")
    open(101,file="cnts.out")
  else
    write(10,*)
    write(10,*)"Total time=",tim_tot
    close(10);close(11);close(12);close(13);close(14);close(15)
    close(100);close(101)
  endif

end subroutine files
!-----------------------------------------------------------------  

subroutine initialize_averages
  implicit none

  cnt_frust=0.d0
  cnt_collapse=0.d0
  cnt_init=0.d0
  cnt_term=0.d0
  pop=0.d0
  pop_surf=0.d0
  pop_amp=0.d0

end subroutine initialize_averages
!-----------------------------------------------------------------  

subroutine init_cond
  implicit none
  integer i
  real*8 sig_x,sig_p,rnd,ak
  real*8 energy_0,ci_diab(nbasis)

  do i=1,nclass
    !ak=2/(hbar*omg(i))*dtanh(beta*hbar*omg(i)/2.d0) !! Wigner
    ak=beta    !! Classical
    sig_x=1.d0/dsqrt(ak*mass(i)*omg(i)**2)
    sig_p=dsqrt(mass(i)/ak)
    call gaussian_random_number(rnd)
    x(i)=rnd*sig_x + x_min(i)
    call gaussian_random_number(rnd)
    v(i)=(1.d0/mass(i)*(rnd*sig_p))
  enddo

  state=1
  call evaluate_variables(0)

  if(ifriction==0) then
    energy_0=-1000.d0*wave_to_J
    v(1)=dsqrt(2.d0/mass(1)*(energy_0-V_k(1)))
    v(2)=0.d0
  endif

  call evaluate_variables(1)

  ci=0.d0
  if(x(1)<0.d0) state=1
  if(x(1)>0.d0) state=2

  ci(State)=1

  ihop=1
  iaverage=1
  iterminate=0
  flag_terminate=0

  curr_time=0.d0
  call evaluate_variables(0)
  call evaluate_variables(1)
  call compute_mat_diab

  en_avg=0.d0;ensq_avg=0.d0
  nst_av=0

end subroutine init_cond
!-----------------------------------------------------------------  

subroutine equilibrate(if_reactant)
  implicit none
  integer,intent(out)::if_reactant
  integer i

  iaverage=0
  en_avg=0.d0;ensq_avg=0.d0
  call evolve(nsteps_eq)

  call check_reactant(if_reactant)

end subroutine equilibrate
!-----------------------------------------------------------------  

subroutine evolve(nsteps)
  implicit none
  integer,intent(in) :: nsteps
  integer i,j,nstep_sm,iflag_coll,i_do_something
  real*8 t1,t2
  integer iterm

  !call cpu_time(t1)

  call write_output(1,1)
  do i=1,nsteps
    call write_output(i,0)
    call average(i)
    call save_old_state
    call evolve_classical(dtc)
    i_do_something=0
    if(ifriction==0) then
      do while(dabs(energy-energy_old)>energy_cutoff.and.i>1) 
        i_do_something=i_do_something+1
        call do_something(i_do_something)
      enddo
    endif
    if(i_do_something.ne.1)call evolve_quantum_small_dtq
    !call evolve_quantum_small_dtq
    if(ihop==1)call hop
    if(icollapse==1)call collapse(dtc,iflag_coll)
    if(flag_terminate==1) call traj_terminate(iterm)
      if(iterm==1)exit

    curr_time=curr_time+dtc
  enddo
  call write_output(1,1)

  !call cpu_time(t2)
  !tim_evolve=tim_evolve+t2-t1

end subroutine evolve
!-----------------------------------------------------------------  

subroutine do_something(i_do_something)
  implicit none
  integer,intent(in)::i_do_something
  real*8 acc_tent(nclass),dt,dtm(3)
  integer i,nstep

  if(i_do_something==1) then
    call evolve_quantum_small_dtq
    if(flag_hop==1) then
      state=state_tentative
      call evaluate_variables(0)
      v=v_old+0.5*(acc_old+acc_tent)*dtc
      call evaluate_variables(1)
    endif
  else
    dtm=1.d0 !! some large value
    dtm(1)=0.1d0/maxval(vdotd)
    dtm(2)=0.5*dtc*dsqrt(energy_cutoff/dabs(energy-energy_old))
    dtm(3)=dtc

    dt=minval(dtm)

    dt=dt/real(i_do_something)
    nstep=nint(dtc/dt)
    dt=dtc/real(nstep)
    call revert_state
    do i=1,nstep
      call evolve_classical(dt)
    enddo
    !if(dabs(energy-energy_old)<energy_cutoff) write(6,*)curr_time*1.d15,dt*1.d15
  endif


end subroutine do_something
!-----------------------------------------------------------------  

subroutine average(i)
  implicit none
  integer,intent(in) :: i
  integer j
  complex*16 ci_diab(nquant)
  real*8 r_avg
  integer if_reactant
  real*8 t1,t2

  !call cpu_time(t1)

  if(iwrite==1) then
    en_avg=en_avg+energy
    ensq_avg=ensq_avg+energy*energy
    nst_av=nst_av+1
  endif

  if(iaverage==1.and.(mod(i,nstep_avg)==1.or.nstep_avg==1)) then
    if(nstep_avg==1) then
      j=i
    else
      j=i/nstep_avg+1
    endif

!    pop(:,j)=pop(:,j)+si_adiab(:,state)**2
!    pop_surf(:,j)=pop_surf(:,j)+si_adiab(:,state)**2
!    pop(:,j)=pop(:,j)+2*realpart(ci(1)*dconjg(ci(2)))*si_adiab(:,1)*si_adiab(:,2)
!    pop_amp(1,j)=pop_amp(1,j)+cdabs(sum(si_adiab(1,:)*ci))**2
!    pop_amp(2,j)=pop_amp(2,j)+cdabs(sum(si_adiab(1,:)*ci))**2

    call check_reactant(if_Reactant)
    if(if_reactant==1)pop(1,j)=pop(1,j)+1

  endif

  !call cpu_time(t2)
  !tim_coll=tim_coll+t2-t1

end subroutine average
!-----------------------------------------------------------------  

subroutine average_end
  implicit none

end subroutine average_end
!-----------------------------------------------------------------  

subroutine check_reactant(if_reactant)
  implicit none
  integer,intent(out)::if_reactant

  if_reactant=0
  if(r_exp>r_cutoff) if_reactant=1

end subroutine check_reactant
!-----------------------------------------------------------------  

subroutine save_old_state
  implicit none

  x_old=x
  v_old=v
  acc_old=acc
  ci_old=ci
  state_old=state
  !ci2_old=ci2
  si_adiab_prev=si_adiab
  V_k_old=V_k
  force_old=force
  energy_old=energy
  delr_old=delr
  delp_old=delp

end subroutine save_old_state
!-----------------------------------------------------------------  

subroutine revert_state
  implicit none

  x=x_old
  v=v_old
  state=state_old
  ci=ci_old
  delr=delr_old
  delp=delp_old
  force=force_old
  !ci2=ci2_old
  call evaluate_variables(0)
  call evaluate_variables(1)

end subroutine revert_state
!-----------------------------------------------------------------  

subroutine evolve_quantum_exponential
  !! Basic quantum propagator
  !! CAUTION AMBER !!
  !! delr and delp not implmented!!
  implicit none
  integer i,j
  real*8 W_overlap(nquant,nquant),si_adiab_old(nbasis,nquant),V_k_old(nquant)

  si_adiab_old=si_adiab
  V_k_old=V_k

  call evolve_classical(dtc)

  do i=1,nquant
    do j=1,nquant
      W_overlap(i,j)=sum(si_adiab(:,i)*si_adiab_old(:,j))
    enddo
  enddo

  ci=ci*cdexp(-iota*V_k_old*dtc/hbar)
  ci=matmul(W_overlap,ci)

end subroutine evolve_quantum_exponential
!-----------------------------------------------------------------  

subroutine evolve_quantum_small_dtq
  implicit none
  integer i,nstep_qm
  real*8 dtq,dtq1,dtq2
  real*8 V_k_hold(nquant),dVk_dt(nquant)
  real*8 dforce_dt(nquant,nclass)
  complex*16 ci_prev(nquant),dci_dt(nquant)

  call compute_vdotd
  dVk_dt=(V_k-V_k_old)/dtc
  if(icollapse==1) then
    call compute_delH_dels_ad
  endif

  dtq1=0.02/maxval(vdotd)
  dtq2=0.02*hbar/maxval(V_k-sum(V_k)/real(nquant))
  dtq=dtq1
  if(dtq>dtq2)dtq=dtq2

  if(dtq>dtc)dtq=dtc
  nstep_qm=nint(dtc/dtq)
  dtq=dtc/real(nstep_qm)
  hop_prob=0.d0
  hop_prob_net=0.d0
  V_k_hold=V_k
  V_k=V_k_old
  call compute_mat_adiab

  flag_hop=0
  do i=1,nstep_qm
    call compute_hop_prob(dtq)
    if(flag_hop==0)call check_hop(i*dtq)
    ci_prev=ci
    call rk4(ci,dtq,dVk_dt)
    if((cdabs(ci_prev(state))**2-0.5d0)*(cdabs(ci(state))**2-0.5d0)<0.d0)tim_hop=i*dtq
  enddo

  if(icollapse==1) then
    call verlet_decoherence(dtc,W_overlap,V_k_old,dvk_dt)
  endif

  do i=1,nquant
    if(hop_prob_net(i)<0.d0)hop_prob_net=0.d0
    hop_prob_net(i)=1.d0-dexp(-hop_prob_net(i))
  enddo

end subroutine evolve_quantum_small_dtq
!-----------------------------------------------------------------  

subroutine compute_hop_prob(dtq)
  implicit none
  real*8,intent(in)::dtq
  integer i
  real*8 pr

  do i=1,nquant
    if(i.ne.state) then
      pr=-2*real(ci(i)*dconjg(ci(state)))*vdotd(i,state)
      pr=pr*dtq/cdabs(ci(state))**2
      if(pr<0.d0)pr=0.d0     !!!! CAUTION AMBER CHECK !!!!
      hop_prob(i)=pr
      hop_prob_net(i)=hop_prob_net(i)+pr
    endif
  enddo

end subroutine compute_hop_prob
!-----------------------------------------------------------------  

subroutine check_hop(tim)
  implicit none
  real*8,intent(in)::tim
  integer i
  real*8 rnd,pr

  call random_number(rnd)
  pr=0.d0
  flag_hop=0
  do i=1,nquant
    if(i.ne.state) then
      pr=pr+hop_prob(i)
      if(rnd<pr) then
        state_tentative=i
        flag_hop=1
        exit
      endif
    endif
  enddo

end subroutine check_hop
!-----------------------------------------------------------------  

subroutine rk4(ci,dtq,dVk_dt)
  implicit none
  complex*16,intent(inout)::ci(nquant)
  real*8,intent(in) :: dtq,dVk_dt(nquant)
  complex*16,dimension(1:nquant):: k1,k2,k3,k4

  k1=matmul(mat_adiab,ci)

  V_k=V_k+dVk_dt*dtq/2.d0
  call compute_mat_adiab

  k2=matmul(mat_adiab,ci+0.5*dtq*k1)
  k3=matmul(mat_adiab,ci+0.5*dtq*k2)

  V_k=V_k+dVk_dt*dtq/2.d0
  call compute_mat_adiab

  k4=matmul(mat_adiab,ci+dtq*k3)

  ci=ci+dtq/6.d0*(k1+2*k2+2*k3+k4)

end subroutine rk4
!-----------------------------------------------------------------  

subroutine verlet_decoherence(dt,W_mat,V_k0,dvk_dt)
  implicit none
  real*8,intent(in):: dt,W_mat(nquant,nquant),V_k0(nquant),dvk_dt(nquant)
  real*8 acc_dec(nquant,nclass),delf(nquant,nclass),temp(nclass)
  complex*16 temp_delr(nquant,nclass),temp_delp(nquant,nclass)
  !complex*16 ci_diab(nquant)
  integer i,j,k

  delF=force_old
  temp=delF(state,:)
  do i=1,nquant
    delF(i,:)=delF(i,:)-temp
    acc_dec(i,:)=delF(i,:)*cdabs(ci_old(i))**2/mass
  enddo

  do i=1,nquant
    delr(i,:)=delr(i,:)+delp(i,:)/mass*dt+0.5*acc_dec(i,:)*dt**2
    delp(i,:)=delp(i,:)+0.5*mass*acc_dec(i,:)*dt
  enddo

  !ci_diab=cdexp(iota*V_k0*dt/hbar)*cdexp(0.5*iota*dvk_dt*dt**2/hbar)*ci
  !ci_diab=matmul_lap(W_mat,ci_diab)
  delF=0.d0
  do j=1,nquant
    do k=1,nquant
      delF(j,:)=delF(j,:)+dabs(W_mat(j,k)**2)*(force(k,:)-force(state,:))
    enddo
  enddo
  !temp=delF(state,:)
  do i=1,nquant
  !  delF(i,:)=delF(i,:)-temp
  !  !acc_dec(i,:)=delF(i,:)*cdabs(ci_diab(i))**2/mass
    acc_dec(i,:)=delF(i,:)*cdabs(ci_old(i))**2/mass
  enddo

  do i=1,nquant
    delp(i,:)=delp(i,:)+0.5*mass*acc_dec(i,:)*dt
  enddo

  temp_delr=0.d0;temp_delp=0.d0
  do j=1,nquant
    do k=1,nquant
      temp_delr(j,:)=temp_delr(j,:)+dabs(W_mat(k,j)**2)*delr(k,:)
      temp_delp(j,:)=temp_delp(j,:)+dabs(W_mat(k,j)**2)*delp(k,:)
    enddo
  enddo
  delr=temp_delr
  delp=temp_delp

  !do i=1,nclass
  !  delr(:,i)=delr(:,i)-delr(state,i)
  !  delp(:,i)=delp(:,i)-delp(state,i)
  !enddo

end subroutine verlet_decoherence
!-----------------------------------------------------------------  

subroutine evolve_quantum_adiabatic
  implicit none
  complex*16,dimension(1:nquant):: k1,k2,k3,k4

  k1=matmul(mat_adiab,ci)

  call evolve_classical(dtc/2.d0)
  call compute_mat_adiab

  k2=matmul(mat_adiab,ci+0.5*dtc*k1)
  k3=matmul(mat_adiab,ci+0.5*dtc*k2)

  call evolve_classical(dtc/2.d0)
  call compute_mat_adiab

  k4=matmul(mat_adiab,ci+dtc*k3)

  ci=ci+dtc/6.d0*(k1+2*k2+2*k3+k4)

end subroutine evolve_quantum_adiabatic
!-----------------------------------------------------------------  

subroutine evolve_classical(dt)
  !! Velocity Verlet
  implicit none
  real*8,intent(in) :: dt
  real*8 gama_dt,c0,c1,c2
  real*8 delta_r(nclass),delta_v(nclass),acc_sav(nclass)
  real*8 t1,t2

  !call cpu_time(t1)

  if(ifriction==0) then
    x=x+v*dt+0.5*acc*dt*dt
    acc_sav=acc
    call evaluate_variables(0)
    v=v+0.5*(acc+acc_sav)*dt
    call evaluate_variables(1)
  endif

  if(ifriction==1) then
    gama_dt=gamma_B*dt
    c0=dexp(-gama_dt)
    c1=1.d0/gama_dt*(1.d0-c0)
    c2=1.d0/gama_dt*(1.d0-c1)
   
     call stochastic_force(delta_r,delta_v,dt)
     x=x+c1*dt*v+c2*dt*dt*acc+delta_r
     acc_sav=acc
     call evaluate_variables(0)
     v=c0*v+(c1-c2)*dt*acc_sav+c2*dt*acc+delta_v
     call evaluate_variables(1)
  
  endif

  !call cpu_time(t2);tim_ev_cl=tim_ev_cl+t2-t1

end subroutine evolve_classical
!-----------------------------------------------------------------  

subroutine rk4_classical(dtc)
  implicit none
  real*8,intent(in)::dtc
  real*8,dimension(2*nclass):: k1,k2,k3,k4,vec,vec0

  vec(1:nclass)=x;vec(nclass+1:2*nclass)=v
  vec0=vec

  call deriv_xv(vec,acc,k1)

  vec=vec0+dtc/2.d0*k1
  x=vec(1:nclass)
  call evaluate_variables(0)
  call deriv_xv(vec,acc,k2)

  vec=vec0+dtc/2.d0*k2
  x=vec(1:nclass)
  call evaluate_variables(0)
  call deriv_xv(vec,acc,k3)

  vec=vec0+dtc*k3
  x=vec(1:nclass)
  call evaluate_variables(0)
  call deriv_xv(vec,acc,k4)

  vec=vec0+dtc/6.d0*(k1+2*k2+2*k3+k4)
  x=vec(1:nclass)
  v=vec(nclass+1:2*nclass)
  call evaluate_variables(0)
  call evaluate_variables(1)

end subroutine rk4_classical
!-----------------------------------------------------------------  

subroutine deriv_xv(vec,acc,kk)
  implicit none
  real*8,intent(in)::vec(2*nclass),acc(nclass)
  real*8,intent(out)::kk(2*nclass)

  kk(1:nclass)=vec(nclass+1:2*nclass)
  kk(nclass+1:2*nclass)=acc

end subroutine deriv_xv
!-----------------------------------------------------------------  

subroutine traj_terminate(iterm)
  implicit none
  integer,intent(out) :: iterm

  iterm=0

end subroutine traj_terminate
!-----------------------------------------------------------------  

subroutine compute_mat_diab
  implicit none
  integer i,j
  real*8 t1,t2

  !call cpu_time(t1)

  mat=0.d0
  do i=1,nbasis
    do j=1,nbasis
      mat(i,j)=-iota/hbar*sum(si_adiab(i,:)*si_adiab(j,:)*V_k(1:nquant))
    enddo
  enddo

  !call cpu_time(t2)
  !tim_mat=tim_mat+t2-t1

end subroutine compute_mat_diab
!-----------------------------------------------------------------  

subroutine compute_mat_adiab
  implicit none
  integer i,j
  real*8 t1,t2
  real*8 V_avg
  
  !call cpu_time(t1)

  mat_adiab=-vdotd
  V_avg=sum(V_k)/real(nquant)
  do i=1,nquant
    mat_adiab(i,i)=mat_adiab(i,i)-iota/hbar*(V_k(i)-V_avg)
  enddo
      
  !call cpu_time(t2)
  !tim_mat=tim_mat+t2-t1
  
end subroutine compute_mat_adiab
!-----------------------------------------------------------------  

subroutine hop
  implicit none
  integer ifrust

  if(flag_hop==1) then
    call velocity_adjust(state_tentative,ifrust)
  endif

end subroutine hop
!-----------------------------------------------------------------  

subroutine velocity_adjust(state_tentative,ifrust)
  implicit none
  integer,intent(in)::state_tentative
  integer,intent(out)::ifrust
  real*8 gij,gama,aa,bb,cc,discr,dp(nclass),vd,f1,f2
  integer i,j,k,kp

  k=state;kp=state_tentative
  cc=V_k(state)-V_k(state_tentative)

  !if(dabs(W_overlap(k,kp))<0.4d0.or.dabs(W_overlap(k,kp))>0.6d0) then
  if(dabs(W_overlap(k,kp))<0.3d0) then
    call compute_dij_2state(x,k,kp,dp)
  else if(dabs(W_overlap(k,kp))>0.7d0) then
    call compute_delH_dels_ad
    dp=delH_dels_ad(k,k,:)-delH_dels_ad(kp,kp,:)
  else
    call compute_dij_2state(x,k,kp,dp)
    if(dabs(sum(v*dp))<0.1*dabs(vdotd(k,kp))) then
      call compute_delH_dels_ad
      dp=delH_dels_ad(k,k,:)-delH_dels_ad(kp,kp,:)
    endif
  endif
  dp=dp/dsqrt(sum(dp*dp))

  aa=0.d0
  bb=0.d0
  do i=1,nclass

    aa=aa+0.5/mass(i)*(dp(i)*dp(i))
    bb=bb+(v(i)*dp(i))

  enddo

  discr=bb**2+4*aa*cc
  if(discr<0.d0) then
    ifrust=1
    cnt_frust=cnt_frust+1.d0
    if(flag_frust==0)then
      gama=0.d0
      call compute_delH_dels_ad
      f1=sum(force(k,:)*dp)
      f2=sum(force(kp,:)*dp)
      vd=sum(v*dp)
      if(f1*f2<0.d0.and.vd*f2<0.d0) then
      !if(f1*f2<0.d0) then
        gama=bb/aa
      endif
    endif
    if(flag_frust>0)gama=0.d0
  else
    ifrust=0
    if(bb>=0.d0) gama=(bb-dsqrt(discr))/(2*aa)
    if(bb<0.d0)  gama=(bb+dsqrt(discr))/(2*aa)
    state=state_tentative
    delr=0.d0
    delp=0.d0
  endif

  do i=1,nclass
    v(i)=v(i)-gama*dp(i)/mass(i)
  enddo

!write(20,*)curr_time*1.d15,dp/dsqrt(sum(dp*dp)),x(1),ifrust
!write(21,*)curr_time*1.d15,k,kp,gama

  call evaluate_variables(0)
  call evaluate_variables(1)

end subroutine velocity_adjust
!-----------------------------------------------------------------  

subroutine reverse_velocity
  implicit none
  

end subroutine reverse_velocity
!-----------------------------------------------------------------  

subroutine collapse(dt,iflag_coll)
  implicit none
  real*8,intent(in) :: dt
  integer,intent(out) :: iflag_coll
  real*8 rnd,gama_collapse,gama_reset
  complex*16 su1
  integer n,i,j

  i=state

  if(icollapse==1) then

    iflag_coll=0
    do n=1,nquant
      if(n.ne.state) then
        gama_reset=sum((force(n,:)-force(i,:))*dble(delr(n,:)-delr(state,:)))/(2*hbar)
        !! CAUTION !! !! Assumes delr(n,n,:) is in direction of v(:) !!
        su1=cdabs((V_k(i)-V_k(n))*vdotd(i,n)*sum((delr(n,:)-delr(state,:))*v))/sum(v*v)
        gama_collapse=gama_reset-2/hbar*cdabs(su1)
        gama_collapse=gama_collapse*dt
        gama_reset=-gama_reset*dt
        call random_number(rnd)

        if(rnd<gama_collapse) then
          iflag_coll=1
          cnt_collapse=cnt_collapse+1
          if(icollapse==1) then
            !do j=1,nquant
            !  if(j.ne.n) ci(j)=ci(j)/dsqrt(1-cdabs(ci(n)**2))
            !enddo
            !! Erratum: Landry, Subotnik JCP 137, 229901 (2012)
            ci(i)=ci(i)/cdabs(ci(i))*dsqrt(cdabs(ci(i))**2+cdabs(ci(n))**2)
            ci(n)=0.d0

          endif
        endif
        if(rnd<gama_collapse.or.rnd<gama_reset) then
          if(icollapse==1) then
            do j=1,nquant
              delr(n,:)=0.d0;delr(j,:)=0.d0
              delp(n,:)=0.d0;delp(j,:)=0.d0
            enddo
          endif
        endif
      endif
    enddo

  endif

end subroutine collapse
!-----------------------------------------------------------------  

subroutine write_output(n,nflag)
  implicit none
  integer,intent(in)::nflag,n
  integer i
  real*8 t1,t2
  real*8 phase

  !call cpu_time(t1)

  if(nflag==0) then
    if(iwrite==1) then
      if(mod(n,nstep_write)==1.or.nstep_write==1) then
        write(10,'(4es17.7,i5)')curr_time*1.d15,energy/wave_to_J,sum(cdabs(ci)**2),temperature,state
        write(11,'(es15.5$)')curr_time*1.d15
        write(12,'(5f15.5)')curr_time*1.d15,cdabs(ci(1:2))**2,datan2(dimag(ci(1:2)),real(ci(1:2)))*180/pi
        write(13,'(5es15.5)')curr_time*1.d15,vdotd(1,2),dasin(W_overlap(1,2))/dtc,hop_prob_net(3-state),state*1.d0
        write(14,'(6f15.5)')curr_time*1.d15,W_overlap(1,1:2),W_overlap(2,1:2),determinant(W_overlap,nquant)
        write(15,'(6es15.5)')curr_time*1.d15,delr(1,1)*1.d10,delr(2,1)*1.d10
        do i=1,nclass
          write(11,'(2es15.5$)')x(i)*1.d10,v(i)
        enddo
        write(11,*)
      endif
    endif
  endif

  if(nflag==1) then
    if(iwrite==0)write(10,'(5es15.5)')traj_num,energy/wave_to_J,sum(cdabs(ci)**2),temperature
    if(iwrite==1) then
      write(10,*)"traj num=",traj_num
      write(10,*)"standard deviation=",dsqrt((ensq_avg-en_avg**2/dfloat(nst_av))/dfloat(nst_av))/wave_to_J
      write(10,*)"ci**2=",sum(cdabs(ci)**2)
      write(10,*);write(10,*)
      write(11,*);write(11,*)
      write(12,*);write(12,*)
      write(13,*);write(13,*)
      write(14,*);write(14,*)
      write(15,*);write(15,*)
    endif
  endif

  !call cpu_time(t2)
  !tim_wr_out=tim_wr_out+t2-t1

end subroutine write_output
!-----------------------------------------------------------------  

subroutine write_average
  implicit none
  integer i,j
  real*8 nf

  nf=dfloat(n_traj)
  cnt_frust=cnt_frust/nf
  cnt_collapse=cnt_collapse/nf

  !pop=pop/nf
  !pop_surf=pop_surf/nf
  !pop_amp=pop_amp/nf

  do i=1,nsteps/nstep_avg+1
    write(100,'(6f15.7)')(i-1)*nstep_avg*dtc*1.d15,pop(1,i),cnt_init
  enddo

  write(101,*) gamma_B/omg(1),cnt_frust,cnt_collapse

end subroutine write_average
!-----------------------------------------------------------------  

subroutine evaluate_variables(flag)
  implicit none
  integer,intent(in):: flag
  integer i,j

  if(flag==0) then
    call tise
    r_exp=sum(si_adiab(:,state)**2*r_grid)
    !call check_reactant
  endif

  if(flag==1) then
    KE_en=0.d0
    do i=1,nclass
      KE_en=KE_en+0.5*mass(i)*v(i)*v(i)
    enddo

    energy=pot_en+KE_en
    !temperature=2*KE_en/(nclass*kb)

    !vdotd=0.d0
    !do i=1,nclass
    !  vdotd=vdotd+v(i)*d_ij(:,:,i)
    !enddo
    !call compute_vdotd
    
  endif

end subroutine evaluate_variables
!-----------------------------------------------------------------  

subroutine tise
  !! Output - pot_en,acc
  !! Output - V_k,d_ij
  implicit none
  integer i,j,k
  real*8 Hamil(nbasis,nbasis),ens(nbasis),vect(nbasis,nquant)
  real*8 pot_cl,acc_cl(nclass),acc_qm(nclass),dpotcl_dx(nclass)
  real*8 si_adiab_old(nquant,nbasis)
  real*8 t1,t2

  !call cpu_time(t1)

  call compute_potential(Hamil,delH_dels)
  Hamil_diab=Hamil
  call diag(Hamil,nbasis,ens,vect,nquant)

!  si_adiab_old=si_adiab
  do i=1,nquant
    si_adiab(:,i)=vect(:,i)
    if(sum(si_adiab(:,i)*si_adiab_prev(:,i))<0.d0)si_adiab(:,i)=-si_adiab(:,i)
  enddo

  do i=1,nclass
    delH_dels_ad(state,state,i)=sum(si_adiab(:,state)*matmul(delH_dels(:,:,i),si_adiab(:,state)))
  enddo
  !call cpu_time(t2);tim_diag=tim_diag+(t2-t1)

  !call cpu_time(t1)

  call potential_classical(pot_cl,dpotcl_dx)
  acc_qm=-1.d0/mass*delH_dels_ad(state,state,:)
  acc_cl=-1.d0/mass*dpotcl_dx

  pot_en=pot_cl+ens(state)
  V_k=pot_cl+ens(1:nquant)
  acc=acc_cl+acc_qm

  !call cpu_time(t2);tim_cl=tim_cl+(t2-t1)

end subroutine tise
!-----------------------------------------------------------------  

subroutine compute_delH_dels_ad
  implicit none
  integer i,k

  force=0.d0
  do k=1,nquant
    do i=1,nclass
      delH_dels_ad(k,k,i)=sum(si_adiab(:,k)*matmul(delH_dels(:,:,i),si_adiab(:,k)))
    enddo
    force(k,:)=-delH_dels_ad(k,k,:)
  enddo

end subroutine compute_delH_dels_ad
!-----------------------------------------------------------------  

subroutine compute_dij
  implicit none
  integer i,k,kp

  do k=1,nquant-1
    do kp=k+1,nquant
      do i=1,nclass
        d_ij(k,kp,i)=sum(si_adiab(:,k)*matmul(delH_dels(:,:,i),si_adiab(:,kp)))
      enddo
      d_ij(k,kp,:)=d_ij(k,kp,:)/(V_k(kp)-V_k(k))
      d_ij(kp,k,:)=-d_ij(k,kp,:)
    enddo
  enddo

end subroutine compute_dij
!-----------------------------------------------------------------  

subroutine compute_dij_2state(x_hop,k,kp,dp)
  implicit none
  integer,intent(in):: k,kp
  real*8,intent(in):: x_hop(nclass)
  real*8,intent(out):: dp(nclass)
  real*8 x_sav(nclass)
  integer i

  x_sav=x
  x=x_hop
  call evaluate_variables(0)

  do i=1,nclass
    dp(i)=sum(si_adiab(:,k)*matmul(delH_dels(:,:,i),si_adiab(:,kp)))
  enddo
  dp=dp/(V_k(kp)-V_k(k))

  x=x_sav
  call evaluate_variables(0)

end subroutine compute_dij_2state
!-----------------------------------------------------------------  

subroutine compute_vdotd
  ! Meek, Levine, JPCL 5, 2351 (2014). Look at Supp info.
  implicit none
  integer i,j,k
  real*8,dimension(nquant,nquant) :: W,ci_W,si_W
  real*8 A,B,C,D,E
  real*8 Wlj,Wlk

  !Method 1
  !call compute_dij
  !vdotd=0.d0
  !do i=1,nclass
  !  vdotd=vdotd+v(i)*d_ij(:,:,i)
  !enddo

  !Method 2
!  do j=1,nquant
!    do k=1,nquant
!      W(j,k)=sum(si_adiab_prev(:,j)*si_adiab(:,k))
!      ci_W(j,k)=dacos(W(j,k))
!      si_W(j,k)=dasin(W(j,k))
!    enddo
!  enddo
!
!  vdotd=0.d0
!  do k=1,nquant-1
!    do j=k+1,nquant
!      A=-sinx_x(ci_W(j,j)-si_W(j,k))
!      B=sinx_x(ci_W(j,j)+si_W(j,k))
!      C=sinx_x(ci_W(k,k)-si_W(k,j))
!      D=sinx_x(ci_W(k,k)+si_W(k,j))
!      Wlj=dsqrt(1.d0-W(j,j)**2-W(k,j)**2)
!      if(Wlj==0.d0.or.nquant==2) then
!        E=0.d0
!      else
!        Wlk=(-W(j,k)*W(j,j)-W(k,k)*W(k,j))/Wlj
!        E=2*dasin(Wlj)/(dasin(Wlj)**2-dasin(Wlk)**2)
!        E=E*(Wlj*Wlk*dasin(Wlj)+dasin(Wlk)*(dsqrt((1-Wlj**2)*(1-Wlk**2))-1.d0))
!      endif
!      vdotd(k,j)=0.5/dtc*(ci_W(j,j)*(A+B)+si_W(k,j)*(C+D)+E)
!      vdotd(j,k)=-vdotd(k,j)
!    enddo
!  enddo

  !Method 3
  do i=1,nquant
    do j=1,nquant
      W_overlap(i,j)=sum(si_adiab_prev(:,i)*si_adiab(:,j))
    enddo
  enddo

  call orthoganalize(W_overlap,nquant)
  call logm(W_overlap,vdotd,nquant)
  vdotd=vdotd/dtc

end subroutine compute_vdotd
!-----------------------------------------------------------------  

subroutine orthoganalize(mat,n)
  integer,intent(in)::n
  real*8,intent(inout)::mat(n,n)
  real*8 S_mat(n,n)

  S_mat=matmul(transpose(mat),mat)
  call inverse_squareroot(S_mat,n)
  mat=matmul(mat,S_mat)

end subroutine orthoganalize
!-----------------------------------------------------------------  

subroutine setup_parameters
  implicit none
  integer i
  real*8 si_diab(nbasis,2),Vb

  r_min=-2.d-10
  r_max=2.d-10

  !-----------------------------------------------------------------  
  r_del=(r_max-r_min)/real(nbasis-1)
  mass_qm=1.d0*amu2kg
  do i=1,nbasis
    r_grid(i)=r_min+(i-1)*r_del
  enddo
  call compute_KE_matrix_dvr(KE_DVR,nbasis,r_del,mass_qm)
  !-----------------------------------------------------------------  

  omg_b=1500.d0*2*pi*clight
  Vb=2085.d0*wave_to_J
  A_pot=-0.5*mass_qm*omg_b**2
  B_pot=mass_qm**2*omg_b**4/(16*Vb)

  mass=20.d0*amu2kg
  omg(1)=200.d0*2*pi*clight
  omg(2)=287.d0*2*pi*clight

  k_coup(1)=1.5d4*wave_to_J*1.d20
  !k_coup(2)=-5.d4*wave_to_J*1.d30
  k_coup(2)=k2_val*wave_to_J*1.d30
  k_coup(3)=-7.d3*wave_to_J*1.d20

  x_min(1)=-0.2d-10
  x_min(2)=0.17d-10
  !x_min=0.d0

  x(1)=0.d-10
  x(2)=0.d0
  state=1
  call evaluate_variables(0)
  si_diab(:,1)=1.d0/dsqrt(2.d0)*(si_adiab(:,1)+si_adiab(:,2))
  si_diab(:,2)=1.d0/dsqrt(2.d0)*(si_adiab(:,1)-si_adiab(:,2))
  r_cutoff=dabs(0.5*sum(si_diab(:,1)**2*r_grid))

  gamma_B=gamma_B*omg(1)

write(6,*) k_coup(2)/wave_to_J/1.d30,r_cutoff*1.d10

end subroutine setup_parameters
!-----------------------------------------------------------------  

subroutine compute_potential(H_diab,delV_dels)
  implicit none
  real*8,intent(out) :: H_diab(nbasis,nbasis),delV_dels(nbasis,nbasis,nclass)
  real*8 V1,dv1_dx(nclass)
  real*8 V2,dv2_dx(nclass)
  integer i

  H_diab=KE_DVR
  delV_dels=0.d0
  do i=1,nbasis
    call pot_q(V1,dv1_dx,r_grid(i))
    call pot_coup(V2,dv2_dx,r_grid(i))
    H_diab(i,i)=H_diab(i,i)+V1+V2
    delV_dels(i,i,:)=dv1_dx+dv2_dx
  enddo

end subroutine compute_potential
!-----------------------------------------------------------------  

subroutine potential_classical(pot_cl,dv_dx)
  implicit none
  real*8,intent(out) :: pot_cl,dv_dx(nclass)

  pot_cl=sum(0.5*mass*omg**2*x**2)+k_coup(3)*x(1)*x(2)
  dv_dx=mass*omg**2*x
  dv_dx(1)=dv_dx(1)+k_coup(3)*x(2)
  dv_dx(2)=dv_dx(2)+k_coup(3)*x(1)

end subroutine potential_classical
!-----------------------------------------------------------------  

subroutine pot_q(V,dv_dx,q)
  implicit none
  real*8,intent(in)::q
  real*8,intent(out)::V,dv_dx(nclass)

  V=A_pot*q**2+B_pot*q**4
  dv_dx=0.d0

end subroutine pot_q
!-----------------------------------------------------------------  

subroutine pot_coup(V,dv_dx,q)
  implicit none
  real*8,intent(in)::q
  real*8,intent(out)::V,dv_dx(nclass)

  V=k_coup(1)*q*x(1)+k_coup(2)*q**2*x(2)
  dv_dx(1)=k_coup(1)*q
  dv_dx(2)=k_coup(2)*q**2

end subroutine pot_coup
!-----------------------------------------------------------------  

subroutine check_acceleration
  implicit none
  integer i,nflag
  real*8 delx,en_old,acc_sav(nclass)
  real*8 q0,rnd

  delx=1.d-17
  state=1

  do i=1,nclass
    call random_number(rnd)
    x(i)=(rnd*2-1.d0)*1.d-10
  enddo

  call evaluate_variables(0)
  en_old=pot_en;acc_sav=acc

  write(6,*) "delx=",delx
  write(6,*)

  do i=1,nclass
      x(i)=x(i)+delx
      call evaluate_variables(0)
      acc(i)=-(pot_en-en_old)/delx/mass(i)
      write(6,*)"Analytical acceleration =",acc_sav(i)
      write(6,*)"Numerical acceleration  =",acc(i)
      write(6,*)"Error =",(acc(i)-acc_sav(i))/acc(i)*100.d0
      write(6,*)
      x(i)=x(i)-delx
  enddo

  stop

end subroutine check_acceleration
!---------------------------------------------------------- 

function inverse(mat)
  complex*16,intent(in) :: mat(2,2)
  complex*16 inverse(2,2),det
  complex*16 a,b,c,d

  a=mat(1,1);b=mat(1,2)
  c=mat(2,1);d=mat(2,2)

  det=a*d-b*c
  inverse(1,1)=d
  inverse(1,2)=-b
  inverse(2,1)=-c
  inverse(2,2)=a

  inverse=inverse/det

end function inverse
!-----------------------------------------------------------------  

subroutine stochastic_force(delr,delv,dt)
  implicit none
  real*8,intent(in)::dt
  real*8,intent(out) :: delr(nclass),delv(nclass)!f(nclass)
  integer i
  real*8 rnd1,rnd2,sig_r,sig_v,sig_rv,gdt

  gdt=gamma_B*dt

  do i=1,nclass

    sig_r=dt*dsqrt(kb*temperature/mass(i) *1.d0/gdt*(2-1.d0/gdt*(3-4*dexp(-gdt)+dexp(-2*gdt))))
    sig_v=dsqrt(kb*temperature/mass(i)*(1-dexp(-2*gdt)))
    sig_rv=(dt*kb*temperature/mass(i)* 1.d0/gdt *(1-dexp(-gdt))**2)/(sig_r*sig_v)  !! correlation coeffecient

    call gaussian_random_number(rnd1)
    call gaussian_random_number(rnd2)
    delr(i)=sig_r*rnd1
    delv(i)=sig_v*(sig_rv*rnd1+dsqrt(1-sig_rv**2)*rnd2)
  enddo

!  delr=delr-sum(delr)/dfloat(nclass)
!  delv=delv-sum(delv)/dfloat(nclass)

end subroutine stochastic_force
!-----------------------------------------------------------------  

function commute(A,B,iflag)
  integer,intent(in) :: iflag
  real*8,intent(in) :: A(:,:)
  complex*16,intent(in) :: B(:,:)
  complex*16 commute(nquant,nquant)
  real*8 delA(nquant,nquant)
  complex*16 tmp
  integer j,k

  if(iflag==0) commute=matmul(A,B)-matmul(B,A)

  if(iflag==1) then
    !! Assume A is diagonal
    do j=1,nquant
      do k=1,nquant
        commute(j,k)=B(j,k)*(A(j,j)-A(k,k))
      enddo
    enddo
  endif

  if(iflag==2) then
    !! Assume A is tridiagonal, with a_ii=0, and a_ij=-a_ji (a is assumed to be d_ij)
    do j=1,nquant
      do k=1,nquant
        tmp=0.d0
        if(j<nquant) tmp=tmp+A(j,j+1)*B(j+1,k)
        if(j>1) tmp=tmp-A(j-1,j)*B(j-1,k)
        if(k>1) tmp=tmp-A(k-1,k)*B(j,k-1)
        if(k<nquant) tmp=tmp+A(k,k+1)*B(j,k+1)
      enddo
    enddo
  endif

end function commute
!-----------------------------------------------------------------  

function anti_commute(A,B,iflag)
  integer,intent(in) :: iflag
  real*8,intent(in) :: A(:,:)
  complex*16,intent(in) :: B(:,:)
  complex*16 anti_commute(nquant,nquant)
  real*8 delA(nquant,nquant)
  integer i,j

  if(iflag==0) anti_commute=matmul(A,B)+matmul(B,A)

  if(iflag==1) then
    !! Assume A is diagonal
    do i=1,nquant
      do j=1,nquant
       anti_commute(i,j)=B(i,j)*(A(i,i)+A(j,j))
      enddo
    enddo
  endif

end function anti_commute
!-----------------------------------------------------------------  

subroutine draw_pes
  implicit none
  integer i,j,n
  real*8 H(2,2),dH(2,2)
  real*8 del,dd
  real*8 vm

  state=1
  vm=0.d0
  open(10,file="pes.out")
  open(11,file="dij.out")
  do i=1,50
    x(1)=-0.5d-10+1.0d-10*i/50.d0
    do j=1,50
      x(2)=-0.1d-10+0.6d-10*j/50.d0
      !call compute_potential(H,dH)
      call evaluate_variables(0)
      call compute_dij
      write(10,'(4f19.7)')x*1.d10,V_k(1:2)/wave_to_J
      write(11,'(4f19.7)')x*1.d10,dsqrt(sum(d_ij(1,2,:)*d_ij(1,2,:)))/1.d10
      if(x(1)<0.d0.and.V_k(1)<vm)vm=V_k(1)
    enddo
    write(10,'(4f15.7)')
    write(11,'(4f15.7)')
  enddo
  close(10)
  close(11)

  x(2)=0.18d-10
  n=1000
  do j=1,n
    x(1)=-0.5d-10+1.d-10*(j-1)/real(n-1)
    call evaluate_variables(0)
    write(20,'(9es15.7)')x(1)*1.d10,V_k/wave_to_J!,dabs(d_ij(1,2,1:2)/1.d10)
  enddo

  x(1)=0.d0
  n=1000
  do j=1,n
    x(2)=-0.1d-10+0.6d-10*(j-1)/real(n-1)
    call evaluate_variables(0)
    call compute_dij
    write(21,'(9f15.7)')x(2)*1.d10,dabs(d_ij(1,2,1:2)/1.d10),(V_k(2)-V_k(1))/wave_to_J
    write(22,'(9f15.7)')x(2)*1.d10,V_k(1:2)/wave_to_J
  enddo
write(6,*) vm/wave_to_J


  stop

end subroutine draw_pes
!-----------------------------------------------------------------  

subroutine draw_wavefn
  implicit none
  integer i
  real*8 V1,dv1_dx(nclass)
  real*8 V2,dv2_dx(nclass)

  !write(23,*)x(2)*1.d10
  do i=1,nbasis
    !call pot_q(V1,dv1_dx,r_grid(i))
    !call pot_coup(V2,dv2_dx,r_grid(i))
    !write(22,*)r_grid(i),(V1+V2)/wave_to_J
    write(23,*)r_grid(i),si_adiab(i,1:2)
  enddo
  !write(22,*);write(22,*)
  write(23,*);write(23,*)

end subroutine draw_wavefn
!-----------------------------------------------------------------  

function gamma_fn(z)
  !http://www1.uprh.edu/rbaretti/GammaFunction7dic2010.htm
  implicit none
  integer nstep,i
  real*8 ti,tf,dt,t
  complex*16 sum,gamma_fn,z,f

  data ti,tf,nstep /0.d0,30.d0,20001/

  f(t)=t**(z)*dexp(-t)

  dt=(tf-ti)/dfloat(nstep)
!  sum=dt**z/z -dt**(z+1.d0)/(z+1.d0)
  sum=0.d0

!  do i=2,nstep,2
  do i=1,nstep
    t=ti+dt*dfloat(i)
!    sum=sum+(dt/3.d0)*(f(t-dt)+ 4.d0*f(t) +f(t+dt))
    sum=sum+f(t)
  enddo

  gamma_fn=sum*dt/z

end function gamma_fn
!---------------------------------------------------------- 

pure function arg(z)
  implicit none
  real*8 arg
  complex*16,intent(in) :: z
  real*8 zr,zi

  zr=dble(z)
  zi=dimag(z)
  arg=datan2(zi,zr)
!  arg=datan(zi/zr)

end function arg
!----------------------------------------------------------   

subroutine diag(mat,n,eigen_value,eigen_vect,m_values)
  !! Diaganalizing matrix using dsyevr. First m_values eigen values and eigenvectors computed.
  !! The module's common variables should contain:

  !! Initialize nold=0 

  !! nold makes sure that everytime value of n changes, work and iwork are re-allocated for optimal performance.
  !! mat is destroyed after use.

  implicit none
  integer,intent(in) :: n,m_values
  real*8,intent(out) :: eigen_value(n),eigen_vect(n,m_values)
  real*8,intent(inout) :: mat(n,n)
  real*8 vl,vu,abstol
  integer il,iu,info,m,AllocateStatus
  integer lwork,liwork

  vl=0.d0;vu=0.d0   !! not referenced
  il=1;iu=m_values
  abstol=0.d0
  info=0

  if(nold.ne.n .or. .not.allocated(work) .or. .not.allocated(iwork) .or. .not.allocated(isuppz)) then
  !if(nold.ne.n) then
    lwork=-1;liwork=-1
    if(allocated(isuppz))deallocate(isuppz)
    if(allocated(work))deallocate(work)
    if(allocated(iwork))deallocate(iwork)
    allocate(isuppz(2*m_values),work(n),iwork(n))
    call dsyevr('V','I','U',n,mat,n,vl,vu,il,iu,abstol,m,eigen_value,eigen_vect,n,isuppz,work,lwork,iwork,liwork,info)
    lwork=nint(work(1)); liwork=iwork(1)
    deallocate(work,iwork)
    allocate(work(lwork),STAT=AllocateStatus)
    if(allocatestatus.ne.0) write(6,*)"problem in diag, allocation"
    allocate(iwork(liwork),STAT=AllocateStatus)
    if(allocatestatus.ne.0) write(6,*)"problem in diag, allocation"
    nold=n
  endif

  lwork=size(work)
  liwork=size(iwork)

  call dsyevr('V','I','U',n,mat,n,vl,vu,il,iu,abstol,m,eigen_value,eigen_vect,n,isuppz,work,lwork,iwork,liwork,info)
  if(info.ne.0) then
    write(6,*) "problem in diagonalization",info
    stop
  endif

end subroutine diag
!---------------------------------------------------------- 

subroutine print_matrix(mat)
  implicit none
  real*8,intent(in)::mat(:,:)
  integer i,j,n

  write(24,*)x(2)*1.d10

  n=size(mat(1,:))
  do i=1,n
    do j=1,n
      write(24,'(es12.2$)')mat(i,j)
    enddo
    write(24,*)
  enddo
  write(24,*)

end subroutine print_matrix
!-----------------------------------------------------------------  

subroutine logm(mat,log_mat,n)
  !! http://arxiv.org/pdf/1203.6151v4.pdf
  implicit none
  integer,intent(in):: n
  real*8,intent(in):: mat(n,n)
  real*8,intent(out):: log_mat(n,n)
  integer i
  complex*16 T(n,n),en(n),vect(n,n)
  complex*16 dd(n,n)

  call schur(mat,T,n,en,vect,nold,cwork)

  dd=0.d0
  do i=1,n
    dd(i,i)=cdlog(t(i,i)/cdabs(t(i,i)))
  enddo

  log_mat=matmul(vect,matmul(dd,conjg(transpose(vect))))

end subroutine logm
!-----------------------------------------------------------------  

subroutine inverse_squareroot(mat,n)
  !! http://arxiv.org/pdf/1203.6151v4.pdf
  implicit none
  integer,intent(in):: n
  real*8,intent(inout):: mat(n,n)
  integer i
  complex*16 T(n,n),en(n),vect(n,n)
  complex*16 dd(n,n)

  call schur(mat,T,n,en,vect,nold,cwork)

  dd=0.d0
  do i=1,n
    dd(i,i)=1.d0/t(i,i)**0.5d0
  enddo

  mat=matmul(vect,matmul(dd,conjg(transpose(vect))))

end subroutine inverse_squareroot
!-----------------------------------------------------------------  

subroutine schur(mat,T,n,eigen_value,eigen_vect,nold,cwork)
  !! Diaganalizing matrix using dsyevr. First m_values eigen values and eigenvectors computed.
  !! The module's common variables should contain:

  !! Initialize nold=0 

  !! nold makes sure that everytime value of n changes, work and iwork are re-allocated for optimal performance.
  !! mat is destroyed after use.

  implicit none
  integer,intent(in) :: n
  integer,intent(inout) :: nold
  complex*16,intent(out) :: eigen_value(n),eigen_vect(n,n)
  real*8,intent(in) :: mat(n,n)
  complex*16,intent(out) :: T(n,n)
  complex*16,allocatable,intent(inout):: cwork(:)
  real*8 rwork(n)
  complex*16 mat_c(n,n)

  integer lwork
  logical:: select
  logical bwork(n)
  integer sdim,info,AllocateStatus

  T=mat

  info=0
  sdim=0

  if(nold.ne.n .or. .not.allocated(cwork)) then
  !if(nold.ne.n) then
    lwork=-1
    if(allocated(cwork))deallocate(cwork)
    allocate(cwork(n))
    call zgees('V','N',SELECT,N,T,n,SDIM,eigen_value,eigen_vect,n,cWORK,LWORK,rwork,BWORK,INFO)
    lwork=int(cwork(1))
    deallocate(cwork)
    allocate(cwork(lwork),STAT=AllocateStatus)
    if(allocatestatus.ne.0) write(6,*)"problem in diag, allocation"
    nold=n
  endif

  lwork=size(cwork)
  call zgees('V','N',SELECT,N,T,n,SDIM,eigen_value,eigen_vect,n,cWORK,LWORK,rwork,BWORK,INFO)
  if(info.ne.0) then
    write(6,*) "problem in diagonalization",info
    stop
  endif

end subroutine schur
!---------------------------------------------------------- 

REAL FUNCTION determinant(matrix, n)
    !!http://web.hku.hk/~gdli/UsefulFiles/Example-Fortran-program.html
    IMPLICIT NONE
    REAL*8, DIMENSION(n,n) :: matrix
    INTEGER, INTENT(IN) :: n
    REAL*8 :: m, temp
    INTEGER :: i, j, k, l
    LOGICAL :: DetExists = .TRUE.
    l = 1
    !Convert to upper triangular form
    DO k = 1, n-1
        IF (matrix(k,k) == 0) THEN
            DetExists = .FALSE.
            DO i = k+1, n
                IF (matrix(i,k) /= 0) THEN
                    DO j = 1, n
                        temp = matrix(i,j)
                        matrix(i,j)= matrix(k,j)
                        matrix(k,j) = temp
                    END DO
                    DetExists = .TRUE.
                    l=-l
                    EXIT
                ENDIF
            END DO
            IF (DetExists .EQV. .FALSE.) THEN
                determinant= 0
                return
            END IF
        ENDIF
        DO j = k+1, n
            m = matrix(j,k)/matrix(k,k)
            DO i = k+1, n
                matrix(j,i) = matrix(j,i) - m*matrix(k,i)
            END DO
        END DO
    END DO
    
    !Calculate determinant by finding product of diagonal elements
    determinant= l
    DO i = 1, n
        determinant= determinant* matrix(i,i)
    END DO
    
END FUNCTION determinant
!-----------------------------------------------------------------  

subroutine gaussian_random_number(rnd)
  !! generates gaussian distribution with center 0, sigma 1
  !! q0+sig*rnd gives center=q0, sigma=sig
  implicit none
  real*8,intent(out)::rnd
  real*8 rnd1,rnd2,pi

  pi=dacos(-1.d0)

  call random_number(rnd1)
  call random_number(rnd2)
  rnd = dsqrt(-2*log(rnd1))*dcos(2*pi*rnd2)

end subroutine gaussian_random_number
!---------------------------------------------------------- 

subroutine compute_KE_matrix_dvr(KE_matrix,ngrid,delq,mass)
  !! computes KE matrix in DVR basis
  !! Appendix A of JCP 96, 1982 (1991)
  implicit none
  integer,intent(in) :: ngrid       !! size of DVR grid
  real*8,intent(inout) :: KE_matrix(ngrid,ngrid)
  real*8,intent(in) :: delq         !! step size in position of DVR basis
  real*8,intent(in) :: mass         !! mass
  integer i,j
  real*8 pi,hbar

  pi=dacos(-1.d0);hbar=1.05457266D-34

  KE_matrix=hbar**2/(2*mass*delq**2)
  do i=1,ngrid
    do j=1,ngrid
      KE_matrix(i,j)=KE_matrix(i,j)*(-1.d0)**(i-j)
      if(i==j) then
        KE_matrix(i,j)=KE_matrix(i,j)*pi**2/3.d0
      else
        KE_matrix(i,j)=KE_matrix(i,j)*2.d0/real(i-j)**2
      endif
    enddo
  enddo
end subroutine compute_KE_matrix_dvr
!---------------------------------------------------------- 

End Module mod_afssh

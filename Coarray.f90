!  Coarray.f90 
!
!  FUNCTIONS:
!  Coarray - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: Coarray
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

program Coarray
!use results
!use constants
!use start
!use interactions
!use sim_step
implicit none
    !real(8)::l[*],pressure[*],temp[*]
    !real(8)::th_eta[*],mtk_nib[*]
    !real(8)::co_press_tensor(3,3)[*],co_virial_tensor(3,3)[*]
    !real(8),allocatable :: x(:,:,:)[:],v(:,:,:)[:],a(:,:,:)[:],coms(:,:,:)[:] !x-položaji atomov, v- hitrosti atomov, a-vsota sil na particles 
    !real(8),allocatable :: particles(:,:)  !particles(i,(m,q))- masa,naboj delca
    !real(8)::co_potencial_el[*],co_potencial_lj[*],co_virial[*]
    !real(8)::dipole(3)[*],tests(10)
    !integer::step,i,imag
    !integer::file_unit1[*],file_unit2[*],file_unit3[*],file_unitOO[*],file_unitOH[*],file_unitHH[*],file_unitDIP[*],file_unitXCON
    !integer:: beginning,ending,rate
    !!real(8),allocatable::RDF_n(:),rdf_oh(:),rdf_hh(:)

    
        !if (this_image()==1) then 
        !    open(newunit=file_unit1,status='replace',file='./results/basic.txt')
        !    open(newunit=file_unit2,status='replace',file='./results/energy.txt')
        !    open(newunit=file_unit3,status='replace',file='./results/rho.txt')
        !    open(newunit=file_unitOO,status='replace',file='./results/rdfOO.txt')
        !    open(newunit=file_unitOH,status='replace',file='./results/rdfOH.txt')
        !    open(newunit=file_unitHH,status='replace',file='./results/rdfHH.txt')
        !    open(newunit=file_unitDIP,status='replace',file='./results/dipole.txt')
        !    open(newunit=file_unitXCON,status='old',action='read',file='./starting_configs/256.txt')
        !end if
        !file_unit1=file_unit1[1]
        !file_unit2=file_unit2[1]
        !file_unit3=file_unit3[1]
        !file_unitOO=file_unitOO[1]
        !file_unitOH=file_unitOH[1]
        !file_unitHH=file_unitHH[1]
        !file_unitDIP=file_unitDIP[1]
        !
        !model=5 
        !!ensamble=2 !1-NVT,2_NPT
        !!TBstat=2 !1-berendesen,2-NH/MTK
        !!starting_conditions=2 !1-file,2-random
        !
        !N=256
        !samplingst=1
        !time_step=2.d-15
        !temperature0=298
        !rho0=0.03327
        !pressure0=1
        !TB_strength=1.d-12
        !
        !ekv_time_steps=100
        !nvt_time_steps=100
        !npt_time_steps=100
        !
        !!RDF           del_r,rdf_max
        !rdf_parameters=(/0.02,8.0/)
        !!ewaldova vsota
        !K=10
        !!LJ cutoff
        !r_lower=9.0
        !r_upper=9.5
        !
        !
        !l=(N/rho0)**(1./3.)
        !!l=2*12.43864445
        !call ref_units(l)
        !MP=int(models(model,1))
        !
        !allocate(particles(mp*N,2))
        !call set_particles_values(particles)
        !
        !allocate(x(num_images(),mp*N/num_images(),3)[*],v(num_images(),3*N/num_images(),3)[*],a(num_images(),3*N/num_images(),3)[*],coms(num_images(),mp*N/num_images(),3)[*])
        !
        !Sync All
        !if (this_image()==1) then
        !    call stconx_from_file(x,v,a,l,particles,file_unitXCON)
        !    !call random_starting_conditions(x,v,particles,a,l)
        !end if
        !Sync All
        !call CO_BROADCAST(x,1)
        !call CO_BROADCAST(v,1)
        !call CO_BROADCAST(a,1)
        !
        !
        !call calc_self_int(particles)
        !call calc_v_mat(v_mat)
        !!call calc_lj_integrals(l)
        !
        !co_virial=0
        !co_virial_tensor=0
        !co_potencial_el=0
        !co_potencial_lj=0
        !if (this_image()==1) then
        !    print*,'start'
        !end if
        !call system_clock(beginning, rate)
        !call calc_coms(x,particles,l,coms)
        !call ewald_sum(x,particles,a,l,co_virial,co_potencial_el,co_virial_tensor,coms)
        !call lj_del_pot(x,particles,a,l,co_potencial_lj,co_virial,co_virial_tensor,coms)
        !call sym_tensor(co_virial_tensor)
        !sync all
        !call CO_SUM(co_potencial_el)
        !call CO_SUM(co_potencial_lj)
        !call CO_SUM(co_virial_tensor)
        !call CO_SUM(a)
        !sync all
        !call system_clock(ending)
        !if (this_image()==1) then
        !    print *, "elapsed time: ", real(ending - beginning) / real(rate)
        !    !print*,(co_potencial_el)/4184!(co_potencial_lj!(
        !    print*,(co_virial_tensor(1,1)+co_virial_tensor(2,2)+co_virial_tensor(3,3))/4184
        !    print*,(co_virial_tensor)/4184
        !end if
        !do step=0,ekv_time_steps
        !    call equilibration_step(x,v,a,particles,l,co_potencial_el,co_potencial_lj,step,pressure,temp,co_press_tensor,co_virial_tensor,coms,co_virial)
        !    !call calc_center_of_mass_motion(x,v,particles,file_unitDIP,l,step)
        !end do
        !
        !do step=0,nvt_time_steps
        !    call nvt_nh_step(x,v,a,particles,l,co_virial,co_potencial_el,co_potencial_lj,step,pressure,temp,th_eta,co_press_tensor,coms,co_virial_tensor)
        !    if (this_image()==1) then
        !        call write_basic(temp,pressure,step,file_unit1)
        !    !end if
        !    !if (this_image()==2) then
        !        call write_energy(co_potencial_lj,co_potencial_el,step,file_unit2)
        !    !end if
        !    !if (this_image()==3) then
        !        call write_rho(l,step,file_unit3)
        !    end if
        !    !call calc_dipole(x,particles,dipole,file_unitDIP,step,l)
        !end do
        !!
        !do step=0,npt_time_steps
        !    call npt_mtk_step(x,v,a,particles,l,co_virial,co_potencial_el,co_potencial_lj,step,pressure,temp,co_press_tensor,co_virial_tensor,th_eta,mtk_nib,coms)
        !    if (this_image()==1) then
        !        call write_basic(temp,pressure,step,file_unit1)
        !    !end if
        !    !if (this_image()==2) then
        !        call write_energy(co_potencial_lj,co_potencial_el,step,file_unit2)
        !    !end if
        !    !if (this_image()==3) then
        !        call write_rho(l,step,file_unit3)
        !    end if
        !    call calc_dipole(x,particles,dipole,file_unitDIP,step,l)
        !end do

        print*,'Hello'
        if (this_image()==1) then
        read*
        end if
end program Coarray


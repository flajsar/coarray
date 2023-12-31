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
use results
use constants
use start
use interactions
use sim_step
implicit none
    real(8)::l[*],pressure[*],temp[*]
    real(8)::th_eta[*],mtk_nib[*]
    real(8)::co_press_tensor(3,3)[*],co_virial_tensor(3,3)[*]
    real(8),allocatable :: x(:,:,:)[:],v(:,:,:)[:],a(:,:,:)[:],coms(:,:,:)[:] !x-polozaji atomov, v- hitrosti atomov, a-vsota sil na particles 
    real(8),allocatable :: particles(:,:),LJ_spheres(:,:),LJ_spheres_force(:,:)  !particles(i,(m,q))- masa,naboj delca
    real(8)::co_potencial_el[*],co_potencial_lj[*],co_virial[*]
    real(8)::dipole(3)[*]
    integer::step,i,imag
    integer::file_unit1,file_unit2,file_unit3,file_unitOO,file_unitOH,file_unitHH,file_unitDIP,file_unitXCON,file_unit_LJSF
    integer:: beginning,ending,rate
    !real(8),allocatable::RDF_n(:),rdf_oh(:),rdf_hh(:)

    
        if (this_image()==1) then 
            open(newunit=file_unit1,status='replace',file='./results/basic.txt')
            open(newunit=file_unitXCON,status='old',action='read',file='./starting_configs/256.txt')
            open(newunit=file_unit2,status='replace',file='./results/energy.txt')
            open(newunit=file_unit3,status='replace',file='./results/rho.txt')
            !open(newunit=file_unitOO,status='replace',file='./results/rdfOO.txt')
            !open(newunit=file_unitOH,status='replace',file='./results/rdfOH.txt')
            !open(newunit=file_unitHH,status='replace',file='./results/rdfHH.txt')
            open(newunit=file_unit_LJSF,status='replace',file='./results/LJSF.txt')
            open(newunit=file_unitDIP,status='replace',file='./results/dipole.txt')
        end if

        
        model=5 
        !ensamble=2 !1-NVT,2_NPT
        !TBstat=2 !1-berendesen,2-NH/MTK
        !starting_conditions=2 !1-file,2-random
        lj_atom=1  ! 1 Ne, 2 Ar, 3 Kr, 4 Xe
        N=256
        N_LJ_spheres=2
        samplingst=1
        time_step=2.d-15
        temperature0=298
        rho0=0.03327
        pressure0=1
        TB_strength=1.d-12
    
        ekv_time_steps=1500
        nvt_time_steps=-1
        npt_time_steps=10000

        !RDF           del_r,rdf_max
        rdf_parameters=(/0.02,8.0/)
        !ewaldova vsota
        K=7
        !LJ cutoff
        
        
        l=(N/rho0+N_LJ_spheres*4/3*pi*LJ_atom_sigma**3)**(1./3.)
        !l=2*12.43864445
        call ref_units(l)
        MP=int(models(model,1))
        
        allocate(particles(mp*N,2))
        call set_particles_values(particles)

        allocate(x(num_images(),mp*N/num_images(),3)[*],v(num_images(),3*N/num_images(),3)[*],a(num_images(),3*N/num_images(),3)[*],coms(num_images(),mp*N/num_images(),3)[*])
        allocate(LJ_spheres(N_LJ_spheres,3),LJ_spheres_force(N_LJ_spheres,3))
        LJ_spheres(1,1)=l/2+0.25
        LJ_spheres(2,1)=l/2-0.25

        LJ_spheres(1,2)=l/2
        LJ_spheres(2,2)=l/2
        LJ_spheres(1,3)=l/2
        LJ_spheres(2,3)=l/2
        
        Sync All
        if (this_image()==1) then
            call stconx_from_file(x,v,a,l,particles,file_unitXCON)
            !call random_starting_conditions(x,v,particles,a,l)
        end if
        Sync All
        call CO_BROADCAST(x,1)
        call CO_BROADCAST(v,1)
        call CO_BROADCAST(a,1)
        
        
        call calc_self_int(particles)
        call calc_v_mat(v_mat)
        !call calc_lj_integrals(l)
        
        co_virial=0
        co_virial_tensor=0
        co_potencial_el=0
        co_potencial_lj=0
        if (this_image()==1) then
            print*,'start'
        end if
        call system_clock(beginning, rate)
        call calc_coms(x,particles,l,coms)
        call ewald_sum(x,particles,a,l,co_virial,co_potencial_el,co_virial_tensor,coms)
        call lj_del_pot(x,particles,a,l,co_potencial_lj,co_virial,co_virial_tensor,coms)
        call sym_tensor(co_virial_tensor)
        sync all
        call CO_SUM(co_potencial_el)
        call CO_SUM(co_potencial_lj)
        call CO_SUM(co_virial_tensor)
        call CO_SUM(a)
        sync all
        call system_clock(ending)
        if (this_image()==1) then
            print *, "elapsed time: ", real(ending - beginning) / real(rate)
        end if
        
        do step=0,1500
            call equilibration_step(x,v,a,particles,l,co_potencial_el,co_potencial_lj,step,pressure,temp,co_press_tensor,co_virial_tensor,coms,co_virial,LJ_spheres,LJ_spheres_force)
        end do
        


        do step=0,nvt_time_steps
            call nvt_nh_step(x,v,a,particles,l,co_virial,co_potencial_el,co_potencial_lj,step,pressure,temp,th_eta,co_press_tensor,coms,co_virial_tensor,LJ_spheres,LJ_spheres_force)
            if (this_image()==1) then
                call write_basic(temp,pressure,step,file_unit1)
                call write_energy(co_potencial_lj,co_potencial_el,step,file_unit2)
                call write_rho(l,step,file_unit3)
                call write_LJSF(LJ_spheres_force,step,file_unit_LJSF)
            end if
            call calc_dipole(x,particles,dipole,file_unitDIP,step,l)
        end do
        !
    do rate=0,25
        LJ_spheres(1,1)=LJ_spheres(1,1)+0.25*(1-kronecker(1,rate+1))
        LJ_spheres(2,1)=LJ_spheres(2,1)-0.25*(1-kronecker(1,rate+1))
        if (this_image()==1) then
            write(file_unit_LJSF,*)sqrt(len_sq(LJ_spheres(1,:)-LJ_spheres(2,:)))
        end if
        do step=0,ekv_time_steps
            call equilibration_step(x,v,a,particles,l,co_potencial_el,co_potencial_lj,step,pressure,temp,co_press_tensor,co_virial_tensor,coms,co_virial,LJ_spheres,LJ_spheres_force)
        end do
        do step=0,npt_time_steps
            call npt_mtk_step(x,v,a,particles,l,co_virial,co_potencial_el,co_potencial_lj,step,pressure,temp,co_press_tensor,co_virial_tensor,th_eta,mtk_nib,coms,LJ_spheres,LJ_spheres_force)
            if (this_image()==1) then
                call write_basic(temp,pressure,step,file_unit1)
                call write_energy(co_potencial_lj,co_potencial_el,step,file_unit2)
                call write_rho(l,step,file_unit3)
                call write_LJSF(LJ_spheres_force,step,file_unit_LJSF)
            end if
            call calc_dipole(x,particles,dipole,file_unitDIP,step,l)
        end do
        
    end do
end program Coarray


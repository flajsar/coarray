module termo_barostat
    use basic_fnc
    use constants
    implicit none
    
    
    contains
    
    !subroutine berendsen_termostat(v,temp0,temp,particles,tau)
    !    implicit none
    !    real(8),intent(in)::temp0,particles(mp*N,2),tau,temp
    !    real(8),intent(inout)::v(num_images(),3*n/num_images(),3)[*]
    !    real(8)::ksi
    !    integer::i,image,npi
    !    image=this_image()
    !    npi=mp*N/num_images()
    !                                        !prejšna vsota se deli z 3*(N-1)-N_c kjer je N_c število constraints-v tem primeru so 3 na molekulo
    !    ksi=sqrt(1-dt/tau *(1-temp0/temp)) !berendsenov faktor
    !    v(:,:,:)=ksi*v(:,:,:)
    !    sync all
    !    call synh_arrays(v)
    !end subroutine
    !!
    !!
    !!
    !!
    !!
    !!
    !!subroutine berendsen_barostat(pressure,pressure0,ni_bar,tau_bar,l,izoterm_st,dt)
    !!    implicit none
    !!    real(8),intent(in)::pressure,pressure0,tau_bar,izoterm_st,dt
    !!    real(8),intent(inout)::ni_bar,l
    !!    ni_bar=(1-izoterm_st * dt/tau_bar *(pressure0 - pressure))**(1./3.)
    !!    l=l*ni_bar
    !!end subroutine
    !!
    ! subroutine remove_comm(x,v,particles)
    !    implicit none
    !    real(8),intent(in)::x(num_images(),mp*n/num_images(),3)[*],particles(mp*N,2)
    !    real(8),intent(inout)::v(num_images(),3*n/num_images(),3)[*]
    !    
    !    integer::i,image,npi
    !    real(8)::comm(3)[*],mass[*]
    !    SAVE comm,mass
    !    comm=0
    !    npi=mp*N/num_images()
    !    mass=0
    !    image=this_image()
    !    do i=1,3*N/num_images()
    !        comm=comm+v(image,i,:)*particles(i+(mp-3)*(i-1)/3+(image-1)*npi,1)
    !        mass=mass+particles(i+(mp-3)*(i-1)/3+(image-1)*npi,1)
    !    end do
    !    
    !    sync all
    !    call CO_SUM(comm)
    !    call CO_SUM(mass)
    !    sync all
    !    comm=comm/mass
    !    do i=1,3*N/num_images()
    !        v(image,i,:)=v(image,i,:)-comm
    !    end do
    !    SYNC ALL
    !    call synh_arrays(v)
    !    SYNC ALL
    !    
    ! end subroutine
    !! 
    !! 
    !! 
    !! 
    ! subroutine calc_th_eta(th_eta,v,particles)
    ! real(8),intent(in)::particles(mp*N,2),v(num_images(),3*n/num_images(),3)[*]
    ! real(8),intent(inout)::th_eta[*]
    ! integer::i,npi,image
    ! real(8)::sum[*]
    ! save sum
    ! sum=0
    ! npi =mp*N/num_images()
    ! image=this_image()
    ! do i=1,3*N/num_images()
    !        sum=sum+particles(i+(mp-3)*(i-1)/3+(image-1)*npi,1)*len_sq(v(this_image(),i,:))
    ! end do
    ! call CO_SUM(sum)
    ! th_eta=th_eta + dt* (sum-(6*N-3)*temp0) /(4*th_q)
    ! end subroutine
    ! 
    ! subroutine nh_v_adjust(v,th_eta)
    ! real(8)::th_eta[*]
    ! real(8),intent(inout):: v(num_images(),3*n/num_images(),3)[*]
    ! 
    !    v=v-v*dt*th_eta/2
    ! end subroutine
    ! 
    ! 
    ! subroutine calc_temp(v,particles,temp)
    ! real(8),intent(in)::particles(mp*N,2),v(num_images(),3*n/num_images(),3)[*]
    ! real(8),intent(out)::temp[*]
    ! integer::i,npi,image
    ! temp=0
    ! npi =mp*N/num_images()
    ! image=this_image()
    ! do i=1,3*N/num_images()
    !        temp=temp+particles(i+(mp-3)*(i-1)/3+(image-1)*npi,1)*len_sq(v(this_image(),i,:))
    !end do
    !
    !sync all
    !call CO_SUM(temp)
    !temp=temp/(6*N-3)
    !end subroutine
    !! 
    !! 
    ! subroutine mtk_calc_th_eta(th_eta,particles,mtk_nib,temp)
    ! implicit none
    ! real(8),intent(in)::particles(mp*N,2),mtk_nib,temp
    ! real(8),intent(inout)::th_eta
    ! real(8)::vs,nf
    ! integer::i
    ! nf=6*n-3
    !
    ! 
    ! th_eta=th_eta+ dt/(4*th_q) * (  temp*nf - nf*temp0 + mtk_w * mtk_nib**2-temp0)
    ! 
    ! end subroutine
    !! 
    !! 
    ! subroutine calc_mtk_nib(th_eta,l,mtk_nib,p,t)
    ! implicit none
    ! real(8),intent(inout)::mtk_nib
    ! real(8),intent(in):: l,th_eta,p,t   
    ! 
    !
    ! mtk_nib=mtk_nib + dt/4 * ( 3/mtk_w * l**3 * (p-pressure0) + 3/mtk_w * t - th_eta* mtk_nib)
    ! end subroutine
    !! 
    !! 
    ! subroutine mtk_v_adjust(v,th_eta,mtk_nib)
    ! implicit none
    ! real(8),intent(in)::th_eta,mtk_nib
    ! real(8),intent(inout)::v(num_images(),3*N/num_images(),3)
    ! integer::i
    ! 
    ! do i=1,3*N/num_images()
    !     v(this_image(),i,:)=v(this_image(),i,:) - dt*(th_eta +(1+3/(6*N-3))* mtk_nib)*v(this_image(),i,:)/2
    ! end do
    ! call synh_arrays(v)
    ! 
    ! end subroutine
    !! 
    !! 
    ! subroutine mtk_box(l,mtk_nib)
    ! implicit none
    ! real(8),intent(in)::mtk_nib
    ! real(8),intent(inout)::l
    ! 
    ! l=l * dexp(dt*mtk_nib)
    ! 
    ! 
    ! 
    ! end subroutine
    !! 
    !! 
    !! 
    ! subroutine shift_coms(x,mtk_nib,particles,l,coms)
    ! real(8),intent(in)::mtk_nib,particles(mp*N,2),l
    ! real(8),intent(inout)::x(num_images(),mp*n/num_images(),3)[*],coms(num_images(),mp*n/num_images(),3)[*]
    ! integer::i
    ! real(8)::mol(3,3),cntr_mass(3),d_i(3,3)
    ! 
    ! do i=1,mp*n/num_images(),mp
    !     mol(1,:)=x(this_image(),i,:)
    !     mol(2,:)=x(this_image(),i+1,:)
    !     mol(3,:)=x(this_image(),i+2,:)
    !     mol(2,:)=mol(2,:)+l*int((mol(1,:)-mol(2,:))*2/l)
    !     mol(3,:)=mol(3,:)+l*int((mol(1,:)-mol(3,:))*2/l)
    !     cntr_mass=coms(this_image(),i,:)
    !     d_i(1,:)=mol(1,:)-cntr_mass
    !     d_i(2,:)=mol(2,:)-cntr_mass
    !     d_i(3,:)=mol(3,:)-cntr_mass
    !     cntr_mass=cntr_mass*(1+dt*mtk_nib)
    !     coms(this_image(),i,:)=cntr_mass
    !     coms(this_image(),i+1,:)=cntr_mass
    !     coms(this_image(),i+2,:)=cntr_mass
    !     x(this_image(),i,:)=cntr_mass+ d_i(1,:)
    !     x(this_image(),i+1,:)=cntr_mass+ d_i(2,:)
    !     x(this_image(),i+2,:)=cntr_mass+ d_i(3,:)
    !     x(this_image(),i,:)=periodic_bound(x(this_image(),i,:),l)
    !     x(this_image(),i+1,:)=periodic_bound(x(this_image(),i+1,:),l)
    !     x(this_image(),i+2,:)=periodic_bound(x(this_image(),i+2,:),l)
    ! end do
    ! end subroutine
    !! 
    !! 
    !subroutine calc_pressure(co_virial_tensor,v,l,p,particles,x,a,co_press_tensor,temp)
    !    implicit none
    !    real(8),intent(in)::l[*],co_virial_tensor(3,3)[*],particles(mp*N,2),x(num_images(),mp*N/num_images(),3)[*],v(num_images(),3*n/num_images(),3)[*],a(num_images(),3*n/num_images(),3)[*],temp
    !    real(8),intent(inout)::p[*],co_press_tensor(3,3)[*]
    !    
    !    integer::i,j1,j2,image,npi
    !    real(8)::mass,v_cm(3),T,v_i(3),virial,tempm(3)
    !    npi=mp*N/num_images()
    !    mass=0
    !    v_cm=0
    !    t=0
    !    !do i=1,3*n
    !    !    v_cm=v_cm+v(i,:)*particles(i,1)
    !    !    mass=mass+particles(i,1)
    !    !    t=t+particles(i,1)*len_sq(v(i,:))
    !    !end do
    !    !v_cm=v_cm/mass
    !    !t=t/(6*N-3)
    !    image=this_image()
    !    p=0
    !    co_press_tensor=0
    !    do i=1,3*N/num_images(),3
    !        tempm(1)=particles(i+(mp-3)*(i-1)/3+(image-1)*npi,1)
    !        tempm(2)=particles(i+(mp-3)*(i)/3+(image-1)*npi+1,1)
    !        tempm(3)=particles(i+(mp-3)*(i+1)/3+(image-1)*npi+2,1)
    !        v_i=(v(image,i,:)*tempm(1)+v(image,i+1,:)*tempm(2)+v(image,i+2,:)*tempm(3))/(tempm(1)+tempm(2)+tempm(3))
    !        !v_i=v_i-v_cm
    !        do j1=1,3
    !            do j2=j1,3
    !                co_press_tensor(j1,j2)=co_press_tensor(j1,j2)+v_i(j1)*v_i(j2)*(tempm(1)+tempm(2)+tempm(3))
    !            end do
    !        end do
    !        !p=p+len_sq(v_i(:)-v_cm(:))*(particles(i,1)+particles(i+n,1)+particles(i+2*n,1))
    !    end do
    !    !p=p+3*t
    !    SYNC ALL
    !    call CO_SUM(co_press_tensor)
    !    !virial=(virial_tensor(1,1)+virial_tensor(2,2)+virial_tensor(3,3))
    !    !p=(p-virial)/(3*l**3)
    !    if (image==1) then
    !        co_press_tensor=(co_press_tensor-co_virial_tensor)/(l**3)
    !        p=(co_press_tensor(1,1)+co_press_tensor(2,2)+co_press_tensor(3,3))/3
    !        p=p+temp/l**3
    !    end if
    !    call CO_BROADCAST(co_press_tensor,1)
    !    call CO_BROADCAST(p,1)
    ! end subroutine
     
    end module
        
        
        
        
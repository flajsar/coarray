module results
    use basic_fnc
    use constants
    implicit none
    contains
    
    subroutine write_energy(potencial_lj,potencial_el,step,file_unit2)
        real(8),intent(in)::potencial_lj,potencial_el
        integer,intent(in)::step,file_unit2
        
        write(file_unit2,*)step,potencial_lj*ener_ref,potencial_el*ener_ref
    end subroutine
    
    subroutine write_basic(temp,pressure,step,file_unit1)
        real(8),intent(in)::temp,pressure
        integer,intent(in)::step,file_unit1
        
        write(file_unit1,*)step,temp*temp_ref, pressure*press_ref/10**5
    end subroutine
    
    subroutine write_rho(l,step,file_unit3)
        real(8),intent(in)::l
        integer,intent(in)::step,file_unit3
        
        write(file_unit3,*)step,l, N/(l*angstrom)**3 * 18.0153*1d-3/avog_num
    end subroutine
    !
    !
    !subroutine calc_rdf_pict(x,rdf_n,rdf_oh,rdf_hh,l)
    !    implicit none
    !    real(8),intent(in)::x(mp*N,3),l
    !    real(8),intent(inout)::rdf_n(n_points),rdf_oh(n_points),rdf_hh(n_points)
    !    
    !    integer::i,j,point
    !    real(8)::razd,rho,c1,c2
    !    
    !
    !    rho=N/l**3
    !    count_rdf=count_rdf+1
    !    do i=1,N
    !        do j=1,N
    !            razd=sqrt(len_sq(diff12(x(i,:),x(j,:),l)))
    !            if (razd<rdf_max) then 
    !                point=1+floor(mod(razd,rdf_max)/del_r)
    !                rdf_n(point)=rdf_n(point)+(1)*kronecker(i,j)/(mp*N*pi*rho*point**2*del_r**3)
    !            end if
    !                
    !        end do
    !    end do 
    !    do i=N+1,3*N
    !        do j=N+1,3*N
    !            razd=sqrt(len_sq(diff12(x(i,:),x(j,:),l)))
    !            if ((razd<rdf_max) .and. ((i-j) .ne. N)) then 
    !                point=1+floor(mod(razd,rdf_max)/del_r)
    !                rdf_hh(point)=rdf_hh(point)+(1)*kronecker(i,j)/(mp*N*pi*rho*point**2*del_r**3)
    !            end if
    !                
    !        end do
    !    end do 
    !    do i=1,N
    !        do j=N+1,3*N
    !            razd=sqrt(len_sq(diff12(x(i,:),x(j,:),l)))
    !            if ((razd<rdf_max) .and. ((i-j) .ne. N)) then  
    !                point=1+floor(mod(razd,rdf_max)/del_r)
    !                rdf_oh(point)=rdf_oh(point)+(1)*kronecker(i,j)/(mp*N*pi*rho*point**2*del_r**3)
    !            end if
    !                
    !        end do
    !    end do 
    !    
    !end subroutine 
    !
    !subroutine calc_rdfs(rdf_n,rdf_oh,rdf_hh,file_unitOO,file_unitOH,file_unitHH)
    !    integer,intent(in)::file_unitOO,file_unitOH,file_unitHH
    !    real(8),intent(inout)::rdf_n(n_points),rdf_oh(n_points),rdf_hh(n_points)
    !    integer::i,ask
    !    rdf_n=rdf_n/count_Rdf
    !    rdf_oh=rdf_oh/count_Rdf
    !    rdf_hh=rdf_hh/count_Rdf
    !    ask=n_points
    !    do i=1,n_points
    !        write(file_unitOO,*)i*del_r,rdf_n(i)
    !        !if (i*del_r .lt. 1.5) then
    !        !     write(file_unitOH,*)i*del_r,0.0
    !        !     write(file_unitHH,*)i*del_r,0.0
    !        !else
    !            write(file_unitOH,*)i*del_r,rdf_oh(i)
    !            write(file_unitHH,*)i*del_r,rdf_hh(i)
    !        !end if
    !    end do
    !end subroutine
    !
    !
    !subroutine calc_dipole(x,particles,dipole,file_unitDIP,step,l)
    !    real(8),intent(in)::particles(mp*N,2),x(mp*N,3),l
    !    integer,intent(in):: file_unitDIP,step
    !    real(8),intent(out)::dipole(3)
    !    integer::i,j
    !    real(8)::mol(mp,3)
    !    dipole=0
    !    if (mp==4) then
    !        do i=1,N
    !            mol(1,:)=x(i,:)
    !            mol(2,:)=x(i+n,:)
    !            mol(3,:)=x(i+2*n,:)
    !            mol(4,:)=x(i+3*n,:)
    !            mol(2,:)=mol(2,:)+l*int((mol(1,:)-mol(2,:))*2/l)
    !            mol(3,:)=mol(3,:)+l*int((mol(1,:)-mol(3,:))*2/l)
    !            mol(4,:)=mol(4,:)+l*int((mol(1,:)-mol(4,:))*2/l)
    !
    !            do j=1,mp
    !                dipole=dipole+mol(j,:)*particles(i+(j-1)*N,2)
    !            end do
    !        end do
    !    else
    !        do i=1,N
    !            mol(1,:)=x(i,:)
    !            mol(2,:)=x(i+n,:)
    !            mol(3,:)=x(i+2*n,:)
    !            mol(2,:)=mol(2,:)+l*int((mol(1,:)-mol(2,:))*2/l)
    !            mol(3,:)=mol(3,:)+l*int((mol(1,:)-mol(3,:))*2/l)
    !            do j=1,mp
    !                dipole=dipole+mol(j,:)*particles(i+(j-1)*N,2)
    !            end do
    !        end do
    !    end if
    !    
    !    
    !
    !    write(file_unitDIP,*)step,dipole
    !end subroutine
    
    subroutine calc_center_of_mass_motion(x,v,particles,file_unitDIP,l,step)
        real(8),intent(in)::x(num_images(),3*n/num_images(),3)[*],particles(mp*N,2),v(num_images(),3*n/num_images(),3)[*],l[*]
        integer,intent(in)::file_unitDIP,step
        integer::i,image,npi
        real(8)::comm(3)[*],mass
        SAVE comm
        comm=0
        npi=mp*N/num_images()
        mass=0
        image=this_image()
        do i=1,3*N/num_images()
            comm=comm+v(image,i,:)*particles(i+(mp-3)*(i-1)/3+(image-1)*npi,1)
            mass=mass+particles(i+(mp-3)*(i-1)/3+(image-1)*npi,1)
        end do
        comm=comm/mass
        sync all
        call CO_SUM(comm)
        sync all
        if (this_image()==1) then
            write(file_unitDIP,*)step,comm
        end if
    end subroutine
end module
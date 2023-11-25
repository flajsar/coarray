module motion
    use basic_fnc
    use constants
    !use omp_lib
    implicit none
    !rattle_položaji -> zračuna halfstep hiutrosti, unconstrained položaje, ki jih potem popravi
    !rattle_hitrosti -> ko z novimi položaji zračunamo sile, zračunamo hitrosti in jih popravimo
    !verlet1/2 -> osnovni velocity verlet algoritem

    contains 

        !subroutine verlet1(N,x,a,v,particles,dt,l,ni_bar)
        !    implicit none
        !    integer,intent(in):: N
        !    real(8),intent(in):: dt,l,particles(mp*N,2),ni_bar
        !    real(8),intent(inout)::x(mp*N,3),a(3*N,3),v(3*N,3)
        !    integer::i
        !    !zračunamo halfstep hitrosti in uporabimo halfstep hitrosti in položaje da zračuna nove položaje
        !    do i=1,3*N
        !        v(i,:)=v(i,:)+a(i,:)*dt/(2*particles(i,1))           
        !        x(i,:)=x(i,:)*ni_bar +dt * v(i,:)
        !        x(i,:)=periodic_bound(x(i,:),l)
        !    end do
        !
        !end subroutine
        !
        !subroutine verlet2(N,x,a,v,particles,dt) 
        !        implicit none
        !        integer,intent(in):: N
        !        real(8),intent(in):: dt,particles(mp*N,2)
        !        real(8),intent(inout)::x(mp*N,3),a(3*N,3),v(3*N,3)
        !        integer::i
        !        do i=1,3*N 
        !            v(i,:)=v(i,:) +dt * a(i,:) / (2.0 * particles(i,1) )
        !        end do
        !
        !end subroutine

        subroutine rattle_polozaji(x,a,v,particles,l)
            implicit none
            real(8),intent(in):: l[*],particles(mp*N,2),a(num_images(),3*n/num_images(),3)[*]
            real(8),intent(inout)::x(num_images(),mp*N/num_images(),3)[*],v(num_images(),3*n/num_images(),3)[*]

            real(8):: tempv(3,3),tempx(3,3),g,tempx_old(3,3),popravek(3,3),tempm(3)
            real(8):: skalarni,d_old(3),d_new(3),d_hoh(3,3),c_mat(3),a_mat(3,3),lg_mat(3)  
            real(8):: dtkv_inv,redm,odlocitev,g_vs(3)
            integer::i,j,k,image,npi
        
            image=this_image()
            npi=mp*N/num_images()
            
            !zračunamo halfstep hitrosti
            do i=1,3*N/num_images()
                v(image,i,:)=v(image,i,:)+ a(image,i,:) *dt / (2*particles(i+(mp-3)*(i-1)/3+(image-1)*npi,1) )
            end do
            
            dtkv_inv=dt**(-2) !inverz kvadrata dt
            do i=1,npi,mp !uporabimo halfstep hitrosti in položaje da zračuna unconstrained položaje za eno molekulo
                odlocitev=1
                tempx_old(1,:)=x(image,i,:)
                tempx_old(2,:)=x(image,i+1,:)
                tempx_old(3,:)=x(image,i+2,:)!shranimo stare položaje, in zračunamo nove unconstrained položaje
                tempx(1,:)=x(image,i,:)+dt * v(image,i-(mp-3)*(i/mp),:)
                tempx(1,:)=periodic_bound(tempx(1,:),l)
                tempx(2,:)=x(image,i+1,:) + dt * v(image,i-(mp-3)*((i+1)/mp)+1,:)
                tempx(2,:)=periodic_bound(tempx(2,:),l)
                tempx(3,:)=x(image,i+2,:) + dt * v(image,i-(mp-3)*((i+2)/mp)+2,:)
                tempx(3,:)=periodic_bound(tempx(3,:),l)
                
                tempm(1)=particles(i,1) !shranimo mase atomov molekule
                tempm(2)=particles(i+1,1)
                tempm(3)=particles(i+2,1)
                tempv(1,:)=v(image,i-(mp-3)*(i/mp),:)
                tempv(2,:)=v(image,i-(mp-3)*((i+1)/mp)+1,:)
                tempv(3,:)=v(image,i-(mp-3)*((i+2)/mp)+2,:)
                d_hoh=r_oh
                d_hoh(2,3)=r_hh    
                d_hoh(3,2)=r_hh

                !zračunamo faktorje g za te tri atome

                !g_vs=0
                giteracije: do while (odlocitev > 0.00003)
                    
                    c_mat(1)=(len_sq(diff12(tempx(1,:),tempx(2,:),l))-r_oh**2)/4
                    c_mat(2)=(len_sq(diff12(tempx(1,:),tempx(3,:),l))-r_oh**2)/4
                    c_mat(3)=(len_sq(diff12(tempx(2,:),tempx(3,:),l))-r_hh**2)/4
                    
                    a_mat(1,1)=(tempm(1)**-1+tempm(2)**-1)*dot_product(diff12(tempx_old(1,:),tempx_old(2,:),l),diff12(tempx(1,:),tempx(2,:),l))
                    a_mat(1,2)=(tempm(1)**-1)*dot_product(diff12(tempx_old(1,:),tempx_old(3,:),l),diff12(tempx(1,:),tempx(2,:),l))
                    a_mat(1,3)=(-tempm(2)**-1)*dot_product(diff12(tempx_old(2,:),tempx_old(3,:),l),diff12(tempx(1,:),tempx(2,:),l))
                    
                    a_mat(2,1)=(tempm(1)**-1)*dot_product(diff12(tempx_old(1,:),tempx_old(2,:),l),diff12(tempx(1,:),tempx(3,:),l))
                    a_mat(2,2)=(tempm(1)**-1+tempm(3)**-1)*dot_product(diff12(tempx_old(1,:),tempx_old(3,:),l),diff12(tempx(1,:),tempx(3,:),l))
                    a_mat(2,3)=(tempm(3)**-1)*dot_product(diff12(tempx_old(2,:),tempx_old(3,:),l),diff12(tempx(1,:),tempx(3,:),l))
                    
                    a_mat(3,1)=(-tempm(2)**-1)*dot_product(diff12(tempx_old(1,:),tempx_old(2,:),l),diff12(tempx(2,:),tempx(3,:),l))
                    a_mat(3,2)=(tempm(3)**-1)*dot_product(diff12(tempx_old(1,:),tempx_old(3,:),l),diff12(tempx(2,:),tempx(3,:),l))
                    a_mat(3,3)=(tempm(2)**-1+tempm(3)**-1)* dot_product(diff12(tempx_old(2,:),tempx_old(3,:),l),diff12(tempx(2,:),tempx(3,:),l))
                    
                    
                    
                    call invert_33_matrix(a_mat)
                    
                    lg_mat=matmul(a_mat,c_mat)
                    
                    tempx(1,:)=tempx(1,:)+(lg_mat(1)*2*diff12(tempx_old(1,:),tempx_old(2,:),l)+lg_mat(2)*2*diff12(tempx_old(1,:),tempx_old(3,:),l))/(2*tempm(1))
                    tempx(2,:)=tempx(2,:)+(lg_mat(1)*2*diff12(tempx_old(2,:),tempx_old(1,:),l)+lg_mat(3)*2*diff12(tempx_old(2,:),tempx_old(3,:),l))/(2*tempm(2))
                    tempx(3,:)=tempx(3,:)+(lg_mat(2)*2*diff12(tempx_old(3,:),tempx_old(1,:),l)+lg_mat(3)*2*diff12(tempx_old(3,:),tempx_old(2,:),l))/(2*tempm(3))
                    tempx(1,:)=periodic_bound(tempx(1,:),l)
                    tempx(2,:)=periodic_bound(tempx(2,:),l)
                    tempx(3,:)=periodic_bound(tempx(3,:),l)
                    
                    tempv(1,:)=tempv(1,:)+(lg_mat(1)*2*diff12(tempx_old(1,:),tempx_old(2,:),l)+lg_mat(2)*2*diff12(tempx_old(1,:),tempx_old(3,:),l))/(2*tempm(1)*dt)
                    tempv(2,:)=tempv(2,:)+(lg_mat(1)*2*diff12(tempx_old(2,:),tempx_old(1,:),l)+lg_mat(3)*2*diff12(tempx_old(2,:),tempx_old(3,:),l))/(2*tempm(2)*dt)
                    tempv(3,:)=tempv(3,:)+(lg_mat(2)*2*diff12(tempx_old(3,:),tempx_old(1,:),l)+lg_mat(3)*2*diff12(tempx_old(3,:),tempx_old(2,:),l))/(2*tempm(3)*dt)

                    !zračunamo mamo vse g faktorje, lahko popravimmo pozicije
                    odlocitev=abs(len_sq(diff12(tempx(2,:),tempx(1,:),l)) -r_oh**2)+ abs(len_sq(diff12(tempx(1,:),tempx(2,:),l))-r_oh**2) +abs(len_sq(diff12(tempx(2,:),tempx(3,:),l))-r_hh**2)
                    !print*,i,image,odlocitev!zračunamo parameter odločitve                
                end do giteracije
                !sprejmemo nove položaje in hitrosti
                x(image,i,:)=tempx(1,:)
                x(image,i+1,:)=tempx(2,:)
                x(image,i+2,:)=tempx(3,:)
                v(image,i-(mp-3)*(i/mp),:)=tempv(1,:)
                v(image,i-(mp-3)*((i+1)/mp)+1,:)=tempv(2,:)
                v(image,i-(mp-3)*((i+2)/mp)+2,:)=tempv(3,:)
                !preveri in upošteva periodne meje
                x(image,i,:)=periodic_bound(x(image,i,:),l)
                x(image,i+1,:)=periodic_bound(x(image,i+1,:),l)
                x(image,i+2,:)=periodic_bound(x(image,i+2,:),l)
            end do
        end subroutine
        



        subroutine rattle_hit(x,a,v,particles,l) 
            implicit none
            real(8),intent(in):: l[*],particles(mp*N,2)
            real(8),intent(inout)::x(num_images(),mp*N/num_images(),3)[*],v(num_images(),3*n/num_images(),3)[*],a(num_images(),3*n/num_images(),3)[*]
        
        
            real(8):: tempv(3,3),h,tempx(3,3),popravek(3,3),tempm(3)
            real(8):: skalarni1,skalarni2,skalarni3
            real(8):: skalarni,d_old(3),v_new(3),d_hoh(3,3),h_vs(3),g_mat(3),vr_mat(3)  
            real(8):: redm,d_kv,odlocitev
            integer::i,j,k,image,npi
            
            image=this_image()
            npi=mp*N/num_images()
            
            do i=1,npi,mp 
                tempx(1,:)=x(image,i,:)
                tempv(1,:)=v(image,i-(mp-3)*(i/mp),:) +dt * a(image,i-(mp-3)*(i/mp),:) / (2.0 * particles(i,1) )
                tempx(2,:)=x(image,i+1,:)
                tempv(2,:)=v(image,i-(mp-3)*((i+1)/mp)+1,:)+dt * a(image,i-(mp-3)*((i+1)/mp)+1,:) /(2.0 * particles(i+1,1))
                tempx(3,:)=x(image,i+2,:)
                tempv(3,:)=v(image,i-(mp-3)*((i+2)/mp)+2,:)+dt * a(image,i-(mp-3)*((i+2)/mp)+2,:) /(2.0 * particles(i+2,1))
                tempm(1)=particles(i,1) !shranimo mase atomov molekule
                tempm(2)=particles(i+1,1)
                tempm(3)=particles(i+2,1)
                d_hoh=r_oh
                d_hoh(2,3)=r_hh    
                d_hoh(3,2)=r_hh
                
                vr_mat(1)=dot_product(tempv(2,:)-tempv(1,:),diff12(tempx(2,:),tempx(1,:),l))
                vr_mat(2)=dot_product(tempv(3,:)-tempv(1,:),diff12(tempx(3,:),tempx(1,:),l))
                vr_mat(3)=dot_product(tempv(3,:)-tempv(2,:),diff12(tempx(3,:),tempx(2,:),l))
                
                g_mat=matmul(v_mat,vr_mat)
                v(image,i-(mp-3)*(i/mp),:)=tempv(1,:)-(diff12(tempx(1,:),tempx(2,:),l)*g_mat(1)+diff12(tempx(1,:),tempx(3,:),l)*g_mat(2))/tempm(1)
                v(image,i-(mp-3)*((i+1)/mp)+1,:)=tempv(2,:)-(diff12(tempx(2,:),tempx(1,:),l)*g_mat(1)+diff12(tempx(2,:),tempx(3,:),l)*g_mat(3))/tempm(2)
                v(image,i-(mp-3)*((i+2)/mp)+2,:)=tempv(3,:)-(diff12(tempx(3,:),tempx(1,:),l)*g_mat(2)+diff12(tempx(3,:),tempx(2,:),l)*g_mat(3))/tempm(3)
                
            end do
        end subroutine
 
        subroutine postavitev_m(x,l)
            implicit none
            real(8),intent(in)::l
            real(8),intent(inout)::x(num_images(),mp*N/num_images(),3)
            real(8)::mol(3,3),vec(3),gama
            integer::i,image
            image=this_image()
            do i=1,N*mp/num_images(),mp
                mol(1,:)=x(image,i,:)
                mol(2,:)=x(image,i+1,:)
                mol(3,:)=x(image,i+2,:)
                mol(2,:)=mol(2,:)+l*int((mol(1,:)-mol(2,:))*2/l)
                mol(3,:)=mol(3,:)+l*int((mol(1,:)-mol(3,:))*2/l)
                
                vec=0.5*(mol(2,:)+mol(3,:))-mol(1,:)
                gama=r_om/(dsqrt(len_sq(vec)))
                x(image,i+3,:)=mol(1,:)+gama*vec
                x(image,i+3,:)=periodic_bound(x(image,i+3,:),l)
            end do
        end subroutine
        
end module
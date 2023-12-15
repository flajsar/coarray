module constants
    implicit none
    SAVE
    
    real(8)::atomic_mass=1.66053873d-27,Angstrom=1.d-10,unit_charge=1.602176462d-19,avog_num=6.02214199d23,k_boltz=1.380649d-23,epsilon0=8.85418781d-12
    real(8)::pi,time_ref,temp_ref,press_ref,velocity_ref,force_ref,ener_ref,eps0,kb
    real(8)::mtk_w,th_q,dt,g,tau
    real(8)::del_r,rdf_max
    integer::n_points
    integer::count_rdf
    integer::N,K
    integer::lj_atom
    real(8)::temp0,temperature0,pressure0,press0,rho0
    real(8)::mass_o, charge_o, mass_h,charge_h,mass_m,charge_m,  r_oh,r_hh ,angle0,r_om,a1,a2 
    real(8),dimension(5,13)::models
    integer::model,ensamble,TBstat,starting_conditions,MP
    integer::samplingst,N_Lj_spheres
    real(8)::time_step,LJ_sphere_parameters(4,2)
    real(8)::rdf_parameters(2),TB_strength,self_int,r_lower,r_upper,lj_integral1,lj_integral2
    integer::ekv_time_steps,nvt_time_steps,npt_time_steps
    real(8)::switch_d,switch_a,switch_b,switch_c
    real(8)::v_mat(3,3)
    real(8)::LJ_atom_sigma,LJ_atom_eps
    
   
end module
    
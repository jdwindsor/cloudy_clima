      subroutine relhum_vect(vapor_mmr,T,P,mw_atmos,relhum_vec) 
      
      !implicit none 
      INCLUDE 'CLIMA/INCLUDE/header.inc'
      COMMON/EBLOK/PG,TG,PG0,IMW,RSURF,OMEGA,POCEAN,IMOIST,
     2  BETA1,BETA2,FVDRY,PDRY

      !common/molec_weight/mw_atmos
      double precision vapor_mmr(ND),T(ND)!,relhum_vec(ND)
      real relhum_vec(ND)
      double precision mw_atmos,dtdlnp,dlnp
      double precision P(ND),pe(ND-1),te(ND-1)
      integer I,J,K,itop,ibot,iz

      do iz=ND-1,1,-1
            itop=iz
            ibot=iz+1
            dlnp=log( P(ibot)/P(itop))
            pe(iz) = 0.5*( p(itop) + p(ibot))
            dtdlnp = (t(itop) -t(ibot) )/dlnp
            te(iz) = t(ibot) + log( p(ibot)/pe(iz))*dtdlnp
      enddo



      do I=1,ND-1
      call satrat(te(I),psat)
      IF(PSAT.GT.PG) psat = POCEAN
      !print*,'mw_atmos',mw_atmos
      relhum_vec(I)=(vapor_mmr(I)*(mw_atmos/18.01528)*pe(I))/(psat)!vapor_mmr(I)/(psat*(18./mw_atmos)/P(I))
      IF (relhum_vec(I) .ge. 1.0)then
             relhum_vec(I)=1.0
      endif
      ! relhum_vec(I)=1.0
      enddo
      !print*,'psat',(P)
      relhum_vec=RSURF!(ND)=RSURF

      return 

      end





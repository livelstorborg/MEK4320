*******************************************************************
      subroutine setbhkorr(pluss)
      real pluss
      real add
      common/solhas/add
c     add er en forel|pig korreksjonsfaktor for aa leke med
c     solitonhastigheter
*--------------------------------------------------------
      add=pluss

      return
      end 

***************************************************************
*                                                             *
*  BEREGNER BOLGEHASTIGHETEN TIL SOLITONEN                    *
*  PARAMETER:  A =                                            *
*                                                             *
***************************************************************

      FUNCTION BHAST(A)
      REAL A,Bhast
ccc      real add
ccc      common/solhas/add
c     add er en forel|pig korreksjonsfaktor for aa leke med
c     solitonhastigheter, den er siden fjernet
*--------------------------------------------------------------
      REAL EPA

      if(a.gt.0.001) then
        EPA   = 1+A 
        BHAST = EPA*SQRT((EPA*LOG(EPA)-A)/
     %          (A*A*(A/3.0+0.5))) 
      else
        bhast=1.0+0.5*a
      end if

      RETURN
      END



************************************************************************
      function halvb(a,eps)
      real a,eps,halvb
*----------------------------------------------------------------------

      halvb=-log(0.25*eps)/sqrt(3.0*a)
      return
      end 


*****************************************************************************
*              B O U N
*     Beregner randverdier for en inkommende soliton.
*     parametere:
*             b - array for randverdier                                  O
*             ds - gitteravstand langs rand                              I
*             dn - gitteravstand normalt rand                            I
*             dt - tidssteg                                              I
*             k2 - det fylles verdier i b(0:k2)                          I
*             snut - posisjon av topp-punkt. snut = ds gir topp i b(1)   I
*                     etc. Topp kan godt ligge utenfor beregningsintervall.
*                     Dersom psi er lik 90 grader har snut tolkning
*                     som avstand fra rand til boelgetopp.
*             psi -  Vinkel mellom boelgetallsvektor og rand.             I
*                    psi er positiv naar btall peker i pos.
*                    normalretning og beregnes i grader
*             a -    amplitude
*             eps -  trunkeringskriterium, likevekt settes naar          I
*                    overfl. hevn. < eps*a
*             islag - type randbetingelse                                I
*                 =0  df/dn
*                 =1  f
*                 =2  dy/dn
*                 =3  y
*                  der f er potensial og y overflatehelning
*             ik - verdi lik 1 angir korreksjon                          I
*****************************************************************************
      subroutine boun(b,ds,dn,dt,k2,snut,psi,a,eps,islag,ik)   
      integer islag,ik,k2   
      real b(0:k2),ds,dn,dt,snut,psi,a,eps
*-----------------------------------------------------------------------------
      integer irand,j
      real cpsi,spsi,pi,fc,topp,c,ua,dinc

      pi=4.0*atan(1.0)
      cpsi=cos(pi*psi/180.0)
      spsi=sin(pi*psi/180.0)
   

      if(islag.eq.0) then
        fc=ik*(dt*dt-ds*ds*cpsi*cpsi)/24.0
        irand=2
      else
        if(islag.eq.1) then
          fc=ik*(dt*dt-dn*dn*spsi*spsi)/24.0
          irand=3
        else
          fc=0.0
          if(islag.eq.2) then
            irand=5
          else
            irand=1
          end if
        end if
      end if

      if(abs(abs(psi)-90.0).lt.0.001) then
        call soliprgen(b,0,0,ds,A,snut,C,UA,irand,fc,eps)   
        do 100 j=1,k2
        b(j)=b(0)
 100    continue
        return
      end if

      dinc=ds*cpsi
      topp=snut*cpsi
      call soliprgen(b,0,k2,Dinc,A,TOPP,C,UA,irand,fc,eps)   

      if(islag.eq.0 .or. islag.eq.2) then
        do 200 j=0,k2
        b(j)=b(j)*spsi
 200    continue
      end if

      return
      end

      
*****************************************************************************
*
*               S O L I P R G E N
*
*
*     Tabulerer verdier svarende til en soliton-form ved kall p} 'SOLITGEN',
*     men gj|r tre ting i tillegg:
*                 1. Trunkerer solitonen der h|yden er mindre  enn amp*eps
*                 2. fyller y med nuller utenfor trunkert omr}de.
*                 3. S|rger for  at soliton kalles med increment (dx) som
*                    aldri overstiger 1.0. Dette kan v}re viktig for
*                    konvergens av underliggende likningsl|sere.
*     Parameterene svarer til de for 'solitgen' untatt eps som alts} gir
*     trunkerings-kriteriet.
**************************************************************************
      subroutine soliprgen(Y,n0,N,DX,AMP,TOPP,C,UA,ik,fc,eps)   
      INTEGER n0,N,ik
      REAL Y(N0:N),DX,AMP,topp,C,UA,fc,eps
*---------------------------------------------------------------------------
      integer nf,i,na,nb,nrmax,nr
      parameter(nrmax=15000)
      real dxred,yr(nrmax),blen,txr,halvb

      do 50 i=n0,n
         y(i)=0.0
 50   continue

      if(amp.le.0.0) return

      blen=halvb(amp,eps)

      na= (topp-blen)/dx -1
      na=max(na,n0)
      nb=(topp+blen)/dx +1
      nb=min(nb,n)

      nf=dx+1
      dxred= dx/nf

      if(na.gt.nb) then
        if(ik.eq.3 .or. ik.eq.4) then
          txr=blen+2.0*dxred
          nr=2*txr/dxred+2
          call solitgen(Yr,1,Nr,DXred,AMP,TXR,C,UA,ik,fc,eps)   
          do 60 i=n0,n
          y(i)=yr(1)
 60       continue
        end if
        return
      end if

      nr=(nb-na)*nf+1
      txr=topp-na*dx+dxred

      do 70 i=1,nr
      yr(i)=0.0
 70   continue

      if(nr.gt.nrmax) then
        nr=nrmax
        nb=(nr-1)/nf + na
        write(0,*)'advarsel, for mange punkter i soliprod'
      end if


      call solitgen(Yr,1,Nr,DXred,AMP,TXR,C,UA,ik,fc,eps)   

      do 100 i=na,nb
       y(i)=yr((i-na)*nf+1)
 100   continue

      if( (ik.eq.3 .or. ik.eq.4) .and. na.gt.n0) then
        do 200 i=n0,na
        y(i)=y(na)
 200    continue
      end if

      return
      end




*****************************************************************************
*
*               S O L I T G E N
*
*
*     Tabulerer verdier svarende til en soliton-form. Det antas at solitonen
*     forplanter seg mot |kende indekser i y. Har den motsatt forplantning
*     beholdes verdiene for overflatehevninger mens foretegnet byttes for
*     hastighet og potensial. Den aktuelle solitonl|sningen er beskrevet
*     i Pedersen og Rygg (1987) og Pedersen (1988).

*     parametere:
*             y  - array som inneholder verdiene                          O
*             n0,n - arraygrenser for y i det kallende programmet         I
*                    maksimalt antall punkter er 1000.
*             dx - avstand (regnet i dyp) mellom punktene.                I
*             amp - maks overflatehevning (regnet i dyp)                  I
*             topp - posisjon av makspunkt, merk at topp=0                I
*                    vil plassere makspunkt hos y(0), topp=dx
*                    hos y(1) etc. Toppen kan godt v{re lokalisert
*                    utenfor arrayendene. V{r oppmerksom p{ at y(1) svarer
*                    til ulike fysiske posisjoner for ulike ukjente i ett
*                    alternerende gitter.
*             c - forplantningshastighet for solitonen (skalert med       O
*                 line{r gruntvannshastighet)
*             ua - maksverdi for hastighet, skalert som c                 O
*             ik - styringsparameter                                      I 
*                 =1  :  Y fylles med overflatehevninger
*                 =2  :  Y fylles med hastigheter
*                 =3  :  y fylles med potensialverdier. I dette tilfellet
*                          er verdiene satt slik at midtpunktdifferens gir
*                          eksakte hastigheter, og y(n) er valgt lik null 
*                          (potensialet er jo bestemt bare p} en konstant n{r)
*                 =4  :  Som 3, bortsett fra at verdiene finnes ved integrasjon
*                        av en spline-interpolant for hastigheten.
*                 =5  :  Y fylles med rom-deriverte av overflatehevningen
*
*             fc - faktor for modifikasjon av hast. og potensial i hht.   I
*                  hydrostatisk korreksjon av regneskjemaer.         
*                  Den mod. hastighet blir satt ved:
*
*                    U_mod = U +fc*U''
*             eps - trunkeringsgrense                                     I
*
**************************************************************************
      SUBROUTINE SOLITGEN(Y,n0,N,DX,AMP,Topp,C,UA,ik,fc,eps)
      INTEGER n0,N,ik
      REAL Y(N0:N),DX,AMP,topp,C,UA,fc,eps
*---------------------------------------------------------------------------
      INTEGER J,nant,ier,nmax,i1
      parameter(nmax=15000)
      real bhast, toppu,f(nmax),der1,der2,wk(5*nmax+10),uder
      real amplim


      amplim=0.035

      nant=n-n0+1

      if(nant.gt.nmax) then
        write(0,*)'for mange punkter i kall p} soliton'
        return
      end if
      
      c=bhast(amp)     
      ua=c*amp/(1.0+amp)            

      if(ik.eq.3) then
c      
c       hastigheter skal beregnes og legges i y slik at toppunkt er
c       forskj|vet -0.5*dx
c
c       pos av pot:     y(n-2)     y(n-1)      y(n)    ------------->x
c                                                                   
c       pos av hast:          u(n-2)     u(n-1)     
c
        toppu = topp -0.5*dx
      else
        toppu = topp
      end if

      if(amp.gt.amplim) then
         CALL STappg(y,n0,n,AMP,TOPPU,Dx,eps)
C       denne returnerer verdier for hastigheten, vi m} derfor regne
c       om n}r ik=1.
      else
c       for sm} amplituder er denne tryggere.
        call  solipert(Y,n0,N,Dx,0.0,0.0,0.0,0.0,0,0,AMP,TOPPU,C,eps)
        c=bhast(amp)   
        do 80 j=n0,n
         y(j)=c*y(j)/(1.0+y(j))
  80    continue
      end if

      if(ik.eq.1) then
       do 100 j=n0,n
         y(j)=y(j)/(c-y(j))
  100  continue
       return
      end if

      if(ik.eq.5) then
       i1=topp/dx
       i1=min(i1,n)
       do 150 j=n0,i1
         y(j)=c*uder(y(j),amp,ier)/(c-y(j))**2
  150  continue
       i1=max(n0-1,i1)
       do 160 j=i1+1,n
         y(j)=-c*uder(y(j),amp,ier)/(c-y(j))**2
  160  continue
       return

      end if   

      if(ik.eq.4) then

        do 300 j=1,nant
        f(j)=y(j+n0-1)
 300    continue

        der1=uder(f(1),amp,ier)
        der2=uder(f(nant),amp,ier)
        call ingrer(y(n0),f,nant,dx,der1,der2,nant,wk)
        return
      end if

c
c    her modifiseres hastigheter i hht. verdi p} fc
c
     
      do 250 j=n0,n
      y(j)=y(j)+3.0*fc*y(j)*( c-0.5*y(j)-1.0/(c-y(j)) )
 250  continue

      if(ik.eq.3) then
c       vi integrerer fram hastighetspotensialet
c       hastpot. settes lik null i punkt n

        y(n)=0.0

        do 200 j=n-1,n0,-1
        y(j)=y(j+1)-dx*y(j)
  200   continue
      end if
              


      RETURN
      END





*****************************************************************************
*
*               S O L I P R O D
*
*
*     Tabulerer verdier svarende til en soliton-form ved kall p} 'SOLITON',
*     men gj|r tre ting i tillegg:
*                 1. Trunkerer solitonen der h|yden er mindre  enn amp*eps
*                 2. fyller y med nuller utenfor trunkert omr}de.
*                 3. S|rger for  at soliton kalles med increment (dx) som
*                    aldri overstiger 1.0. Dette kan v}re viktig for
*                    konvergens av underliggende likniongsl|sere.
*     Parameterene svarer til de for 'soliton' untatt eps som alts} gir
*     trunkerings-kriteriet.
**************************************************************************
      subroutine soliprod(Y,n0,N,DX,AMP,TOPP,C,UA,ik,eps)   
      INTEGER n0,N,ik
      REAL Y(N0:N),DX,AMP,topp,C,UA,eps
*---------------------------------------------------------------------------
      real fdumm

      fdumm=0.0

      call soliprgen(Y,n0,N,DX,AMP,TOPP,C,UA,ik,fdumm,eps)   

      return
      end




*****************************************************************************
*
*               S O L I T O N
*
*
*     Tabulerer verdier svarende til en soliton-form. Det antas at solitonen
*     forplanter seg mot |kende indekser i y. Har den motsatt forplantning
*     beholdes verdiene for overflatehevninger mens foretegnet byttes for
*     hastighet og potensial. Den aktuelle solitonl|sningen er beskrevet
*     i Pedersen og Rygg (1987) og Pedersen (1988).

*     parametere:
*             y  - array som inneholder verdiene                          O
*             n0,n - arraygrenser for y i det kallende programmet         I
*                    maksimalt antall punkter er 1000.
*             dx - avstand (regnet i dyp) mellom punktene.                I
*             amp - maks overflatehevning (regnet i dyp)                  I
*             topp - posisjon av makspunkt, merk at topp=0                I
*                    vil plassere makspunkt hos y(0), topp=dx
*                    hos y(1) etc. Toppen kan godt v{re lokalisert
*                    utenfor arrayendene. V{r oppmerksom p{ at y(1) svarer
*                    til ulike fysiske posisjoner for ulike ukjente i ett
*                    alternerende gitter.
*             c - forplantningshastighet for solitonen (skalert med       O
*                 line{r gruntvannshastighet)
*             ua - maksverdi for hastighet, skalert som c                 O
*             ik - styringsparameter                                      I 
*                 =1  :  Y fylles med overflatehevninger
*                 =2  :  Y fylles med hastigheter
*                 =3  :  y fylles med potensialverdier. I dette tilfellet
*                          er verdiene satt slik at midtpunktdifferens gir
*                          eksakte hastigheter, og y(n) er valgt lik null 
*                          (potensialet er jo bestemt bare p} en konstant n{r)
*                 =4  :  Som 3, bortsett fra at verdiene finnes ved integrasjon
*                        av en spline-interpolant for hastigheten.
**************************************************************************
      SUBROUTINE SOLITON(Y,n0,N,DX,AMP,Topp,C,UA,ik)
      INTEGER n0,N,ik
      REAL Y(N0:N),DX,AMP,topp,C,UA
*---------------------------------------------------------------------------
      real fcdum,eps
      fcdum=0.0
      eps=0.0001
      call SOLITGEN(Y,n0,N,DX,AMP,Topp,C,UA,ik,fcdum,eps)

      RETURN
      END


*********************************************************************
*     u    -   verdi av hastighet i soliton                       I
*     a    -   amplitude av soliton                               I
*     ier  -   feil-parameter ; ier=1 svarer til ulovlig          O
*              verdi p} u
*     uder -   verdi for du/dx                                    O
*
**********************************************************************
      function uder(u,a,ier)
      integer ier
      real u,a,uder
*--------------------------------------------------------------------
      real bhast,hj,c,ua

      c=bhast(a)

      ua=c*a/(1.0+a)

      if(u.ge.ua .and. (u-ua).lt.0.001*a) then
        ier=0
        uder=0.0
        return
      end if
      if(u.gt.ua .or. u.lt.0.0) then
        ier=1
        uder=0.0
        return
      else
        ier=0
      end if

      hj=c*log(1.0-u/c)+u+0.5*c*u*u-u*u*u/6.0
      if(hj.lt.0.0) hj=0.0
      uder=sqrt(6.0*hj/c)

      return
      end




***************************************************************************
      SUBROUTINE STAPPG(VER,N0,N,A,TOPP,DKSI,eps)
      INTEGER N0,N
      REAL VER(N0:N),A,TOPP,DKSI,eps
*----------------------------------------------------------------------
      REAL UV(0:15000),XTOPP,BL,halvb,tmp,rest,ds,xuv,xver,w(50000)
      INTEGER NADD,NTP,J,nl,nmin,nmax,n1,n2,nb,nuv,iflag
                                  
      BL=halvb(a,eps)
      nl=bl/dksi+2
c     totallet er lagt til for sikring
c    
c     vi finner avstand fra topp til n{rmeste venstre (mindre) punkt      
      NTP=TOPP/DKSI
      tmp=topp +abs((ntp+1)*dksi)
c     sikrer at tmp blir positiv, samtid som avstand beholdes.
      ntp=tmp/dksi
      rest=tmp-ntp*dksi
      
      xtopp = nl*dksi +rest
      if(xtopp.gt.topp) then
        nadd=(xtopp-topp+0.01*dksi)/dksi
      else
        nadd=-(topp-xtopp+0.01*dksi)/dksi
      end if

c     vi maa naa finne hvilket intervall n1,n2 som ligger innenfor
c     trunkeringsender.

      nmax=(topp+bl)/dksi
      nmin=(topp-bl)/dksi

      n1=max(n0,nmin)
      n2=min(n,nmax)      

      DO 5 J=N0,N
      ver(J)=0.0
 5    CONTINUE                             

      if(n1.gt.n2) return


      XTOPP=TOPP+NADD*DKSI                           

      if(2*xtopp/dksi.lt.800) then
         DO 10 J=N1,N2
         UV(J+NADD)=0.0
  10     CONTINUE                             
         CALL STUFF(UV,A,XTOPP,DKSI)
         DO 100 J=N1,N2
         VER(J)=UV(J+NADD)
  100    CONTINUE                        
      else
c       for } unng} overskrivning i stuff m} vi g} via en spline interpolant
        ds=bl/100.0
        xtopp= 101*ds
        call stuff(uv,a,xtopp,ds)
        nuv=202
        uv(0)=0.0
        uv(nuv-1)=0.0
        xuv=topp-xtopp
        xver=n0*dksi
        nb=n-n0+1
        call stag(uv(0),nuv,ds,xuv,ver,nb,dksi,xver,iflag,w)
        if(iflag.gt.0) then
          write(0,*)' problemer med stag i stappg, iflag=',iflag
        end if
      end if
      RETURN
      END





******************************************************************************

      SUBROUTINE STAPP(VER,N0,N,A,TOPP,DKSI)
      INTEGER N0,N
      REAL VER(N0:N),A,TOPP,DKSI
*----------------------------------------------------------------------
      REAL UV(0:10000),XTOPP,HLP,BL
      INTEGER NADD,NTP,J,IBL
                                  
      BL=-LOG(0.0005)/SQRT(3.0*A)
      IBL=BL/DKSI
      NTP=TOPP/DKSI
      NADD=MAX(2+N+N0-2*NTP,2-N0)
      NADD=MAX(0,NADD)
      NADD=MIN(IBL,NADD)
      DO 10 J=N0,N
      UV(J+NADD)=0.0
  10  CONTINUE                             
      XTOPP=TOPP+NADD*DKSI                           
      HLP=XTOPP/DKSI
      CALL STUFF(UV,A,XTOPP,DKSI)
      DO 100 J=N0,N
      VER(J)=UV(J+NADD)
  100 CONTINUE                             
      RETURN
      END


***************************************************************
*                                                             *
* FINNER FORMEN TIL SOLITONEN VED "ANALYTISKE" BEREGNINGER.   *
* PARAMETERE:  UV   =                                         *
*              TOPP =                                         *
*              DKSI =                                         *
* -KALLER 'SIIE,SING,UEN'-                                    *
*                                                             *
***************************************************************

      SUBROUTINE STUFF(UV,AA,TOPP,DKSI)
      REAL UV(0:10000),TOPP,DKSI,AA
      EXTERNAL SING,UEN
*-------------------------------------------------------------
      REAL EPS1,EPS,DARG,S1,S0,HJ(10000),UA,C ,A
      REAL BHAST,sdum
      INTEGER NV,NH,NP,IFLAG,NF,NL,J
      COMMON/SOL/A,UA,C

 11   format(1x,'advarsel,full konvergens i soliton-beregning')
 12   format(1x,'er ikke oppnaadd ved pos. ',i1,' iflag=',i4)     
      A=AA
      C=BHAST(A)
      UA=C*A/(1.0+A)
      EPS  = 1.0E-4 
      EPS1=0.01
      NV   = TOPP/DKSI
      NH   = NV+1
      DARG = 0.5*(TOPP-NV*DKSI)
      NP   = 1
      S1   = SQRT(UA*0.7)

      sdum=0.0
      CALL SIIE(SING,HJ,NP,DARG,sdum,S1,20,EPS,IFLAG,50)    

      if(iflag.gt.0) then
        write(0,11)
        write(0,12) 1,iflag
      end if                     

      S0   = HJ(1)
      UV(NV) = UA-S0*S0
      NP   = NV-1
      DARG = DKSI*0.5

      CALL SIIE(SING,HJ,NP,DARG,S0,S1,30,EPS,IFLAG,50)
      if(iflag.gt.0) then
        write(0,11)
        write(0,12) 2,iflag
           write(0,*)'np=',np
           do 32 j=1,np
           write(0,*) j,hj(j)
 32        continue
      end if                     

      DO 100 J=1,NP
         UV(NV-J) = UA-HJ(J)*HJ(J)
  100 CONTINUE
      NF   = NV-1-NP
      NP   = NF
      IF(NP.GT.0)THEN
         S0 =-LOG(UV(NF+1))
         S1 = 20

         CALL SIIE(UEN,HJ,NP,DKSI,S0,S1,20,EPS1,IFLAG,50)
         if(iflag.gt.0) then
           write(0,11)
           write(0,12) 3 ,iflag
         end if                     

         DO 85 J=1,NP
            UV(NF+1-J) = EXP(-HJ(J))
   85    CONTINUE
      END IF
      NL   = NF+1-NP
      DARG = 0.5*(TOPP-NV*DKSI)
      NP   = 1
      S1   = SQRT(UA*0.7)
      DARG = 0.5*(NH*DKSI-TOPP)

      sdum=0.0
      CALL SIIE(SING,HJ,NP,DARG,sdum,S1,20,EPS,IFLAG,50)
      if(iflag.gt.0) then
        write(0,11)
        write(0,12) 4,iflag
      end if                     

      S0   = HJ(1)
      UV(NH) = UA-S0*S0
      NP   = NV-1
      DARG = DKSI*0.5

      CALL SIIE(SING,HJ,NP,DARG,S0,S1,20,EPS,IFLAG,50)
      if(iflag.gt.0) then
        write(0,11)
        write(0,12) 5,iflag
      end if                     

      DO 200 J=1,NP
         UV(NH+J) = UA-HJ(J)*HJ(J)
  200 CONTINUE
      NF   = NH+NP
      NP   = NV-1-NP
      IF(NP.GT.0)THEN
         S0 =-LOG(UV(NF))
         S1 = 20

         CALL SIIE(UEN,HJ,NP,DKSI,S0,S1,20,EPS1,IFLAG,50)
         if(iflag.gt.0) then
           write(0,11)
           write(0,12) 6,iflag
         end if                     

         DO 185 J=1,NP
            UV(NF+J) = EXP(-HJ(J))
  185    CONTINUE
      END IF

      RETURN
      END

***************************************************************
*                                                             *
* INTEGRANDEN I INTERGRALET SOM TABULERES AV SIIE.            *
* KALLES AV DEN EKSTERNE PROSEDYREN SIIE.                     *
* PARAMETER:  R =                                             *
*                                                             *
***************************************************************

      FUNCTION SING(R)
      REAL R
      REAL A,UA,C
      COMMON/SOL/A,UA,C
*--------------------------------------------------------------
      REAL S,F,TEM1,SF,TEM2,SING
      INTEGER I

      S    = R*R
      F    = 1.0-UA/C
      TEM1 =-6.0/C+3.0*(C*F*F-C+F*S)+S*S/C
      SF   = S/(C-UA)
      IF(SF.GT.0.01) THEN
         TEM2 = 6.0*LOG(1.0+SF)/S
      ELSE
         TEM2 = 0.0
         DO 100 I=0,4
            TEM2 = 1.0/(5.0-I)-SF*TEM2
  100    CONTINUE
         TEM2 = 6.0*TEM2/(C-UA)
      END IF

      SING = 1.0/SQRT(TEM1+TEM2)

      RETURN
      END



***************************************************************
*                                                             *
* INTEGRANDEN I INTERGRALET SOM TABULERES AV SIIE.            *
* KALLES AV DEN EKSTERNE PROSEDYREN SIIE.                     *
* PARAMETER:  T =                                             *
*                                                             *
***************************************************************

      FUNCTION UEN(T)
      REAL T
      REAL A,UA,C
      COMMON/SOL/A,UA,C
*--------------------------------------------------------------
      REAL S,SC,TEM,UEN
      INTEGER I

      S   = EXP(-T)
      SC  = S/C
      IF(S.LT.0.01) THEN
         TEM = 0.0
         DO 70 I=0,5
            TEM = 1.0/(7.0-I)+SC*TEM
   70    CONTINUE
         TEM =-6.0*TEM/(C*C)
      ELSE
         TEM = 6.0*(1.0/(C*S)+LOG(1-SC)/(S*S))
      END IF

      UEN = 1.0/SQRT(3.0-SC+TEM)

      RETURN
      END


*****************************************************************
*
*                        GKDVSOL
*      amp - amplitude( to depth)                         I
*      x - distance (length/depth) from top             I 
*      kdv - kdv sol. if .true., Serre otherwise        I
*      eta - surface elevation                          O
*      etax - x derivative of surface elevation
****************************************************************
      subroutine gkdvsol(amp,x,kdv,eta,etax)
      real amp,x,eta,etax
      logical kdv
*--------------------------------------------------------------
      integer isgn
      real fak,c,arg,e,cc,ss

      if(kdv) then
cc      standard kdv equation with ..... Yxxx=0
        c=1.0+0.5*amp
        fak=sqrt(3.0*amp)*0.5
      else
        c=sqrt(1.0+amp)
        fak=sqrt(3.0*amp)*0.5/c
      end if

      arg=fak*x
      if(arg.gt.0.0) then
        isgn=1
      else
        isgn=-1
        arg=-arg
      end if

      if(arg.gt.20.0) then
        eta=0.0
        etax=0.0
      else
        e=exp(arg)
        cc=0.5*(e+1.0/e)
        ss=0.5*(e-1.0/e)
        eta=amp/(cc*cc)
        etax=-2.0*isgn*amp*fak*(ss/cc)/(cc*cc)
      end if


      return
      end

***************************************************************************
* The subroutine computes strain parameter (alf) and Froude number (f) from
* amplitude/depth (a) according to solution of Grimshaw (1971) and Fenton (1972).
*
* Leading order surface elevation then becomes (propagation toward pos. x)
*         eta=A*sech^2(alf*(x-f*t)) 
*
***************************************************************************
      subroutine parfent(a,alf,c)
      real a,alf,c
*-------------------------------------------------
      real f2,uinf
      alf=sqrt(0.75*a)*(1.0 +a*(-5.0/8.0+a*71.0/128.0))
      f2=1.0+a*(1.0 -a*(0.05+3.0*a/70.0))
      c=sqrt(f2)
cc    This value is consistent with u=0 at infinity


cc    half wavelength -log(0.25*eps)/(2.0*alf)
      return
      end 


***************************************************************************
* The subroutine computes the surface elevation, velocities and potential
*  according to the solution of Fenton (1972) or Grimshaw (1971)
*
* Scaling: same scaling vertically and horizontally (d -equilibrium depth
*          ) , char velocity is sqrt(g*d). Zero velocity at infinity.
*
*  Velocity components are
*    u=u0+y*y*u1+y*y*y*y*u2
*    v=y*(v0+y*y*v1+y*y*y*y*v2)
*  where y is height over bottom.
*  The velocity potential is defined as zero at x=-infinity and is computed at
*  the bottom
*    fb=phi(x,0)
*  The potential at other locations are then given as
*    phi(x,y)=fb(x)+y^2*v0/2+y^4*v1/4+y^6*v2/6
*  fu is the value of fb at x=+infty
*
*  Parameters:
*   a - amplitude (max eta)                                         I
*   alf - strain parameter                                          I
*   x - value of argument (x-c*t) or (x+c*t)                        I
*  eta - surface elevation                                          O
*  etax - derivative of eta with respect to x                        O
* u0,u1,u2,v0,v1,v2,fb,fu - velocty and potential factors as explained above O 
***************************************************************************
      subroutine fent3(a,alf,x,eta,etax,u0,u1,u2,v0,v1,v2,fb,fu)
      real a,alf,x,eta,etax,u0,u1,u2,v0,v1,v2,fb,fu
*-------------------------------------------------
      real s,t,ss,tt,ssx,ttx,e1,e2,e3,fak,i1,i2,i3
      real iu1,iu2,iu3

      s=1.0/cosh(alf*x)
      t=sinh(alf*x)*s
      ss=s*s
      tt=t*t
      ssx=-2.0*alf*t*ss
      ttx=-ssx

cc    i1 is the integral of ss from x=-\infty, i2 is the integral of ss*ss etc
cc    recursion formula is used to compute values
cc    in=(t*s^(2*n-2)+ (2*n-2)*i(n-1))/(2*n-1)
      i1=t+1
      i2=(t*ss+2.0*i1)/3.0
      i3=(t*ss*ss+4.0*i2)/5.0
cc    iu1 etc is value at infinity
      iu1=2.0
      iu2=2.0*i1/3.0
      iu3=4.0*i2/5.0

      e1=ss
      e2=-0.75*ss*tt
      e3=ss*tt*(5.0/8.0 - ss*101.0/80.0)
      eta=a*(e1+a*(e2+a*e3))
      e1=ssx
      e2=-0.75*(ss*ttx+ssx*tt)
      e3=(ss*ttx+ssx*tt)*5.0/8.0 - ss*(ss*ttx+2.0*ssx*tt)*101.0/80.0
      etax=a*(e1+a*(e2+a*e3))

cc    Fenton has propagation toward decreasing x
cc    Constant term in x in solution of Fenton vanish due to coordinate system
      e1=ss
      e2=-ss*(ss-0.25)
      e3=-ss*(19.0/40.0+ss*(0.2-ss*1.2))       
      u0=a*(e1+a*(e2+a*e3))

      e2=-ss*(1.5-ss*2.25)
      e3=-ss*(-1.5+ss*(-3.75+ss*7.5))
      u1=a*a*(e2+a*e3)
    
      u2=-a*a*a*ss*(-3.0/8.0+ss*(1.0-ss)*45.0/16.0)
      
cc    Calculation of bottom potental, formula for u0 with ss**n replaced by in
cc      e1=ss
cc      e2=-ss*ss+0.25*ss
cc      e3=-19.0/40.0*ss -ss*ss*0.2+ss*ss*ss*1.2    
      e1=i1
      e2=-i2+0.25*i1
      e3=-i1*19.0/40.0 -i2*0.2 + i3*1.2    
      fb=a*(e1+a*(e2+a*e3))/alf

      e1=iu1
      e2=-iu2+0.25*iu1
      e3=-iu1*19.0/40.0 -iu2*0.2 + iu3*1.2    
      fu=a*(e1+a*(e2+a*e3))/alf

      fak=sqrt(3.0*a)*t
      e1=-ss
      e2=ss*(3.0/8.0+2.0*ss)
      e3=ss*(49.0/640.0+ss*(-17.0/20.0-ss*3.6))       
      v0=-fak*a*(e1+a*(e2+a*e3))

      e2=ss*(0.5-ss*1.5)
      e3=ss*(-13.0/16.0+ss*(-25.0/16+ss*7.5))
      v1=-fak*a*a*(e2+a*e3)
    
      v2=-fak*a*a*a*ss*(-3.0/40.0+ss*(1.125-ss*27.0/16.0))


      return
      end 


***************************************************************************
* The subroutine computes the surface elevation, velocities and potential
*  according to the solution of Fenton (1972) or Grimshaw (1971)
*
* Scaling: same scaling vertically and horizontally (d -equilibrium depth
*          ) , char velocity is sqrt(g*d). Zero velocity at infinity.
*
*  Velocity components are
*    u=u0+y*y*u1+y*y*y*y*u2
*    v=y*(v0+y*y*v1+y*y*y*y*v2)
*  where y is height over bottom.
*  The velocity potential is defined as zero at x=-infinity and is computed at
*  the bottom
*    fb=phi(x,0)
*  The potential at other locations are then given as
*    phi(x,y)=fb(x)+y^2*v0/2+y^4*v1/4+y^6*v2/6
*  and an alternative expression for u is defined as
*    phi_x=u0+y^2*fx1+y^4*fx2+y^6*fx3
*  where fx1=v0'/2 etc.
*  fu is the value of fb at x=+infty
*
*  Parameters:
*   a - amplitude (max eta)                                         I
*   alf - strain parameter                                          I
*   x - value of argument (x-c*t) or (x+c*t)                        I
*  eta - surface elevation                                          O
*  etax - derivative of eta with respect to x                        O
* u0,u1,u2,v0,v1,v2,fb,fu,fx1,fx2,fx3 - velocty and potential factors as explained above O 
***************************************************************************
      subroutine Xfent3(a,alf,x,eta,etax,u0,u1,u2,v0,v1,v2,fb,fu,fx1,fx2
     % ,fx3)
      real a,alf,x,eta,etax,u0,u1,u2,v0,v1,v2,fb,fu,fx1,fx2,fx3
*-------------------------------------------------
      real s,t,ss,tt,ssx,ttx,e1,e2,e3,fak,i1,i2,i3
      real iu1,iu2,iu3,p0,p1,p2,p0d,p1d,p2d,a2,a3

      s=1.0/cosh(alf*x)
      t=sinh(alf*x)*s
      ss=s*s
      tt=t*t
      ssx=-2.0*alf*t*ss
      ttx=-ssx
      a2=a*a
      a3=a2*a

cc    i1 is the integral of ss from x=-\infty, i2 is the integral of ss*ss etc
cc    recursion formula is used to compute values
cc    in=(t*s^(2*n-2)+ (2*n-2)*i(n-1))/(2*n-1)
      i1=t+1
      i2=(t*ss+2.0*i1)/3.0
      i3=(t*ss*ss+4.0*i2)/5.0
cc    iu1 etc is value at infinity
      iu1=2.0
      iu2=2.0*i1/3.0
      iu3=4.0*i2/5.0

      e1=ss
      e2=-0.75*ss*tt
      e3=ss*tt*(5.0/8.0 - ss*101.0/80.0)
      eta=a*(e1+a*(e2+a*e3))
      e1=ssx
      e2=-0.75*(ss*ttx+ssx*tt)
      e3=(ss*ttx+ssx*tt)*5.0/8.0 - ss*(ss*ttx+2.0*ssx*tt)*101.0/80.0
      etax=a*(e1+a*(e2+a*e3))

cc    Fenton has propagation toward decreasing x
cc    Constant term in x in solution of Fenton vanish due to coordinate system
      e1=ss
      e2=-ss*(ss-0.25)
      e3=-ss*(19.0/40.0+ss*(0.2-ss*1.2))       
      u0=a*(e1+a*(e2+a*e3))

      e2=-ss*(1.5-ss*2.25)
      e3=-ss*(-1.5+ss*(-3.75+ss*7.5))
      u1=a*a*(e2+a*e3)
    
      u2=-a*a*a*ss*(-3.0/8.0+ss*(1.0-ss)*45.0/16.0)
      
cc    Calculation of bottom potental, formula for u0 with ss**n replaced by in
cc      e1=ss
cc      e2=-ss*ss+0.25*ss
cc      e3=-19.0/40.0*ss -ss*ss*0.2+ss*ss*ss*1.2    
      e1=i1
      e2=-i2+0.25*i1
      e3=-i1*19.0/40.0 -i2*0.2 + i3*1.2    
      fb=a*(e1+a*(e2+a*e3))/alf

      e1=iu1
      e2=-iu2+0.25*iu1
      e3=-iu1*19.0/40.0 -iu2*0.2 + iu3*1.2    
      fu=a*(e1+a*(e2+a*e3))/alf

cc      fak=sqrt(3.0*a)*t
      fak=sqrt(3.0*a)
      e1=-ss
      e2=ss*(3.0/8.0+2.0*ss)
      e3=ss*(49.0/640.0+ss*(-17.0/20.0-ss*3.6))
      p0=a*(e1+a*(e2+a*e3)) 
      p0d=-a+a2*(3.0/8.0+4.0*ss)+a3*(49.0/640.0+ss*(-17.0/10.0-ss*10.8))      
      v0=-fak*t*p0
      fx1=-0.5*alf*fak*ss*(p0-2*t*t*p0d)

      e2=ss*(0.5-ss*1.5)
      e3=ss*(-13.0/16.0+ss*(-25.0/16+ss*7.5))
      p1=a*a*(e2+a*e3)
      p1d=a2*(0.5-3.0*ss)+a3*(-13.0/16.0+ss*(-25.0/8.0+ss*22.5))
      v1=-fak*t*p1
      fx2=-0.25*alf*fak*ss*(p1-2*t*t*p1d)
    
      p2=a*a*a*ss*(-3.0/40.0+ss*(1.125-ss*27.0/16.0))
      p2d=a3*(-3.0/40.0+ss*(2.25-ss*81.0/16.0))
      v2=-fak*t*p2
      fx3=-alf*fak*ss*(p2-2*t*t*p2d)/6.0

      return
      end 

*****************************************************************************
*
*               S O L F E N T O N
*
*
*     Tabulerer verdier svarende til en soliton-form. Det antas at solitonen
*     forplanter seg mot |kende indekser i y. Har den motsatt forplantning
*     beholdes verdiene for overflatehevninger mens foretegnet byttes for
*     hastighet og potensial. 
*     parametere:
*             eta,etax,u,v,phi  - arrayer som inneholder verdiene av   O
*                    feltene i overflaten 
*             um -  umiddel                                                O
*             ub,phib - felter ved bunn                                 O
*             n - arraygrenser for y i det kallende programmet         I
*                    maksimalt antall punkter er 10000. (strat er 0)
*             dx - avstand (regnet i dyp) mellom punktene.                I
*             amp - maks overflatehevning (regnet i dyp)                  I
*             topp - posisjon av makspunkt, merk at topp=x0                I
*                    vil plassere makspunkt hos y(0), topp=x0+dx
*                    hos y(1) etc. Toppen kan godt v{re lokalisert
*                    utenfor arrayendene. V{r oppmerksom p{ at y(1) svarer
*                    til ulike fysiske posisjoner for ulike ukjente i ett
*                    alternerende gitter.
*             c - forplantningshastighet for solitonen (skalert med       O
*                 line{r gruntvannshastighet)
*             fu - potential at + infinity
*             eps - trunkeringsgrense                                     I
*
**************************************************************************
      subroutine solfenton(eta,etax,u,v,phi,um,ub,phib,n,x0,dx,amp,
     %           topp,c,fu,eps)
      integer n
      real eta(0:n),etax(0:n),u(0:n),v(0:n),phi(0:n)
      real um(0:n),ub(0:n),phib(0:n),x0,dx,amp,topp,c,fu,eps
*---------------------------------------------------------------------------
      integer i
      real alf,halvb,xeff,u0,u1,u2,v0,v1,v2,d,d2

      call parfent(amp,alf,c)
      halvb=-log(0.25*eps)/(2.0*alf)

      do 100 i=0,n
        xeff=x0+i*dx-topp
        if(xeff.gt.halvb.or.xeff.lt.-halvb) then
         eta(i)=0.0
         etax(i)=0.0
         ub(i)=0.0
         u(i)=0.0
         um(i)=0.0
         v(i)=0.0
        else
          call fent3(amp,alf,xeff,eta(i),etax(i),u0,u1,u2,v0,v1,v2
     %     ,phib(i),fu)
          d=1.0+eta(i)
          d2=d*d
          ub(i)=u0
          u(i)=u0+d*d*(u1+d*d*u2)
          um(i)=u0+d*d*(u1/3.0+d*d*u2/5.0)
          v(i)=d*(v0+d*d*(v1+d*d*v2))
          phi(i)=phib(i)+d2*(0.5*v0+d2*(0.25*v1+d2*v2/6.0))
        end if
 100  continue
      return
      end





*****************************************************************************
*
*               S O L F E N T O N
*
*
*     Tabulerer verdier svarende til en soliton-form. Det antas at solitonen
*     forplanter seg mot |kende indekser i y. Har den motsatt forplantning
*     beholdes verdiene for overflatehevninger mens foretegnet byttes for
*     hastighet og potensial. 
*     parametere:
*             eta,etax,u,v,phi  - arrayer som inneholder verdiene av   O
*                    feltene i overflaten 
*             um -  umiddel                                                O
*             ub,phib - felter ved bunn                                 O
*             phix - d(phi)/dx ved overflate. Alternativ til u.         O 
*             n - arraygrenser for y i det kallende programmet         I
*                    maksimalt antall punkter er 10000. (strat er 0)
*             dx - avstand (regnet i dyp) mellom punktene.                I
*             amp - maks overflatehevning (regnet i dyp)                  I
*             topp - posisjon av makspunkt, merk at topp=x0                I
*                    vil plassere makspunkt hos y(0), topp=x0+dx
*                    hos y(1) etc. Toppen kan godt v{re lokalisert
*                    utenfor arrayendene. V{r oppmerksom p{ at y(1) svarer
*                    til ulike fysiske posisjoner for ulike ukjente i ett
*                    alternerende gitter.
*             c - forplantningshastighet for solitonen (skalert med       O
*                 line{r gruntvannshastighet)
*             fu - potential at + infinity
*             eps - trunkeringsgrense                                     I
*
**************************************************************************
      subroutine xsolfenton(eta,etax,u,v,phi,um,ub,phib,phix,n,x0,dx,
     %           amp,topp,c,fu,eps)
      integer n
      real eta(0:n),etax(0:n),u(0:n),v(0:n),phi(0:n)
      real um(0:n),ub(0:n),phib(0:n),phix(0:n),x0,dx,amp,topp,c,fu,eps
*---------------------------------------------------------------------------
      integer i
      real alf,halvb,xeff,u0,u1,u2,v0,v1,v2,d,d2,fx1,fx2,fx3

      call parfent(amp,alf,c)
      halvb=-log(0.25*eps)/(2.0*alf)

      do 100 i=0,n
        xeff=x0+i*dx-topp
        if(xeff.gt.halvb.or.xeff.lt.-halvb) then
         eta(i)=0.0
         etax(i)=0.0
         ub(i)=0.0
         u(i)=0.0
         um(i)=0.0
         v(i)=0.0
         phix(i)=0
        else
cc          call fent3(amp,alf,xeff,eta(i),etax(i),u0,u1,u2,v0,v1,v2
cc     %     ,phib(i),fu)
          call xfent3(amp,alf,xeff,eta(i),etax(i),u0,u1,u2,v0,v1,v2
     %     ,phib(i),fu,fx1,fx2,fx3)
          d=1.0+eta(i)
          d2=d*d
          ub(i)=u0
          u(i)=u0+d*d*(u1+d*d*u2)
          um(i)=u0+d*d*(u1/3.0+d*d*u2/5.0)
          v(i)=d*(v0+d*d*(v1+d*d*v2))
          phi(i)=phib(i)+d2*(0.5*v0+d2*(0.25*v1+d2*v2/6.0))
          phix(i)=u0+d2*(fx1+d2*(fx2+d2*fx3))
        end if
 100  continue
      return
      end



*****************************************************************************
*
*               ssolip
*
*
*     Tabulerer verdier svarende til en soliton-form. Det antas at solitonen
*     forplanter seg mot |kende indekser i y. Har den motsatt forplantning
*     beholdes verdiene for overflatehevninger mens fortegnet byttes for
*     hastighet og potensial. 

*     parametere:
*             y  - array som inneholder verdiene                          O
*             n0,n - arraygrenser for y i det kallende programmet         I
*             dx - avstand (regnet i dyp) mellom punktene.                I
*             amp - maks overflatehevning (regnet i dyp)                  I
*             topp - posisjon av makspunkt, merk at topp=0                I
*                    vil plassere makspunkt hos y(0), topp=dx
*                    hos y(1) etc. Toppen kan godt v{re lokalisert
*                    utenfor arrayendene. V{r oppmerksom p{ at y(1) svarer
*                    til ulike fysiske posisjoner for ulike ukjente i ett
*                    alternerende gitter.
*             c - wave celerity                                           O
*             ivis - determine type of solitary wave                      I
*                   = 0 KdV
*                     1 Serre with no correction
*                     2 Serre with correction of Nwogu type
*                     3 Serre which reproduces linear dispersion incl. O(k^4)
*                     4 Serre with a specified value on kappa 
*             gamma - parameter in equations, selected automatically       I
*                     in case ivis <=3                      
*
**************************************************************************
      subroutine ssolip(Y,n0,N,DX,AMP,TOPP,c,ivis,gamma)   
      INTEGER n0,N,ivis
      REAL Y(N0:N),DX,AMP,topp,c,gamma
*---------------------------------------------------------------------------
      integer i
      real fak,arg,e,rgamma,z0,cc

      if(ivis.le.1) rgamma=0.0
      if(ivis.eq.2) rgamma=-0.0566862
      if(ivis.eq.3) rgamma=-1.0/15.0
      if(ivis.eq.4) rgamma=gamma

      
      if(ivis.eq.0) then
cc      standard kdv equation with ..... Yxxx=0
        c=1.0+0.5*amp
        fak=sqrt(3.0*amp)*0.5
      else
        c=sqrt(1.0+amp+0.75*rgamma*amp*amp)
        fak=sqrt(3.0*amp)*0.5*sqrt(1.0+3.75*rgamma*amp)/c
      end if


      do 100 i=n0,n
        arg=abs(fak*(i*dx-topp))
        if(arg.gt.20.0) then
          z0=0.0
        else
          e=exp(arg)
          cc=(e+1.0/e)*0.5
          z0=1.0/(cc*cc)
        end if
        y(i)=amp*z0*(1.0+amp*11.25*rgamma*(1.0-z0))
  100 continue

      return
      end


*****************************************************************
*
*                        ZKDVSOL
*      amp - amplitude( to depth)                         I
*      x - distance (length/depth) from top             I 
*      ivis -  determine type of solitary wave        I
*                   = 0 KdV
*                     1 Serre with no correction
*                     2 Serre with correction of Nwogu type
*                     3 Serre which reproduces linear dispersion incl. O(k^4)
*                     4 Serre with a specified value on kappa 
*             gamma - parameter in equations, selected automatically       I
*                     in case ivis <=3                      
*      eta - surface elevation                          O
*      etax - x derivative of surface elevation         O
****************************************************************
      subroutine zkdvsol(amp,x,ivis,gamma,eta,etax)
      integer ivis
      real amp,x,gamma,eta,etax
*--------------------------------------------------------------
      real fak,c,arg,e,cc,ss,z0,z0x,rgamma

      if(ivis.le.1) rgamma=0.0
      if(ivis.eq.2) rgamma=-0.0566862
      if(ivis.eq.3) rgamma=-1.0/15.0
      if(ivis.eq.4) rgamma=gamma

      if(ivis.eq.0) then
cc      standard kdv equation with ..... Yxxx=0
        c=1.0+0.5*amp
        fak=sqrt(3.0*amp)*0.5
      else
        c=sqrt(1.0+amp+0.75*rgamma*amp*amp)
        fak=sqrt(3.0*amp)*0.5*sqrt(1.0+3.75*rgamma*amp)/c
      end if

      arg=fak*x

      if(arg.gt.20.0) then
        eta=0.0
        etax=0.0
      else
        e=exp(arg)
        cc=0.5*(e+1.0/e)
        ss=0.5*(e-1.0/e)
        z0=1.0/(cc*cc)

        z0x=-2.0*fak*ss/(cc*cc*cc)
        eta=amp*z0*(1.0+amp*11.25*rgamma*(1.0-z0))
        etax=z0x*amp*(1.0+amp*11.25*rgamma*(1.0-2.0*z0) )
      end if


      return
      end




**************************************************************
      SUBROUTINE KVAIN(YM,Y0,YP,EM,XM,IER) 
      REAL YM,Y0,YP,EM,XM               
      INTEGER IER
*--------------------------------------------------------------
      REAL A,B
               
      IER=0
      A     = 0.5*(YP-YM)
      B     = 0.5*(YP+YM-2.0*Y0)
      IF(B.EQ. 0.0) THEN
         IER=1
      ELSE
         XM =-0.5*A/B
         EM   = Y0+A*XM+B*XM*XM
         IF( (EM-Y0) .GT. 1.5*MAX(Y0-YM,Y0-YP)) IER=2
      END IF

      RETURN
      END




******************************************************************************
      SUBROUTINE GLMAX(Y,YINP,NS,N,K,GR,AMP,XA,NA)
      INTEGER NS,N,K,NA
      REAL Y(0:NS,0:K),YINP(0:N+1),XA(20),AMP(20),GR
*----------------------------------------------------------------------------
      INTEGER NP2,J

      NP2=N+1

      CALL LIMAX(Y(0,K),YINP,NP2,GR,AMP,XA,NA)                          
                                                                       
      DO 100 J=1,NA
  100 XA(J)=XA(J)-1.0

      RETURN
      END




******************************************************************************
*                                                                            *
*             L I M A X                                                      *
*                                                                            *
*     MAKSIMALPUNKTER FOR SPLINE-INTERPOLANTEN TIL ET TALLSETT FINNES.       *
*     PARAMETERE:                                                            *
*            YR - YR(1)...YR(N) UTGJOER TALLSETTET.            I             *
*            YINP - YINP(0)...YINP(N+1) ER SPLINE KOEFF.       O             *
*            N - ANTALL PUNKTER.                               I             *
*            GR - MAKSP. LOKALISERES I ET INTERVALL (K,K+1)                  *
*                 BARE NAAR YR(K) EL. YR(K+1)>GR               I             *
*            AMP ,XA - MAKSVERDIER AMP(1)...AMP(NA) ER FUNNET                *
*                 I PUNKTER XA(1)<XA(2)<...XA(NA). KOORDINAT-                *
*                 AKSEN ER SLIK AT YR(K) ER VERDIEN I K.       O             *
*            NA - ANTALL LOKALISERTE MAKSIMA. MAKSIMALT                      *
*                   ANTALL ER 20.                              O             *
*                                                                            *
*     RUTINEN BENYTTER "BSPL","SPB","KVAIN","TRI"                            *
******************************************************************************

      SUBROUTINE LIMAX(YR,YINP,N,GR,AMP,XA,NA)                               
      IMPLICIT LOGICAL (A-Z)
      INTEGER N,NA
      REAL YR(N),YINP(0:N+1),XA(20),AMP(20),GR
*-------------------------------------------------
      REAL DER1,DER2,QPP,QP,QM,QMM,Y0,YM,YP,XM
      REAL YMAX
      INTEGER I,NM ,IER,IER1,IM

      NA=0
      NM=N-1

      DER1=0.5*(4.0*YR(2)-3.0*YR(1)-YR(3))
      DER2=-0.5*(4.0*YR(N-1)-3.0*YR(N) -YR(N-2))
      CALL BSPL(YR,YINP,N,1.0,DER1,DER2)

      Y0=YR(1)
      YP=YR(2)

      DO 200 I=2,NM
      YM=Y0
      Y0=YP
      YP=YR(I+1)
      IF( Y0.GE.YP .AND. Y0.GT.YM .AND. Y0.GT.GR) THEN
        NA=NA+1
        IF(YM.GT.YP) THEN
          IM=I-1
        ELSE
          IM=I
        END IF
        QMM=YINP(IM-1)
        QM=YINP(IM)
        QP=YINP(IM+1)
        QPP=YINP(IM+2)
        CALL SPB(QMM,QM,QP,QPP,XM,YMAX,IER)
        IF(IER.EQ.0) THEN
          XA(NA)= IM+XM
          AMP(NA)= YMAX
        else
          IM=2*I-1-IM
          QMM=YINP(IM-1)
          QM=YINP(IM)
          QP=YINP(IM+1)
          QPP=YINP(IM+2)
          IER1=IER
          CALL SPB(QMM,QM,QP,QPP,XM,YMAX,IER)
          IF(IER.EQ.0) THEN
            XA(NA)= IM+XM
            AMP(NA)= YMAX
          else                             
            CALL KVAIN(YM,Y0,YP,YMAX,XM,IER)    
            IF(IER.EQ.0) THEN
               XA(NA)=I+XM
               AMP(NA)=YMAX
            ELSE
               XA(NA)=I
               AMP(NA)=Y0
            end if
          END IF
        END IF
C
COMMENT: MAKSIMALT TILLATT ANTALL MAXIMA ER SATT LIK 20
C
      IF(NA.EQ.20)RETURN
      END IF
  200 CONTINUE
      RETURN
      END



******************************************************************************
*                                                                            *
*             N E W M A X                                                    *
*                                                                            *
*     MAKS(MIN)IMALPUNKTER FOR SPLINE-INTERPOLANTEN TIL ET TALLSETT FINNES.
*
*     PARAMETERE:                                                            *
*            YR - YR(1)...YR(N) UTGJOER TALLSETTET.            I             *
*            YINP - YINP(0)...YINP(N+1) ER SPLINE KOEFF.       O             *
*            N - ANTALL PUNKTER.                               I             *
*            x0 - posisjon av punkt 1
*            dx - gitteravstand
*            ivis - type maks                                  I
*                 <=1    fra venstre, posetive ekstrema
*                  =2    fra hoyre, posetive ekstrema
*                  =3    fra venstre, negative ekstrema
*                  =4    fra hoyre, negative ekstrema
*            GR - MAKSP. LOKALISERES I ET INTERVALL (K,K+1)                  *
*                 BARE NAAR YR(K) EL. YR(K+1)>GR               I
*                 (dersom ivis >2 YR(k) <-GR)                                *
*            AMP ,XA - MAKSVERDIER AMP(1)...AMP(NA) ER FUNNET                *
*                 I PUNKTER XA(1)<XA(2)<...XA(NA). KOORDINAT-                *
*                 AKSEN ER SLIK AT YR(K) ER VERDIEN I K.       O             *
*            namax - maximum number of max-points
*            NA - ANTALL LOKALISERTE MAKSIMA. MAKSIMALT                      *
*                   ANTALL ER 20.                              O            
*            w -arbeidsarray                                                 *
*                                                                            *
*     RUTINEN BENYTTER "BSPL","SPB","KVAIN","TRI"                            *
******************************************************************************

      SUBROUTINE NEWMAX(YR,YINP,N,x0,dx,ivis,GR,AMP,XA,namax,NA,w) 
      INTEGER N,ivis,namax,NA
      REAL YR(N),YINP(0:N+1),x0,dx,XA(namax),AMP(namax),GR,w(5*(n+2))
*-------------------------------------------------
      REAL DER1,DER2,QPP,QP,QM,QMM,Y0,YM,YP,XM
      REAL YMAX,xend,steg
      INTEGER I,NM ,IER,IER1,IM,nw,isgn

      
      if(ivis.le.2) then
        isgn=1
      else
        isgn=-1
      end if

      if(ivis.eq.1 .or. ivis.eq.3) then
        xend=x0
        steg=dx
        do 100 i=1,n
          w(i)=isgn*yr(i)
 100    continue
      else
        xend=x0+(n-1)*dx
        steg=-dx
        do 150 i=1,n
          w(i)=isgn*yr(n+1-i)
 150    continue
      end if
   

      nw=n+1
      NA=0
      NM=N-1

      DER1=0.5*(4.0*W(2)-3.0*W(1)-W(3))
      DER2=-0.5*(4.0*W(N-1)-3.0*W(N) -W(N-2))
      CALL bsplw(W,YINP,N,1.0,DER1,DER2,w(nw))

      Y0=W(1)
      YP=W(2)

      DO 200 I=2,NM
      YM=Y0
      Y0=YP
      YP=W(I+1)
      IF( Y0.GE.YP .AND. Y0.GT.YM .AND. Y0.GT.GR) THEN
        NA=NA+1
        IF(YM.GT.YP) THEN
          IM=I-1
        ELSE
          IM=I
        END IF
        QMM=YINP(IM-1)
        QM=YINP(IM)
        QP=YINP(IM+1)
        QPP=YINP(IM+2)
        CALL SPB(QMM,QM,QP,QPP,XM,YMAX,IER)
        IF(IER.EQ.0) THEN
          XA(NA)= IM+XM
          AMP(NA)= YMAX
        else
          IM=2*I-1-IM
          QMM=YINP(IM-1)
          QM=YINP(IM)
          QP=YINP(IM+1)
          QPP=YINP(IM+2)
          IER1=IER
          CALL SPB(QMM,QM,QP,QPP,XM,YMAX,IER)
          IF(IER.EQ.0) THEN
            XA(NA)= IM+XM
            AMP(NA)= YMAX
          else                             
            CALL KVAIN(YM,Y0,YP,YMAX,XM,IER)    
            IF(IER.EQ.0) THEN
               XA(NA)=I+XM
               AMP(NA)=YMAX
            ELSE
               XA(NA)=I
               AMP(NA)=Y0
            end if
          END IF
        END IF
        xa(na)=xend+(xa(na)-1)*steg
        amp(na)=isgn*amp(na)
C
COMMENT: MAKSIMALT TILLATT ANTALL MAXIMA ER SATT LIK namax
C
        IF(NA.EQ.namax)return
      END IF
  200 CONTINUE
      
      RETURN
      END




*****************************************************************
      SUBROUTINE SPB(QMM,QM,QP,QPP,XM,YMAX,IER)
      IMPLICIT LOGICAL(A-Z)
      REAL QMM,QM,QP,QPP,XM,YMAX
      INTEGER IER
*--------------------------------------------------------------
      REAL A,B,C,D,DISK,S1,S2,GR,AS1,AS2,DTOL

      IER=0
      GR=(ABS(QP)+ABS(QMM))*1.0E-9

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     UTREGNING AV TOLERANSE FOR DISKRIMAND I DEN
C     DERIVERTE BASERER SEG PAA DIFFERANSE FOR
C     ANNEN DERIVERT.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      DTOL=QP-2.0*QM+QMM
      DTOL=0.001*DTOL*DTOL

      A=( QMM +4.0*QM +QP)/6.0
      B=0.5*( QP -QMM )
      C=0.5*( QMM -2*QM +QP )
      D=( QPP -3*QP +3*QM -QMM )/6.0

      IF(ABS(D).LT.GR) THEN
        IF(ABS(C).LT.GR) THEN
          IER=3
          RETURN
        END IF
        XM=-B*0.5/C
        IF(ABS(XM-0.5).GT.(0.5+GR)) THEN
          IER=4
          RETURN
        END IF
        YMAX=A+B*XM+C*XM*XM
        RETURN
      END IF

      DISK=C*C-3.0*B*D

      IF(DISK.LT.0) THEN
        IER=1
        RETURN
      END IF

      IF(DISK.LT.DTOL) THEN
        IER=2
        RETURN
      END IF

      S1=(-C+SQRT(DISK))/(3.0*D)
      S2=(-C-SQRT(DISK))/(3.0*D)
      AS1=ABS(S1-0.5)
      AS2=ABS(S2-0.5)

      IF(AS1.LT.AS2) THEN
        IF(AS1.GT.(0.5+GR)) THEN
          IER=-1
          RETURN
        END IF
        IF(AS2.LT.(0.5-GR)) THEN
          IER=-2
          RETURN
        END IF
        XM=S1
      ELSE
        IF(AS2.GT.(0.5+GR)) THEN
          IER=-1
          RETURN
        END IF
        IF(AS1.LT.(0.5-GR)) THEN
          IER=-2
          RETURN
        END IF
        XM=S2
      END IF

      YMAX=A +B*XM +C*XM*XM +D*XM*XM*XM

      RETURN
      END




******************************************************************************
*                                                                            *
*             L I M A W                                                      *
*                                                                            *
*     MAKSIMALPUNKTER FOR SPLINE-INTERPOLANTEN TIL ET TALLSETT FINNES.       *
*     PARAMETERE:                                                            *
*            YR - YR(1)...YR(N) UTGJOER TALLSETTET.            I             *
*            YINP - YINP(0)...YINP(N+1) ER SPLINE KOEFF.       O             *
*            N - ANTALL PUNKTER.                               I             *
*            GR - MAKSP. LOKALISERES I ET INTERVALL (K,K+1)                  *
*                 BARE NAAR YR(K) EL. YR(K+1)>GR               I             *
*            AMP ,XA - MAKSVERDIER AMP(1)...AMP(NA) ER FUNNET                *
*                 I PUNKTER XA(1)<XA(2)<...XA(NA). KOORDINAT-                *
*                 AKSEN ER SLIK AT YR(K) ER VERDIEN I K.       O             *
*            NA - ANTALL LOKALISERTE MAKSIMA. MAKSIMALT                      *
*                   ANTALL ER 20.                              O             *
*            WRK - arbeisdomraade                              D
*                                                                            *
*     RUTINEN BENYTTER "BSPL","SPB","KVAIN","TRI"                            *
******************************************************************************

      SUBROUTINE LIMAW(YR,YINP,N,GR,AMP,XA,NA,WRK)             
      INTEGER N,NA
      REAL YR(N),YINP(0:N+1),XA(20),AMP(20),GR,wrk(4*N)
*-------------------------------------------------
      REAL DER1,DER2,QPP,QP,QM,QMM,Y0,YM,YP,XM
      REAL YMAX,ein
      INTEGER I,NM ,IER,IER1,IM

      NA=0
      NM=N-1
      ein=1.0

      DER1=0.5*(4.0*YR(2)-3.0*YR(1)-YR(3))
      DER2=-0.5*(4.0*YR(N-1)-3.0*YR(N) -YR(N-2))
      CALL BSPLW(YR,YINP,N,ein,DER1,DER2,WRK)

      Y0=YR(1)
      YP=YR(2)

      DO 200 I=2,NM
      YM=Y0
      Y0=YP
      YP=YR(I+1)
      IF( Y0.GE.YP .AND. Y0.GT.YM .AND. Y0.GT.GR) THEN
        NA=NA+1
        IF(YM.GT.YP) THEN
          IM=I-1
        ELSE
          IM=I
        END IF
        QMM=YINP(IM-1)
        QM=YINP(IM)
        QP=YINP(IM+1)
        QPP=YINP(IM+2)
        CALL SPB(QMM,QM,QP,QPP,XM,YMAX,IER)
        IF(IER.EQ.0) THEN
          XA(NA)= IM+XM
          AMP(NA)= YMAX
        else
          IM=2*I-1-IM
          QMM=YINP(IM-1)
          QM=YINP(IM)
          QP=YINP(IM+1)
          QPP=YINP(IM+2)
          IER1=IER
          CALL SPB(QMM,QM,QP,QPP,XM,YMAX,IER)
          IF(IER.EQ.0) THEN
            XA(NA)= IM+XM
            AMP(NA)= YMAX
          else                             
            CALL KVAIN(YM,Y0,YP,YMAX,XM,IER)    
            IF(IER.EQ.0) THEN
               XA(NA)=I+XM
               AMP(NA)=YMAX
            ELSE
               XA(NA)=I
               AMP(NA)=Y0
            end if
          END IF
        END IF
C
COMMENT: MAKSIMALT TILLATT ANTALL MAXIMA ER SATT LIK 20
C
      IF(NA.EQ.20)RETURN
      END IF
  200 CONTINUE
      RETURN
      END




******************************************************************************
      SUBROUTINE GLMAW(Y,YINP,NS,N,K,GR,AMP,XA,NA,WRK)
      INTEGER NS,N,K,NA
      REAL Y(0:NS,0:K),YINP(0:N+1),XA(20),AMP(20),GR,wrk(4*n+4)
*----------------------------------------------------------------------------
      INTEGER NP2,J

      NP2=N+1

      CALL LIMAW(Y(0,K),YINP,NP2,GR,AMP,XA,NA,wrk)                          
                                                                       
      DO 100 J=1,NA
  100 XA(J)=XA(J)-1.0

      RETURN
      END



**************************************************************************
*
*                     B L I N 2 G
*   
*     Beregner funksjonsverdier av en bilinear interpolant. 
*      parametere:
*              X,y - koordinater                                        I
*              X0,y0 - posisjon tilh|rende B-spline med koeff. bs(1,1)  I
*              xd,dy - punktavstand                                      I
*              B - array med funksjonsverd (bi-lin koeff).         I
*              N,m - antall interpolasjonspunkter                       I
*              eps - justerings-grense, faller et punkt mindre enn      I
*                    eps*dx etc. utenfor, trekkes det inn.
*              fdef - verdi som settes dersom x,y er utenfor
*              f -  verdi av interpolant                                O
*****************************************************************************
      subroutine blin2g(X,y,X0,y0,Dx,dy,B,N,m,eps,fdef,f)   
      INTEGER N,m
      REAL B(N,m),X,Y,X0,Y0,Dx,dy,eps,fdef,f
*------------------------------------------------ 
      INTEGER NX,ny,nn,i,j,nxp,nyp
      REAL Sy,sx,xv,yv

                 

      XV=(X-X0)/DX+1.0                          
      NX=XV         
      SX=XV-NX      

      if(nx.le.0 .or. nx.gt.n) then
        if(nx.eq.0 .and. sx .gt.(1.0-eps)) then
           nx=1
           sx=0.0
        else
          if(nx.eq.n .and. sx .lt.eps) then
           nx=n-1
           sx=1
          else 
            write(0,*)'x-punkt utenfor interval i blin2g'
            f=fdef
            return
          end if
        end if
      end if

      YV=(Y-Y0)/DY+1.0                          
      NY=YV         
      SY=YV-NY       


      if(ny.le.0 .or. ny.gt.m) then
        if(ny.eq.0 .and. sy .gt.(1.0-eps)) then
           ny=1
           sy=0.0
        else
          if(ny.eq.m .and. sy .lt.eps) then
           ny=m-1
           sy=1
          else 
            write(0,*)'y-punkt utenfor interval i blin2g'
            f=0
            return
          end if
        end if
      end if



      nxp=nx+1
      nyp=ny+1
      f=(1.0-sx)*(1.0-sy)*b(nx,ny)+sx*(1.0-sy)*b(nxp,ny)
     %  +sx*sy*b(nxp,nyp)+(1.0-sx)*sy*b(nx,nyp)


      RETURN        
      END           

*******************************************************
*             BASJOIN
*
*     Computes the value of a normalized spline S(x) defined
*     on -1<x<1 according to
*       S(1)=f; S'(1)=fd; S''(1)=fdd;  S(-1)=S'(-1)=S''(-1)=0 
***************************************************************
      function basjoin(f,fd,fdd,x)
      real f,fd,fdd,x,basjoin
*-------------------------------------------------------
      real a0,a1,a2,xm,xp

      xm=x-1.0
      xp=x+1.0

      a0=f*0.125
      a1=fd*0.125-1.5*a0
      a2=fdd/16.0-1.5*a1-0.75*a0

      basjoin=xp*xp*xp*(a0+xm*(a1+xm*a2))
      return
      end


*******************************************************
*             SPJOIN
*
*     Computes the value of a spline P(x) defined
*     on x0-h<x<x0+h according to
*       P(x+h)=r; P'(x+h)=rd; P''(x+h)=rdd;  
*       P(x-h)=l; P'(x-h)=ld; P''(x-h)=ldd;  
***************************************************************
      function spjoin(l,ld,ldd,r,rd,rdd,h,x0,x)
      real r,rd,rdd,l,ld,ldd,h,x0,x,spjoin
*-------------------------------------------------------
      real f,fd,fdd,arg,left,right,rh,basjoin
      external basjoin

      rh=1.0/h
      f=l
      fd=-ld*h
      fdd=ldd*h*h
      arg=(x0-x)*rh
      left=basjoin(f,fd,fdd,arg)

      f=r
      fd=rd*h
      fdd=rdd*h*h
      arg=(x-x0)*rh
      right=basjoin(f,fd,fdd,arg)
 
      spjoin=left+right
      return
      end


*******************************************************
*             glstep
*
*     Computes the value of two shelves joined by a spline P(x), giving
*     a continuous second derivative
*        ya,yb - left and right function values
*        xa,xb - transition interval
***************************************************************
      function glstep(xa,xb,ya,yb,x)
      real ya,yb,xa,xb,x,glstep
*-------------------------------------------------------
      real xm,h,dn,spjoin
      external spjoin

      xm=0.5*(xa+xb)
      h=0.5*(xb-xa)
      dn=0.0

      if(x.lt.xa) then
        glstep=ya
      else
        if(x.lt.xb) then
          glstep=spjoin(ya,dn,dn,yb,dn,dn,h,xm,x)
        else
          glstep=yb
        end if
      end if

      return
      end

*******************************************************
*             smpolyg
*
*     Computes a smoothed polygon with sharp angles replaced by
*     splines with a continuous second derivative
*        xv,yv - x, y values at corners
*        lv - transition intervals for each point
*        n - number of points
*             The interval must not be overlapping
*        x -coordinate at which value is to be computed
*        ia,ib - interval to search for x 
*        ileft - 1...n-1  : OK. it is the position left of x 
*                <=0  : x to small
*                >n  : x to large
*        eps - toleranse for x being outside domain
***************************************************************
      function smpolyg(xv,yv,lv,n,x,ia,ib,ileft,eps)
      integer n,ia,ib,iflag
      real xv(n),yv(n),lv(n),x,smpolyg,eps
*-------------------------------------------------------
      integer i,ipos,ileft
      real spjoin
      real xa,xb,l,ld,ldd,r,rd,rdd,h,x0
      external spjoin


      if(x.le.xv(ia)) then
        if(xv(ia)-x.lt.eps) then
          ileft=ia
          smpolyg=yv(ia)
        else
          ileft=0
        end if
        return
      end if

      if(x.ge.xv(ib)) then
        if(x-xv(ib).lt.eps) then
          ileft=ib
          smpolyg=yv(ib)
        else
          ileft=2*ib
        end if
        return
      end if

      ipos=ia

 100  continue

      if(ipos.ge.n) then
        ileft=ipos
        return
      end if



      xa=xv(ipos)
      xb=xv(ipos+1)
      if(x.ge.xa .and. x.le.xb) go to 200
      ipos=ipos+1
      go to 100

 200  continue

      ileft=ipos
      
      if( (x.lt.xa+lv(ipos)).and.ipos.gt.1) then
       ipos=ipos     
      else
         if(x.gt.xb-lv(ipos+1).and.ipos.lt.(n-1)) then
           ipos=ipos+1
         else
          ld=(yv(ipos+1)-yv(ipos))/(xb-xa)
          smpolyg= yv(ipos)+ld*(x-xa)
          return
         end if
      end if

       
       h=lv(ipos)
       x0=xv(ipos)
       ld=(yv(ipos)-yv(ipos-1))/(xv(ipos)-xv(ipos-1))
       l=yv(ipos)-h*ld
       rd=(yv(ipos+1)-yv(ipos))/(xv(ipos+1)-xv(ipos))
       r=yv(ipos)+h*rd
       ldd=0.0
       rdd=0.0
       smpolyg=spjoin(l,ld,ldd,r,rd,rdd,h,x0,x)

      return
      end

****************************************************************************
*
*     Fill an array with values from smoothed polygon
*         y(n)  - array for  values                                      O
*         n  - array-grense                                                 I
*         dx - gitterinkrement (betydning bare nar noun=.false.)          I
*         x - posisjoner  av y -verdier, naar noun=.false.                  I
*             benyttes bare x(1) som er pos av y(1)
*         noun - .false. angir uniformt gitter                            I
*         undef - value set otside polygon
*         eps - tolerance
*******************************************************************
      subroutine smarrp(xv,yv,lv,nv,y,n,dx,x,noun,undef,eps)
      integer nv,n
      real xv(nv),yv(nv),lv(nv),y(n),x(n),dx,undef,eps
      logical noun
*--------------------------------------------------
      integer ia,ib,i,ileft,ifak
      real xx,dxeff,val,smpolyg
      external smpolyg

      if(noun) then
        dxeff=0.0
        ifak=1
      else
        dxeff=dx
        ifak=0
      end if

      ia=1

      do 100 i=1,n
       xx=x(1+(i-1)*ifak)+(i-1)*dxeff
       val=smpolyg(xv,yv,lv,nv,xx,ia,nv,ileft,eps)
       if(ileft.ge.1 .and.ileft.le.nv) then
        ia=ileft
        y(i)=val
       else
        y(i)=undef
       end if
 100  continue

      return
      end 
************************************************************************
      subroutine noreg(y,n,a,b,c,d,dev1,dev2,dev3)
      integer n
      real a,b,c,d,dev1,dev2,dev3,y(n)
*----------------------------------------------------------------------
      integer i
      real sig0,sig1,sig2,sig3,sig4,rn,an,bn,syy,s1,s2,s0
      real p0,p1,p2,yi,ci,hj,g1,g2,g3,cn,s3,p3

c
c   definisjoner
c      sig2 = sum( i*i/n*n)  i=1.....n   etc.
c      ci = (i/n-an)   der an=(n+1)/(2n)
c   utviklingsfunksjoner er gitt ved:
c      g0(i)=0
c      g1(i) = ci
c      g2(i) = ci*ci -bn
c      g3(i) = ci*ci*ci -cn*ci
c   der
c      bn = sum(ci**3)/sum(ci)
c      cn = sum(ci**4)/sum(ci**2)
c   kvadatnormer av utv. funk er gitt ved
c      p0 = sum(g0*g0)
c      p1 =sum(g1*g1)  etc.
c
  
      if(n.le.0) then
        write(*,*)'ulovlig verdi: n=',n
        return
      end if

      rn=1.0/n
      an=0.5*(1.0+rn)
      bn=(1.0-rn*rn)/12.0
      sig0=n
      sig1=0.5*(n+1.0)
      sig2=(n+1)*(2.0+rn)/6.0
      sig3=0.25*(n+1.0)*(1.0+rn)
      sig4=sig2*(3.0*(1.0+rn)-rn*rn)*0.2

      p0=sig0
      p1=sig2-2*an*sig1+an*an*sig0
      cn=sig4-4*an*sig3+6*an*an*sig2-4*an*an*an*sig1+an*an*an*an*sig0
      p2=cn-2.0*bn*p1+bn*bn*sig0
      cn=cn/p1

      syy=0
      s0=0
      s1=0
      s2=0
      s3=0
      p3=0

      do 100 i=1,n
      yi=y(i)
      ci=i*rn-an
      g1=ci
      g2=ci*ci-bn
      g3=ci*(ci*ci-cn)
      s0=s0+yi
      syy=syy+yi*yi
      s1=s1+g1*yi
      s2=s2+g2*yi
      s3=s3+g3*yi
      p3=p3+g3*g3

 100  continue



      a=s0/p0
      b=s1/p1
      c=s2/p2
      d=s3/p3

      dev1=0.0
      dev2=0.0
      dev3=0.0 

      do 200 i=1,n
      yi=y(i)
      ci=i*rn-an
      g1=ci
      g2=ci*ci-bn
      g3=ci*(ci*ci-cn)
      hj=yi-a-b*g1
      dev1=dev1+hj*hj
      hj=hj-c*g2
      dev2=dev2+hj*hj
      hj=hj-d*g3
      dev3=dev3+hj*hj
 200  continue

ctt      dev=syy-a*s0-b*s1-c*s2
ctt
ctt      write(*,*)'a,b,c,dev1,dev2 i nyreg:'
ctt      write(*,*)a,b,c,dev1,dev2
ctt
      dev1=sqrt(dev1/n)
      dev2=sqrt(dev2/n)
      dev3=sqrt(dev3/n)   

      return
      end



*******************************************************************
      subroutine tra(n,dy,y0,a,b,c,d,id)
      integer n,id
      real dy,y0,a,b,c,d
*-----------------------------------------------------------------
      real hj,bn,rn,cn,an,sig0,sig1,sig2,sig3,sig4,p1
c
c      regner om til x=aa+bb*y+cc*y*y +dd*y*y*y der y=(i-1)*dy +y0
c
      rn=1.0/n
      an=0.5*(1.0+rn)
      hj=0.5*(1.0-rn)+y0/(n*dy)
      bn=(1.0-rn*rn)/12.0

      sig0=n
      sig1=0.5*(n+1.0)
      sig2=(n+1)*(2.0+rn)/6.0
      sig3=0.25*(n+1.0)*(1.0+rn)
      sig4=sig2*(3.0*(1.0+rn)-rn*rn)*0.2
      p1=sig2-2*an*sig1+an*an*sig0
      cn=sig4-4*an*sig3+6*an*an*sig2-4*an*an*an*sig1+an*an*an*an*sig0
      cn=cn/p1

      if(id.le.2) then
        d=0.0
        if(id.le.1) c=0.0
      end if

      a=a-b*hj+c*(hj*hj-bn) +d*hj*(cn-hj*hj)
      b=(b-2*c*hj +d*(3*hj*hj-cn) )*rn/dy
      c=(c -d*3*hj )*rn*rn/(dy*dy)
      d=d*rn*rn*rn/(dy*dy*dy)

ctt
ctt      write(*,*)'a,b,c,d omregnet i nyreg:'
ctt      write(*,*)a,b,c,d
ctt

      return
      end


****************************************************************************
      subroutine angleold(x,v,n,dy,ibs,id)
      integer n,ibs,id
      real x(n),v(n),dy
*-------------------------------------------------------------------------
      integer i,i1,i2,nint
      real x1,dev1,dev2,dev3,a,b,c,d,pi,y0,dydum

      pi=4.0*atan(1.0)
      x1=x(1)

      do 100 i=1,n
      x(i)=x(i)-x1
 100  continue


      do 200 i=1,n
      i1=i-ibs
      i2=i+ibs
      i1=max(1,i1)
      i2=min(n,i2)
      nint=i2-i1+1
      call noreg(x(i1),nint,a,b,c,d,dev1,dev2,dev3)
      y0=i1-i
      dydum=1.0
      call tra(nint,dydum,y0,a,b,c,d,id)
      v(i)=-180.0*atan(b/dy)/pi
 200  continue


      do 300 i=1,n
      x(i)=x(i)+x1
 300  continue

      return
      end


****************************************************************************
      subroutine angle(x,v,n,dy,ibs,id)
      integer n,ibs,id
      real x(n),v(n),dy
*-------------------------------------------------------------------------
      integer i
      real pi

      pi=4.0*atan(1.0)

      call grad(x,v,n,dy,ibs,id)

      do 200 i=1,n
      v(i)=-180.0*atan(v(i))/pi
 200  continue

      return
      end

****************************************************************************
      subroutine grad(x,der,n,dy,ibs,id)
      integer n,ibs,id
      real x(n),der(n),dy
*-------------------------------------------------------------------------
      integer i,i1,i2,nint
      real x1,dev1,dev2,dev3,a,b,c,d,y0,dydum


      x1=x(1)

      do 100 i=1,n
      x(i)=x(i)-x1
 100  continue


      do 200 i=1,n
      i1=i-ibs
      i2=i+ibs
      i1=max(1,i1)
      i2=min(n,i2)
      nint=i2-i1+1
      call noreg(x(i1),nint,a,b,c,d,dev1,dev2,dev3)
      y0=i1-i
      dydum=1.0
      call tra(nint,dydum,y0,a,b,c,d,id)
      der(i)=b/dy
 200  continue


      do 300 i=1,n
      x(i)=x(i)+x1
 300  continue

      return
      end

****************************************************************************
      subroutine gradg(x,der,n,dy,ibs,id,ird)
      integer n,ibs,id,ird
      real x(n),der(n),dy
*-------------------------------------------------------------------------
      integer i,i1,i2,nint
      real x1,dev1,dev2,dev3,a,b,c,d,y0,f0,f1,f2,ff,rdy,rdy2,dydum

      rdy=1.0/dy
      rdy2=rdy*rdy

      if(ird.eq.0) then
          f0=1.0
          f1=0.0
          f2=0.0
          ff=1.0
      else
        if(ird.eq.1) then
          f0=0.0
          f1=1.0
          f2=0.0
          ff=0.0
        else
          f0=0.0
          f1=0.0
          f2=1.0
          ff=0.0
         end if
      end if

      x1=x(1)

      do 100 i=1,n
      x(i)=x(i)-x1
 100  continue


      do 200 i=1,n
      i1=i-ibs
      i2=i+ibs
      i1=max(1,i1)
      i2=min(n,i2)
      nint=i2-i1+1
      call noreg(x(i1),nint,a,b,c,d,dev1,dev2,dev3)
      y0=i1-i
      dydum=1.0
      call tra(nint,dydum,y0,a,b,c,d,id)
      der(i)=ff*x1+f0*a+f1*b*rdy +f2*2.0*c*rdy2
 200  continue


      do 300 i=1,n
      x(i)=x(i)+x1
 300  continue

      return
      end



****************************************************************************
      subroutine gradpoin(x,der,n,dy,ipkt,ibs,id)
      integer n,ibs,id,ipkt
      real x(n),der(0:2),dy
*-------------------------------------------------------------------------
      integer i,i1,i2,nint
      real x1,dev1,dev2,dev3,a,b,c,d,y0,rdy,rdy2,dydum

      rdy=1.0/dy
      rdy2=rdy*rdy




cc      i=ipkt
      i1=ipkt-ibs
      i2=ipkt+ibs
      i1=max(1,i1)
      i2=min(n,i2)

      x1=x(i1)

      do 100 i=i1,i2
      x(i)=x(i)-x1
 100  continue


      nint=i2-i1+1
      call noreg(x(i1),nint,a,b,c,d,dev1,dev2,dev3)
      y0=i1-ipkt
      dydum=1.0
      call tra(nint,dydum,y0,a,b,c,d,id)
      der(0)=x1+a
      der(1)=b*rdy
      der(2)=2.0*c*rdy2


      do 300 i=i1,i2
      x(i)=x(i)+x1
 300  continue

      return
      end

****************************************************************************
      subroutine ggdstretch(x,der,n,y0,dy,yp,span,ibmin,id,a,b,c,d)
      integer n,ibmin,id
      real x(n),der(0:2),y0,dy,yp,span,a,b,c,d
*-------------------------------------------------------------------------
      integer ibs


      ibs=span/dy+1
      if(ibs.lt.ibmin) ibs=ibmin
      call ggdpoin(x,der,n,y0,dy,yp,ibs,id,a,b,c,d)

      return
      end 


****************************************************************************
      subroutine ggdpoin(x,der,n,y0,dy,yp,ibs,id,a,b,c,d)
      integer n,ibs,id
      real x(n),der(0:2),y0,dy,yp,a,b,c,d
*-------------------------------------------------------------------------
      integer i,i1,i2,nint,ipkt
      real x1,dev1,dev2,dev3,rdy,rdy2,dydum,ynp,ynorm

      rdy=1.0/dy
      rdy2=rdy*rdy

      ynp=(yp-y0)/dy+1
      ipkt=ynp

cc      i=ipkt
      i1=ipkt-ibs
      i2=ipkt+ibs
      i1=max(1,i1)
      i2=min(n,i2)

      x1=x(i1)

      do 100 i=i1,i2
      x(i)=x(i)-x1
 100  continue


      nint=i2-i1+1
      call noreg(x(i1),nint,a,b,c,d,dev1,dev2,dev3)
cc      ynorm=i1-ynp
      ynorm=y0
cc      dydum=1.0
      dydum=dy
      call tra(nint,dydum,ynorm,a,b,c,d,id)
cc      der(0)=x1+a
cc      der(1)=b*rdy
cc      der(2)=2.0*c*rdy2
      a=a+x1
      der(0)=a +yp*(b+yp*(c+yp*d))
      der(1)=b+yp*(2.0*c+yp*3.0*d)
      der(2)=2.0*c+ 6.0*yp*d


      do 300 i=i1,i2
      x(i)=x(i)+x1
 300  continue

      return
      end


************************************************************************
*
*               RUTINER FOR SKRIVING AV MATRISER P] FIL.
*
*     Eksempler p} bruk:
*     (1):
*        I et program er datafiler deklarert ved : 
*
*      REAL A(2:100,50),B(2:100,50),C(2:100,50)
*
*        Etter at beregningene er utf|rt |nsker vi } skrive ut
*        delene (10:60,1:30) uformattert for plotting. 
*        Dette kan gj|res ved:
*                         
*             .
*             .  (vi antar dx etc. har f}tt verdier her)
*             .
*      OPEN(UNIT=ITAPE ... ,form='unformatted' .... )
*      CALL SKHODEG(0,51,30,DX,DY,X0,Y0,3,ITAPE,.false.)
*      CALL SKRIVG(A,2,1,100,50,10,1,60,30,ITAPE,.false.)
*      CALL SKRIVG(B,2,1,100,50,10,1,60,30,ITAPE,.false.)
*      CALL SKRIVG(C,2,1,100,50,10,1,60,30,ITAPE,.false.)
*      CLOSE(ITAPE)
*
*      (2):
*         Samme som i (i) men vi |nsker } plotte alt som om det l} i en
*         matrise med 51 ganger 90 punkter og denne gangen formattert:
*
*             .
*             .  (som f|r)
*      OPEN(UNIT=ITAPE ... ,form='formatted' .... )
*      CALL SKHODEG(0,51,90,DX,DY,X0,Y0,1,ITAPE,.true.)
*             .
*             .  (som f|r bortsett fra at .false. byttes med .true.)
*             .
*
***********************************************************************


***********************************************************************
*
*               S K I P D A T
*
*     Rutinenen posisjonerer innp.fila ett datasett framover.
*     parametere:
*           N,M - antall verdier i ett datasett er N ganger M        (I)
*           ITAPE - tapenummer                                       (I)
***********************************************************************   
      SUBROUTINE SKIPDAT(N,M,ITAPE)
      INTEGER N,M,ITAPE
*---------------------------------------------------------------------
      INTEGER I,IREC

      IREC=N/10
      IF(N.GT.IREC*10) IREC=IREC+1
      IREC=IREC*M

      DO 100 I=1,IREC
      READ(ITAPE,50)
   50 FORMAT(1X)
  100 CONTINUE
      RETURN
      END

***********************************************************************
*
*               S K I P D A T G
*
*     Rutinenen posisjonerer innp.fila ett datasett framover.
*     parametere:
*           N,M - antall verdier i ett datasett er N ganger M        (I)
*           ITAPE - tapenummer                                       (I)
*           asc - angir formatering                                  (I)
***********************************************************************   
      SUBROUTINE SKIPDATG(N,M,ITAPE,asc)
      INTEGER N,M,ITAPE
      logical asc
*---------------------------------------------------------------------
      INTEGER I,IREC
      character c

      IREC=N/10
      IF(N.GT.IREC*10) IREC=IREC+1
      IREC=IREC*M

      if(asc) then
        DO 100 I=1,IREC
        READ(ITAPE,50)
   50   FORMAT(1X)
  100   CONTINUE
      else
        DO 200 I=1,IREC
        READ(ITAPE)c
  200   CONTINUE
      end if
      RETURN
      END




*********************************************************************
*
*            S K I P H O D E
*
*     Rutinen posisjonerer forbi hodedeklarasjonen p} innp.fila.
*     Parametere:
*           ITAPE  - unitnummer for innputfil
*********************************************************************
      SUBROUTINE SKIPHODE(ITAPE)
      INTEGER ITAPE
*--------------------------------------------------------------------
      READ(ITAPE,100)
  100 FORMAT(1X)
      RETURN
      END 

*********************************************************************
*
*            S K I P H O D E G
*
*     Rutinen posisjonerer forbi hodedeklarasjonen p} innp.fila.
*     Parametere:
*           ITAPE  - unitnummer for innputfil
*********************************************************************
      SUBROUTINE SKIPHODEG(ITAPE,ASC)
      INTEGER ITAPE
      logical asc
*--------------------------------------------------------------------
      integer NC,N,M,NPLOT
      real DX,DY,X0,Y0

      call headles(NC,N,M,DX,DY,X0,Y0,NPLOT,ITAPE,asc)

      RETURN
      END 



**********************************************************************
*
*                S K R I V 
*
*      skriver ut ett datasett p} outputfila.
*      Parametere:
*             Y - array som inneholder dataene.                     (I)
*             NS0,MS0,NS,MS - fysiske arraygrenser for den          (I) 
*                             todimensjonale nummerering av Y.
*             N0,M0 - nedre venstre element som skrives ut.
*             N,M   - angir |vre h|yre element som skal skrives ut.   (I) 
*                     Dvs. det utskrevne datasett svarere til :
*                     
*                     Y(N0,M0)     ....   Y(N,M0)  
*                        :                   :
*                        :                   :
*                     Y(N0,M)      ....   Y(N,M) 
*
*                     merk at antall punkter er (n-n0+1)X(m-m0+1)
*        
*             ITAPE - unitnummer for fil                            (I)
****************************************************************************
      SUBROUTINE SKRIV(Y,NS0,MS0,NS,MS,N0,M0,N,M,ITAPE)
      INTEGER N,M,NS,MS,ITAPE,N0,M0,NS0,MS0
      INTEGER K,J,NA,NB,NDEL,IRAD,NREST,L
      REAL Y(NS0:NS,MS0:MS)

      IRAD=N-N0+1
      NDEL=IRAD/10
      NREST=IRAD-NDEL*10
   11 format(10F12.6)

      DO 100 K=M0,M
      DO 50 L=1,NDEL
      NA=N0+(L-1)*10
      NB=NA+9
      WRITE(ITAPE,11) (Y(J,K) , J=NA,NB)
   50 CONTINUE
      IF(NREST.GT.0) THEN
        NA=N-NREST+1
        WRITE(ITAPE,11) (Y(J,K) , J=NA,N)
      END IF
  100 CONTINUE
      RETURN
      END



**********************************************************************
*
*                S K R I V G
*
*      skriver ut ett datasett p} outputfila.
*      Parametere:
*             Y - array som inneholder dataene.                     (I)
*             NS0,MS0,NS,MS - fysiske arraygrenser for den          (I) 
*                             todimensjonale nummerering av Y.
*             N0,M0 - nedre venstre element som skrives ut.
*             N,M   - angir |vre h|yre element som skal skrives ut.   (I) 
*                     Dvs. det utskrevne datasett svarere til :
*                     
*                     Y(N0,M0)     ....   Y(N,M0)  
*                        :                   :
*                        :                   :
*                     Y(N0,M)      ....   Y(N,M) 
*
*                     merk at antall punkter er (n-n0+1)X(m-m0+1)
*        
*             ITAPE - unitnummer for fil                            (I)
*             asc - verdi = .true. angir ascii format, =.false. binaerf.  I 
****************************************************************************
      SUBROUTINE SKRIVG(Y,NS0,MS0,NS,MS,N0,M0,N,M,ITAPE,ASC)
      INTEGER N,M,NS,MS,ITAPE,N0,M0,NS0,MS0
      REAL Y(NS0:NS,MS0:MS)
      logical asc
*--------------------------------------------------------------
      INTEGER K,J,NA,NB,NDEL,IRAD,NREST,L

      IRAD=N-N0+1
      NDEL=IRAD/10
      NREST=IRAD-NDEL*10
   11 format(10F12.6)

      if(asc) then
        DO 200 K=M0,M
        DO 150 L=1,NDEL
        NA=N0+(L-1)*10
        NB=NA+9
        WRITE(ITAPE,11) (Y(J,K) , J=NA,NB)
  150   CONTINUE
        IF(NREST.GT.0) THEN
          NA=N-NREST+1
          WRITE(ITAPE,11) (Y(J,K) , J=NA,N)
        END IF
  200   CONTINUE
      else
        DO 100 K=M0,M
        DO 50 L=1,NDEL
        NA=N0+(L-1)*10
        NB=NA+9
        WRITE(ITAPE) (Y(J,K) , J=NA,NB)
   50   CONTINUE
        IF(NREST.GT.0) THEN
          NA=N-NREST+1
          WRITE(ITAPE) (Y(J,K) , J=NA,N)
        END IF
  100   CONTINUE
      end if
      RETURN
      END




**********************************************************************      
*
*               L E S V
*
*      leser et datasett skrevet ut av 'skriv'. Parameterene har 
*      tilsvarende betydninger. Differensene (N -N0)og (M-M0) m} 
*      ha samme verdi som i kall p} 'skriv', men dette gjelder ikke
*      ns0,ms0,ns,ms
***********************************************************************
      SUBROUTINE LESV(Y,NS0,MS0,NS,MS,N0,M0,N,M,ITAPE)
      INTEGER N,M,NS,MS,ITAPE,N0,M0,NS0,MS0
      INTEGER K,J,NA,NB,NDEL,IRAD,NREST,L
      REAL Y(NS0:NS,MS0:MS)

      IRAD=N-N0+1
      NDEL=IRAD/10
      NREST=IRAD-NDEL*10
   11 FORMAT(10F12.6)


      DO 100 K=M0,M
      DO 50 L=1,NDEL
      NA=N0+(L-1)*10
      NB=NA+9
      READ(ITAPE,11) (Y(J,K) , J=NA,NB)
   50 CONTINUE
      IF(NREST.GT.0) THEN
        NA=N-NREST+1
        READ(ITAPE,11) (Y(J,K) , J=NA,N)
      END IF
  100 CONTINUE
      RETURN
      END

**********************************************************************      
*
*               L E S G
*
*      leser et datasett skrevet ut av 'skriv'/'skrivg'. Parameterene har 
*      tilsvarende betydninger. Differensene (N -N0)og (M-M0) m} 
*      ha samme verdi som i kall p} 'skriv', men dette gjelder ikke
*      ns0,ms0,ns,ms
***********************************************************************
      SUBROUTINE LESG(Y,NS0,MS0,NS,MS,N0,M0,N,M,ITAPE,ASC)
      INTEGER N,M,NS,MS,ITAPE,N0,M0,NS0,MS0
      REAL Y(NS0:NS,MS0:MS)
      logical asc
*------------------------------------------------------------------
      INTEGER K,J,NA,NB,NDEL,IRAD,NREST,L

      IRAD=N-N0+1
      NDEL=IRAD/10
      NREST=IRAD-NDEL*10
   11 FORMAT(10F12.6)

      if(asc) then
        DO 200 K=M0,M
        DO 150 L=1,NDEL
        NA=N0+(L-1)*10
        NB=NA+9
        READ(ITAPE,11) (Y(J,K) , J=NA,NB)
  150   CONTINUE
        IF(NREST.GT.0) THEN
          NA=N-NREST+1
          READ(ITAPE,11) (Y(J,K) , J=NA,N)
        END IF
  200   CONTINUE
      else
        DO 100 K=M0,M
        DO 50 L=1,NDEL
        NA=N0+(L-1)*10
        NB=NA+9
        READ(ITAPE) (Y(J,K) , J=NA,NB)
   50   CONTINUE
        IF(NREST.GT.0) THEN
          NA=N-NREST+1
          READ(ITAPE) (Y(J,K) , J=NA,N)
        END IF
  100   CONTINUE
      end if
      RETURN
      END

**********************************************************************      
*
*               SJEKKLES 
*
*      leser et datasett skrevet ut av 'skriv'/'skrivg' mens det sjekkes på
*      tilgjengelighet av data. 
*      parametere:
*        y - dataarray                                             O
*        nmax - maks antall data i y                                I
*        N,m - antall kolonner og maksimalt antall rader             I
*        nlest - antall tall i siste rad. Dersom nlest=n er det     O
*                lest et helt antall rader
*        mlest - antall rader lest, helt eller påbegynt.            O
*        itape - unitnummer                                         I
*        asc - type datafil                                          I
***********************************************************************
      SUBROUTINE sjekkles(Y,nmax,N,M,nlest,mlest,ITAPE,ASC)
      INTEGER nmax,N,M,nlest,mlest,itape
      REAL Y(nmax)
      logical asc
*------------------------------------------------------------------
      INTEGER K,J,Ntot,NDEL,NREST,L
      real buff(10)
    
      ndel=n/10
      nrest=n-ndel*10

   11 FORMAT(10F12.6)

      ntot=0
      k=1

 200  continue

      l=1

      if(ndel.gt.0) then
 150  continue
        if(asc) then
          READ(ITAPE,11,end=300,err=400) (buff(j) , J=1,10)
        else
          READ(ITAPE,end=300) (buff(j) , J=1,10)
        end if

        if(ntot+10.gt.nmax) go to 300        

        do 120 j=1,10
          y(ntot+j)=buff(j)
 120    continue
        ntot=ntot+10

        l=l+1
        if(l.le.ndel)go to 150
      end if 


      IF(NREST.GT.0) THEN
        if(asc) then
          READ(ITAPE,11,end=300,err=400) (buff(j) , J=1,nrest)
        else
          READ(ITAPE,end=300) (buff(j) , J=1,nrest)
        end if
        if(ntot+nrest.gt.nmax) go to 300        

        do 140 j=1,nrest
          y(ntot+j)=buff(j)
 140    continue
        ntot=ntot+nrest


      END IF

      if(k.lt.m) then
       k=k+1
       go to 200
      end if

c     Here everything is read OK
      nlest=n
      mlest=m
      return

 400  continue
c     Reading has failed because of error. Treated as premature end of file
      write(0,*)'sjekkles: wrong data at', k,l 
 300  continue
c     Here some attempted reading has failed 
      
      if(l.eq.1) then
c      whole row k has failed
       mlest=k-1
       if(mlest.gt.0) then
        nlest=n
       else
        nlest=0
       end if
      else
       mlest=k
       nlest=(l-1)*10
      end if
      

      RETURN
      END

**************************************************************************
*
*                S K H O D E G
*     skriver deklarasjonshode for datasett.
*     Parametere: (alle er innput)
*            NC - kontrollparameter, verdi lik 0 er det vanlige.
*                 Verdi lik 1 viser at f|rste datasett er } oppfatte
*                 som dybdematrise og plotteprogram gphov vil maske
*                 ut automatisk p} grunnlag av dette.
*            N,M - antall punkter i datasettene som forklart i kommentar
*                  til skriv.
*            DX,DY - avstander mellom punkter.
*            X0,Y0 - posisjon av nedre venstre punkt.
*            NPLOT - antall datasett som skal skrives ut.
*            ITAPE - unitnummer for outputfil.
*            asc  - .true. angir ascii, .false. bin{r
**************************************************************************
      SUBROUTINE SKHODEG(NC,N,M,DX,DY,X0,Y0,NPLOT,ITAPE,asc)
      INTEGER NC,N,M,NPLOT,ITAPE
      REAL DX,DY,X0,Y0
      logical asc

      if(asc) then
        WRITE(ITAPE,10) NC,N,M,DX,DY,X0,Y0,NPLOT
      else
        WRITE(ITAPE) NC,N,M,DX,DY,X0,Y0,NPLOT
      end if
   10 FORMAT(3I6,4F12.6,I6)

      RETURN
      END
**************************************************************************
*
*                S K H O D E
*     skriver deklarasjonshode for datasett.
*     Parametere: (alle er innput)
*            NC - kontrollparameter, verdi lik 0 er det vanlige.
*                 Verdi lik 1 viser at f|rste datasett er } oppfatte
*                 som dybdematrise og plotteprogram gphov vil maske
*                 ut automatisk p} grunnlag av dette.
*            N,M - antall punkter i datasettene som forklart i kommentar
*                  til skriv.
*            DX,DY - avstander mellom punkter.
*            X0,Y0 - posisjon av nedre venstre punkt.
*            NPLOT - antall datasett som skal skrives ut.
*            ITAPE - unitnummer for outputfil.
**************************************************************************
      SUBROUTINE SKHODEB(NC,N,M,DX,DY,X0,Y0,NPLOT,ITAPE)
      INTEGER NC,N,M,NPLOT,ITAPE
      REAL DX,DY,X0,Y0

        WRITE(ITAPE,10) NC,N,M,DX,DY,X0,Y0,NPLOT

   10 FORMAT(3I6,4F12.6,I6)

      RETURN
      END

**************************************************************************
*
*            L E H O D E
*
*     leser deklarasjonene skrevet ut med 'skhode'
*************************************************************************
      SUBROUTINE LEHODE(NC,N,M,DX,DY,X0,Y0,NPLOT,ITAPE)
      INTEGER NC,N,M,NPLOT,ITAPE
      REAL DX,DY,X0,Y0

      READ(ITAPE,10) NC,N,M,DX,DY,X0,Y0,NPLOT
   10 FORMAT(3I6,4F12.6,I6)

       RETURN
       END






*****************************************************************
*
*    Rutinen leser inn parametere i hodet som lehode, men hopper
*    over  innledende kommentarer. Som innledende kommentar regnes
*    linjer som er blanke, eller linjer der f|rste ikke-blanke tegn 
*    er '!' . De 10 f|rste kommentarlinjene tas vare p} i common-
*    omraadet hode og kan listes ut ve kall p} listexpl
*
*************************************************************************
      subroutine headles(NC,N,M,DX,DY,X0,Y0,NPLOT,ITAPE,asc)
      INTEGER NC,N,M,NPLOT,ITAPE
      REAL DX,DY,X0,Y0
      logical asc
      integer nlin
      character*80 clin(10)
      common/hode/nlin,clin
*----------------------------------------------
      integer iflag,i0,i1,ib,ie
      character*80 cimage

      nlin=0

      if(.not.asc) then 
        READ(ITAPE) NC,N,M,DX,DY,X0,Y0,NPLOT
cc        write(0,*) NC,N,M,DX,DY,X0,Y0,NPLOT
        return
      end if

 100  continue

      read(itape,10) cimage(1:80)

 10   format(a80)

      ib=1

      call lokord(cimage,ib,80,i0,i1,iflag)

      if(iflag.eq.0) go to 100
      if(cimage(i0:i0).ne.'!')  go to 200

      if(nlin.lt.10) then
        nlin=nlin+1
        ie=80-i0
        clin(nlin)(1:ie)=cimage(I0+1:80)
      end if

      go to 100

 200  continue

      call GETI(Cimage,ib,80,I0,I1,nc,IFLAG)
      if(iflag.ne.1) write(0,*)'t|ys i filhode nc', Cimage
      ib=i1+1

      call GETI(Cimage,ib,80,I0,I1,n,IFLAG)
      if(iflag.ne.1) write(0,*)'t|ys i filhode n'
      ib=i1+1

      call GETI(Cimage,ib,80,I0,I1,m,IFLAG)
      if(iflag.ne.1) write(0,*)'t|ys i filhode m'
      ib=i1+1

      call GETR(Cimage,ib,80,I0,I1,dx,IFLAG)
      if(iflag.ne.1) write(0,*)'t|ys i filhode dx'
      ib=i1+1

      call GETR(Cimage,ib,80,I0,I1,dy,IFLAG)
      if(iflag.ne.1) write(0,*)'t|ys i filhode dy'
      ib=i1+1

      call GETR(Cimage,ib,80,I0,I1,x0,IFLAG)
      if(iflag.ne.1) write(0,*)'t|ys i filhode x0'
      ib=i1+1

      call GETR(Cimage,ib,80,I0,I1,y0,IFLAG)
      if(iflag.ne.1) write(0,*)'t|ys i filhode y0'
      ib=i1+1

      call GETI(Cimage,ib,80,I0,I1,nplot,IFLAG)
      if(iflag.ne.1) write(0,*)'t|ys i filhode nplot'
      ib=i1+1

      return
      end





*********************************************************************
      subroutine listexpl(itape,asc)
      integer itape
      logical asc

      integer nlin
      character*80 clin(10)
      common/hode/nlin,clin
*----------------------------------------------
      integer i

      if(nlin.eq.0) return

      if(asc.and.(itape.eq.0.or.itape.eq.5)) then
        write(itape,*) ' '
        write(itape,*)'FILBESKRIVELSE'
        write(itape,*)'-------------------------------------'
      end if

      do 100 i=1,nlin
      if(asc) then
        write(itape,10) clin(i)(1:79)
      else
        write(itape) clin(i)(1:80)
      end if
 10   format(a)
 100  continue

      if(asc.and.(itape.eq.0.or.itape.eq.5)) then
        write(itape,*) ' '
        write(itape,*) ' '
      end if
      
      return
      end


*****************************************************************
*                        HEADSKIP
*     Skipper forbi hode.
*       itape - unitnummer                                 I
*       asc - verdi=.true./.false. angir ascii/binaert     I
*****************************************************************
      subroutine headskip(ITAPE,asc)
      integer itape
      logical asc
*----------------------------------------------
      integer iflag,i0,i1,ib
      character*80 cimage



 100  continue

      if(asc) then
        read(itape,10) cimage(1:80)
      else
        read(itape) cimage(1:80)
      end if

 10   format(a80)

      ib=1

      call lokord(cimage,ib,80,i0,i1,iflag)

      if(iflag.eq.0) go to 100
      if(cimage(i0:i0).ne.'!')  go to 200



      go to 100

 200  continue


      return
      end





*****************************************************************
*
*                KRYMP
*
*      Leser en del av et datasett inn i u
*      parametere:
*            U - array                                          O
*            n,m - antall punkter i lest grid                   I
*            nt,mt - antall punkter i datasett paa fil          I
*            nb,mb - startpunkt for lesning, dvs. U(1,1)        I
*                    svarer til punkt (nb,mb) i datasett
*            ix,iy - tetthet i hhv. x og y retning              I
*                    ix=3 gir feks at bare hver tredje kolonne
*                    leses.
*            itape - unitnummer                                 I
*            asc - verdi=.true./.false. angir ascii/binaert     I
*            w - arbeidsarray av lengde nt                      D
*****************************************************************            
      SUBROUTINE KRYMP(U,N,M,NT,MT,NB,MB,IX,IY,ITAPE,asc,W)
      INTEGER N,M,NT,MT,NB,MB,IX,IY,ITAPE
      REAL U(N,M),W(NT)
      logical asc
*-----------------------------------------------------------------
      INTEGER L,J,K,MH

      IF(N.EQ.NT .AND. M.EQ.MT) THEN
        CALL LESG(U,1,1,N,M,1,1,N,M,ITAPE,asc)
        RETURN
      END IF

      DO 50 J=1,MB
   50 CALL LESG(W,1,1,NT,1,1,1,NT,1,ITAPE,asc)

      DO 55 L=1,N
   55 U(L,1)=W(NB+(L-1)*IX)

      DO 100 J=2,M
      DO 60 K=1,IY
   60 CALL LESG(W,1,1,NT,1,1,1,NT,1,ITAPE,asc)
      DO 70 L=1,N
   70 U(L,J)=W(NB+(L-1)*IX)
  100 CONTINUE
                      
      MH=MT-(MB+(M-1)*IY)
      DO 200 J=1,MH
  200 CALL LESG(W,1,1,NT,1,1,1,NT,1,ITAPE,asc)
      RETURN
      END




*****************************************************************
* 
*                  REDSKG
*
*       Skriver del av datasett.
*      Parametere:
*             U - array som inneholder dataene.                     (I)
*             NS0,MS0,NS,MS - fysiske arraygrenser for den          (I) 
*                             todimensjonale nummerering av U.
*             N0,M0 - nedre venstre element som skrives ut.
*             N,M   - angir |vre h|yre grense for utskrift. NB: selve   (I)
*                     elementet skrives ut barre dersom det gaar opp med
*                     ix,iy - se nedenfor.
*             ix,iy - tetthet av utskrift. Antall punkter blir da nb,mb   I
*                     med    mb=(m-m0)/iy+1    nb=(n-n0)/ix+1
*                     Dvs. det utskrevne datasett svarere til :
*                     
*                     Y(N0,M0)     Y(n0+ix,m0)    ....      .
*                     Y(n0,m0+iy)                 
*                        :                                  :
*                        :                                  :
*                        :           .    ....   Y(N0+(nb-1)*ix,M0+(mb-1)*ix) 
*
*                     merk at antall punkter er (n-n0+1)X(m-m0+1)
*        
*             ITAPE - unitnummer for fil                            (I)
*             asc - verdi = .true. angir ascii format, =.false. binaerf.  I 
*             w - arbeidsarray min. lengde: nb                            D
****************************************************************************

      SUBROUTINE redskg(U,NS0,MS0,NS,MS,N0,M0,N,m,IX,IY,ITAPE,asc,W)
      INTEGER N,M,NS0,MS0,NS,MS,N0,M0,IX,IY,ITAPE
      REAL U(ns0:ns,ms0:ms),w(n)
      logical asc
*-----------------------------------------------------------------
      INTEGER L,J,K,Mb,nb

      mb=(m-m0)/iy+1
      nb=(n-n0)/ix+1

      IF(ix.EQ.1 .AND. iy.EQ.1) THEN
        CALL skrivg(U,ns0,ms0,Ns,Ms,n0,m0,N,M,ITAPE,asc)
        RETURN
      END IF

      DO 60 J=1,MB
      k=m0+(j-1)*iy
      DO 55 L=1,Nb
   55 w(L)=u(N0+(L-1)*IX,k)

      CALL skrivg(W,1,1,Nb,1,1,1,Nb,1,ITAPE,asc)
 60   continue

      RETURN
      END

*****************************************************************
* 
*                  MREDSKG
*
*       Skriver del av datasett.
*      Parametere:
*             U - array som inneholder dataene.                     (I)
*             NS0,MS0,NS,MS - fysiske arraygrenser for den          (I) 
*                             todimensjonale nummerering av U.
*             N0,M0 - nedre venstre element som skrives ut.
*             NB,MB   - number of values to be written              (I)
*                     OBS OBS: difference from redskg
*             ix,iy - tetthet av utskrift. Antall punkter blir da nb,mb   I
*                     med    mb=(m-m0)/iy+1    nb=(n-n0)/ix+1
*                     Dvs. det utskrevne datasett svarere til :
*                     
*                     Y(N0,M0)     Y(n0+ix,m0)    ....      .
*                     Y(n0,m0+iy)                 
*                        :                                  :
*                        :                                  :
*                        :           .    ....   Y(N0+(nb-1)*ix,M0+(mb-1)*iy) 
*
*                    
*        
*             ITAPE - unitnummer for fil                            (I)
*             sf - scaling factor; sf*data is written                 (I)
*             asc - verdi = .true. angir ascii format, =.false. binaerf.  I 
*             w - arbeidsarray min. lengde: nb                            D
****************************************************************************

      SUBROUTINE mredskg(U,NS0,MS0,NS,MS,N0,M0,Nb,Mb,IX,IY,ITAPE,sf,
     %    asc,W)
      INTEGER Nb,Mb,NS0,MS0,NS,MS,N0,M0,IX,IY,ITAPE
      REAL U(ns0:ns,ms0:ms),w(ns),sf
      logical asc
*-----------------------------------------------------------------
      INTEGER L,J,K


      DO 60 J=1,MB
      k=m0+(j-1)*iy
      DO 55 L=1,Nb
   55 w(L)=sf*u(N0+(L-1)*IX,k)

      CALL skrivg(W,1,1,Nb,1,1,1,Nb,1,ITAPE,asc)
 60   continue

      RETURN
      END



****************************************************************************
      SUBROUTINE SKRIVI(Y,NS0,MS0,NS,MS,N0,M0,N,M,ITAPE,ASC)
      INTEGER N,M,NS,MS,ITAPE,N0,M0,NS0,MS0
      integer Y(NS0:NS,MS0:MS)
      logical asc
*--------------------------------------------------------------
      INTEGER K,J

      if(asc) then

        do 200 k=m0,m
        write(itape,'(10i10)')(y(j,k), j=n0,n)
 200    continue
  
      else

        do 100 k=m0,m
        write(itape)(y(j,k), j=n0,n)
 100    continue

      end if

      RETURN
      END


****************************************************************************
      SUBROUTINE LESI(Y,NS0,MS0,NS,MS,N0,M0,N,M,ITAPE,ASC)
      INTEGER N,M,NS,MS,ITAPE,N0,M0,NS0,MS0
      integer Y(NS0:NS,MS0:MS)
      logical asc
*--------------------------------------------------------------
      INTEGER K,J

      if(asc) then

        do 200 k=m0,m
        read(itape,'(10i10)')(y(j,k), j=n0,n)
 200    continue
  
      else

        do 100 k=m0,m
        read(itape)(y(j,k), j=n0,n)
 100    continue

      end if

      RETURN
      END





*****************************************************************
*
*     GPHOVPRI
*
* printer array i gphov-format p} fil
* parametere:
*         a  - array                          I
*         ns - fysisk array grense                   I
*         nc - kontroll parameter
*         n,m  - array-grenser                  I
*         dx,dy  - gitterinkrementer           I
*         x0,y0  - posisjon av (1,1)           I
*         isyk - nummer paa fil              I
*         stem - navnestamme, feks stem='eta#', isyk=11 gir datafil eta11 I
*         asc - asci eller binaer
***************************************************************************
      subroutine gphovpri(a,ns,n,m,dx,dy,x0,y0,nc,stem,isyk,asc)
      integer ns,n,m,nc,isyk
      real a(ns,m),dx,dy,x0,y0
      logical asc
      character*80 stem
*--------------------------------------------------------------
      integer itape
    
      itape=37

      call noppg(stem,isyk,itape,asc)
      CALL SKHODEG(nc,n,m,dx,dy,x0,y0,1,itape,asc)
      CALL SKRIVG(a,1,1,ns,m,1,1,n,m,itape,asc)
      close(itape)

      return
      end



*****************************************************************
*
*     GPPRSCAL
*
* printer array i gphov-format p} fil
* parametere:
*         a  - array                          I
*         ns - fysisk array grense                   I
*         nc - kontroll parameter
*         n,m  - array-grenser                  I
*         dx,dy  - gitterinkrementer           I
*         x0,y0  - posisjon av (1,1)           I
*         scal - scaling factor                I  
*         isyk - nummer paa fil              I
*         stem - navnestamme, feks stem='eta#', isyk=11 gir datafil eta11 I
*         asc - asci eller binaer                    I
*         wrk - work array                           D
***************************************************************************
      subroutine gpprscal(a,ns,n,m,dx,dy,x0,y0,scal,nc,stem,isyk,
     %         asc,wrk)
      integer ns,n,m,nc,isyk
      real a(ns,m),wrk(n),dx,dy,x0,y0,scal
      logical asc
      character*80 stem
*--------------------------------------------------------------
      integer itape,k,i
    
      itape=37

      call noppg(stem,isyk,itape,asc)
      CALL SKHODEG(nc,n,m,dx,dy,x0,y0,1,itape,asc)
      do 100 k=1,m
        do 50 i=1,n
          wrk(i)=scal*a(i,k)
 50     continue
        CALL SKRIVG(wrk,1,1,n,1,1,1,n,1,itape,asc)
 100  continue  
      close(itape)

      return
      end

*****************************************************************
*
*     spanprint
*
* printer 1D array i gphov-format p} fil
* parametere:
*         a  - array,  only 1D                        I
*         nc - kontroll parameter
*         n,m  - array-grenser                  I
*         dx,dy  - gitterinkrementer           I
*         x0,y0  - posisjon av (1,1)           I
*         isyk - nummer paa fil              I
*         stem - navnestamme, feks stem='eta#', isyk=11 gir datafil eta11 I
*         asc - asci eller binaer
***************************************************************************
      subroutine spanprint(a,n,m,dx,dy,x0,y0,scal,nc,stem,isyk,
     %         asc,wrk)
      integer n,m,nc,isyk
      real a(n),wrk(n),dx,dy,x0,y0,scal
      logical asc
      character*80 stem
*--------------------------------------------------------------
      integer itape,k,i
    
      itape=37

      do 50 i=1,n
          wrk(i)=scal*a(i)
 50   continue

      call noppg(stem,isyk,itape,asc)
      CALL SKHODEG(nc,n,m,dx,dy,x0,y0,1,itape,asc)
      do 100 k=1,m
         CALL SKRIVG(wrk,1,1,n,1,1,1,n,1,itape,asc)
 100  continue  
      close(itape)

      return
      end


*****************************************************************
*
*     spanprint
*
* printer 1D array i gphov-format p} fil
* parametere:
*         a  - array,  only 1D                        I
*         nc - kontroll parameter
*         n,m  - array-grenser                  I
*         dx,dy  - gitterinkrementer           I
*         x0,y0  - posisjon av (1,1)           I
*         isyk - nummer paa fil              I
*         stem - navnestamme, feks stem='eta#', isyk=11 gir datafil eta11 I
*         asc - asci eller binaer
***************************************************************************
      subroutine kanpri(a,n,dx,x0,scal,nc,stem,isyk,asc,wrk)
      integer n,nc,isyk
      real a(n),wrk(n),dx,x0,scal
      logical asc
      character*80 stem
*--------------------------------------------------------------
      integer itape,k,i,m
      real y0,dy
    
      dy=dx
      y0=-dy
      m=3
      itape=37


      call noppg(stem,isyk,itape,asc)
      CALL SKHODEG(nc,n,m,dx,dy,x0,y0,1,itape,asc)
 
      do 40 i=1,n
          wrk(i)=-1.0
 40   continue
      CALL SKRIVG(wrk,1,1,n,1,1,1,n,1,itape,asc)

      do 50 i=1,n
          wrk(i)=scal*a(i)
 50   continue
      CALL SKRIVG(wrk,1,1,n,1,1,1,n,1,itape,asc)
 
      do 60 i=1,n
          wrk(i)=-1.0
 60   continue
      CALL SKRIVG(wrk,1,1,n,1,1,1,n,1,itape,asc)

      close(itape)

      return
      end




*****************************************************************
*
*     GPPRSCAL
*
* inflates a series of rows to fields and write the result in gphov format
* parametere:
*         a  - array                          I
*         ns - fysisk array grense                   I
*         nc - kontroll parameter
*         n,m  - array-grenser                  I
*         nz - width of fields                 I
*         dx,dz  - grid-increments           I
*         x0,z0  - position av (1,1) in each field          I
*         isyk - nummer paa fil              I
*         stem - navnestamme, feks stem='eta#', isyk=11 gir datafil eta11 I
*         asc - asci eller binaer
***************************************************************************
      subroutine gpprinflate(a,ns,n,m,nz,dx,dz,x0,z0,scal,nc,stem,isyk,
     %         asc,wrk)
      integer ns,n,m,nz,nc,isyk
      real a(ns,m),wrk(n),dx,dz,x0,z0,scal
      logical asc
      character*80 stem
*--------------------------------------------------------------
      integer itape,k,i,iz
    
      itape=37

      call noppg(stem,isyk,itape,asc)
      CALL SKHODEG(nc,n,nz,dx,dz,x0,z0,m,itape,asc)
      do 100 k=1,m
        do 50 i=1,n
          wrk(i)=scal*a(i,k)
 50     continue
        do 60 iz=1,nz
        CALL SKRIVG(wrk,1,1,n,1,1,1,n,1,itape,asc)
 60     continue
 100  continue  
      close(itape)

      return
      end



******************************************************************
*
*  The routine prompts for file-name, reads heading and 
*  reads first adatset, without checking
*
*  Identical to beolesfelt from bouss/trappe
*
*
************************************************************
      subroutine blesfelt(h,n,m,dx,dy,x0,y0,nc,ns,spor,
     %                def,navn,kf,ierr)
      integer n,m,nc,ns,ierr,kf
      real h(ns),dx,dy,x0,y0
      character*80 spor,def,navn
*---------------------------------------------------------
      integer itape,ndat,ntot
      logical asc

      ierr=0
      itape=55
      kf=79
      call filgen(itape,spor,def,navn,kf,asc,.true.)
      if(itape.lt.0 ) then
        write(0,*)'ERROR:beolesfelt --> filgen'
        ierr=2
        return
      end if  

      call headles(nc,n,m,dx,dy,x0,y0,ndat,itape,asc)

      ntot=n*m
      if(ntot.gt.ns) then
          write(0,*)'ERROR:beolesfelt  to many points'
          ierr=1
          return
      end if

      call lesg(h,1,1,n,m,1,1,n,m,itape,asc)
      close(itape)

      return
      end

******************************************************************
*      beglesf
*  The routine prompts for file-name, reads heading and 
*  reads first dataset, without checking.
*  The file is lest open
*    nc,n,m,dx,dy,x0,y0,ndat - header data              O
*    h - array for field                                O
*    ns - size of h                                     I
*    itape - unit number                                I 
*    spor,def,navn - question, default name and file name I,I,O
*    kf - length of name                                O
*    asc - binary or ascii                                O
*    lefir - =.true. implies first field is to be read    I
*    ierr - error flag, value=0 all is OK               O
*
************************************************************
      subroutine beglesf(nc,n,m,dx,dy,x0,y0,ndat,h,ns,itape,spor,
     %                def,navn,kf,asc,lefir,ierr)
      integer n,m,nc,ndat,ns,itape,ierr,kf
      real h(ns),dx,dy,x0,y0
      character*80 spor,def,navn
      logical asc,lefir
*---------------------------------------------------------
      integer ntot
  

      ierr=0
      kf=79
      call filgen(itape,spor,def,navn,kf,asc,.true.)
      if(itape.lt.0 ) then
        write(0,*)'ERROR:beglesf --> filgen'
        ierr=2
        return
      end if  

      call headles(nc,n,m,dx,dy,x0,y0,ndat,itape,asc)

      ntot=n*m
      if(ntot.gt.ns) then
          write(0,*)'ERROR:beglesf  to many points'
          ierr=1
          return
      end if

      if(lefir) then
      call lesg(h,1,1,n,m,1,1,n,m,itape,asc)
      end if

      return
      end


******************************************************************
*
*  The routine prompts for file-name, reads heading and 
*  reads first adatset, with checking of data
*
*        h - dataarray, filled as n times n                         O
*        nmax - maks antall data i y                                I
*        N,m - antall kolonner og maksimalt antall rader             I
*        dx,dy,x0,y0 - grid and lower left corner                   I
*        nlest - antall tall i siste rad. Dersom nlest=n er det     O
*                lest et helt antall rader
*        mlest - antall rader lest, helt eller påbegynt.            O
*        ierr= 0  Ok                                                O
*           = 1  datafile to large, as many rows as possible are read        
*           = 2  some data read, but not the entire array  
*           = 3  no data read
*           = 4  unable to open the file
*
************************************************************
      subroutine blsjekk(h,n,m,nlest,mlest,dx,dy,x0,y0,nc,ns,spor,
     %                def,navn,kf,ierr)
      integer n,m,nlest,mlest,nc,ns,ierr,kf
      real h(ns),dx,dy,x0,y0
      character*80 spor,def,navn
*---------------------------------------------------------
      integer itape,ndat,ntot,meff
      logical asc

      ierr=0
      itape=55
      kf=79
      call filgen(itape,spor,def,navn,kf,asc,.true.)
      if(itape.lt.0 ) then
        write(0,*)'ERROR:blsjekk --> filgen'
        ierr=4
        return
      end if  

      call headles(nc,n,m,dx,dy,x0,y0,ndat,itape,asc)

      ntot=n*m
      meff=m
      if(ntot.gt.ns) then
          write(0,*)'ERROR:blsjekk  too many points'
          ierr=1
          meff=ntot/n
          if(meff.eq.0) then
            ierr=3
          return
          end if
      end if

      call sjekkles(h,ns,N,meff,nlest,mlest,ITAPE,ASC)

      if(mlest.eq.0) ierr=3
      if(mlest.lt.meff.or.nlest.lt.n) ierr=2
      
      close(itape)

      return
      end

*************************************************************************
      SUBROUTINE SIIE(G,FI,N,DARG,FI0,FI1,M,EPS,IFLAG,NIT)
      INTEGER N,M,IFLAG,ANTP,NIT
      REAL FI(N),G,DARG,FI0,FI1,EPS
*--------------------------------------------------------
      integer nhelp,i,iit
      real dfi,arg,fih,argh,gh,fiv,fim,argv,gv,a,fix,feil
      real gx,fimid,argx

      IFLAG=0
      ANTP=N
      N=0
      DFI=(FI1-FI0)/M
      ARG=DARG
      FIH=FI0
      ARGH=0
      GH=G(FI0)
      DO 100 I=1,M
      FIV=FIH
      FIH=FIH+DFI
      FIM=0.5*(FIH+FIV)
      ARGV=ARGH
      GV=GH
      GH=G(FIH)
      ARGH=ARGV+DFI*(GH+GV+4.0*G(FIM))/6.0
      A=(ARGH-ARGV)/DFI
   50 CONTINUE
      IF(ARG.GT.ARGH) GO TO 70
      FIX=FIV+(ARG-ARGV)/A
      FEIL=10*EPS
      IIT=0
   60 CONTINUE
      IF(ABS(FEIL).LT.EPS.OR.IIT.GT.NIT) GO TO 65
      GX=G(FIX)
      IIT=IIT+1
      FIMID=(FIX+FIV)*0.5
      ARGX=ARGV+(FIX-FIV)*(GV+GX+4.0*G(FIMID))/6.0
      FEIL=(ARG-ARGX)/GX
      FIX=FIX+FEIL
      GO TO 60
   65 CONTINUE
      IF(ABS(FEIL).GT.EPS) THEN
        IFLAG=IFLAG+1                        
        nhelp=n+1
        write(0,*)'konvergenssvikt i siie'
        write(0,*)'n= ',nhelp,' feil=',feil,'  fix=',fix
      END IF
      N=N+1      
      FI(N)=FIX
      IF(N.GE.ANTP) RETURN
      ARG=ARG+DARG
      GO TO 50
   70 CONTINUE
  100 CONTINUE
      RETURN
      END



*************************************************************

      SUBROUTINE SECANT(F,X0,X1,XTOL,FTOL,NTOL,IFLAG)
      integer ntol,iflag
      real f,x0,x1,xtol,ftol
*-------------------------------------------------------
      integer n
      real f0,f1,deltax,deltaf

      IFLAG=0
      F0=F(X0)
      DELTAX=X1-X0
      DO 20 N=1,NTOL
      F1=F(X1)
      IF(ABS(F1).LE.FTOL)  GO TO 30
      DELTAF=F0-F1
      IF(DELTAF.EQ.0) GO TO 999
      DELTAX=F1/DELTAF*DELTAX
      X0=X1
      X1=X1+DELTAX
      IF(ABS(DELTAX).LE.XTOL) RETURN
   20 F0=F1
  999 IFLAG=2
      RETURN
   30 IFLAG=1
      RETURN
      END



******************************************************************
      SUBROUTINE SECAUT(F,X0,X1,XTOL,FTOL,NTOL,IFLAG)
      integer ntol,iflag
      real f,x0,x1,xtol,ftol
*-------------------------------------------------------
      integer n
      real f0,f1,deltax,deltaf

      write(0,*) 'N,X0,X1,F0,F1'
      IFLAG=0
      F0=F(X0)
      DELTAX=X1-X0
      DO 20 N=1,NTOL
      F1=F(X1)
      write(0,*) N,X0,X1,F0,F1

      IF(ABS(F1).LE.FTOL)  GO TO 30
      DELTAF=F0-F1
      IF(DELTAF.EQ.0) GO TO 999
      DELTAX=F1/DELTAF*DELTAX
      X0=X1
      X1=X1+DELTAX
      IF(ABS(DELTAX).LE.XTOL) RETURN
   20 F0=F1
  999 IFLAG=2
      RETURN
   30 IFLAG=1
      RETURN
      END


*************************************************************
      SUBROUTINE SECA2(H,X0,y0,X1,y1,XTOL,ytol,FTOL,Gtol,
     %         NTOL,ipr,IFLAG)
      integer ntol,ipr,iflag
      real x0,y0,x1,y1,xtol,ytol,ftol,gtol
      external H
*-------------------------------------------------------
      integer n,itry
      real f1,g1,deltax,deltay,fm,gm,fx,gx,fy,gy,det

      IFLAG=0
      itry=1
 
      deltay=y1-y0
      DELTAX=X1-X0

      if(deltax.eq.0.0 .or. deltay.eq.0.0) then
        iflag=5
        return
      end if

      DO 20 N=1,NTOL
      call h(x1,y1,f1,g1)
      IF(ABS(F1).LE.FTOL .and. abs(g1).lt.gtol) then
        if(ipr.eq.1) then
           write(0,'(a,2e14.6)')'utg, f,g=',f1,g1
        end if
        iflag=1
        return
      end if

 70   continue

      call h(x0,y1,fm,gm)
      fx=(f1-fm)/deltax
      gx=(g1-gm)/deltax   

      call h(x1,y0,fm,gm)
      fy=(f1-fm)/deltay
      gy=(g1-gm)/deltay   
      
      det=fx*gy-fy*gx

      IF(DET.EQ.0.0) then
        write(0,*)'seca2, det=0, fx,fy,gx,gy='
        write(0,*)fx,fy,gx,gy
        if(itry.ge.5) then
          iflag=2
          return
        else
          itry=itry+1
          deltax=10*xtol
          deltay=10*ytol
          x0=x1-deltax
          y0=y1-deltay
          go to 70
        end if
      end if

      deltax=-(f1*gy-g1*fy)/det
      deltay=-(g1*fx-f1*gx)/det
 
      if(deltax.ne.0.0) then
        X0=X1
        X1=X1+DELTAX
      else
        deltax=x1-x0
      end if

      if(deltay.ne.0.0) then
        y0=y1
        y1=y1+deltay
      else
        deltay=y1-y0
      end if

      if(ipr.eq.1) then
        write(0,'(i3,7e12.4)')n,x1,y1,deltax,deltay,f1,g1,det
      end if
      IF(ABS(DELTAX).LE.XTOL .and. abs(deltay).le.ytol) RETURN
 20   continue
 
      iflag=3
      RETURN
      END


*************************************************************                  
*                                                           *
*                 T R I                                     * 
*                                                           *  
*  'TRI' solves a tri-diagonal system of equations
*                                                           *                  
*    N = number of equations                             I  *                  
*    A = sub-diagonal elements in matrix              I     *                  
*    D = diagonal elements in matrix                  I     *                  
*    C = super-diagonal elements in matrix            I     *                 
*    B = right hand side at input,                    I/O   *                 
*        the solution at output                             *             
*                                                           *         
*************************************************************                
      SUBROUTINE TRI(A,D,C,B,N)                   
                    
      INTEGER N     
      REAL A(N),D(N),C(N),B(N)                    
*--------------------------------------------------------                     
      REAL XMULT    
      INTEGER I,NM,NMI                            
                    
      NM = N - 1    
                    
      DO 10 I=2,N   
         XMULT = A(I)/D(I-1)                      
         D(I)  = D(I) - XMULT*C(I-1)              
         B(I)  = B(I) - XMULT*B(I-1)              
                    
 10   CONTINUE      
                    
      B(N) = B(N)/D(N)                            
                    
      DO 20 I=1,NM  
                    
         NMI = N - I
         B(NMI) = ( B(NMI) - C(NMI)*B(NMI+1) )/D(NMI)                          
                    
 20   CONTINUE      
                    
      RETURN        
      END           




*************************************************************
*                                                           *
*                 T R I C                                     *
*                                                           *
*  'TRIC' LOESER ET TRIDIAGONALT KOMPLEKST LIKNINGSETT.     *
*                                                           *
*    N = ANTALL LIKNINGER                             I     *
*    A = ELEMENTER TIL VENSTRE FOR DIAGONALEN         I     *
*    D = ELEMENTER PAA DIAGONALEN                     I     *
*    C = ELEMENTER TIL HOEYRE FOR DIAGONALEN          I     *
*    B = HOEYRESIDEN AV LIKNINGEN (INPUT)             I/O   *
*      = LOESNINGEN X  (OUTPUT)                             *
*                                                           *
*************************************************************

      SUBROUTINE TRIC(A,D,C,B,N)
      INTEGER N
      complex  A(N),D(N),C(N),B(N)
*--------------------------------------------------------
      complex XMULT
      INTEGER I,NM,NMI
      NM = N - 1

      DO 10 I=2,N
         XMULT = A(I)/D(I-1)
         D(I)  = D(I) - XMULT*C(I-1)
         B(I)  = B(I) - XMULT*B(I-1)

 10   CONTINUE

      B(N) = B(N)/D(N)

      DO 20 I=1,NM

         NMI = N - I
         B(NMI) = ( B(NMI) - C(NMI)*B(NMI+1) )/D(NMI)

 20   CONTINUE

      RETURN
      END








**************************************************************************
*            Rutinen loeser en andregradslikning p} formen
*             
*               Y'' +F(x,Y,Y')=0
*    
*     vha. midtpunktdiskretisering og Newton iterasjon p} ikkelineariteten.
*     Parameterere
*              y - den ukjente. Ve kall inneholder den et gjett       I/O
*                  p} loesningen.
*              n -antall indre gitterpunkter.                         I
*              dx - gitteravstand                                     I
*              fver - subrutine som definerer F ovenfor.              I
*                     kallet m} se ut som
*                       call fver(x,yv,yder,f,fy,fyd)
*                     der yv og yd er hhv. verdi av y og midtpunktdiff
*                     for y, mens f,fy og fyd er output og svarer til hhv.
*                     F, dF/dY og dF/dY'
*              al,bl,ar,br - definerer randbetingelser ihht:             I
*                     y(0)+al*y(1)=bl
*                     y(n+1)+ar*y(n)=br
*              iit - maksimalt antall newton iterasjoner                 I
*              r   - relaksasjonsfaktor.                                 I
*              epy,epr - konvergenskrit for hhv. den ukjente og residu.  I
*              skr - verdi=.true. angir mellomutskrifter fra iterasjon   I
*                    denne kommer p} standard output.
*              iflag - feilparameter                                     O
*                    =0 konvergerte utfra epy
*                     1 konvergerte utfra epr
*                     2 konvergens ikke oppnaadd
*                     3 for mange punkter
*
***************************************************************************
      subroutine aglikn(y,n,dx,fver,al,bl,ar,br,iit,r,epy,epr
     %,skr,iflag)
      integer n,iit,iflag
      real y(0:n+1),dx,al,bl,ar,br,r,epy,epr
      logical skr
*-----------------------------------------------------------------------
      integer i,nmax,kit,np
      parameter (nmax=1002)
      real vdig(nmax),hdig(nmax),hs(nmax),dig(nmax),rdx2,rdx
      real x,y1,yder,f,fy,fyd,rmax,umax
      external fver

      if(n+2.gt.nmax) then
        write(0,*)'for mange punkter i aglikn'
        iflag=3
        return
      end if
  
      iflag=0
      rdx=1.0/dx
      rdx2=rdx*rdx
      np=n+1

      kit=0

 100  continue

      kit=kit+1
      if(kit.gt.iit) then
        iflag=2
        return
      end if

      rmax=0.0

      do 150 i=1,n
      x=i*dx
      y1=y(i)      
      yder=0.5*rdx*(y(i+1)-y(i-1))
      call fver(x,y1,yder,f,fy,fyd)
      hs(i)=-rdx2*(y(i+1)-2*y(i)+y(i-1))-f
      if(abs(hs(i)).gt.rmax) rmax=abs(hs(i))
      vdig(i)= rdx2 - 0.5*rdx*fyd
      dig(i)= -2*rdx2 + fy
      hdig(i)=rdx2 + 0.5*rdx*fyd
 150  continue

      if(skr) write(6,*)'kit,rmax=',kit,rmax

      if(rmax.le.epr) then
        iflag=1
        return
      end if

      dig(1)=dig(1)-al*vdig(1)
      hs(1)=hs(1)-vdig(1)*(bl -y(0)-al*y(1))
      dig(np)=1.0
      vdig(np)=ar
      hs(np)=br -y(np)-ar*y(n)


      call tri(vdig,dig,hdig,hs,np)

      umax=abs( bl-y(0)-al*y(1)-al*hs(1))
   
      do 160 i=1,np
      if(abs(hs(i)).gt.umax) umax=abs(hs(i))
      y(i)=y(i)+r*hs(i)
 160  continue

      if(skr) write(6,*)'umax=',umax
      y(0)=bl-al*y(1)
      if(umax.le.epy) return
     
      go to 100

      end
      










**************************************************************
      SUBROUTINE SMOX(A,NS,MS,NG,MG)
      integer ns,ms,ng,mg
      REAL A(0:NS,0:MS)
*------------------------------------------------------
      REAL AM,A0,AP
      INTEGER I,K ,NGM

      NGM=NG-1

      DO 100 K=0,MG
      A0=A(0,K)
      AP=A(1,K)
      A(0,K)=0.75*A0+0.5*AP-0.25*A(2,K)
      DO 50 I=1,NGM
      AM=A0
      A0=AP
      AP=A(I+1,K)
      A(I,K)=0.5*A0+0.25*(AM+AP)
  50  CONTINUE
      A(NG,K)=0.75*AP+0.5*A0-0.25*AM
 100  CONTINUE
      RETURN
      END






**************************************************************
      SUBROUTINE SMOY(A,NS,MS,NG,MG)
      INTEGER NS,MS,NG,MG
      REAL A(0:NS,0:MS)
*------------------------------------------------------
      REAL AM,A0,AP
      INTEGER I,K ,MGM

      MGM=MG-1

      DO 100 I=0,NG
      A0=A(I,0)
      AP=A(I,1)
      A(I,0)=0.75*A0+0.5*AP-0.25*A(I,2)
      DO 50 K=1,MGM
      AM=A0
      A0=AP
      AP=A(I,K+1)
      A(I,K)=0.5*A0+0.25*(AM+AP)
  50  CONTINUE
      A(I,MG)=0.75*AP+0.5*A0-0.25*AM
 100  CONTINUE
      RETURN
      END





*******************************************************************
      SUBROUTINE REG(XR,KY,IP,A,B,VAR)
      integer ip
      REAL XR(IP),A,B,VAR
      INTEGER KY(IP)
*-------------------------------------------------------------------
      REAL SX,SA,SXX,SAA,SXA,DET
      INTEGER I

      IF(IP.LT.2) THEN
        VAR=-1.0
        RETURN
      END IF

      SX=0.0
      SA=0.0
      SXA=0.0
      SAA=0.0
      SXX=0.0

      DO 100 I=1,IP
      SX=SX+XR(I)
      SXX=SX+XR(I)*XR(I)
      SA=SA+KY(I)
      SAA=SAA+KY(I)*KY(I)
      SXA=SXA+KY(I)*XR(I)
  100 CONTINUE

      DET=IP*SAA-SA*SA
      A=(SX*SAA-SXA*SA)/DET
      B=(IP*SXA-SA*SX)/DET
      VAR=SXX+IP*A*A+SAA*B*B-2.0*(SX*A+SXA*B+SA*A*B)

      RETURN
      END



**************************************************************
      SUBROUTINE NULL(A,N0,N1,M0,M1)
      INTEGER N0,N1,M0,M1,I,K
      REAL A(N0:N1,M0:M1)
*--------------------------------------------------------------
      DO 100 K=M0,M1
      DO 100 I=N0,N1
  100 A(I,K)=0.0
      RETURN
      END

**************************************************************
      SUBROUTINE skonst(A,N0,N1,M0,M1,val)
      INTEGER N0,N1,M0,M1,I,K
      REAL A(N0:N1,M0:M1),val
*--------------------------------------------------------------
      DO 100 K=M0,M1
      DO 100 I=N0,N1
  100 A(I,K)=val
      RETURN
      END

******************
      SUBROUTINE NOISE(A,ns,N,M,Mode,styrke)
      INTEGER Ns,n,m,mode
      REAL A(0:ns,0:m+1),styrke
*--------------------------------------------------------------
      integer k,i
      real pi,btall,ver

      pi=4.0*atan(1.0)

      btall = mode*pi/m

      DO 100 K=0,M+1
      ver=styrke*cos(btall*(k-0.5))
      DO 100 I=0,N+1
  100 A(I,K)=A(I,K)+ver


      btall = mode*pi/n

      DO 200 I=0,N+1
      ver=styrke*cos(btall*(I-0.5))
      DO 200 k=0,m+1
  200 A(I,K)=a(i,k)+ver

      RETURN
      END


*******************************************      
      SUBROUTINE BSPL(HV,BS,N,DELT,DER1,DER2)    
      INTEGER N
      REAL HV(N),BS(0:N+1),DELT,DER1,DER2                      
*--------------------------------------------------                            
      integer nmz
      parameter(nmz=3000)
      REAL WK(4*nmz)

      IF(N.GT.nmz) THEN
        WRITE(*,*)'for mange punkter i BSPL'
      ELSE  
        CALL  BSPLW(HV,BS,N,DELT,DER1,DER2,WK)
      END IF

      RETURN
      END 


*******************************************      
      SUBROUTINE BSPLOLD(HV,BS,N,DELT,DER1,DER2)    
      INTEGER N
      REAL HV(N),BS(0:N+1),DELT,DER1,DER2                      
*--------------------------------------------------                            
      REAL SUB(1500),SUP(1500),DIAG(1500),HS(1500)                             
      INTEGER S
                        
      DO 100 S=1,N      
      SUB(S)=1.0       
      SUP(S)=1.0    
      DIAG(S)=4.0   
      HS(S)=6.0*HV(S)                             
  100 CONTINUE      
      SUP(1)=2.0    
      HS(1)=HS(1)+2.0*DELT*DER1                   
      SUB(N)=2.0    
      HS(N)=HS(N)-2.0*DELT*DER2                   
      CALL TRI(SUB,DIAG,SUP,HS,N)                 
      DO 200 S=1,N  
      BS(S)=HS(S)   
  200 CONTINUE      
      BS(0)=BS(2)-2.0*DELT*DER1                   
      BS(N+1)=BS(N-1)+2.0*DELT*DER2               
      RETURN        
      END           
                    


*******************************************    
*
*       Automatisk setting av ende-deriverete i 
*       splineinterpolasjon.
*****************************************************
      SUBROUTINE BSPLAUT(HV,BS,N,W)    
      INTEGER N
      REAL HV(N),BS(0:N+1),W(4*N)                      
*--------------------------------------------------
      REAL DELT,DER1,DER2

      DELT=1.0
      DER1= 2.0*hv(2)-1.5*hv(1)-0.5*hv(3)
      DER2= -( 2.0*hv(n-1)-1.5*hv(n)-0.5*hv(n-2))

      call BSPLW(HV,BS,N,DELT,DER1,DER2,W)    

      return
      end

                    

*******************************************      
      SUBROUTINE BSPLW(HV,BS,N,DELT,DER1,DER2,W)    
      INTEGER N
      REAL HV(N),BS(0:N+1),DELT,DER1,DER2,W(4*N)                      
*--------------------------------------------------
      INTEGER S,ISUB,ISUP,IDIAG,IHS
      INTEGER ISUBA,ISUPA,IDIAGA,IHSA

      ISUB=1
      IDIAG=N+1
      ISUP=2*N+1
      IHS=3*N+1
      ISUBA=ISUB-1
      IDIAGA=IDIAG-1
      ISUPA=ISUP-1
      IHSA=IHS-1
                        
      DO 100 S=1,N      
      W(ISUBA+S)=1.0       
      W(ISUPA+S)=1.0    
      W(IDIAGA+S)=4.0   
      W(IHSA+S)=6.0*HV(S)                             
  100 CONTINUE      
      W(ISUPA+1)=2.0    
      W(IHSA+1)=W(IHSA+1)+2.0*DELT*DER1                   
      W(ISUBA+N)=2.0    
      W(IHSA+N)=W(IHSA+N)-2.0*DELT*DER2                   
      CALL TRI(W(ISUB),W(IDIAG),W(ISUP),W(IHS),N)                 
      DO 200 S=1,N  
      BS(S)=W(IHSA+S)   
  200 CONTINUE      
      BS(0)=BS(2)-2.0*DELT*DER1                   
      BS(N+1)=BS(N-1)+2.0*DELT*DER2               
      RETURN        
      END           
                    

************************************************************
*                 BSPLGEN
*
*     HV - array of values to be interpolated                         I
*     BS - (0:N+1) array with spline coefiicients                     O
*     N - nymber of points                                            I
*     DELT - grid increment (constant)                                I
*     ik1,ik2 type of boundary condition for left and right boundary   I
*         respectively.
*         value=0 : derivative adapted to value d1,d2
*         value=1 : one sided differences used for derivative
*         value=2 : second derivative set to d1 and d2.
*                   if hv(0)=0 and d1=0 we then obtain an antisymmetric
*                   condition at left boundary etc. 
*     d1,d2 - values of derivativea t left and right boundary          I
*     w  - work array, at least 4*N long                                D 
*
****************************************************************************
      SUBROUTINE BSPLGEN(HV,BS,N,DELT,ik1,D1,ik2,D2,W)    
      INTEGER N,ik1,ik2
      REAL HV(N),BS(0:N+1),DELT,D1,D2,W(4*N)                      
*--------------------------------------------------
      INTEGER S,ISUB,ISUP,IDIAG,IHS
      INTEGER ISUBA,ISUPA,IDIAGA,IHSA
      real der1,der2,p1,p2,r1,r2,q1,q2,h1,h2

c
c     The left hand boundary condition, say, is formulated as
c        p1*bs(0)+q1*bs(1)+r1*bs(2)=h1
c     where the coefficients depend on the actual condition. To retain
c     the tri-diagonal structure this equation is then used to eliminate
c     bs(0) from the interpolation equation for hv(1): 
c         bs(0)+4*bs(1)+bs(2)=6*hv(1)
c

      if(ik1.lt.2) then
        if(ik1.eq.0) then
          der1=d1*delt
       else
          if(n.ge.3) then
             der1=2.0*hv(2)-1.5*hv(1)-0.5*hv(3)
          else
             der1=hv(2)-hv(1)
          end if
          
        end if
        p1=-0.5
        q1=0.0
        r1=0.5
        h1=der1
      else
        p1=1.0
        q1=-2.0
        r1=1.0
        h1=d1
      end if
        
      if(ik2.lt.2) then
        if(ik2.eq.0) then
          der2=d2*delt
        else
          if(n.ge.3) then
             der2=-( 2.0*hv(n-1)-1.5*hv(n)-0.5*hv(n-2))
          else
             der2=hv(n)-hv(n-1)
          end if
          
        end if
        p2=-0.5
        q2=0.0
        r2=0.5
        h2=der2
      else
        p2=1.0
        q2=-2.0
        r2=1.0
        h2=d2
      end if
        


      ISUB=1
      IDIAG=N+1
      ISUP=2*N+1
      IHS=3*N+1
      ISUBA=ISUB-1
      IDIAGA=IDIAG-1
      ISUPA=ISUP-1
      IHSA=IHS-1
                        
      DO 100 S=1,N      
      W(ISUBA+S)=1.0       
      W(ISUPA+S)=1.0    
      W(IDIAGA+S)=4.0   
      W(IHSA+S)=6.0*HV(S)                             
  100 CONTINUE 
      if(p1.ne.0.0) then     
        W(ISUP)=W(isup)-r1/p1
        W(idiag)=w(idiag)-q1/p1    
        W(IHS)=W(IHS)-h1/p1
      else
        W(ISUP)=r1
        W(idiag)=q1    
        W(IHS)=h1
      end if   

      if(r2.ne.0.0) then     
        W(ISUBA+N)=W(ISUBA+N)-p2/r2   
        W(IHSA+N)=W(IHSA+N)-h2/r2
        W(idiaga+n)=w(idiaga+n)-q2/r2    
      else
        W(ISUBA+N)=p2
        W(idiaga+n)=q2    
        W(IHSA+N)=h2
      end if
   

      CALL TRI(W(ISUB),W(IDIAG),W(ISUP),W(IHS),N)                 
      DO 200 S=1,N  
      BS(S)=W(IHSA+S)   
  200 CONTINUE      

      if(p1.ne.0.0) then     
        BS(0)=-q1*bs(1)/p1 -r1*BS(2)/p1 +h1/p1
      else
        BS(0)=-4.0*bs(1) -BS(2) +6.0*hv(1)
      end if

      if(r2.ne.0.0) then     
        BS(n+1)=-q2*bs(n)/r2 -p2*BS(n-1)/r2 +h2/r2
      else
        BS(n+1)=-4.0*bs(n) -BS(n-1) +6.0*hv(N)
      end if


      RETURN        
      END           
                    
                    
                    
                    
                    
**************************************************************************
*
*                     S P D E R
*   
*     Beregner funksjonsverdier av en splines-interpolant generert ved
*     bsplw. 
*      parametere:
*              X - koordinat
*              X0 - posisjon tilh|rende B-spline med koeff. bs(1)       I
*              DELT - punktavstand                                      I
*              BS - array med spline-koeff.                             I
*              NP - antall interpolasjonspunkter (np+2 B-splines)       I
*              d0,d1,d2 - hhv. funk.ver 1.ste og 2.dre derivert         O
*****************************************************************************
      subroutine spder(X,X0,DELT,BS,NP,d0,d1,d2)                
      INTEGER NP    
      REAL BS(0:NP+1),X,X0,DELT,d0,d1,d2                
*------------------------------------------------ 
      REAL XV,S,sx,bm,bp,bmm,bpp     
      INTEGER NX    
                  
      sx=1.0/delt  
      XV=(X-X0)/DELT+1.0                          
      NX=XV         
      S=XV-NX       
      XV=S-1.0      
      if(nx.lt.0 .or. nx.gt.np) then
c       write(0,*)'punkt utenfor interval i spder'
        return
      end if
                    
      bm=bs(nx)
      bp=bs(nx+1)

      IF(NX.EQ.0) THEN                            
        bmm=0.0
      ELSE          
        bmm=bs(nx-1)
      END IF        

      IF(NX.EQ.np) THEN                            
        bpp=0.0
      ELSE          
        bpp=bs(nx+2)
      END IF        

      d0=-bmm*XV*XV*XV/6.0+bm*(2.0/3.0-S*S*(1.0-0.5*S))          
     %   +Bp*(2.0/3.0-XV*XV*(1+0.5*XV))+bpp*S*S*S/6.0               

      d1=sx*(-bmm*XV*XV*0.5+bm*(-2.0*S+1.5*s*s)          
     %   -Bp*(2.0*xv+1.5*xv*xv)+bpp*S*S*0.5)               

      d2=sx*sx*(-bmm*XV + bm*(-2.0+3.0*s)          
     %   -Bp*(2.0+3.0*xv)+bpp*S)               


      RETURN        
      END           

                    
                    
                    
**************************************************************************
*
*                     G P D E R
*   
*     Beregner funksjonsverdier av en splines-interpolant generert ved
*     bsplw. 
*      parametere:
*              X - koordinat
*              X0 - posisjon tilh|rende B-spline med koeff. bs(1)       I
*              DELT - punktavstand                                      I
*              BS - array med spline-koeff.                             I
*              NP - antall interpolasjonspunkter (np+2 B-splines)       I
*              d0,d1,d2,d3 - hhv. funk.ver 1.ste og 2.dre derivert      O
*                           Mark: third derivative is not continuous
*****************************************************************************
      subroutine gpder(X,X0,DELT,BS,NP,d0,d1,d2,d3)                
      INTEGER NP    
      REAL BS(0:NP+1),X,X0,DELT,d0,d1,d2,d3        
*------------------------------------------------ 
      REAL XV,S,sx,bm,bp,bmm,bpp     
      INTEGER NX    
                  
      sx=1.0/delt  
      XV=(X-X0)/DELT+1.0                          
      NX=XV         
      S=XV-NX       
      XV=S-1.0      
      if(nx.lt.0 .or. nx.gt.np) then
c       write(0,*)'punkt utenfor interval i spder'
        return
      end if
                    
      bm=bs(nx)
      bp=bs(nx+1)

      IF(NX.EQ.0) THEN                            
        bmm=0.0
      ELSE          
        bmm=bs(nx-1)
      END IF        

      IF(NX.EQ.np) THEN                            
        bpp=0.0
      ELSE          
        bpp=bs(nx+2)
      END IF        

      d0=-bmm*XV*XV*XV/6.0+bm*(2.0/3.0-S*S*(1.0-0.5*S))          
     %   +Bp*(2.0/3.0-XV*XV*(1+0.5*XV))+bpp*S*S*S/6.0               

      d1=sx*(-bmm*XV*XV*0.5+bm*(-2.0*S+1.5*s*s)          
     %   -Bp*(2.0*xv+1.5*xv*xv)+bpp*S*S*0.5)               

      d2=sx*sx*(-bmm*XV + bm*(-2.0+3.0*s)          
     %   -Bp*(2.0+3.0*xv)+bpp*S)               


      d3=sx*sx*sx*(-bmm + 3.0*bm          
     %   -Bp*3.0+bpp)               


      RETURN        
      END           

                    
                    
*********************************************************                       
      FUNCTION HG(X,X0,DELT,BS,NP)                
      INTEGER NP    
      REAL BS(0:NP+1),X,X0,DELT,HG                
*------------------------------------------------ 
      REAL XV,S     
      INTEGER NX    
                    
      XV=(X-X0)/DELT+1.0                          
      NX=XV         
      S=XV-NX       
      XV=S-1.0      
                    
      IF(NX.EQ.0) THEN                            
      HG=BS(NX)*(2.0/3.0-S*S*(1.0-0.5*S))         
     %   +BS(NX+1)*(2.0/3.0-XV*XV*(1+0.5*XV))+BS(NX+2)*S*S*S/6.0                
      ELSE          
      IF(NX.GE.NP) THEN                           
      HG=-BS(NX-1)*XV*XV*XV/6.0+BS(NX)*(2.0/3.0-S*S*(1.0-0.5*S))          
     %   +BS(NX+1)*(2.0/3.0-XV*XV*(1+0.5*XV))     
      ELSE          
      HG=-BS(NX-1)*XV*XV*XV/6.0+BS(NX)*(2.0/3.0-S*S*(1.0-0.5*S))          
     %   +BS(NX+1)*(2.0/3.0-XV*XV*(1+0.5*XV))+BS(NX+2)*S*S*S/6.0               
      END IF        
      END IF        
      RETURN        
      END           




******************************************************************
*
*                   S T A G
*
*     Innholdet i et en-D gitter overf|res til et annet med forskj|vet
*     startpunkt og annen gitteravstand. Det benyttes spline interpolasjon.
*     parametere:
*             a -  array for kilde-gitter                                 I
*             na,dxa,xa - hhv antall punkter, gitteravstand og pos av    I
*                         f|rste punkt i kilde-gitter 
*             b - array for kopigitter                                   O
*             nb,dxb,xb - hhv antall punkter, gitteravstand og pos av    I
*                         f|rste punkt i kopi-gitter. Dersom dette gitteret 
*                         g}r utenfor  kilde-gitteret etterfylles med 
*                         hhv. a(1) og a(na).
*             iflag - feilparameter                                      O
*             w  - arbeidsarray                                          D
*
***************************************************************************
      subroutine stag(a,na,dxa,xa,b,nb,dxb,xb,iflag,w)
      integer na,nb,iflag
      real a(na),dxa,xa,b(nb),dxb,xb,w(5*na+15)
*------------------------------------------------------------------
      integer k,k1,k2,na3,k1m
      real xver,hg,x1,x2,der1,der2

      na3=na+3
      iflag=0

      x1=max(xa,xb)
      x2=min(xa+(na-1)*dxa,xb+(nb-1)*dxb)

      if(x1.ge.x2) then
        iflag=1
        return
      end if

      k1=(x1-xb)/dxb+1
      if( (xb+(k1-1)*dxb).lt.x1) k1=k1+1
      k2=(x2-xb)/dxb+1

      DER1=0.5*(4.0*a(2)-3.0*A(1)-A(3))/dxa
      DER2=-0.5*(4.0*A(Na-1)-3.0*A(Na) -A(Na-2))/dxa

      call bsplw(a,w(1),na,dxa,der1,der2,w(na3))

      do 100 k=k1,k2
        xver=(k-1)*dxb+xb
        b(k)=hg(xver,xa,dxa,w(1),na)
 100    continue

      k1m=k1-1

      do 200 k=1,k1m
        b(k)=a(1)
 200    continue

      do 300 k=k2+1,nb
        b(k)=a(na)
 300    continue

      
      return

      end 

**************************************************************
*     OLD variant, impelemented as wrap of staggeg
***********************************************************
      subroutine stagge(a,na,dxa,xa,b,nb,dxb,xb,eks,bleft,bright,
     %           iflag,w)
      integer na,nb,iflag
      real a(na),dxa,xa,b(nb),dxb,xb,bleft,bright,w(5*na+15)
      logical eks
*------------------------------------------------------------------
      real lddum,hddum
      call staggeg(a,na,dxa,xa,b,nb,dxb,xb,eks,bleft,bright,
     %           .false.,lddum,.false.,hddum,iflag,w)
      return
      end


******************************************************************
*
*                   S T A G G E G
*
*     Innholdet i et en-D gitter overf|res til et annet med forskj|vet
*     startpunkt og annen gitteravstand. Det benyttes spline interpolasjon.
*     parametere:
*             a -  array for kilde-gitter                                 I
*             na,dxa,xa - hhv antall punkter, gitteravstand og pos av    I
*                         f|rste punkt i kilde-gitter 
*             b - array for kopigitter                                   O
*             nb,dxb,xb - hhv antall punkter, gitteravstand og pos av    I
*                         f|rste punkt i kopi-gitter. 
*             eks       -Dersom dette gitteret                          I
*                         g}r utenfor  kilde-gitteret etterfylles med 
*                         hhv. a(1) og a(na) hvis eks er sann, med
*                         bleft (utenfor a(1)) og bright ellers
*             bleft, bright - verdier som settes paa b der den gaar      I
*                             utenfor a
*             lse,hse  - mark if boundary derivatives in spline approximant I
*                        are set or are to be  
*                        found from one-sided differences.
*                        lse=.true. indicates that a derivative to the
*                        left is provided etc.
*             lder,hder - Values for boundary derivatives in spline      I
*                         approximant. Used if lse or hse are true. 
*             iflag - feilparameter                                      O
*             w  - arbeidsarray                                          D
*
***************************************************************************
      subroutine staggeg(a,na,dxa,xa,b,nb,dxb,xb,eks,bleft,bright,
     %           lse,lder,hse,hder,iflag,w)
      integer na,nb,iflag
      real a(na),dxa,xa,b(nb),dxb,xb,lder,hder,bleft,bright,w(5*na+15)
      logical eks,lse,hse
*------------------------------------------------------------------
      integer k,k1,k2,na3,k1m
      real xver,hg,x1,x2,der1,der2,lval,rval

      na3=na+3
      iflag=0

      x1=max(xa,xb)
      x2=min(xa+(na-1)*dxa,xb+(nb-1)*dxb)

      if(x1.ge.x2) then
        iflag=1
        return
      end if

      k1=(x1-xb)/dxb+1
      if( (xb+(k1-1)*dxb).lt.x1-0.001*dxb) k1=k1+1
      k2=(x2-xb)/dxb+1
      if( (xb+(k2)*dxb).lt.x2+0.001*dxb) k2=k2+1

      if(lse) then
        der1=lder
      else
        DER1=0.5*(4.0*a(2)-3.0*A(1)-A(3))/dxa
      end if

      if(hse) then
        der2=hder
      else
        DER2=-0.5*(4.0*A(Na-1)-3.0*A(Na) -A(Na-2))/dxa
      end if

      call bsplw(a,w(1),na,dxa,der1,der2,w(na3))

      do 100 k=k1,k2
        xver=(k-1)*dxb+xb
        b(k)=hg(xver,xa,dxa,w(1),na)
 100    continue

      k1m=k1-1

      if(eks) then
       rval=a(na)
       lval=a(1)
      else
       rval=bright
       lval=bleft
      end if

      do 200 k=1,k1m
        b(k)=lval
 200    continue

      do 300 k=k2+1,nb
        b(k)=rval
 300    continue

      
      return

      end 

******************************************************************
*
*                   S T A G G E G
*
*     Innholdet i et en-D gitter overf|res til et annet med forskj|vet
*     startpunkt og annen gitteravstand. Det benyttes spline interpolasjon.
*     parametere:
*             a -  array for kilde-gitter                                 I
*             na,dxa,xa - hhv antall punkter, gitteravstand og pos av    I
*                         f|rste punkt i kilde-gitter 
*             b - array for kopigitter                                   O
*             nb,dxb,xb - hhv antall punkter, gitteravstand og pos av    I
*                         f|rste punkt i kopi-gitter. 
*             irvalg     -Dersom dette gitteret                          I
*                         g}r utenfor  kilde-gitteret etterfylles med: 
*                         irvalg=0 hhv. a(1) og a(na)  
*                         irvalg=1 bleft (utenfor a(1)) og bright ellers
*                         irvalg=2 intet settes utenfor
*             spline - .true. for spline, linear otherwise
*             bleft, bright - verdier som settes paa b der den gaar      I
*                             utenfor a
*             lse,hse  - mark if boundary derivatives in spline approximant I
*                        are set or are to be  
*                        found from one-sided differences.
*                        lse=.true. indicates that a derivative to the
*                        left is provided etc.
*             lder,hder - Values for boundary derivatives in spline      I
*                         approximant. Used if lse or hse are true. 
*             iflag - feilparameter                                      O
*             eps - toleranse for innenfor..                             I
*             w  - arbeidsarray                                          D
*
***************************************************************************
      subroutine gstaggeg(a,na,dxa,xa,b,nb,dxb,xb,irvalg,spline,
     %           bleft,bright,lse,lder,hse,hder,iflag,eps,w)
      integer na,nb,irvalg,iflag
      real a(na),dxa,xa,b(nb),dxb,xb,lder,hder,bleft,bright,eps
      real w(5*na+15)
      logical spline,lse,hse
*------------------------------------------------------------------
      integer k,k1,k2,na3,k1m,nk
      real xver,hg,x1,x2,der1,der2,lval,rval,xbk,umdef

      umdef=-666.0
      na3=na+3
      iflag=0

      x1=max(xa,xb)
      x2=min(xa+(na-1)*dxa,xb+(nb-1)*dxb)

      if(x1.gt.x2) then
        iflag=1
        return
      end if

      k1=(x1-xb)/dxb+1
      if( (xb+(k1-1)*dxb).lt.x1-eps*dxb) k1=k1+1
      k2=(x2-xb)/dxb+1
      if( (xb+(k2)*dxb).lt.x2+eps*dxb) k2=k2+1

      if(spline) then
        if(lse) then
          der1=lder
        else
          DER1=0.5*(4.0*a(2)-3.0*A(1)-A(3))/dxa
        end if

        if(hse) then
          der2=hder
        else
          DER2=-0.5*(4.0*A(Na-1)-3.0*A(Na) -A(Na-2))/dxa
        end if

        call bsplw(a,w(1),na,dxa,der1,der2,w(na3))

        do 100 k=k1,k2
          xver=(k-1)*dxb+xb
          b(k)=hg(xver,xa,dxa,w(1),na)
 100    continue
      else
        do 110 k=1,na
          w(k)=(k-1)*dxa+xa
 110    continue
        nk=(k2-k1)+1
        xbk=xb+(k1-1)*dxb
        call ureglin(w,a,na,b(k1),nk,dxb,xbk,.true.,umdef,umdef,
     %                   eps,iflag)
      end if


      if(irvalg.eq.2) return
      k1m=k1-1

      if(irvalg.eq.0) then
       rval=a(na)
       lval=a(1)
      else
       rval=bright
       lval=bleft
      end if

      do 200 k=1,k1m
        b(k)=lval
 200    continue

      do 300 k=k2+1,nb
        b(k)=rval
 300    continue

      
      return

      end 

******************************************************************
*
*                   K S T A G E G
*
*     Innholdet i et en-D gitter overf|res til et annet med ikke-uniform
*     oppdeling. Det benyttes spline interpolasjon.
*     parametere:
*             a -  array for kilde-gitter                                 I
*             na,dxa,xa - hhv antall punkter, gitteravstand og pos av    I
*                         f|rste punkt i kilde-gitter 
*             b - array for kopigitter                                    O
*             xb - array for b's gitterposisjoner                         I
*                     Dersom dette gitteret 
*                     g}r utenfor  kilde-gitteret etterfylles med 
*                     hhv. a(1) og a(na).
*             nb - antall punkter i b                                    I
*             eks       -Dersom kopi-gitteret                          I
*                         g}r utenfor  kilde-gitteret etterfylles med 
*                         hhv. a(1) og a(na) hvis eks er sann, med
*                         bleft (utenfor a(1)) og bright ellers
*             bleft, bright - verdier som settes paa b der den gaar      I
*                             utenfor a
*             lse,hse  - mark if boundary derivatives in spline approximant I
*                        are set or are to be  
*                        found from one-sided differences.
*                        lse=.true. indicates that a derivative to the
*                        left is provided etc.
*             lder,hder - Values for boundary derivatives in spline      I
*                         approximant. Used if lse or hse are true. 
*             iflag - feilparameter                                      O
*             w  - arbeidsarray                                          D
*
***************************************************************************
      subroutine kstageg(a,na,dxa,xa,b,xb,nb,eks,bleft,bright,
     %           lse,lder,hse,hder,iflag,w) 
      integer na,nb,iflag
      real a(na),dxa,xa,b(nb),xb(nb),lder,hder,bleft,bright,w(5*na+15)
      logical eks,lse,hse
*------------------------------------------------------------------
      integer k,k1,k2,na3,k1m
      real hg,x1,x2,der1,der2,lval,rval

      na3=na+3
      iflag=0

      x1=max(xa,xb(1))
      x2=min(xa+(na-1)*dxa,xb(nb))

      if(x1.ge.x2) then
        iflag=1
        return
      end if

      k1=1
 50   continue
      if(xb(k1).ge.x1) go to 60
      k1=k1+1
      go to 50

 60   continue

      k2=nb
 70   continue
      if(xb(k2).le.x2) go to 80
      k2=k2-1
      go to 70

 80   continue


      if(lse) then
        der1=lder
      else
        DER1=0.5*(4.0*a(2)-3.0*A(1)-A(3))/dxa
      end if

      if(hse) then
        der2=hder
      else
        DER2=-0.5*(4.0*A(Na-1)-3.0*A(Na) -A(Na-2))/dxa
      end if

      call bsplw(a,w(1),na,dxa,der1,der2,w(na3))

      do 100 k=k1,k2
        b(k)=hg(xb(k),xa,dxa,w(1),na)
 100    continue

      k1m=k1-1

      if(eks) then
       rval=a(na)
       lval=a(1)
      else
       rval=bright
       lval=bleft
      end if

      do 200 k=1,k1m
        b(k)=lval
 200    continue

      do 300 k=k2+1,nb
        b(k)=rval
 300    continue

      
      return

      end 


******************************************************************
*
*                   Q S T G
*
*     Innholdet i et en-D gitter overf|res til et annet med ikke-uniform
*     oppdeling. Det benyttes spline interpolasjon. More robust than
*     kstageg, in the sense that xb does not need to be monotoneous 
*     parametere:
*             a -  array for kilde-gitter                                 I
*             na,dxa,xa - hhv antall punkter, gitteravstand og pos av    I
*                         f|rste punkt i kilde-gitter 
*             b - array for kopigitter                                    O
*             xb - array for b's gitterposisjoner                         I
*                     Dersom dette gitteret 
*                     g}r utenfor  kilde-gitteret etterfylles med 
*                     hhv. a(1) og a(na).
*             nb - antall punkter i b                                    I
*             eks       -Dersom kopi-gitteret                          I
*                         g}r utenfor  kilde-gitteret etterfylles med 
*                         hhv. a(1) og a(na) hvis eks er sann, med
*                         bleft (utenfor a(1)) og bright ellers
*             bleft, bright - verdier som settes paa b der den gaar      I
*                             utenfor a
*             lse,hse  - mark if boundary derivatives in spline approximant I
*                        are set or are to be  
*                        found from one-sided differences.
*                        lse=.true. indicates that a derivative to the
*                        left is provided etc.
*             lder,hder - Values for boundary derivatives in spline      I
*                         approximant. Used if lse or hse are true. 
*             iflag - feilparameter                                      O
*             w  - arbeidsarray                                          D
*
***************************************************************************
      subroutine qstg(a,na,dxa,xa,b,xb,nb,eks,bleft,bright,
     %           lse,lder,hse,hder,iflag,w) 
      integer na,nb,iflag
      real a(na),dxa,xa,b(nb),xb(nb),lder,hder,bleft,bright,w(5*na+15)
      logical eks,lse,hse
*------------------------------------------------------------------
      integer k,na3,k1m
      real hg,x1,x2,der1,der2,lval,rval

      na3=na+3
      iflag=0

      x1=xa
      x2=xa+(na-1)*dxa


      if(eks) then
       rval=a(na)
       lval=a(1)
      else
       rval=bright
       lval=bleft
      end if

      if(lse) then
        der1=lder
      else
        DER1=0.5*(4.0*a(2)-3.0*A(1)-A(3))/dxa
      end if

      if(hse) then
        der2=hder
      else
        DER2=-0.5*(4.0*A(Na-1)-3.0*A(Na) -A(Na-2))/dxa
      end if

      call bsplw(a,w(1),na,dxa,der1,der2,w(na3))

      do 100 k=1,nb
        if(xb(k).le.x1) then
          b(k)=lval
        else
          if(xb(k).le.x2) then
            b(k)=hg(xb(k),xa,dxa,w(1),na)
          else
            b(k)=rval
          end if
        end if
 100  continue


      
      return

      end 


******************************************************************
*
*                   UREGLIN
*
*     A one-D dataset given on a non-uniform grid is interpolated onto
*     a uniform grid by linear interpolation. The source grid must 
*     contain coordinates in increasing order.
*     ERRORS:      do not work for nb=1
*     parameters:
*             xs, ys -  arrays for original grid and data               I
*                       The xs must increase monotonously
*             na - number of points in original grid
*             b - array for interpolated values on regular grid          O
*             nb,dxb,xb - number of points, grid increment and position  I
*                         of first point in regular grid 
*             eks       - If the regular grid exeeds the original
*                         grid this parameter defines the extrapolation. I
*                         If eks is .true. a(1) and a(n) is used for
*                         regular nodes that are outside the original grid.
*                         Otherwise the values  bleft ( left of xs(1)) and
*                         bright are used 
*             bleft, bright - see above                                  I
*             eps  - tolerance for defining a point inside original grid   I
*             iflag - feilparameter                                      O
***************************************************************************
      subroutine ureglin(xs,ys,na,b,nb,dxb,xb,eks,bleft,bright,
     %                   eps,iflag)
      integer na,nb,iflag
      real xs(na),ys(na),b(nb),dxb,xb,bleft,bright,eps
      logical eks
*------------------------------------------------------------------
      integer k,k1,k2,k1m,ipos
      real xver,x1,x2,lval,rval,xm,xp,wg

 
      iflag=0

      x1=max(xs(1),xb)
      x2=min(xs(na),xb+(nb-1)*dxb)


      if(na.eq.1) then
        iflag=2
        return
      end if

      if(x1.ge.x2) then
        if(xb.lt.xs(1)) then
          k1=nb+1
          k2=nb
        else
          k1=1
          k2=0
        end if
      else
        k1=(x1-xb)/dxb+1
        if( (xb+(k1-1)*dxb).lt.x1-eps*dxb) k1=k1+1
        k2=(x2-xb)/dxb+1
        if( (xb+k2*dxb).lt.x2+eps*dxb) k2=k2+1
      end if


      ipos=1

      do 100 k=k1,k2
        xver=(k-1)*dxb+xb
 60     continue
        if(xver.le.xs(ipos+1) .or. ipos.eq.(na-1)) go to 80
        ipos=ipos+1
        go to 60
 80     continue

        xp=xs(ipos+1)
        xm=xs(ipos)
        if(xp.le.xm) then
          iflag=3
          return
        end if

        wg=(xver-xm)/(xp-xm)
       
        b(k)=(1.0-wg)*ys(ipos)+wg*ys(ipos+1)

 100  continue

      k1m=k1-1

      if(eks) then
       rval=ys(na)
       lval=ys(1)
      else
       rval=bright
       lval=bleft
      end if

      do 200 k=1,k1m
        b(k)=lval
 200    continue

      do 300 k=k2+1,nb
        b(k)=rval
 300    continue

      
      return

      end 

******************************************************************
*
*                   GUREGLIN
*
*     A one-D dataset given on a non-uniform grid is interpolated onto
*     another non-uniform grid by linear interpolation. Both grids must 
*     contain coordinates in increasing order.
*     ERRORS:      do not work for nb=1
*     parameters:
*             xs, ys -  arrays for original grid and data               I
*                       The xs must increase monotonously
*             na - number of points in original grid
*             xb - array for interpolated grid                            I
*             yb - array for interpolated values                         O
*             nb  - number of interpolated points                        I
*                         of first point in regular grid 
*             eks       - If the regular grid exeeds the original
*                         grid this parameter defines the extrapolation. I
*                         If eks is .true. a(1) and a(n) is used for
*                         regular nodes that are outside the original grid.
*                         Otherwise the values  bleft ( left of xs(1)) and
*                         bright are used 
*             bleft, bright - see above                                  I
*             eps  - tolerance for defining a point inside original grid   I
*             iflag - error parameter                O
***************************************************************************
      subroutine gureglin(xs,ys,na,xb,yb,nb,eks,bleft,bright,
     %                   eps,iflag)
      integer na,nb,iflag
      real xs(na),ys(na),yb(nb),xb(nb),bleft,bright,eps
      logical eks
*------------------------------------------------------------------
      integer k,k1,k2,k1m,ipos
      real xver,x1,x2,lval,rval,xm,xp,wg

 
      iflag=0

      x1=xs(1)
      x2=xs(na)


      if(na.eq.1) then
        iflag=2
        return
      end if


      if(eks) then
       rval=ys(na)
       lval=ys(1)
      else
       rval=bright
       lval=bleft
      end if

      ipos=1

      do 100 k=1,nb
        xver=xb(k)
 60     continue
        if(xver.le.xs(ipos+1) .or. ipos.eq.(na-1)) go to 80
        ipos=ipos+1
        go to 60
 80     continue

        xp=xs(ipos+1)
        xm=xs(ipos)
        if(xp.le.xm) then
          iflag=3
          return
        end if

        if(xver.lt.x1-eps) then
           yb(k)=lval
        else
         if(xver.gt.x2+eps) then
           yb(k)=rval
         else
           wg=(xver-xm)/(xp-xm)
       
           yb(k)=(1.0-wg)*ys(ipos)+wg*ys(ipos+1)
         end if
        end if

 100  continue


      
      return

      end 

******************************************************************
*
*                   NUREGLIN
*
*     A one-D dataset given on a non-uniform grid is interpolated onto
*     a uniform grid by linear interpolation. The source grid must 
*     contain coordinates in increasing order.
*     parameters:
*             xs, ys -  arrays for original grid and data               I
*                       The xs must increase monotonously
*             na - number of points in original grid
*             yb - array for interpolated values on regular grid          O
*             nb,dxb,xb - number of points, grid increment and position  I
*                         of first point in regular grid 
*             eks       - If the regular grid exeeds the original
*                         grid this parameter defines the extrapolation. I
*                         If eks is .true. a(1) and a(n) is used for
*                         regular nodes that are outside the original grid.
*                         Otherwise the values  bleft ( left of xs(1)) and
*                         bright are used 
*             bleft, bright - see above                                  I
*             eps  - tolerance for defining a point inside original grid   I
*             iflag - error parameter                                      O
***************************************************************************
      subroutine nureglin(xs,ys,na,yb,nb,dxb,xb,eks,bleft,bright,
     %                   eps,iflag)
      integer na,nb,iflag
      real xs(na),ys(na),yb(nb),xb,dxb,bleft,bright,eps
      logical eks
*------------------------------------------------------------------
      integer k,k1,k2,k1m,ipos
      real xver,x1,x2,lval,rval,xm,xp,wg

 
      iflag=0

      x1=xs(1)
      x2=xs(na)


      if(na.eq.1) then
        iflag=2
        return
      end if


      if(eks) then
       rval=ys(na)
       lval=ys(1)
      else
       rval=bright
       lval=bleft
      end if

      ipos=1

      do 100 k=1,nb
        xver=xb+(k-1)*dxb
 60     continue
        if(xver.le.xs(ipos+1) .or. ipos.eq.(na-1)) go to 80
        ipos=ipos+1
        go to 60
 80     continue

        xp=xs(ipos+1)
        xm=xs(ipos)
        if(xp.le.xm) then
          iflag=3
          return
        end if

        if(xver.lt.x1-eps) then
           yb(k)=lval
        else
         if(xver.gt.x2+eps) then
           yb(k)=rval
         else
           wg=(xver-xm)/(xp-xm)
       
           yb(k)=(1.0-wg)*ys(ipos)+wg*ys(ipos+1)
         end if
        end if

 100  continue


      
      return

      end 


******************************************************************
*
*                   FGUREGLIN
*
*     A one-D dataset given on a non-uniform grid is interpolated onto
*     another non-uniform grid by linear interpolation. Both grids must 
*     contain coordinates in increasing order.
*     ERRORS:      do not work for nb=1
*     parameters:
*             xs, ys -  arrays for original grid and data               I
*                       The xs must increase monotonously
*             na - number of points in original grid
*             xb - array for interpolated grid                            I
*             yb - array for interpolated values                         O
*             nb  - number of interpolated points                        I
*                         of first point in regular grid 
*             ieks       - If the regular grid exeeds the original
*                         grid this parameter defines the extrapolation. I
*                         ieks=0 : values of yb are left unchanged
*                         ieks=1 :  a(1) and a(n) is used for
*                            regular nodes that are outside the original grid.
*                         ieks >=2 : the values  bleft ( left of xs(1)) and
*                              bright are used 
*                         If eks is .true. a(1) and a(n) is used for
*                         regular nodes that are outside the original grid.
*                         Otherwise the values  bleft ( left of xs(1)) and
*                         bright are used 
*             bleft, bright - see above                                  I
*             eps  - tolerance for defining a point inside original grid   I
*             iflag - error parameter                O
***************************************************************************
      subroutine fgureglin(xs,ys,na,xb,yb,nb,ieks,bleft,bright,
     %                   eps,iflag)
      integer na,nb,ieks,iflag
      real xs(na),ys(na),yb(nb),xb(nb),bleft,bright,eps
*------------------------------------------------------------------
      integer k,k1,k2,k1m,ipos
      real xver,x1,x2,lval,rval,xm,xp,wg

 
      iflag=0

      x1=xs(1)
      x2=xs(na)


      if(na.eq.1) then
        iflag=2
        return
      end if


      if(ieks.eq.1) then
       rval=ys(na)
       lval=ys(1)
      else
       rval=bright
       lval=bleft
      end if

      ipos=1

      do 100 k=1,nb
        xver=xb(k)
 60     continue
        if(xver.le.xs(ipos+1) .or. ipos.eq.(na-1)) go to 80
        ipos=ipos+1
        go to 60
 80     continue

        xp=xs(ipos+1)
        xm=xs(ipos)
        if(xp.le.xm) then
          iflag=3
          return
        end if

        if(xver.lt.x1-eps) then
         if(ieks.gt.0)yb(k)=lval
        else
         if(xver.gt.x2+eps) then
           if(ieks.gt.0)yb(k)=rval
         else
           wg=(xver-xm)/(xp-xm)
       
           yb(k)=(1.0-wg)*ys(ipos)+wg*ys(ipos+1)
         end if
        end if

 100  continue


      
      return

      end 


******************************************************************
*
*                   punktUREG
*
*     A one-D dataset given on a non-uniform grid is interpolated onto
*     a single point by linear interpolation. 
*     parameters:
*             xs, ys -  arrays for original grid and data               I
*                       The xs must increase monotonously
*             na - number of points in original grid
*             xb - interpolation point                            I
*             yb - interpolated values                         O
*             eps  - tolerance for defining a point inside original grid   I
*             iflag - error parameter                O
***************************************************************************
      subroutine punktureg(xs,ys,na,xb,yb,eps,iflag)
      integer na,iflag
      real xs(na),ys(na),yb,xb,eps
*------------------------------------------------------------------
      integer k,k1,k2,k1m,ipos
      real x1,x2,xm,xp,wg

 
      iflag=0

      x1=xs(1)
      x2=xs(na)
  

      if(na.eq.1) then
        iflag=2
        return
      end if

      if(xb.lt.x1) then
        if(xb+eps.ge.x1) then
          yb=ys(1)
        else
          iflag=1
        end if

        return
      end if

      if(xb.gt.x2) then
        if(xb-eps.le.x2) then
          yb=ys(na)
        else
          iflag=1
        end if

        return
      end if



      ipos=1

 60   continue
        if(xb.le.xs(ipos+1) .or. ipos.eq.(na-1)) go to 80
        ipos=ipos+1
      go to 60

 80   continue

      xp=xs(ipos+1)
      xm=xs(ipos)
      if(xp.le.xm) then
        iflag=3
        return
      end if

      wg=(xb-xm)/(xp-xm) 
      yb=(1.0-wg)*ys(ipos)+wg*ys(ipos+1)
      
      return

      end 

******************************************************************
*
*                   HREGuLIN
*
*     A one-D dataset given on a uniform grid is interpolated onto
*     another non-uniform grid by linear interpolation.  
*     xb must contain coordinates in increasing order.
*     parameters:
*             ys -  array for original  data                              I
*             na - number of points in original grid
*             dx - increment in source grid                               I
*             x0 - position of ys(1)
*             xb - array for interpolated grid                            I
*             yb - array for interpolated values                         O
*             nb  - number of interpolated points                        I
*                         of first point in regular grid 
*             ieks       - If the regular grid exeeds the original
*                         grid this parameter defines the extrapolation. I
*                         ieks=0 : values of yb are left unchanged
*                         ieks=1 :  a(1) and a(n) is used for
*                            regular nodes that are outside the original grid.
*                         ieks >=2 : the values  bleft ( left of xs(1)) and
*                              bright are used 
*             bleft, bright - see above                                  I
*             eps  - tolerance for defining a point inside original grid   I
*             iflag - error parameter                O
***************************************************************************
      subroutine hregulin(ys,na,dx,x0,xb,yb,nb,ieks,bleft,bright,
     %                   eps,iflag)
      integer na,nb,ieks,iflag
      real ys(na),dx,x0,yb(nb),xb(nb),bleft,bright,eps
*------------------------------------------------------------------
      integer k,k1,k2,k1m,ipos
      real xver,x1,x2,lval,rval,xm,xp,wg

 
      iflag=0

      x1=x0
      x2=x0+(na-1)*dx


      if(na.eq.1) then
        iflag=2
        return
      end if


      if(ieks.eq.1) then
       rval=ys(na)
       lval=ys(1)
      else
       rval=bright
       lval=bleft
      end if

      ipos=1

      do 100 k=1,nb
        xver=xb(k)
 60     continue
        if(xver.le.x0+ipos*dx .or. ipos.eq.(na-1)) go to 80
        ipos=ipos+1
        go to 60
 80     continue

        xp=x0+ipos*dx
        xm=xp-dx
 

        if(xver.lt.x1-eps) then
         if(ieks.gt.0)yb(k)=lval
        else
         if(xver.gt.x2+eps) then
           if(ieks.gt.0)yb(k)=rval
         else
           wg=(xver-xm)/(xp-xm)
       
           yb(k)=(1.0-wg)*ys(ipos)+wg*ys(ipos+1)
         end if
        end if

 100  continue


      
      return

      end 

******************************************************************
*
*                   regran
*
*     A one-D dataset given on a uniform grid is interpolated onto
*     another non-uniform grid by linear interpolation. Unlike hregulin xb
*     may contain coordinates in non-monotonous order.
*     parameters:
*             ys -  array for original  data                              I
*             na - number of points in original grid
*             dx - increment in source grid                               I
*             x0 - position of ys(1)
*             xb - array for interpolated grid                            I
*             yb - array for interpolated values                         O
*             nb  - number of interpolated points                        I
*                         of first point in regular grid 
*             ieks       - If the regular grid exeeds the original
*                         grid this parameter defines the extrapolation. I
*                         ieks=0 : values of yb are left unchanged
*                         ieks=1 :  a(1) and a(n) is used for
*                            regular nodes that are outside the original grid.
*                         ieks >=2 : the values  bleft ( left of xs(1)) and
*                              bright are used 
*             bleft, bright - see above                                  I
*             eps  - tolerance for defining a point inside original grid   I
*             iflag - error parameter                O
***************************************************************************
      subroutine regran(ys,na,dx,x0,xb,yb,nb,ieks,bleft,bright,
     %                   eps,iflag)
      integer na,nb,ieks,iflag
      real ys(na),dx,x0,yb(nb),xb(nb),bleft,bright,eps
*------------------------------------------------------------------
      integer k,k1,k2,k1m,ipos
      real xver,x1,x2,lval,rval,xm,xp,wg

 
      iflag=0

      x1=x0
      x2=x0+(na-1)*dx


      if(na.eq.1) then
        iflag=2
        return
      end if


      if(ieks.eq.1) then
       rval=ys(na)
       lval=ys(1)
      else
       rval=bright
       lval=bleft
      end if


      do 100 k=1,nb
        ipos=1
        xver=xb(k)
 60     continue
        if(xver.le.x0+ipos*dx .or. ipos.eq.(na-1)) go to 80
        ipos=ipos+1
        go to 60
 80     continue

        xp=x0+ipos*dx
        xm=xp-dx
 

        if(xver.lt.x1-eps) then
         if(ieks.gt.0)yb(k)=lval
        else
         if(xver.gt.x2+eps) then
           if(ieks.gt.0)yb(k)=rval
         else
           wg=(xver-xm)/(xp-xm)
       
           yb(k)=(1.0-wg)*ys(ipos)+wg*ys(ipos+1)
         end if
        end if

 100  continue


      
      return

      end 

******************************************************************
*
*                   FNUREGLIN
*
*     A one-D dataset given on a non-uniform grid is interpolated onto
*     a uniform grid by linear interpolation. The source grid must 
*     contain coordinates in increasing order.
*     parameters:
*             xs, ys -  arrays for original grid and data               I
*                       The xs must increase monotonously
*             na - number of points in original grid
*             yb - array for interpolated values on regular grid          O
*             nb,dxb,xb - number of points, grid increment and position  I
*                         of first point in regular grid 
*             ieks       - If the regular grid exeeds the original
*                         grid this parameter defines the extrapolation. I
*                         ieks=0 : values of yb are left unchanged
*                         ieks=1 :  a(1) and a(n) is used for
*                            regular nodes that are outside the original grid.
*                         ieks >=2 : the values  bleft ( left of xs(1)) and
*                              bright are used 
*             bleft, bright - see above                                  I
*             eps  - tolerance for defining a point inside original grid   I
*             iflag - error parameter                                      O
***************************************************************************
      subroutine fnureglin(xs,ys,na,yb,nb,dxb,xb,ieks,bleft,bright,
     %                   eps,iflag)
      integer na,nb,ieks,iflag
      real xs(na),ys(na),yb(nb),xb,dxb,bleft,bright,eps
*------------------------------------------------------------------
      integer k,k1,k2,k1m,ipos
      real xver,x1,x2,lval,rval,xm,xp,wg

 
      iflag=0

      x1=xs(1)
      x2=xs(na)


      if(na.eq.1) then
        iflag=2
        return
      end if


      if(ieks.eq.1) then
       rval=ys(na)
       lval=ys(1)
      else
       rval=bright
       lval=bleft
      end if

      ipos=1

      do 100 k=1,nb
        xver=xb+(k-1)*dxb
 60     continue
        if(xver.le.xs(ipos+1) .or. ipos.eq.(na-1)) go to 80
        ipos=ipos+1
        go to 60
 80     continue

        xp=xs(ipos+1)
        xm=xs(ipos)
        if(xp.le.xm) then
          iflag=3
          return
        end if

        if(xver.lt.x1-eps) then
           if(ieks.gt.0)yb(k)=lval
        else
         if(xver.gt.x2+eps) then
           if(ieks.gt.0)yb(k)=rval
         else
           wg=(xver-xm)/(xp-xm)
       
           yb(k)=(1.0-wg)*ys(ipos)+wg*ys(ipos+1)
         end if
        end if

 100  continue


      
      return

      end 

********************************************************************
      subroutine lestwoc(itape,x,y,nmax,n,iflag)
      integer itape,nmax,n,iflag
      real x(nmax),y(nmax)
*------------------------------------------------------------------
      integer ierr,nskip
      real a,b
 
      iflag=0
        
      call skipkom(nskip,itape,ierr)
      if(ierr.gt.0) then
      write(0,*)'skipkom: nskip,ierr=',nskip,ierr
         iflag=2
         close(itape)
         return
      end if

      n=0

 100  continue

      read(itape,*,end=200)a,b
      n=n+1
      if(n.gt.nmax) then
        iflag=5
        close(itape)
        return
      end if
      x(n)=a
      y(n)=b

      go to 100

 200  continue

      return


      end

**********************************************************************
*
*     Checks if a grid is uniform
*     x -grid      I
*     n - no of points
*     dmin,dmax - min and max increment
*     eps - rel tolerence for uniformity
*     ir - =0 uniform, 1 nonuniform                              O 
************************************************************************
      subroutine unifgrid(x,n,dmin,dmax,eps,ir)
      integer n,ir
      real x(n),dmin,dmax,eps
*----------------------------------------------------
      integer i
      real dx

      if(n.eq.1) then
        ir=0
         return
      end if

      dmin=x(2)-x(1)
      dmax=dmin

      do 100 i=3,n
       dx=x(i)-x(i-1)
       if(dx.gt.dmax)dmax=dx      
       if(dx.lt.dmin)dmin=dx      
 100  continue

      if( (dmax-dmin).gt.eps*dmin ) then
        ir=1
      else
        ir=0
      end if

      return
      end


********************************************************************
*
*                REGPFIL
*                
*     Prompts for name and reads data from two-column file, then 
*     interpolates to a uniform grid
*
*     spor   - question                                                I
*     fnavn,defn  - Name on file and def value                         O
*     kl  -   lengtn of filename. If kl<=0 the name is assumed to be   O
*             ended by '#' or '!' and the length is calculated
*     x,y  -  arrays for data from file                                O
*     nmax  - maximum number of data                                   I
*     n  -    number of data read  from file                           O
*     nskip  -  number of leading comment lines in file                O
*     val - array for regular ouput-data                               O
*     x0,dx,nv - left value, step and num,ber of points in val          I
*     spline - if true then splines are used, when x is uniform         I
*     eks       - If the regular grid exeeds the original
*                         grid this parameter defines the extrapolation. I
*                         If eks is .true. a(1) and a(n) is used for
*                         regular nodes that are outside the original grid.
*                         Otherwise the values  bleft ( left of xs(1)) and
*                         bright are used 
*     bleft, bright - see above                                  I
*     eps  - tolerance for defining a point inside original grid   I
*     iflag  - error flag                                              O
*                0     : OK
*                1     : maximum number of data reached before end of file.
*                        The data that is read is returned in x and y
*               11     : file not opened
*                3     : error in comment processing ( no data ?)
*                30    : inappropriate filename 
*                50    : cannot open file
*     w - working array
********************************************************************
      subroutine regpfil(spor,fnavn,defn,kl,x,y,nmax,n,nskip,
     % val,x0,nv,dx,spline,eps,eks,bleft,bright,iflag,w)
      integer kl,nmax,n,nskip,nv,iflag
      real x(nmax),y(nmax),val(nmax),x0,dx,bleft,bright,eps
      real w(5*nmax)
      logical eks,spline
      character*80 spor,fnavn,defn
*------------------------------------------------------------------
      integer ir,irv,ierr
      real dmax,dmin,dxi
      logical spl

      call promptfil(spor,fnavn,defn,kl,x,y,nmax,n,nskip,iflag)
      if(iflag.gt.1) then
        return
      end if
 
      call unifgrid(x,n,dmin,dmax,eps,ir)
      spl=spline.and.ir.eq.0
      if(n.gt.1) then
       dxi=(x(n)-x(1))/(n-1)
      else
       dxi=1.0
      end if

      if(eks) then
       irv=0
      else
       irv=1
      end if

      if(spl) then 
        call gstaggeg(y,n,dxi,x(1),val,nv,x0,dx,irv,spline,
     %           bleft,bright,.false.,0.0,.false.,0.0,ierr,eps,w)
      else
         call ureglin(x,y,n,val,nv,dx,x0,eks,bleft,bright,
     %                   eps,ierr)
      end if

      if(ierr.gt.0) then
       iflag=100+ierr
      end if



      return
      end


********************************************************************
*
*                REGPAUT
*                
*     Prompts for name and reads data from two-column file, then 
*     interpolates to a uniform grid; if necessary
*
*     spor   - question                                                I
*     fnavn,defn  - Name on file and def value                         O
*     kl  -   lengtn of filename. If kl<=0 the name is assumed to be   O
*             ended by '#' or '!' and the length is calculated
*     x,y  -  arrays for data from file                                O
*     nmax  - maximum number of data                                   I
*     n  -    number of data read  from file                           O
*     nskip  -  number of leading comment lines in file                O
*     val - array for regular ouput-data                               O
*     x0,dx - left value and  step  in val                              O
*     bleft, bright - see above                                  
*     iflag  - error flag                                              O
*                0     : OK
*                1     : maximum number of data reached before end of file.
*                        The data that is read is returned in x and y
*               11     : file not opened
*                3     : error in comment processing ( no data ?)
*                30    : inappropriate filename 
*                50    : cannot open file
*     w - working array
********************************************************************
      subroutine regpaut(spor,fnavn,defn,kl,x,y,nmax,n,
     % val,x0,dx,iflag,w)
      integer kl,nmax,n,iflag
      real x(nmax),y(nmax),val(nmax),x0,dx
      real w(5*nmax)
      character*80 spor,fnavn,defn
*------------------------------------------------------------------
      integer ierr,nskip,ir
      real dmax,dmin,dxi,eps,bleft,bright
      logical eks,spline,spl

      eks=.true.
      spline=.false.
      call promptfil(spor,fnavn,defn,kl,x,y,nmax,n,nskip,iflag)
      if(iflag.gt.1) then
        return
      end if
 
      call unifgrid(x,n,dmin,dmax,eps,ir)
      spl=spline.and.ir.eq.0
      if(n.gt.1) then
       dx=(x(n)-x(1))/(n-1)
      else
       dx=1.0
      end if
      x0=x(1)
      eps=0.0001*dx


      call ureglin(x,y,n,val,n,dx,x0,eks,bleft,bright,
     %                   eps,ierr)

      if(ierr.gt.0) then
       iflag=100+ierr
      end if



      return
      end




********************************************************************
*               I N G R E R
*
*     Rutinen integrer en spline interpolant for et ekvidistant punktsett.
*   parametere:
*       y  -  y(j) inneholder integralet av f's spline interpolant.      O
*             Integrasjonskonstanten er tilpasset slik at y(inull)=0.0
*             Dette betyr at y(j)= integral av G(x) fra x=inull*dx til
*             x=j*dx, der G(x) er interpolanten. Merk  at j=1 er 
*             identifisert med x=dx og ikke x=0.
*       f  -  Punktsett som skal integreres                              I
*       n  -  antall punkter                                             I
*       dx -  avstand mellom punkter                                     I
*       der1,der2 - deriverte av interpolant i hhv x=dx og x=n*dx        I
*       inull - Definerer integrasjonskonstant som angitt ovenfor        I
*       wk    -   Arbeidsomr}de, minimum st\rrelse er angitt i rutinehode.
*
*****************************************************************************
      subroutine ingrer(y,f,n,dx,der1,der2,inull,wk)      
      integer n,inull
      real y(n),f(n),wk(5*n+2),dx,der1,der2
*------------------------------------------------
      integer i
      real  r1,r2,diff

      call bsplw(f,wk(1),n,dx,der1,der2,wk(n+3))

      r1=dx/24.0
      r2=dx*11.0/24.0

      y(1)=0.0

      do 100 i=2,n
      y(i)=y(i-1)+r1*(wk(i-1)+wk(i+2))+r2*(wk(i)+wk(i+1))
 100  continue

      diff=y(inull)

      do 200 i=1,n
      y(i)=y(i)-diff
 200  continue

      return
      end 



********************************************************************
*               G I N G R E R
*
*     Rutinen integrer en spline interpolant for et ekvidistant punktsett.
*   parametere:
*       y  -  y(j) inneholder integralet av f's spline interpolant.      O
*             Integrasjonskonstanten er tilpasset slik at y(inull)=0.0
*             Dette betyr at y(j)= integral av G(x) fra x=inull*dx til
*             x=j*dx, der G(x) er interpolanten. Merk  at j=1 er 
*             identifisert med x=dx og ikke x=0.
*       f  -  Punktsett som skal integreres                              I
*       n  -  antall punkter                                             I
*       dx -  avstand mellom punkter                                     I
*       der1,der2 - deriverte av interpolant i hhv x=dx og x=n*dx        I
*      ik1,ik2 type of boundary condition for left and right boundary   I
*         respectively.
*         value=0 : derivative adapted to value d1,d2
*         value=1 : one sided differences used for derivative
*         value=2 : second derivative set to d1 and d2.
*                   if hv(0)=0 and d1=0 we then obtain an antisymmetric
*                   condition at left boundary etc. 
*       inull - Definerer integrasjonskonstant som angitt ovenfor        I
*       wk    -   Arbeidsomr}de, minimum st\rrelse er angitt i rutinehode.
*
*****************************************************************************
      subroutine gingrer(y,f,n,dx,ik1,der1,ik2,der2,inull,wk)      
      integer n,ik1,ik2,inull
      real y(n),f(n),wk(5*n+2),dx,der1,der2
*------------------------------------------------
      integer i
      real  r1,r2,diff

      call bsplgen(f,wk(1),n,dx,ik1,der1,ik2,der2,wk(n+3))

      r1=dx/24.0
      r2=dx*11.0/24.0

      y(1)=0.0

      do 100 i=2,n
      y(i)=y(i-1)+r1*(wk(i-1)+wk(i+2))+r2*(wk(i)+wk(i+1))
 100  continue

      diff=y(inull)

      do 200 i=1,n
      y(i)=y(i)-diff
 200  continue

      return
      end 

*********************************************************************
      subroutine trapes(f,n,dx,sum)
      integer n
      real f(n),dx,sum
*--------------------------------------------------------
      integer i,nm

      sum=0.5*(f(1)+f(n))
      nm=n-1
      do 100 i=2,nm
       sum=sum+f(i)
 100  continue

      sum=sum*dx
      return
      end

************************************************************
      function midpu(a,b,n,f)
      integer n
      real a,b,f,midpu
      external f
*-------------------------------------------------------------
      integer i
      real dx,x,sum

      sum=0.0
      dx=(b-a)/n

      do 100 i=1,n
       x=(i-0.5)*dx+a
       sum=sum+f(x)
 100  continue
      midpu=dx*sum

      return
      end


************************************************************
      function simps(a,b,n,f)
      integer n
      real a,b,f,simps
      external f
*-------------------------------------------------------------
      integer i,j,nm
      real dx,x,suml,sumo

      dx=(b-a)/(2.0*n)
      nm=n-1
      suml=0.0
      sumo=0.0

      do 100 i=1,n
       j=2*i-1
       x=a+j*dx
       sumo=sumo+f(x)
 100  continue
      
      do 200 i=1,nm
       j=2*i
       x=a+j*dx
       suml=suml+f(x)
 200  continue

      simps=dx*(f(a)+f(b)+4.0*sumo+2.0*suml)/3.0


      return
      end

************************************************************
      subroutine regitab(res,n,dx,x0,dstol,f,rmet)
      integer n
      real res(n),dx,x0,dstol,f,rmet
      external f,rmet
*-----------------------------------------------------------
      integer i,ni,nm
      real a,b

      nm=n-1
      ni=dx/dstol
      if(ni*dstol.lt.dx) ni=ni+1

      b=x0
      do 100 i=1,nm
        a=b
        b=b+dx
        res(i+1)=res(i)+rmet(a,b,ni,f)
 100  continue


      end 


************************************************************
      subroutine varitab(res,x,n,dstol,f,rmet)
      integer n
      real res(n),x(n),dstol,f,rmet
      external f,rmet
*-----------------------------------------------------------
      integer i,ni,nm
      real a,b,dx

      nm=n-1


      b=x(1)

      do 100 i=1,nm
        a=b
        b=x(i+1)
        dx=b-a
        ni=dx/dstol
        if(ni*dstol.lt.dx) ni=ni+1
        res(i+1)=res(i)+rmet(a,b,ni,f)
 100  continue


      end 

*********************************************************************
      subroutine gtrapes(y,f,n,dx)
      integer n
      real y(n),f(n),dx
*--------------------------------------------------------
      integer i
      real ym,yp,dx05

      dx05=dx*0.5
      
      yp=y(1)
      f(1)=0

      do 100 i=2,n
       ym=yp
       yp=y(i)
       f(i)=f(i-1)+dx05*(yp+ym)
 100  continue

      return
      end


*********************************************************************
      subroutine ktrapes(y,f,x,n)
      integer n
      real y(n),f(n),x(n)
*--------------------------------------------------------
      integer i
      real ym,yp,xm,xp

      yp=y(1)
      xp=x(1)
      f(1)=0

      do 100 i=2,n
       ym=yp
       yp=y(i)
       xm=xp
       xp=x(i)
       f(i)=f(i-1)+(xp-xm)*(yp+ym)
 100  continue

      return
      end


*********************************************************************
      subroutine intmat(y,n,m,dx,dy,vol)
      integer n,m
      real y(n,m),dx,dy,vol
*--------------------------------------------------------
      integer i,k,nm,mm
      real sum

      nm=n-1
      mm=m-1
      sum=0.25*(y(1,1)+y(n,1)+y(1,m)+y(n,m))

      do 200 k=2,mm
        sum=sum+0.5*(y(1,k)+y(n,k))
        do 100 i=2,nm
          sum=sum+y(i,k)
 100    continue
 200  continue

       do 300 i=2,nm
          sum=sum+0.5*(y(i,1)+y(i,m))
 300   continue

      sum=sum*dx*dy
      vol=sum

      return
      end




***********************************************************************
* Finds minum and maximum of a(n)
*
***********************************************************************
      subroutine ekstrem(a,n,amin,amax)
      integer n
      real a(n),amax,amin
*------------------------------------------------------
      integer i

      amax=a(1)
      amin=a(1)

      do 100 i=2,n
      if(a(i).gt.amax) amax=a(i)
      if(a(i).lt.amin) amin=a(i)
 100  continue

      return
      end

***********************************************************************
* Finds minum and maximum of a(n)
*
***********************************************************************
      subroutine gekstrem(a,n,imin,amin,imax,amax)
      integer n,imin,imax
      real a(n),amax,amin
*------------------------------------------------------
      integer i

      amax=a(1)
      amin=a(1)
      imin=1
      imax=1
      do 100 i=2,n
      if(a(i).gt.amax) then
         amax=a(i)
         imax=i
      end if

      if(a(i).lt.amin) then
        amin=a(i)
        imin=i
      end if
 100  continue

      return
      end


*************************************************************
      subroutine recol(v,nc,nmax,n,itape,iflag)
      integer nc,nmax,n,itape,iflag
      real v(nmax,nc)
*-----------------------------------------------------------
      integer j
      real buff(50)

      if(nc.gt.50) then
        iflag=2
        return
      end if

      iflag=0
      n=0

 100  continue
      n=n+1
      read(itape,*,end=200) (buff(j), j=1,nc)
      if(n.gt.nmax) then
        iflag=1
        n=n-1
        return
      end if

      do 120 j=1,nc
      v(n,j)=buff(j)
 120  continue

      go to 100

 200  continue
      n=n-1
      return
      

      end

************************************************************************
      function av(y,n,na)
      integer n,na
      real y(n),av
*---------------------------------------------------------------------
      integer nn,i
      real sum

      nn=min(n,na)
      sum=0.0

      do 100 i=1,nn
      sum=sum+y(i)
 100  continue

      av=sum/nn

      return
      end


************************************************************************
      SUBROUTINE DEL(A,B,N,M,NB,MB,N0,M0,N1,M1,IX,IY)
      INTEGER N,M,NB,MB,N0,M0,N1,M1,IX,IY
      REAL A(n*m),B(n*m)
*----------------------------------------------------------
      INTEGER JSTART,KST,K,I

      NB=(N1-N0)/IX +1
      MB=(M1-M0)/IY +1
      JSTART=(M0-1)*N+N0

      DO 100 K=1,MB
      KST=(K-1)*NB
      DO 50 I=1,NB
   50 B(KST+I)=A(JSTART+(I-1)*IX)
      JSTART=JSTART+IY*N
  100 CONTINUE
      RETURN
      END




*******************************************************************
*
*              PART
*
*         Klipper ut en del av en array
*      parametere:
*             U - matrise det skal velges fra                          I
*             Y - resultatmatrise                                      O
*             NU,MU,DXU,DYU,XU,YU - gitterpunkter, avstander           I
*                      og startpunkt for U
*             NY,MY,DXY,DYY,X1,Y1 - gitterpunkter, avstander           O
*                      og startpunkt for Y
***************************************************************************
      SUBROUTINE PART(U,Y,NU,MU,DXU,DYU,XU,YU,NY,MY,DXY,DYY,X1,Y1)
      INTEGER NU,MU,NY,MY
      REAL U(1),Y(1),XU,YU,X1,Y1,DXU,DYU,DXY,DYY
*----------------------------------------------------------------
      INTEGER N0,M0,N1,M1,IX,IY
      CHARACTER C

  50  CONTINUE
      WRITE(0,*)'N=',NU,'    M=',MU
      call lesi4(1,1,nu,mu,'GI N. VEN. OG O. HO. HJOERNE#',N0,M0,N1,M1)

      IF(N0.Gt.N1 .OR. M0.Gt.M1) THEN
        WRITE(0,*)'ULOVLIGE DATA:',n0,m0,n1,m1
        GO TO 50
      END IF

      IF(N0.LE.0) N0=1
      IF(N1.GT.NU) N1=NU
      IF(M0.LE.0) M0=1
      IF(M1.GT.MU) M1=MU

      WRITE(0,*)'N0,M0=',N0,' , ',M0,'   N1,M1=',N1,' , ',M1
      call lesi2(1,1,' GI TETTHET I X OG Y-RETNING#',ix,iy)
      IF(IX.LT.1) IX=1
      IF(IY.LT.1) IY=1

      CALL DEL(U,Y,NU,MU,NY,MY,N0,M0,N1,M1,IX,IY)
      X1=XU+(N0-1)*DXU
      Y1=YU+(M0-1)*DYU
      DXY=IX*DXU
      DYY=IY*DYU

      RETURN
      END


*******************************************************************
*
*              DPART
*
*         Leser klipp-data for en array
*      parametere:
*             NU,MU - gitterpunkter, avstander           I
*                      og startpunkt for oppr. matrise
*             NY,MY,IX,IY,NO,MO - gitterpunkter, avstander           O
*                      og startpunkt for ny matrise
***************************************************************************
      SUBROUTINE DPART(NU,MU,NY,MY,IX,IY,N0,M0)
      INTEGER NU,MU,NY,MY,N0,M0,IX,IY
*----------------------------------------------------------------
      INTEGER N1,M1
      CHARACTER C

  50  CONTINUE
      WRITE(0,*)'N=',NU,'    M=',MU
      call lesi4(1,1,nu,mu,'GI N. VEN. OG O. HO. HJOERNE#',N0,M0,N1,M1)

      IF(N0.GT.N1 .OR. M0.GT.M1) THEN
        WRITE(0,*)'ULOVLIGE DATA:',n0,m0,n1,m1
        GO TO 50
      END IF

      IF(N0.LE.0) N0=1
      IF(N1.GT.NU) N1=NU
      IF(M0.LE.0) M0=1
      IF(M1.GT.MU) M1=MU

      WRITE(0,*)'N0,M0=',N0,' , ',M0,'   N1,M1=',N1,' , ',M1
      call lesi2(1,1,' GI TETTHET I X OG Y-RETNING#',ix,iy)
      IF(IX.LT.1) IX=1
      IF(IY.LT.1) IY=1

      ny= 1+(n1-n0)/ix
      my= 1 +(m1-m0)/iy


      RETURN
      END

          






*************************************************************************
*
*                 TRANPK
*
*       speiler og transponerer i denne rekkefoelge.
*      parametere:
*          Y - datarray     ved input en tett nxm, ved output      I/O
*                     tett mxn dersom transponering
*          n,m - arraystorrelsere ( endres ikke pg trans)          I
*          xspeil,yspeil,transp - angir operasjoner                I
*          W - arbeidsarray av samme storrelse som Y               D
*
************************************************************************
      SUBROUTINE TRANPK(Y,N,M,xspeil,yspeil,transp,W)
      INTEGER N,M
      REAL Y(N*M),W(N*M)
      logical xspeil,yspeil,transp
*---------------------------------------------
      REAL TEMP
      INTEGER I,K,NGM,NTEMP,MH,NH

      NGM=N*M

      IF(xspeil)THEN
        MH=M/2
        DO 200 I=1,N
        DO 200 K=1,MH
        TEMP=Y(I+(K-1)*N)
        Y(I+(K-1)*N)=Y(I+(M-K)*N)
        Y(I+(M-K)*N)=TEMP
  200   CONTINUE
      end if

      IF(yspeil)THEN
        NH=N/2
        DO 300 K=1,M
        DO 300 I=1,NH
        TEMP=Y(I+(K-1)*N)
        Y(I+(K-1)*N)=Y(N+1-I +(K-1)*N)
        Y(N+1-I+(K-1)*N)=TEMP
  300   CONTINUE
      end if 

      if(transp) then
        DO 500 K=1,M
        DO 500 I=1,N
         w(K+(I-1)*M)=Y(I+(K-1)*N)
  500   CONTINUE
        DO 600 I=1,NGM
  600   Y(I)=W(I)
      end if

      END

********************************************************************
*          SMGEN
*
*     Performs 3 or 5 points smoothing of dataset
*
*     parameters
*       f   --  array with field tobe smoothened           I/O
*       n   --  number of points                           I
*       iord -- value =1: 3 point, else 5 point             I
*       sgn1, sgn2  -- govern boundary treatment                    I
*                    = -1 antisymmetric boundary assumed
*                    = 0  one sided smoothing formulas used
*                    = 1  symmetric boundary assumed   
*       bet  -- relaxation parameter  (value=1 full smoothing)      I
*       ierr -- error parameter (0: OK, 1: to few points)           O
***************************************************************************
      subroutine smgen(f,n,iord,sgn1,sgn2,bet,ierr)
      integer n,iord,sgn1,sgn2,ierr
      real f(n),bet
*--------------------------------------------------------------------

      if(iord.eq.1)  call smgen3(f,n,sgn1,sgn2,bet,ierr)
      if(iord.eq.2) call  smgen5(f,n,sgn1,sgn2,bet,ierr)
    
      return
      end

      

********************************************************************
*          SMGEN5
*
*     Performs 5 points smoothing of dataset
*
*     parameters
*       f   --  array with field to be smoothened           I/O
*       n   --  number of points                           I
*       sgn1, sgn2  -- govern boundary treatment                    I
*                    = -1 antisymmetric boundary assumed
*                    = 0  one sided smoothing formulas used
*                    = 1  symmetric boundary assumed   
*                    = 2  symmetric boundary at intermediate grid-cite 
*       bet  -- relaxation parameter  (value=1 full smoothing)      I
*       ierr -- error parameter (0: OK, 1: too few points)           O
***************************************************************************
      subroutine smgen5(f,n,sgn1,sgn2,bet,ierr)
      integer n,sgn1,sgn2,ierr
      real f(n),bet
*--------------------------------------------------------------------
      integer j,nm2
      real a0,ap,app,f1,f2,fn,fnm,fmm,fm,f0,fp,fpp

      if(n.lt.5) then
        ierr=1
        return
      else
        ierr=0
      end if

      nm2=n-2


      app=-1.0/16.0
      ap=0.25
      a0=5.0/8.0




      if(sgn2.eq.0) then
        fn=15.0*f(n)/16.0+0.25*f(n-1)-3.0*f(n-2)/8.0
     %     +0.25*f(n-3)-f(n-4)/16.0
        fnm=f(n)/16.0+0.75*f(n-1)+3.0*f(n-2)/8.0
     %     -0.25*f(n-3)+f(n-4)/16.0
      else
        if(sgn2.eq.2) then
          fp=f(n)
          fpp=f(n-1)
        else
          fp=sgn2*f(n-1)+(1-sgn2)*f(n)
          fpp=sgn2*f(n-2)+(1-sgn2)*f(n)
        end if
        fn=a0*f(n)+ap*(f(n-1)+fp)+app*(f(n-2)+fpp)
        fnm=a0*f(n-1)+ap*(f(n-2)+f(n))+app*(f(n-3)+fp)
      end if
 
     
      if(sgn1.eq.0) then
        f1=15.0*f(1)/16.0+0.25*f(2)-3.0*f(3)/8.0+
     %     0.25*f(4)-f(5)/16.0
        f2=f(1)/16.0+0.75*f(2)+3.0*f(3)/8.0
     %     -0.25*f(4)+f(5)/16.0
      else
        if(sgn1.eq.2) then
          fm=f(1)
          fmm=f(2)
        else
          fm=sgn1*f(2)+(1-sgn1)*f(1)
          fmm=sgn1*f(3)+(1-sgn1)*f(1)
        end if
        f1=a0*f(1)+ap*(f(2)+fm)+app*(f(3)+fmm)
        f2=a0*f(2)+ap*(f(1)+f(3))+app*(f(4)+fm)
      end if
 

      fm=f(1)
      f0=f(2)
      fp=f(3)
      fpp=f(4)

      do 200 j=3,nm2
        fmm=fm
        fm=f0
        f0=fp
        fp=fpp
        fpp=f(j+2)

        f(j)=(1.0-bet)*f(j)+bet*(a0*f0+ap*(fm+fp)+app*(fpp+fmm))
 200  continue

      f(1)=(1.0-bet)*f(1)+bet*f1
      f(2)=(1.0-bet)*f(2)+bet*f2
      f(n-1)=(1.0-bet)*f(n-1)+bet*fnm
      f(n)=(1.0-bet)*f(n)+bet*fn
      return
      end


********************************************************************
*          SMGEN3
*
*     Performs 3 points smoothing of dataset
*
*     parameters
*       f   --  array with field tobe smoothened           I/O
*       n   --  number of points                           I
*       sgn1, sgn2  -- govern boundary treatment                    I
*                    = -1 antisymmetric boundary assumed
*                    = 0  one sided smoothing formulas used
*                    = 1  symmetric boundary assumed   
*       bet  -- relaxation parameter  (value=1 full smoothing)      I
*       ierr -- error parameter (0: OK, 1: to few points)           O
***************************************************************************
      subroutine smgen3(f,n,sgn1,sgn2,bet,ierr)
      integer n,sgn1,sgn2,ierr
      real f(n),bet
*--------------------------------------------------------------------
      integer j,nm
      real a0,ap,f1,fn,fm,f0,fp

      if(n.lt.3) then
        ierr=1
        return
      else
        ierr=0
      end if

      nm=n-1



      ap=0.25
      a0=0.5

      if(sgn2.eq.0) then
        fn=0.75*f(n)+0.5*f(n-1)-0.25*f(n-2)
      else
        if(sgn2.eq.2) then
          fp=f(n)
        else
          fp=sgn2*f(n-1)+(1-sgn2)*f(n)
        end if
        fn=a0*f(n)+ap*(f(n-1)+fp)
      end if
 
     
      if(sgn1.eq.0) then
        f1=0.75*f(1)+0.5*f(2)-0.25*f(3)
      else
        if(sgn1.eq.2) then
          fm=f(1)
        else
          fm=sgn1*f(2)+(1-sgn1)*f(1)
        end if
        f1=a0*f(1)+ap*(f(2)+fm)
       end if
 

  
      f0=f(1)
      fp=f(2)

      do 200 j=2,nm

        fm=f0
        f0=fp
        fp=f(j+1)

        f(j)=(1.0-bet)*f(j)+bet*(a0*f0+ap*(fm+fp))
 200  continue

      f(1)=(1.0-bet)*f(1)+bet*f1
      f(n)=(1.0-bet)*f(n)+bet*fn
      return
      end




********************************************************************
*          NSMO
*
*     Performs variable order  smoothing of dataset
*
*     parameters
*       f   --  array with field to be smoothened           I/O
*     n   --  number of points                           I
*     iord - order of smoothing (2*iord+1 points used)   I
*            viable values 1, 2, 3, 4, 5      
*       sgn1, sgn2  -- govern boundary treatment                    I
*                    = -1 antisymmetric boundary assumed
*                    = 1  symmetric boundary assumed   
*                    = 2  symmetric boundary at intermediate grid-cite 
*       bet  -- relaxation parameter  (value=1 full smoothing)      I
*     ierr -- error parameter (0: OK, 1: too few points)           O
*     wrk -- working array                                         D      
***************************************************************************
      subroutine nsmo(f,n,iord,sgn1,sgn2,bet,ierr,wrk)
      integer n,iord,sgn1,sgn2,ierr
      real f(n),bet,wrk(1-iord:n+iord)
*--------------------------------------------------------------------
      integer j,sgnl,sgnh,k,loff,hoff
      real a(0:20),sum,ff

      if(n.lt.(2*iord+1)) then
        ierr=1
        return
      else
        ierr=0
      end if

cc Weights from Dold 1991      
      selectcase(iord)
      case(1)
         a(0)=0.5
         a(1)=0.25
      case(2)
         ff=1.0/16.0
         a(0)=1.0-6.0*ff
         a(1)=4.0*ff
         a(2)=-ff
      case(3)
         ff=1.0/64.0
         a(0)=1.0-20.0*ff
         a(1)=15.0*ff
         a(2)=-6.0*ff
         a(3)=ff
      case(4)
         ff=1.0/256.0
         a(0)=1.0-70.0*ff
         a(1)=56.0*ff
         a(2)=-28.0*ff
         a(3)=8.0*ff
         a(4)=-ff
      case(5)
         ff=1.0/1024.0
         a(0)=1.0-252.0*ff
         a(1)=210.0*ff
         a(2)=-120.0*ff
         a(3)=45.0*ff
         a(4)=-10.0*ff
         a(5)=ff
       end select

      
      do 100 j=1,n
        wrk(j)=f(j) 
 100  continue

      if(sgn1.ge.0) then
         sgnl=1
      else
         sgnl=-1
         if(sgn1.eq.-1) wrk(0)=0.0
      end if

      if(sgn1.eq.2.or.sgn1.eq.-2) then
         loff=-1
      else
         loff=0
      end if

      if(sgn2.ge.0) then
         sgnh=1
      else
         sgnh=-1
         if(sgn1.eq.-1) wrk(n)=0.0
      end if

      if(sgn2.eq.2.or.sgn2.eq.-2) then
         hoff=-1
      else
         hoff=0
      end if


      
      do 150 j=1,iord
         wrk(1-j)=sgnl*(f(1+j+loff) -f(1))+f(1)
         wrk(n+j)=sgnh*(f(n-j+hoff) -f(n))+f(n)
 150  continue
      
      


      do 200 j=1,n
         sum=a(0)*wrk(j)
         do 220 k=1,iord
            sum=sum+a(k)*(wrk(j-k)+wrk(j+k))
 220     continue
         f(j)=(1.0-bet)*wrk(j)+bet*sum
 200  continue

      return
      end


********************************************************************
*          GlatGEN
*
*     Performs 3 or 5 points smoothing of dataset
*
*     parameters
*       f   --  array with field tobe smoothened           I/O
*       n   --  number of points                           I
*     iord -- value =1: 3 points, possible one sided        I
*                   >=2  5 points, possible one sided
*                    <0   (1-2*iord) points only sym/antisym. at boundaries
*                        iord >=-5 is available
*     sgn1, sgn2  -- govern boundary treatment                    I
*                    = -2 antisymmetric half an interval outside end point      
*                    = -1 antisymmetric boundary assumed
*                    = 0  one sided smoothing formulas used (iord>0)
*                         otherwise same acton as value=1
*                    = 1  symmetric boundary assumed
*                    = 2 antisymmetric half an interval outside end point      
*       bet  -- relaxation parameter  (value=1 full smoothing)      I
*     ierr -- error parameter (0: OK, 1: to few points)           O
*     wrk - work array                                            D      
***************************************************************************
      subroutine glatgen(f,n,iord,sgn1,sgn2,bet,ierr,wrk)
      integer n,iord,sgn1,sgn2,ierr
      real f(n),bet,wrk(*)
*--------------------------------------------------------------------
      integer ieff

      if(iord.gt.0) then
         if(iord.eq.1)  then
            call smgen3(f,n,sgn1,sgn2,bet,ierr)
         else    
            call  smgen5(f,n,sgn1,sgn2,bet,ierr)
         end if
      else
         if(iord.lt.0) then
           ieff=-iord
           call  nsmo(f,n,ieff,sgn1,sgn2,bet,ierr,wrk)
         end if
        
      end if
      
    
      return
      end

      


**********************************************************************
      subroutine lagkopi(a,b,nz)
      integer nz
      real a(nz),b(nz)
*---------------------------------------------------------------------
      integer i

      do 100 i=1,nz
       b(i)=a(i)
 100   continue

      return
      end

**********************************************************************
      subroutine affine(x,nz,a,b)
      integer nz
      real x(nz),a,b
*---------------------------------------------------------------------
      integer i

      do 100 i=1,nz
       x(i)=a+b*x(i)
 100   continue

      return
      end


**********************************************************************
      subroutine nullinit(a,na,nz)
      integer na,nz
      real a(na:nz)
*---------------------------------------------------------------------
      integer i

      do 100 i=na,nz
       a(i)=0
 100   continue

      return
      end

**********************************************************************
      subroutine addab(a,b,c,na,nz,wa,wb)
      integer na,nz
      real a(na:nz),b(na:nz),c(na:nz),wa,wb
*---------------------------------------------------------------------
      integer i

      do 100 i=na,nz
       c(i)=wb*b(i)+wa*a(i)
 100   continue

      return
      end

**********************************************************************
      subroutine leggbtoa(a,b,na,nz,wa,wb)
      integer na,nz
      real a(na:nz),b(na:nz),wa,wb
*---------------------------------------------------------------------
      integer i

      do 100 i=na,nz
       a(i)=wb*b(i)+wa*a(i)
 100   continue

      return
      end


******************************************************************
*
*                   RINTLIN
*
*     A one-D dataset on a uniform grid is interpolated at a given point
*             ys -  array  for original    data               I
*             n,dx,xa - number of points, grid increment and position  I
*                         of first point for grid
*             xb -        position  to be interpolated                   I
*             yb -        interpolated value
*             eps  - tolerance for defining a point inside original grid   I
*             iflag - error parameter                                      O
*                      value 0 : OK, 1 xb outside grid                
***************************************************************************
      subroutine rintlin(ys,n,dx,xa,xb,yb,eps,iflag)
      integer n,iflag
      real ys(n),xa,dx,xb,yb,eps
*------------------------------------------------------------------
      integer ipos
      real r,w



 
      iflag=0

cc    find ipos such that xb is between point ipos and ipos+1
      
      ipos=(xb-xa)/dx +1

cc    Special case 1
      if(ipos.eq.0 .and.(xa-xb).gt.eps)ipos=1
cc    Special case 2
      if(ipos.eq.n .and.(xb-xa-(n-1)*dx).lt.eps)ipos=n-1

cc    r is sistance from point ipos
      r=xb-xa-(ipos-1)*dx

      w=r/dx
      yb=w*ys(ipos+1)+(1.0-w)*ys(ipos)

      
      return

      end 


***********************************************************************
      subroutine linspace(xa,xb,n,x)
      integer n
      real xa,xb,x(n)
*----------------------------------------------------------------------
      integer i
      real dx

      if(n.lt.0) return
      if(n.eq.1) then
        x(1)=xa
        return
      end if

      dx=(xb-xa)/(n-1)

      do 100 i=1,n
        x(i)=xa+(i-1.0)*dx
 100  continue
 
      return
      end


***********************************************************************
      subroutine lingitt(x0,dx,n,x)
      integer n
      real x0,dx,x(n)
*----------------------------------------------------------------------
      integer i
 

       do 100 i=1,n
        x(i)=x0+(i-1.0)*dx
 100  continue
 
      return
      end

************************************************************************
*
*     Checks wether r not agrid is uniform
*
*     x - grid array                                        I
*     n - number of points                                  I
*     dx - averaged increment                               O
*     tol - if increments vary less than tol*dx the grid    I
*           is acepted as uniform
*     unif - true if grid is uniform                        O
****************************************************************************
      subroutine klassdat(x,n,dx,tol,unif)
      integer n
      real x(n),dx,tol
      logical unif
*---------------------------------------------------------------
      integer i,imiss
      real dev,dev2,tol2
      if(n.eq.1) then
        dx=0
        unif=.true.
        return
      end if
      tol2=tol*tol*dx*dx
      dx=(x(n)-x(1))/(n-1.0)
      imiss=0

      do 100 i=2,n
       dev=x(i)-x(i-1)-dx
       dev2=dev*dev
       if(dev2.gt.tol2) imiss=imiss+1
 100  continue
      unif=imiss.eq.0
      return
      end

*******************************************************************
*
*                stgrid
*      copied from     laplace/rusource/grigen:linstr
*
*      Defines a map from (0,1) in ksi to (0,c) in x, such that
*
*            x(0)=0     x(ak)=a  x(bk)=b  x(1)=c
*     
*            dx/dksi continuous, constant on  (0,ak), (bk,1) and linear
*                    on (ak,bk)
*      parameters:
*               a,b,c - partition in x                               I
*               ak,bk - partition in ksi                             O
*               f     - ratio between dx/dksi in (bk,1) and (0,ak)   I
*                       defines degree  of local refinement degree      
*               n     - number of points in generated grid in ksi    I  
*               x   -array with ksi-grid                           O  
****************************************************
      subroutine stgrid(a,b,c,f,n,ak,bk,x)
      integer n
      real a,b,c,f,ak,bk,x(n)
*---------------------------------------------------------------------
      integer i,ia,ib
      real kl,kr,nu,dksi,vksi

      kl=2*(b-a)/(1+f)+a+(c-b)/f
      kr=f*kl
      ak=a/kl
      bk=1.0-(c-b)/kr
      nu=(kr-kl)/(bk-ak) 
      write(0,*)kl,kr,ak,bk,nu
      dksi=1.0/(n-1)

      ia=1+ak/dksi

      ib=1+bk/dksi

      do 100 i=1,ia
        x(i)=(i-1)*dksi*kl
 100  continue

      do 200 i=ia+1,ib
        vksi=(i-1)*dksi-ak
        x(i)=a+kl*vksi+0.5*nu*vksi*vksi
 200  continue

      do 300 i=ib+1,n
        vksi=1-(i-1)*dksi
        x(i)=c-kr*vksi
 300  continue

      return
      end
**********************************************************
*
*              Subrutine 'maskintype' kaller p} systemrutiner som
*      er spesielle for dec-stat og setter maskinidentitet i
*      cpuid
*************************************************************
      subroutine maskintype(cpuid)
      character*10 cpuid

      cpuid='decstation'


      return
      end


*******************************************************************
*                                                                 *
*              I R E G                                            *
*                                                                 *
* rett linje tilpassing slik at                                   *
*                                                                 *
*                   var=summasjon(xr(i)-a-b*i)**2                 *
*                       i=1,ip                                    *
* minimaliseres.                                                  *
*                                                                 *
*******************************************************************
      SUBROUTINE IREG(XR,IP,A,B,VAR)
      integer ip
      REAL XR(IP),A,B,VAR
*-------------------------------------------------------------------
      REAL SX,SA,SXX,SAA,SXA,DET
      INTEGER I

      IF(IP.LT.2) THEN
        VAR=-1.0
        RETURN
      END IF

      SX=0.0
      SA=0.0
      SXA=0.0
      SAA=0.0
      SXX=0.0

      DO 100 I=1,IP
      SX=SX+XR(I)
      SXX=SX+XR(I)*XR(I)
      SA=SA+I
      SAA=SAA+I*I
      SXA=SXA+I*XR(I)
  100 CONTINUE

      DET=IP*SAA-SA*SA
      A=(SX*SAA-SXA*SA)/DET
      B=(IP*SXA-SA*SX)/DET
      VAR=SXX+IP*A*A+SAA*B*B-2.0*(SX*A+SXA*B-SA*A*B)

      RETURN
      END





******************************************************************************
*
*
*                 R Y G G E R
*
*     Subroutinene gjenkjenner rygger innenfor omraadet n1<i<n2,m1<k<m2 i
*     arrayen eta(1:n,1:m). Parametere:
*
*       eta   - punktarray
*       n,m   - grenser for eta.
*       xmax,amax  - plassering og styrke av maksimaene i hver linje. 
*              dvs. xmax(1,j) > xmax(2,j)...> xmax( np(j),j)  er maksp.
*              for eta(n1,m1-1+j)...eta(n2,m1-1+j).
*          a  - hver rygg har likning x(y)=a(1,j)+ a(2,j)*y. Koordinatsystem
*               er valgt slik at eta(i,j) hoerer til x=i, y=j
*         ah  - amplitude ved starten av hver rygg.
*        ipt  - antall punkter i hver rygg 
*         np  - np(j) er antall maksima funnet for linje m1-1+j.
*         nr  - antall rygger med start i j=m1 (eller m2 dersom iov=-1)
*               som er funnet. 
*     n1,n2,m1,m2 - omraade ekstrmaene skal soekes i.
*      gtol   -  bare maksima med verdi stoerre enn gtol regnes med.
*      vitol  -  maksimal x forflytning av en rygg fra en linje til den
*                neste.
*      iov    -  verdi lik 1 gir de rygger med dx/dy>0 og start i m1.
*                verdi lik -1 gir rygger med dx/dy<0 og start i m2.
*
*******************************************************************************
      subroutine ryggdob(eta,n,m,xmax,amax,a,ah,ipt,np,nr,n1,n2,
     %m1,m2,gtol,vitol,iov)
      integer n,m,n1,n2,m1,m2,nr,iov
      integer ipt(20),np(m2-m1+1)
      real eta(n,m),xmax(20,m2-m1+1),amax(20,m2-m1+1)
      real a(2,20),ah(20),gtol,vitol,aa,bb,var
*------------------------------------------------------------------------------
      real bs(0:1005),amp(20),xa(20),XVERD(500),XV,AHVER,TEMP
      integer na,npun,mtall,kverd(500),irygg,J,nrygg
      integer istop,itemp,i,l,jp,mrad,irad,isgn,kadd

      npun=n2-n1+1                                 
      mrad=m2-m1+1
      nr=0                                                   
                                          
    
      mtall=mrad
      if(iov.eq.-1) then
        irad=m2+1                                   
        isgn=-1
      else
        irad=m1-1 
        isgn=1
      end if                                        

      kadd=irad

      do 100 i=1,mrad
      irad=irad+iov
      call limax(eta(n1,irad),bs,npun,gtol,amp,xa,na)
      np(i)=na                                
      if(na.eq.0) then 
        mtall=i-1
        go to 60
      end if
      do 50 l=1,na
      xmax(l,i)=xa(na+1-l)+n1 -1
      amax(l,i)=amp(na+1-L)
   50 continue
  100 continue



   60 continue

      if(mtall.lt.2) return
                        


      nrygg=np(1)

      do 200 j=1,nrygg
      xverd(1)=xmax(j,1)
      kverd(1)=j                                                   
      irygg=1
      jp=j           
      ahver=amax(j,1)
      do 150 i=2,mtall
      if(jp.gt.np(i)) jp=np(i)
      xv=xverd(i-1)
  110 continue
      if(xmax(jp,i).gt.(xv-0.1)) go to 120
      jp=jp-1 
c      IF(JP.EQ.0 .AND. I.EQ.2) THEN   
c        WRITE(*,*)'XM0,XM1=',XMAX(1,1),XMAX(1,2)
c      END IF
      if(jp.eq.0) go to 170
      go to 110
  120 continue
      if(xmax(jp,i)-xv.lt.vitol) then
        xverd(i)=xmax(jp,i)
        kverd(i)=jp
        irygg=i
      else                       
        IF(IRYGG.EQ.1) THEN    
        TEMP=XMAX(JP,I)-XV
c        WRITE(*,*)'VITOL<',TEMP 
        END IF
        go to 170
      end if
  150 continue   

  170 continue
      if(irygg.ge.2) then
        nr=nr+1
        ipt(nr)=irygg
        call ireg(xverd,irygg,aa,bb,var) 
        a(1,nr)=aa-isgn*bb*kadd
        a(2,nr)=isgn*bb
        ah(nr)=ahver
      end if
  200 continue    

      if(iov.eq.-1) then
c       I dette tilfellet er rekkefoelgen av radene invertert slik at
c       behandlingen har fulgt tilfellet iov=1. Verdiene p} a er 
c       allerede regnet om, men rader i np,xmax,amax byttes.
c                  

        istop=mrad/2
        do 370 i=1,istop 
          itemp=np(mrad+1-i)
          np(mrad+1-i)=np(i)
          np(i)=itemp
          do 375 l=1,20
            temp=amax(l,i)
            amax(l,i)=amax(l,mrad+1-i)
            amax(l,mrad+1-i)=temp
            temp=xmax(l,i)
            xmax(l,i)=xmax(l,mrad+1-i)
            xmax(l,mrad+1-i)=temp
  375     continue
  370   continue
      end if
      return
      end
         


******************************************************************************
*
*
*                 R Y G G E R
*
*     Subroutinene gjenkjenner rygger innenfor omraadet n1<i<n2,m1<k<m2 i
*     arrayen eta(1:n,1:m). Parametere:
*
*       eta   - punktarray
*       n,m   - grenser for eta.
*       xmax,amax  - plassering og styrke av maksimaene i hver linje. 
*              dvs. xmax(1,j) > xmax(2,j)...> xmax( np(j),j)  er maksp.
*              for eta(n1,m1-1+j)...eta(n2,m1-1+j).
*          a  - hver rygg har likning x(y)=a(1,j)+ a(2,j)*y. Koordinatsystem
*               er valgt slik at eta(i,j) hoerer til x=i, y=j
*         ah  - amplitude ved starten av hver rygg.
*         np  - np(j) er antall maksima funnet for linje m1-1+j.
*         nr  - antall rygger med start i j=m1 (eller m2 dersom iov=-1)
*               som er funnet. 
*     n1,n2,m1,m2 - omraade ekstrmaene skal soekes i.
*      gtol   -  bare maksima med verdi stoerre enn gtol regnes med.
*      vitol  -  maksimal x forflytning av en rygg fra en linje til den
*                neste.
*
*******************************************************************************
      subroutine rygger(eta,n,m,xmax,amax,a,ah,np,nr,n1,n2,m1,m2,
     %gtol,vitol)
      integer n,m,n1,n2,m1,m2,nr
      integer np(m2-m1+1),i,l,jp
      real eta(n,m),xmax(20,m2-m1+1),amax(20,m2-m1+1)
      real a(2,20),ah(20),gtol,vitol,aa,bb,var
*------------------------------------------------------------------------------
      real bs(0:500),amp(20),xa(20),XVERD(500),XV,AHVER,TEMP
      integer na,npun,mtall,kverd(500),irygg,J,nrygg

      npun=n2-n1+1                                 
      mtall=m2-m1+1
      nr=0                                                   


      do 100 i=m1,m2
      call limax(eta(n1,i),bs,npun,gtol,amp,xa,na)
      np(i-m1+1)=na                                
      if(na.eq.0) then 
        mtall=i-m1
        go to 60
      end if
      do 50 l=1,na
      xmax(l,i-m1+1)=xa(na+1-l)+n1 -1
      amax(l,i-m1+1)=amp(na+1-L)
   50 continue
  100 continue

   60 continue
      if(mtall.lt.2) return
      
      nrygg=np(1)

      do 200 j=1,nrygg
      xverd(1)=xmax(j,1)
      kverd(1)=j                                                   
      irygg=1
      jp=j           
      ahver=amax(j,1)
      do 150 i=2,mtall
      if(jp.gt.np(i)) jp=np(i)
      xv=xverd(i-1)
  110 continue
      if(xmax(jp,i).gt.(xv-0.1)) go to 120
      jp=jp-1 
c      IF(JP.EQ.0 .AND. I.EQ.2) THEN   
c        WRITE(*,*)'XM0,XM1=',XMAX(1,1),XMAX(1,2)
c      END IF
      if(jp.eq.0) go to 170
      go to 110
  120 continue
      if(xmax(jp,i)-xv.lt.vitol) then
        xverd(i)=xmax(jp,i)
        kverd(i)=jp
        irygg=i
      else                       
        IF(IRYGG.EQ.1) THEN    
        TEMP=XMAX(JP,I)-XV
c        WRITE(*,*)'VITOL<',TEMP 
        END IF
        go to 170
      end if
  150 continue   

  170 continue
      if(irygg.ge.2) then
        nr=nr+1
        call ireg(xverd,irygg,aa,bb,var) 
        a(1,nr)=aa+bb*(1.0-m1)
        a(2,nr)=bb
        ah(nr)=ahver
      end if
  200 continue    
c      IF(NR.EQ.0) THEN
c        IF(NP(1).EQ.0) THEN                 
c          TEMP=0
c          DO 277 J=N1,N2
c          IF(ETA(J,M1).GT.TEMP) TEMP=ETA(J,M1)
c  277     CONTINUE
c          WRITE(*,*)'FLAT NE. MAX=',TEMP
c        ELSE
c          WRITE(*,*)'NP1,NP2=',NP(1),NP(2)
c        END IF
c      END IF
      return
      end
******************************************************************************
*
*                  b e r p a r
*      Rutinen beregner hoyere ordens korreksjoner til numerisk soliton-
*  fasong. Uttrykket for overflatehevningen er:
*         y= alph*y0 -0.5*d*(alph*y0)**2       ; y0 =(cosh(k*(s-ct)))**-2
*  der s=x*cos(psi)+y*sin(psi), k=1 +alph*ak1+.. 
*  og c= 1+0.5*alph+c2*alph*alph+... Koordinatene
*  x,y og t er skalert med dyp/sqrt(alph) som karakteristisk lengde.
*
*
*       parametere:
*                dx,dy,dt - gitteravstander. Disse er regnet  i             I
*                        normalskalering, dvs. at dx er regnet i antall
*                        dyp, dt i tidsintervall det tar  aa tilbakelegge
*                        ett dyp med linear grundtvannshastighet.
*                psi   - vinkel med x-aksen, regnet i grader                I
*                ik    - verdi lik 1 gir korreksjon                         I
*                c2    - 2-ordens hastighetskorreksjon                      O
*                ak1   - 1-ordens boelgetallskorreksjon                     O
*                d     - koeff. for ikkkelinear formkorr.                   O
*                sig   - koeff i likning for Y0 :                           O
*                        2*sig*Y0''+(1.5*Y0-1)*Y0=0
******************************************************************************
      subroutine berpar(dx,dy,dt,psi,ik,c2,ak1,d,sig)
      integer ik
      real dx,dy,dt,psi,c2,ak1,d,sig
*---------------------------------------------------------------------
      real gamm1,gamm2,gamm3,kapp1,kapp2,kapp3,kapp4,kapp5
      real kapp6,kapp7,kapp8
      real r2sig,a,b,u,d2u,uy,d2y,d4y,yd2y,d2y2,dydy
      real dx2,dt2,aa,bb,dd,f1,f2,s1,s2,q,rfak,sf2,cf2
      real e2,e4,pi,dy2,my

      pi=atan(1.0)*4.0
      sf2=sin(pi*psi/180.0)
      sf2=sf2*sf2
      cf2=cos(pi*psi/180.0)
      cf2=cf2*cf2
      dt2=dt*dt
      dx2=dx*dx
      dy2=dy*dy
      e2=dx2*cf2*cf2+dy2*sf2*sf2
      e4=dx2*dx2*cf2*cf2*cf2 + dy2*dy2*sf2*sf2*sf2

      my=ik*( (dt2-dx2)*cf2 + (dt2-dy2)*sf2 )/12.0
      rfak=1.0/3.0+my

      gamm1=e2/12.0 +ik*(dx2+dy2)*sf2*cf2/12.0
      gamm2=dt2/24.0
      gamm3=gamm2

      kapp1=e4/360.0 + ik*(dx2+dy2)*(dx2*cf2+dy2*sf2)*cf2*sf2/144.0
      kapp2=dt2*dt2/1920.0
      kapp3=e2/24.0
      kapp4=(e2+dt2)/8.0
      kapp5=kapp3
      kapp6=(2*e2+dt2)/72.0-kapp2 + my*dt2/24.0 + ik*(dx2*cf2*cf2
     %       *(dt2-dx2) +dy2*sf2*sf2*(dt2-dy2))/144.0
csep1      kapp4=kapp3/3.0
csep1      kapp5=dt2/16.0
csep1      kapp6=(2*e2+dt2)/72.0
      kapp7=(2*e2+3*dt2)/12.0
      kapp8=dt2/4.0
csep1      kapp8=2*kapp5

      sig=0.5*(rfak+gamm1-gamm2-gamm3)
      r2sig=0.5/sig

      d2u=rfak+gamm1-gamm3
      uy=2.0
      u=-0.5
      d2y=0.5*rfak-1.5*gamm3-gamm2
      d4y=kapp1-kapp2+kapp6
      yd2y=kapp4+kapp5+kapp7
      d2y2=kapp3
      dydy=-kapp8

c      u= y1 +0.5*y -y*y +q*y''

       q=-gamm1+gamm2
       aa = 0.5*u
       bb =-1.0*u
       d2y = d2y +u*q

       bb = bb+ 0.5*uy
       dd = -uy
       yd2y = yd2y +q*uy

       d2y=d2y +d2u*0.5
       d2y2=d2y2 - d2u
       d4y = d4y +d2u*q

c      y''=f1*y+f2*y*y
c      y''''=f1*y''+f2*(y*y)''

       f1=r2sig
       f2=-1.5*r2sig

       d2y=d2y + f1*d4y
       d2y2=d2y2 + f2*d4y

c      (y*y)''=2y*y'' + 2y'*y'

       yd2y=yd2y + 2.0*d2y2
       dydy= dydy+ 2.0*d2y2

       aa=aa +f1*d2y
       bb=bb +f2*d2y + f1*yd2y
       dd= dd+f2*yd2y

c      y'y'= s1*y*y +s2*y*y*y

       s1=r2sig
       s2=-r2sig

       bb= bb + s1*dydy
       dd= dd + s2*dydy


       a=-0.25 +r2sig*rfak + r2sig*(gamm1-1.5*gamm2-2*gamm3)
     %   -r2sig*r2sig*(rfak+gamm1-gamm3)*(gamm1-gamm2)
     %   +r2sig*r2sig*(kapp1-kapp2+kapp6)

       b=1.5-1.5*rfak*r2sig-1.5*r2sig*(gamm1-1.5*gamm2-2*gamm3)
     %  -2.0*r2sig*(gamm1-gamm2)+r2sig*(kapp4+kapp5+kapp7-kapp8)
     %  -4.0*r2sig*(rfak +gamm1-gamm3-kapp3)
     %   +7.5*r2sig*r2sig*(rfak+gamm1-gamm3)*(gamm1-gamm2)
     %   -7.5*r2sig*r2sig*(kapp1-kapp2+kapp6)       

       d=-2 +3.0*r2sig*(gamm1-gamm2) -1.5*r2sig*(kapp4+kapp5+kapp7)
     %  +r2sig*kapp8
     %  +5*r2sig*(rfak +gamm1-gamm3-kapp3)
     %   -7.5*r2sig*r2sig*(rfak+gamm1-gamm3)*(gamm1-gamm2)
     %   +7.5*r2sig*r2sig*(kapp1-kapp2+kapp6)       

c       write(*,*)'a,b,d,sig'
c       write(*,*)a,b,d,sig

c       write(*,*) 'aa,bb,dd'
cmsol       write(*,*) aa,bb,dd

c       fortegn maa byttes om:
       a=-a
       b=-b
       d=-d
c
c    likning ser ut som  L(Y1)= (2*c2-2*ak1 +a)*Y0 +(3*ak1 +b)*Y0**2 +d*y0**3
c
       ak1=-b/3.0-0.5*d
       c2=ak1-0.5*a

       return
       end

*****************************************************************************
*
*     parametre i commonomraade i per.inc settes.
*     variablene betyr det samme som i foorige rutine.
*
*******************************************************************************
       subroutine settkoeff(A,dx,dy,dt,psi,ik)
       integer ik
       real A,dx,dy,dt,psi
       include 'per.inc'
*--------------------------------------------------------------------------
       real c2,ak1,sig

       call berpar(dx,dy,dt,psi,ik,c2,ak1,d,sig)

       if(A*d.gt.0.5) then
         write(*,*)'amplitude for stor'
       else
         alph=2.0*A/(1.0+sqrt(1.0-2.0*A*d))
       end if
       c=1.0 +alph*0.5+alph*alph*c2
       keff0=0.5*sqrt(0.5*A/sig)
       keff=0.5*sqrt(0.5*alph/sig)
       keff=keff*(1.0+ak1*alph)
c
c      keff er slik at argument til sech**2 er keff*(x-ct) der x og t
c      er regnet i normalkoordinater.
c 
       return
       end

****************************************************************************
*
*         Beregner pert. losn. for diskrete solitoner.
*         parameter:
*               x  -  romkoordinat langs forplantningsretning.           I
*                     Denne er dimensjonert med dypet som lengdeskala.
*               app1,app2 - hhv. laveste og neste approx til overfl.hevn O
***************************************************************************
       subroutine solnum(x,app1,app2)
       real app1,app2,x
       include 'per.inc'
*--------------------------------------------------------------------------
       real temp

       temp=cosh(keff0*x)
       temp=1.0/(temp*temp)
       app1= alph*(1.0-0.5*d*alph)*temp
       temp=cosh(keff*x)
       temp=1.0/(temp*temp)
       app2=alph*temp*(1.0-0.5*d*alph*temp)

       return
       end


*************************************************************************
       subroutine listsol(yr,n,ds,x0,xtopp,amp,dx,dy,dt,psi,proj,
     % ik,itape)
       integer n,itape,ik
       real yr(1500),ds,x0,xtopp,amp,dx,dy,dt,psi
       logical proj
*------------------------------------------------------------------------
       include 'per.inc'

       integer i,nsx,ii
       real xver,yfyll(1500),rtopp,cek,ua,app1,app2
       real capp1,err1,err2,cpsi,pi,dsp

       pi=4.0*atan(1.0)
       if(proj) then
         cpsi=cos(pi*psi/180.0)
       else
         cpsi=1.0
       end if

       dsp=ds*cpsi

       nsx=(dx/ds+0.001)
       if(itape.lt.0) return

       rtopp= (xtopp-x0+ds)*cpsi
       call SOLITON(yfyll,1,N,dsp,AMP,rTopp,Cek,UA,1)

       call settkoeff(amp,dx,dy,dt,psi,ik)

       write(itape,*) 4,3
       write(itape,26) amp,dx,dy,dt,psi
 26    format(1x,'#  Pertub. loes. for  A=',f10.5,'   dx=',f10.5,
     % '    dy=',f10.5,'    dt=',f10.5,'   psi=',f10.5)
       capp1= 1.0+0.5*amp
       write(itape,5) capp1,c,cek,keff
 5     format(1x,'#  1.ord. c=',f8.5,'    2.ord. c=',f8.5,
     %'    eks. c=',f8.5,'  btall=',f8.5)
       write(itape,*)'# kolonne 1: x, 2: 1.ord, 3: 2.ord',
     %  ', 4: eks. analyt'


       do 100 i=1,n
       xver=(x0+(i-1)*ds-xtopp)*cpsi
       call solnum(xver,app1,app2)
       xver=xver/cpsi+xtopp
       write(itape,10) xver,app1,app2,yfyll(i)
 10    format(1x,4f12.6)
       if(itape.gt.0) then
         write(itape+10,20) xver,app1
         write(itape+20,20) xver,app2
         write(itape+30,20) xver,yfyll(i)
 20      format(1x,2f12.6)
       end if
       if(nsx.gt.0) then
         ii=1+(i-1)/nsx
         if( (ii-1)*nsx.eq.i-1 ) then
           err1=abs(yr(ii)-app1)
           err2=abs(yr(ii)-app2)
           write(itape+40,20)xver,err1
           write(itape+50,20)xver,err2
         end if
       end if
 100   continue

       return
       end


********************************************************************
       function alamb(a,eps)
       real a,eps,alamb
*--------------------------------------------------
       alamb=-log(0.25*eps)/sqrt(3.0*a)
       return
       end


*********************************************************************
       subroutine listgitt(y,yr,n,m,dx,dt,eps,mv,a,x0,xtopp,ni,ds,ktape)
       integer n,m,mv,ni,ktape
       real y(1:n,1:m),yr(1500),dx,dt,x0,xtopp,ds,a,eps
*------------------------------------------------------------------
       integer na,i,n1,n2,nds
       real yinp(0:1501),xa(20),amp(20),gr,xver
       real xm
       real rtall,alamb
       
       if(ktape.lt.0) return

       do 100 i=1,n
         yr(i)=y(i,mv)
 100     continue

       gr=0.005
       call LIMAX(YR,YINP,N,GR,AMP,XA,NA)
       if(na.lt.1) then
         write(*,*) 'ingen topper funnet, amp. leses'
         a=rtall(0.1,'gi amplitude#')
         xtopp=alamb(a,eps)
         x0=0.0
         ds=0.25
         ni=2.0*xtopp/ds
         return
       end if

       x0=0.0
       nds=dx/0.25+1
       ds=dx/nds
       xm=(xa(1)-1.0)*dx+x0
       a=amp(1)
       n1=(xm-alamb(a,eps))/dx
       n2=(xm+alamb(a,eps)+1.0)/dx
       n1=max(1,n1)
       n2=min(n,n2)
       xtopp=xm-(n1-1)*dx
       ni=(n2-n1+1)*nds

       write(ktape,*) 2,2
       write(ktape,10) dx,dt
 10    format(1x,'#    gitterdata med dx=',f10.5,'   og dt=',f10.5)
       write(ktape,15) a,xtopp
 15    format(1x,'#    ampl=',f10.6,'  funnet ved:',f12.4)

       do 200 i=n1,n2
       xver= x0+(i-n1)*dx
       write(ktape,20) xver,yr(i)
       yr(i+1-n1)=yr(i)
 20    format(2f12.6)
 200   continue


       return
       end

       subroutine konsjekk(v)
       real v(3,2)
       include 'per.inc'
*-----------------------------------------------------------------
       integer i
       real fak2,fak3,app1,app2,y0,xver,yver

       write(*,*)'i,yver,fak2,fak3'
       do 100 i=1,3
       xver=v(i,1)
       yver=v(i,2)
       call solnum(xver,app1,app2)
       y0=1.0/cosh(keff*xver)
       y0=y0*y0
       fak2=(yver-app1)/(alph*alph*(y0*y0-y0))
       fak3=(yver-app2)/(alph*alph*alph*(y0*y0*y0-y0))
       write(*,*)i,yver,fak2,fak3
 100    continue
      return
      end



*****************************************************************************
*
*               S O L I P E R T
*
*
*     Har tisvarende funksjon som soliprgen i solny, men bruker to-ledds
*     perturbasjonsl|sning og beregner alltid overflatehevning. 
*    Endrede/nye parametere:
*                      ds --  increment i array
*                      dx,dy,dt -- gitteravstander,
*                      psi -- vinkel i grader med x-akse
*                      ivis -- =0 det tabuleres normalt kammen
*                              =1 det tabuleres langs x-aksen
*                              =2 det tabuleres langs y-aksen
*                             Topp refererer seg alltid til tabuleringsretning
*                      ik -- =1 angir korreksjon 
**************************************************************************
      subroutine solipert(Y,n0,N,DS,dx,dy,dt,psi,ivis,ik,AMP,TOPP,C,eps)   
      INTEGER n0,N,ivis,ik
      REAL Y(N0:N),DS,AMP,topp,C,UA,eps,dx,dy,dt,psi
*---------------------------------------------------------------------------
      integer i,na,nb
      real dfak,blen,halvb,pi,xver,app1

      pi=4.0*atan(1.0)

      do 50 i=n0,n
         y(i)=0.0
 50   continue

      if(ivis.eq.0) then
        dfak=1.0
      else
        if(ivis.eq.1) then
          dfak=cos(pi*psi/180.0)
        else
          dfak=sin(pi*psi/180.0)
        end if
      end if

      if(amp.le.0.0) return

      blen=halvb(amp,eps)/dfak

      na= (topp-blen)/ds -1
      na=max(na,n0)
      nb=(topp+blen)/ds +1
      nb=min(nb,n)

      call settkoeff(amp,dx,dy,dt,psi,ik)

      do 100 i=na,nb

      xver= (topp-i*ds)*dfak
      call solnum(xver,app1,y(i))

 100  continue

      return

      end
*********************************************************
* 
*                RIM
*
*   Legger til "border" med nuller rundt en array.
*   parametere:
*           a  - array. Ved input er den en tett pakket (n,m) array,   I/O
*                ved output en (n+2*ix,m+2*iy) array der den opprinnelige
*                datamengden finne i  a(i,j) for  ix +1 <= i <= ix+n, 
*                                                 iy +1 <= j <= iy+m
*           n,m - arraygrenser                                         I
*           ix,iy - bordbredder                                        I 
*********************************************************
      subroutine rim(a,n,m,ix,iy)
      integer n,m,ix,iy
      real a((n+2*ix)*(m+2*iy))
*-------------------------------------------------------
      integer i,k,ianew,iaold,ian2,nn

      do 150 k=m,1,-1
       iaold=(k-1)*n
       ianew=(k+iy-1)*(n+2*ix)+ix
       do 100 i=n,1,-1
         a(ianew+i)=a(iaold+i)
 100     continue

       do 110 i=1,ix
         a(ianew-ix+i)=0.0
         a(ianew+n+i)=0.0
 110     continue

 150   continue

      nn=n+2*ix

      do 200 k=1,iy
         ianew=(k-1)*nn
         ian2=(m+iy+k-1)*nn
         do 220 i=1,nn
           a(ianew+i)=0.0
           a(ian2+i)=0.0
 220       continue
 200     continue

       return
       end

**********************************************************
*
*                SETSPLINE
*
*    Bestemmer B-spline koeffisientene i et kubisk bi-spline
*    ved gjentatte kall paa bsplaut. Dette innebaerer at likningene lukkes
*    vha differens-tilnaermelser til den 1ste deriverte paa randa. 
*      Parametere:
*             a  - ved input: (n,m) array med funksjonsverdier           I/O
*                  ved output  (0:n+1,0:m+1) array med spline koeffisienter
*             n,m - antall datapunkter                                   I
*             nms - maksimal storrelse for a                             I
*             wrk - arbeidsarray, min storrelse = (max(n,m)+2)*5         D
*             iflag - returnerer verdi 1 dersom (n+2)*(m+2)> nms         O
*****************************************************************************
      subroutine setspline(a,n,m,nms,wrk,iflag)
      integer n,m,nms,iflag
      real a(nms),wrk((n+2)*5)
*----------------------------------------------------
      integer i,k,ih,ip,n2,m2

      n2=n+2
      m2=m+2

      if(n2*m2.gt.nms) then
        iflag=1
        return
      end if

      iflag=0

      call rim(a,n,m,1,1)

      do 100 k=1,m

       ip=k*n2+1
       call bsplaut(a(ip+1),a(ip),n,wrk)
 100   continue

      ih=m+3

      do 200 i=1,n2

       do 210 k=1,m
        wrk(1+k)= a(k*n2+i)
 210    continue
       call bsplaut(wrk(2),wrk(1),m,wrk(ih))
       do 220 k=1,m2
        a((k-1)*n2+i)=wrk(k)
 220    continue
 200   continue

      return
      end

**************************************************************************
*
*                     S P K O F F                                        *   
*     Beregner splines-koeffisienter i et 1-D spline             
*      parametere:
*     s - koordinat relativt til naermeste venstre punkt                I
*         vi maa alltid ha 0<=s<=1
*     d - koeff array. d(2)  refererer seg til naermeste                O
*                venstre nodepunkt ( s=0 gir dette punkt)
*     ider - orden paa derivert ( 0, 1 eller 2)                         I
*****************************************************************************
      subroutine spkoff(s,d,ider)                
      integer ider
      REAL s,d(4)                
*------------------------------------------------ 
      REAL XV     

                  
      XV=S-1.0      
      if(s.lt.0 .or. s.gt.1.0) then
        write(0,*)'punkt utenfor interval i spkoff'
        return
      end if
                    
      if(ider.le.0) then
       d(1)=-XV*XV*XV/6.0
       d(2)=(2.0/3.0-S*S*(1.0-0.5*S))          
       d(3)=(2.0/3.0-XV*XV*(1+0.5*XV))
       d(4)=S*S*S/6.0               
      else
       if(ider.eq.1) then
        d(1)=-XV*XV*0.5
        d(2)=-2.0*S+1.5*s*s          
        d(3)=-(2.0*xv+1.5*xv*xv)
        d(4)=S*S*0.5               
       else
        d(1)=-XV
        d(2)=-2.0+3.0*s          
        d(3)=-(2.0+3.0*xv)
        d(4)=S               
       end if
      end if

      RETURN        
      END           


**************************************************************************
*
*                     S P D E R 2
*   
*     Beregner funksjonsverdier av en splines-interpolant generert ved
*     setspline. 
*      parametere:
*              X,y - koordinater                                        I
*              X0,y0 - posisjon tilh|rende B-spline med koeff. bs(1,1)  I
*              xd,dy - punktavstand                                      I
*              BS - array med spline-koeff.                             I
*              N,m - antall interpolasjonspunkter                       I
*              idx,idy - orden av derivert i hhv. x og y-retning        I
*              eps - justerings-grense, faller et punkt mindre enn      I
*                    eps*dx etc. utenfor, trekkes det inn.
*              f -  verdi av interpolant                                O
*****************************************************************************
      subroutine spder2(X,y,X0,y0,Dx,dy,BS,N,m,idx,idy,eps,f)   
      INTEGER N,m,idx,idy
      REAL BS(0:N+1,0:m+1),X,Y,X0,Y0,Dx,dy,eps,f
*------------------------------------------------ 
      REAL Sy,sx,xv,yv,bx(4),by(4),w
      INTEGER NX,ny,nn,i,j
                 

      XV=(X-X0)/DX+1.0                          
      NX=XV         
      SX=XV-NX      

      if(nx.le.0 .or. nx.gt.n) then
        if(nx.eq.0 .and. sx .gt.(1.0-eps)) then
           nx=1
           sx=0.0
        else
          if(nx.eq.n .and. sx .lt.eps) then
           nx=n-1
           sx=1
          else 
            write(0,*)'x-punkt utenfor interval i spder2'
            f=0
            return
          end if
        end if
      end if

      YV=(Y-Y0)/DY+1.0                          
      NY=YV         
      SY=YV-NY       


      if(ny.le.0 .or. ny.gt.m) then
        if(ny.eq.0 .and. sy .gt.(1.0-eps)) then
           ny=1
           sy=0.0
        else
          if(ny.eq.m .and. sy .lt.eps) then
           ny=m-1
           sy=1
          else 
            write(0,*)'y-punkt utenfor interval i spder2'
            f=0
            return
          end if
        end if
      end if



      call spkoff(sx,bx,idx)
      call spkoff(sy,by,idy)

      f=0.0

      do 200 i=1,4
       nn=nx-2+i
       w=0.0
       do 100 j=1,4
        w=w+by(j)*bs(nn,ny-2+j)
 100    continue
       f=f+w*bx(i)
 200   continue

      RETURN        
      END           


**************************************************************************
*
*                     S P D E R 2 G
*   
*     Beregner funksjonsverdier av en splines-interpolant generert ved
*     setspline. 
*      parametere:
*              X,y - koordinater                                        I
*              X0,y0 - posisjon tilh|rende B-spline med koeff. bs(1,1)  I
*              xd,dy - punktavstand                                      I
*              BS - array med spline-koeff.                             I
*              N,m - antall interpolasjonspunkter                       I
*              idx,idy - orden av derivert i hhv. x og y-retning        I
*              eps - justerings-grense, faller et punkt mindre enn      I
*                    eps*dx etc. utenfor, trekkes det inn.
*              fdef - verdi som settes dersom x,y er utenfor
*              f -  verdi av interpolant                                O
*****************************************************************************
      subroutine spder2g(X,y,X0,y0,Dx,dy,BS,N,m,idx,idy,eps,fdef,f)   
      INTEGER N,m,idx,idy
      REAL BS(0:N+1,0:m+1),X,Y,X0,Y0,Dx,dy,eps,fdef,f
*------------------------------------------------ 
      REAL Sy,sx,xv,yv,bx(4),by(4),w
      INTEGER NX,ny,nn,i,j
                 

      XV=(X-X0)/DX+1.0                          
      NX=XV         
      SX=XV-NX      

      if(nx.le.0 .or. nx.gt.n) then
        if(nx.eq.0 .and. sx .gt.(1.0-eps)) then
           nx=1
           sx=0.0
        else
          if(nx.eq.n .and. sx .lt.eps) then
           nx=n-1
           sx=1
          else 
            write(0,*)'x-punkt utenfor interval i spder2'
            f=fdef
            return
          end if
        end if
      end if

      YV=(Y-Y0)/DY+1.0                          
      NY=YV         
      SY=YV-NY       


      if(ny.le.0 .or. ny.gt.m) then
        if(ny.eq.0 .and. sy .gt.(1.0-eps)) then
           ny=1
           sy=0.0
        else
          if(ny.eq.m .and. sy .lt.eps) then
           ny=m-1
           sy=1
          else 
            write(0,*)'y-punkt utenfor interval i spder2'
            f=0
            return
          end if
        end if
      end if



      call spkoff(sx,bx,idx)
      call spkoff(sy,by,idy)

      f=0.0

      do 200 i=1,4
       nn=nx-2+i
       w=0.0
       do 100 j=1,4
        w=w+by(j)*bs(nn,ny-2+j)
 100    continue
       f=f+w*bx(i)
 200   continue

      RETURN        
      END           



**************************************************************************
*
*                     G P I N T
*   
*     Beregner funksjonsverdier av en splines (generert ved setspline) 
*     eller bilinear interpolant. 
*      parametere:
*              X,y - koordinater                                        I
*              X0,y0 - posisjon tilh|rende B-spline med koeff. bs(1,1)  I
*              dx,dy - punktavstand                                      I
*              BS - array med spline-koeff.                             I
*              N,m - antall interpolasjonspunkter                       I
*              idx,idy - orden av derivert i hhv. x og y-retning        I
*              eps - justerings-grense, faller et punkt mindre enn      I
*                    eps*dx etc. utenfor, trekkes det inn.
*              fdef - verdi som settes dersom x,y er utenfor
*              f -  verdi av interpolant                                O
*          ityp - verdi=1: bilin, verdi=2: splines                      I
*****************************************************************************
      subroutine gpint(x,y,x0,y0,dx,dy,bs,n,m,idx,idy,eps,fdef,f,ityp)   
      INTEGER N,m,idx,idy,ityp
      REAL BS(2-ityp:n+ityp-1,2-ityp:m+ityp-1)
      real X,Y,X0,Y0,Dx,dy,eps,fdef,f
*------------------------------------------------ 
      REAL Sy,sx,xv,yv,bx(4),by(4),w
      real rdy,rdx,wx1,wx2,wy1,wy2,w1,w2,w3,w4
      INTEGER NX,ny,nn,i,j
                 
      rdx=1.0/dx
      rdy=1.0/dy

      XV=(X-X0)/DX+1.0                          
      NX=XV         
      SX=XV-NX      

      if(nx.le.0 .or. nx.ge.n) then
        if(nx.eq.0 .and. sx .gt.(1.0-eps)) then
           nx=1
           sx=0.0
        else
          if(nx.eq.n .and. sx .lt.eps) then
           nx=n-1
           sx=1
          else 
cc            write(0,*)'x-punkt utenfor interval i gpint',ityp
            f=fdef
            return
          end if
        end if
      end if

      YV=(Y-Y0)/DY+1.0                          
      NY=YV         
      SY=YV-NY       


      if(ny.le.0 .or. ny.ge.m) then
        if(ny.eq.0 .and. sy .gt.(1.0-eps)) then
           ny=1
           sy=0.0
        else
          if(ny.eq.m .and. sy .lt.eps) then
           ny=m-1
           sy=1
          else 
cc            write(0,*)'y-punkt utenfor interval i gpint',ityp
            f=fdef
            return
          end if
        end if
      end if

      if(ityp.eq.1) then
        if(idx.eq.0) then
           wx1=1.0-sx
           wx2=sx
        else
           if(idx.eq.1) then
             wx1=-rdx
             wx2=rdx
           else
             wx1=0.0
             wx2=0.0
           end if
        end if

        if(idy.eq.0) then
           wy1=1.0-sy
           wy2=sy
        else
           if(idy.eq.1) then
             wy1=-rdy
             wy2=rdy
           else
             wy1=0.0
             wy2=0.0
           end if
        end if

        w1=wx1*wy1
        w2=wx2*wy1
        w3=wx1*wy2
        w4=wx2*wy2

        f=w1*bs(nx,ny)+w2*bs(nx+1,ny)+w3*bs(nx,ny+1)+w4*bs(nx+1,ny+1)
      else

        call spkoff(sx,bx,idx)
        call spkoff(sy,by,idy)
  
        f=0.0

        do 200 i=1,4
          nn=nx-2+i
          w=0.0
          do 100 j=1,4
            w=w+by(j)*bs(nn,ny-2+j)
 100      continue
          f=f+w*bx(i)
 200    continue
      end if

      RETURN        
      END           



**************************************************************************
*
*          B S T A G V
*
*   Adderer f*  interpolant av node verdier i matrise a til den delen av 
*   matrise b  som ligger ekte innenfor a. En har valget mellom bilineaer
*   interpolasjon der a inneholder funksjonsverdier og spline-interpolasjon
*   der a inneholder koeff for bi-kubiske B splines.
*     parametere:
*          a - matrise med kilde-data                                   I
*          na,ma - grenser for del av a som inneholder data             I
*                  for biliner int. er a en tett (na,ma) matrise, for
*                  spline int en tett (0:na+1,0:ma+1) matrise. 
*          xa,ya - posisjon av a(1,1)                                   I
*          b - matrise der a's data skal adderes                       I/O
*          nsb - fysisk grense for foerste indeks i b                   I
*          nb,mb - grenser for del av b som inneholder data             I
*          xb,yb - posisjon av b(1,1)                                   I
*          dx,dy - gitteravstander, felles for a og b                   I
*          f - faktor i " b=b+f*a"                                      I
*          ityp - verdi=1: bilin, verdi=2: splines                      I
**************************************************************************
      subroutine bstagv(a,na,ma,xa,ya,b,nsb,nb,mb,xb,yb,dx,dy,f,ityp)
      integer na,ma,nsb,nb,mb,ityp
      real a(2-ityp:na+ityp-1,2-ityp:ma+ityp-1),b(nsb,mb)
      real xa,ya,xb,yb,dx,dy,f
*-------------------------------------------------------------------------
      integer i,j,nadd,madd,ibeg,islutt,jbeg,jslutt,im,ip,jm,jp
      integer jmm,jpp,ipp
      real xpluss,ypluss,w1,w2,w3,w4,q1,q2,q3,q4,r1,r2,r3,r4,r(4),q(4)

c     forklaring av variable
c     nadd,madd: a(i+nadd,j+madd) er naermeste nabo (mot avtakende x og y)
c                for b(i,j)
c     xpluss,ypluss er avstanden mellom dividert med dx,dy
c     (ibeg...islutt)x(jbeg...jslutt) angir del av b som skal gis bidrag 


      if(xb.gt.xa) then
        nadd=(xb-xa)/dx
        ibeg=1
      else
        nadd=-( (xa-xb)/dx+1 )
        ibeg=1-nadd
      end if

      xpluss=( xb+(ibeg-1)*dx - (xa+(ibeg+nadd-1)*dx) )/dx

      if( (xb+(nb-1)*dx).gt.(xa+(na-1)*dx) ) then
        islutt=na-1-nadd
      else
        islutt=nb
      end if



      if(yb.gt.ya) then
        madd=(yb-ya)/dy
        jbeg=1
      else
        madd=-( (ya-yb)/dy+1 )
        jbeg=1-madd
      end if

      ypluss=( yb+(jbeg-1)*dy - (ya+(jbeg+madd-1)*dy) )/dy

      if( (yb+(mb-1)*dy).gt.(ya+(ma-1)*dy) ) then
        jslutt=ma-1-madd
      else
        jslutt=mb
      end if

      if(ityp.eq.1) then

      w1=f*(1-xpluss)*(1-ypluss)
      w2=f*xpluss*(1-ypluss)
      w3=f*(1-xpluss)*ypluss
      w4=f*xpluss*ypluss


      do 100 j=jbeg,jslutt

      jm=j+madd
      jp=jm+1
      ip=ibeg+nadd

      do 100 i=ibeg,islutt
      im=ip
      ip=im+1
      b(i,j)=b(i,j)+w1*a(im,jm)+w2*a(ip,jm)+w3*a(im,jp)+w4*a(ip,jp)
 100  continue

      else

        call spkoff(xpluss,r,0)
        call spkoff(ypluss,q,0)
        r1=f*r(1)        
        r2=f*r(2)        
        r3=f*r(3)        
        r4=f*r(4)        
        q1=q(1)
        q2=q(2)
        q3=q(3)
        q4=q(4)

      do 200 j=jbeg,jslutt

      jm=j+madd
      jp=jm+1
      jpp=jm+2
      jmm=jm-1
      ip=ibeg+nadd
      ipp=ip+1
      im=ip-1
      w2=q1*a(im,jmm)+q2*a(im,jm)+q3*a(im,jp)+q4*a(im,jpp)
      w3=q1*a(ip,jmm)+q2*a(ip,jm)+q3*a(ip,jp)+q4*a(ip,jpp)
      w4=q1*a(ipp,jmm)+q2*a(ipp,jm)+q3*a(ipp,jp)+q4*a(ipp,jpp)

      do 200 i=ibeg,islutt
      im=ip
      ip=im+1
      ipp=im+2
      w1=w2
      w2=w3
      w3=w4
      w4=q1*a(ipp,jmm)+q2*a(ipp,jm)+q3*a(ipp,jp)+q4*a(ipp,jpp)
      b(i,j)=b(i,j)+w1*r1+w2*r2+w3*r3+w4*r4
 200  continue

      end if

      return
      end


**************************************************************************
*
*          G G S T A G V
*
*   Adderer f*  interpolant av node verdier i matrise a til den delen av 
*   matrise b  som ligger ekte innenfor a. En har valget mellom bilineaer
*   interpolasjon der a inneholder funksjonsverdier og spline-interpolasjon
*   der a inneholder koeff for bi-kubiske B splines. Rutinen er ikke spesielt 
*   effektiv.
*     parametere:
*          a - matrise med kilde-data                                   I
*          na,ma - grenser for del av a som inneholder data             I
*                  for biliner int. er a en tett (na,ma) matrise, for
*                  spline int en tett (0:na+1,0:ma+1) matrise.
*          dxa,dya - gitteravstander i a.                               I 
*          xa,ya - posisjon av a(1,1)                                   I
*          b - matrise der a's data skal adderes                       I/O
*          nsb - fysisk grense for foerste indeks i b                   I
*          nb,mb - grenser for del av b som inneholder data             I
*          dxb,dyb - gitteravstander i b                                I
*          xb,yb - posisjon av b(1,1)                                   I
*          fa,fb - faktor i " b=fb*b+fa*a"
*          ityp - verdi=1: bilin, verdi=2: splines                      I
**************************************************************************
      subroutine ggstagv(a,na,ma,dxa,dya,xa,ya,b,nsb,nb,mb,dxb,dyb,xb,yb
     %                  ,fa,fb,ityp,eps,vdef)
      integer na,ma,nsb,nb,mb,ityp
      real a(2-ityp:na+ityp-1,2-ityp:ma+ityp-1),b(nsb,mb)
      real xa,ya,xb,yb,dxa,dya,dxb,dyb,fa,fb,eps,vdef
*-------------------------------------------------------------------------
      integer i,k
      real xv,yv,val


      do 200 k=1,mb
        yv=yb+(k-1)*dyb
        do 100 i=1,nb
          xv=xb+(i-1)*dxb
          call gpint(xv,yv,xa,ya,dxa,dya,a,na,ma,0,0,eps,vdef,val,ityp)
          b(i,k)=fb*b(i,k)+fa*val
 100  continue   
 200  continue   

      return
      end


*************************************************************
*     Finds the grid section ia ib that is within the 
*     interval xa,xb. A margin eps is taken into account
*     If there is no overlap the subroutine returns with
*        ib<ia<0
*************************************************************
      subroutine fsec(n,dx,x0,xa,xb,eps,ia,ib)
      integer n,ia,ib
      real dx,x0,xa,xb,eps
*--------------------------------------------------------------
      

      ia=(xa-eps-x0)/dx+1
      if(xa-eps.gt.x0+(ia-1)*dx) ia=ia+1
      
      if(ia.lt.1) ia=1
      if(ia.gt.n) then
         ia=-1
         ib=ia-1
         return
      end if

      ib=(xb+eps-x0)/dx+1
      if(xb+eps.lt.x0+(ib-1)*dx) ib=ib-1
      
      if(ib.lt.1) then
         ia=-2
         ib=ia-1
         return
      end if

      if(ib.gt.n) ib=n


      return
      end




**************************************************************************
*
*          Q S T A G V
*
*   Adderer f*  interpolant av node verdier i matrise a til den delen av 
*   matrise b  som ligger ekte innenfor a. En har valget mellom bilineaer
*   interpolasjon der a inneholder funksjonsverdier og spline-interpolasjon
*   der a inneholder koeff for bi-kubiske B splines.
*     parametere:
*          a - matrise med kilde-data                                   I
*          na,ma - grenser for del av a som inneholder data             I
*                  for biliner int. er a en tett (na,ma) matrise, for
*                  spline int en tett (0:na+1,0:ma+1) matrise.
*          dxa,dya 
*          xa,ya - posisjon av a(1,1)                                   I
*          b - matrise der a's data skal adderes                       I/O
*          nsb - fysisk grense for foerste indeks i b                   I
*          nb,mb - grenser for del av b som inneholder data             I
*          xb,yb - posisjon av b(1,1)                                   I
*          dxb,dyb - gitteravstander for  b                             I
*          f - faktor i " b=b+f*a"                                      I
*          ityp - verdi=1: bilin, verdi=2: splines                      I
*          eps - tolerense                                              I
**************************************************************************
      subroutine qstagv(a,na,ma,xa,ya,dxa,dya,b,nsb,nb,mb,xb,yb,dxb,
     %    dyb,f,ityp,eps)
      integer na,ma,nsb,nb,mb,ityp
      real a(2-ityp:na+ityp-1,2-ityp:ma+ityp-1),b(nsb,mb)
      real dxa,dya,xa,ya,xb,yb,dxb,dyb,f,eps
*-------------------------------------------------------------------------
      integer i,j,ibeg,islutt,jbeg,jslutt,im,ip,jm,jp
      integer jmm,jpp,ipp,imm
      real xpluss,ypluss,w1,w2,w3,w4,q1,q2,q3,q4,r(4),q(4)
      real xm,ym,xv,yv

c     forklaring av variable
c         a(im,jm) er naermeste nabo (mot avtakende x og y)
c                for b(i,j)
c     xpluss,ypluss er avstanden fra a(im,jm) til b(i,j) dividert med dx,dy
c     (ibeg...islutt)x(jbeg...jslutt) angir del av b som skal gis bidrag 


      xm=xa+(na-1)*dxa
      call fsec(nb,dxb,xb,xa,xm,eps,ibeg,islutt)

      ym=ya+(ma-1)*dya
      call fsec(mb,dyb,yb,ya,ym,eps,jbeg,jslutt)





      if(ityp.eq.1) then



      do 100 j=jbeg,jslutt

      yv=yb+(j-1)*dyb
      jm=(yv-ya)/dya+1
      if(jm.lt.1)jm=1
      ypluss=(yv-ya-(jm-1)*dya)/dya
      jp=jm+1


      do 100 i=ibeg,islutt

      xv=xb+(i-1)*dxb
      im=(xv-xa)/dxa+1
      if(im.lt.1)im=1
      xpluss=(xv-xa-(im-1)*dxa)/dxa
      ip=im+1


      w1=f*(1-xpluss)*(1-ypluss)
      w2=f*xpluss*(1-ypluss)
      w3=f*(1-xpluss)*ypluss
      w4=f*xpluss*ypluss


      b(i,j)=b(i,j)+w1*a(im,jm)+w2*a(ip,jm)+w3*a(im,jp)+w4*a(ip,jp)
 100  continue

      else


      do 200 j=jbeg,jslutt

      yv=yb+(j-1)*dyb
      jm=(yv-ya)/dya+1
      if(jm.lt.1)jm=1
      ypluss=(yv-ya-(jm-1)*dya)/dya
      jp=jm+1
      jpp=jm+2
      jmm=jm-1

      call spkoff(ypluss,q,0)
      q1=q(1)
      q2=q(2)
      q3=q(3)
      q4=q(4)




      do 200 i=ibeg,islutt

      xv=xb+(i-1)*dxb
      im=(xv-xa)/dxa+1
      if(im.lt.1)im=1
      xpluss=(xv-xa-(im-1)*dxa)/dxa
      ip=im+1
      ipp=im+2
      imm=im-1


      call spkoff(xpluss,r,0)


      w1=q1*a(imm,jmm)+q2*a(imm,jm)+q3*a(imm,jp)+q4*a(imm,jpp)
      w2=q1*a(im,jmm)+q2*a(im,jm)+q3*a(im,jp)+q4*a(im,jpp)
      w3=q1*a(ip,jmm)+q2*a(ip,jm)+q3*a(ip,jp)+q4*a(ip,jpp)
      w4=q1*a(ipp,jmm)+q2*a(ipp,jm)+q3*a(ipp,jp)+q4*a(ipp,jpp)
      b(i,j)=b(i,j)+f*(w1*r(1)+w2*r(2)+w3*r(3)+w4*r(4))
 200  continue

      end if

      return
      end


*************************************************************
*
*         AIRYBER
*
*     Computes Airy function by fourth order Runge-Kutta method.
*             ai - array containing function values                I
*             n0,n1 - defines interval:  -n1*dx < x < n0*dx         O
*                     observe: ai(n1+1) corresponds to x=0
*             dx - increment
************************************************************************
      subroutine airyber(ai,n0,n1,dx)
      integer n0,n1
      real ai(1+n0+n1),dx
*----------------------------------------------
      integer j,n,np
      real x0,x,y(2),w(10),ff
      external fa

      n=n0+n1
      np=n+1
      x0=n0*dx
      x=-x0

      y(1)=(x0**(-0.25))*exp(-2*x0**1.5/3)
cc     WARNING: Noe er upresist
cc     y(2)=f*(0.25/x0+sqrt(x0))
      y(2)=0.0
      ai(np)=y(1)

      do 100 j=1,n
       call rkstep(y,w,2,x,dx,fa)
       x=x+dx
       ai(np-j)=y(1)
 100  continue

cc    Ai(0)=0.35502 80538 87817
      ff=0.35502805/ai(n1+1)

      do 200 j=1,np
       ai(j)=ff*ai(j)
 200  continue

      return

      end

*********************************************************************
*     FA defines Airys differential equation 
**********************************************************************
      subroutine fa(n,x,y,dy)
      integer n
      real x,y(2),dy(2)
*-----------------------------------------------------------

      dy(1)=y(2)
      dy(2)=-x*y(1)

      return
      end 



**********************************************************************
*
*            RKSTEP
*    
*     Performs one step in a fourth order Runge-Kutta integration
*     for the set
*             Y'_j=F_j(Y,t)      j=1....n  
*
*      Parameters:
*             y - the unknown function: an array  that on input/output  I/O
*                 holds the old values and new values respectively
*             w - work space; five times the size of array y            D
*             n - dimension of y = number of equations                  I
*             t - time. At input: old time, output: new time = old      I
*             dt - time step                                            I
*             f - subroutine that defines the right hand side of the    I
*                 differential equation. f must be decleared according to
*                   subroutine f(n,t,y,dy)
*                   integer n
*                   real t, y(n), dy(n)
*                 where
*                   n - number of equations  I
*                   t - time                 I
*                   y - function values      I
*                   dy - derivatives         O
*                
*       NOTES:
*                1. f must be decleared external in the calling program.
*                2. additional parameters in f must be introduced
*                through a common block.
************************************************************************
      subroutine rkstep(y,w,n,t,dt,f)
      integer n
      real y(n),w(5*n),t,dt
      external f
*--------------------------------------------------
      integer i,k1,k2,k3,k4,kw
      real tp,tpp

      k1=1
      k2=k1+n
      k3=k2+n
      k4=k3+n
      kw=k4+n

      tp=t+0.5*dt
      tpp=t+dt
      call f(n,t,y,w(k1))

      do 100 i=1,n
       w(kw-1+i)=y(i)+0.5*dt*w(k1-1+i)
 100  continue

      call f(n,tp,w(kw),w(k2))


      do 200 i=1,n
       w(kw-1+i)=y(i)+0.5*dt*w(k2-1+i)
 200  continue

      call f(n,tp,w(kw),w(k3))


      do 300 i=1,n
       w(kw-1+i)=y(i)+dt*w(k3-1+i)
 300  continue

      call f(n,tpp,w(kw),w(k4))



      do 500 i=1,n
       y(i)=y(i)+dt*(w(k1-1+i)+2*w(k2-1+i)+2*w(k3-1+i)+w(k4-1+i) )/6.0
 500  continue

      return
      end







**********************************************************************************
*
*     If dob=.true. a contains the cofficient in a polynomial in x**2
***********************************************************************************
      subroutine modbekoff(a,n,dob)
      integer n
      real a(0:n)
      logical dob
*----------------------------------------------------------------------------------
      integer i,istep,nn

      nn=n/2
      if(dob) then
        istep=1
      else
        istep=2
        do 50 i=1,nn
         a(i*istep-1)=0.0
50      continue

        if(2*nn.ne.n) then
          a(n)=0.0
        end if   
      end if

      a(0)=1
      do 100 i=1,nn
       a(i*istep)=a((i-1)*istep)/(4.0*i*i)
 100  continue   
      return
      end



**********************************************************************************
      subroutine modnekoff(a,n)
      integer n
      real a(0:n)
*----------------------------------------------------------------------------------
      integer i
      real fak,dsum

      a(0)=0
cc    formula in B & Orzag is strange
      do 100 i=1,n
       a(i)=0.0
cc    must be corrected
 100  continue   
      return
      end


**********************************************************************************
      function potber(a,n,s)
      integer n
      real potber,a(0:n),s
*----------------------------------------------------------------------------------
      integer i
      real sum
   
      sum=a(n)
      do 100 i=n-1,0,-1
       sum=a(i)+s*sum
 100  continue   
      potber=sum
      return
      end

**********************************************************************************
      subroutine potder(a,b,n)
      integer n
      real a(0:n),b(0:n-1)
*----------------------------------------------------------------------------------
      integer i
      
      do 100 i=1,n
       b(i-1)=a(i)*i
 100  continue   
   
      return
      end


*********************************************************************
*     MB defines modified Bessel's differential equation 
**********************************************************************
      subroutine mb(n,x,y,dy)
      integer n
      real x,y(2),dy(2)
*-----------------------------------------------------------

      dy(1)=y(2)
      dy(2)=y(1) -y(2)/x

      return
      end 

******************************************************************
      subroutine premobe(bmod,n,dx,nledd,wrk)
      integer n,nledd
      real bmod(0:n),dx,wrk(0:nledd)
*--------------------------------------------------------------
      integer i
      real s,potber
      external potber

      call modbekoff(wrk,nledd,.true.)

      bmod(0)=1
      do 100 i=1,n
       s=i*i*dx*dx
       bmod(i)=potber(wrk,nledd,s)
 100  continue
      return
      end

******************************************************************
      subroutine commobe(bmod,db,n,dx,na,nledd,wrk,w2)
      integer n,na,nledd
      real bmod(0:n),db(0:n),dx,wrk(0:2*nledd),w2(0:2*nledd)
*--------------------------------------------------------------
      integer i,n2
      real x0,x,y(2),w(10),potber
      external mb,potber

      call premobe(bmod,na,dx,nledd,wrk)
      n2=nledd*2
      call modbekoff(wrk,n2,.false.)
      call potder(wrk,w2,n)
      db(0)=w2(0)
      do 50 i=1,n
       x=i*dx
       db(i)=potber(w2,n2,x)
 50   continue

      x=na*dx
      y(1)=bmod(na)
      y(2)=db(na)
cc    compute derivative      
      do 100 i=na+1,n
       call rkstep(y,w,2,x,dx,mb)
       x=x+dx
       bmod(i)=y(1)
 100  continue

      return
      end

      

**********************************************************************************
      function kajgreen(n,r)
      integer n
      real kajgreen,r
*----------------------------------------------------------------------------------
      integer i,sgn,i2p
      real sum,pi,r2,arg

      r2=r*r
      pi=3.14159265358979
      sum=0.0
      sgn=-1
      do 100 i=0,n
       sgn=-sgn
       i2p=2*i+1
       arg=i2p*i2p+r2
       sum=sum+sgn*i2p/(arg*sqrt(arg))
 100  continue   
      kajgreen=sum/pi

      return
      end

**********************************************************************************
      function dkajgreen(n,r2)
      integer n
      real dkajgreen,r2
*----------------------------------------------------------------------------------
      integer i,sgn,i2p
      real sum,pi,arg

     
      pi=3.14159265358979
      sum=0.0
      sgn=-1
      do 100 i=0,n
       sgn=-sgn
       i2p=2*i+1
       arg=i2p*i2p+r2
       sum=sum+sgn*i2p/(arg*sqrt(arg))
 100  continue   
      dkajgreen=sum/pi

      return
      end



**********************************************************************************
*     Computes the Green function and its derivative with respect to r2
*******************************************************************************
      subroutine Gkjgreen(n,r2,g,gd,gdd)
      integer n
      real r2,g,gd,gdd
*---------------------------------------------------------------------------------
      integer i,sgn,i2p
      real sum,sumd,sumdd,ledd,pi,arg

     
      pi=3.14159265358979
      sum=0.0
      sumd=0.0
      sumdd=0.0
      sgn=-1
      do 100 i=0,n
       sgn=-sgn
       i2p=2*i+1
       arg=1.0/(i2p*i2p+r2)
       ledd=sgn*i2p*arg*sqrt(arg)
       sum=sum+ledd
       sumd=sumd-1.5*ledd*arg
       sumdd=sumdd+3.75*ledd*arg*arg
 100  continue   
      g=sum/pi
      gd=sumd/pi
      gdd=sumdd/pi

      return
      end

******************************************************************
      subroutine Kajitab(a,n,ds,nledd,ipol,w)
      integer n,nledd,ipol
      real a(0:n+2),ds,w(5*(n+3))
*----------------------------------------------------------------
      integer i,na,nw,np
      real arg,der1,der2,dkajgreen
      external dkajgreen

      do 100 i=0,n
       arg=i*ds
       a(i)=dkajgreen(nledd,arg)
 100  continue

      if(ipol.le.1) return
      na=n+2
      nw=n+4
      np=n+1       
      der2=-0.5*(4.0*a(n-1)-3.0*a(n) -a(n-2))/ds
      der1=0.5*(4.0*a(1)-3.0*a(0)-a(2))/ds
      call bsplw(a,w(1),np,ds,der1,der2,w(nw))
      do 200 i=0,na
       a(i)=w(i+1)
 200  continue

      return
      end


******************************************************************
      subroutine GKtab(a,ad,n,ds,nledd,ipol,w)
      integer n,nledd,ipol
      real a(0:n+2),ad(0:n+2),ds,w(5*(n+3))
*----------------------------------------------------------------
      integer i,na,nw,np
      real arg,der1,der2,g,gd,gdd,gddiff

      i=0
      arg=i*ds
      call Gkjgreen(nledd,arg,g,gd,gdd)
      a(i)=g
      ad(i)=gd
      gddiff=gdd

      do 100 i=1,n
       arg=i*ds
       call Gkjgreen(nledd,arg,g,gd,gdd)
       a(i)=g
       ad(i)=gd
 100  continue

      if(ipol.le.1) return
      na=n+2
      nw=n+4
      np=n+1       
cc      der2=-0.5*(4.0*a(n-1)-3.0*a(n) -a(n-2))/ds
      der2=ad(n)
cc      der1=0.5*(4.0*a(1)-3.0*a(0)-a(2))/ds
      der1=ad(0)
      call bsplw(a,w(1),np,ds,der1,der2,w(nw))
      do 200 i=0,na
       a(i)=w(i+1)
 200  continue

      der2=gdd
      der1=gddiff
      call bsplw(ad,w(1),np,ds,der1,der2,w(nw))
      do 300 i=0,na
       ad(i)=w(i+1)
 300  continue

      return
      end

******************************************************************
      subroutine WeKaji(wl,ns,nx,ny,dx,dy,h,btru,ipol,a,na,ds,wrk)
      integer ns,nx,ny,ipol,na
      real wl(ns),dx,dy,h,btru,a(0:na),ds,wrk(*)
*----------------------------------------------------------------
      integer nw,i,i0,k,iflag,nap,nam
      real dx2,y2,rh2,eps,bteff,tmp,hg
      external hg

      rh2=1.0/(h*h)
      tmp=h*sqrt(na*ds*0.5)
      if(btru.ge.tmp) then
        bteff=tmp
      else
        bteff=btru
      end if 

      dx2=dx*dx
      eps=dx2*0.0001
      nap=na+1
      nam=na-1
c     nam are number of internal points for splines

      nx=bteff/dx+1
      ny=bteff/dy+1

      if(ny*nx.gt.ns) then
        write(0,*)'nx*ny>ns',nx,ny,ns
        nx=1
        ny=1
      end if

      do 100 k=1,ny
       y2=(k-1)*(k-1)*dy*dy
       do 50 i=1,nx
        wrk(i)=rh2*((i-1)*(i-1)*dx2+y2)
  50  continue

      i0=(k-1)*nx+1

      If(ipol.le.1) then
      call  hregulin(a,nap,ds,0.0,wrk,wl(i0),nx,.false.,0.0,0.0,
     %                   eps,iflag)
      else
       do 60 i=1,nx
         wl(i0+i-1)=hg(wrk(i),0.0,ds,a,nam)
 60    continue
      end if     
 100  continue
      return
      end


******************************************************************
*
*     Computes Kajiuara weights from a point source on an independent grid
*     wl  - weights (Green function, inclusive all h effects)       O
*     ns  - maximum number of weights                     I
*     nx,ny -  size of weight grid                        O
*     ia,ka - lower left point, meaning that wl(1) is     O
*             weight for i=ia, k=ka
*     rf - correction factor to secure volume conservation      O
*     n,m,dx,dy,x0,y0  - specify external grid            I
*     xs, ys - location of source                         I
*     h      -   depth                                    I
*     btru   - extent in x and y of dependence relative to     I
*              h 
*      ipol  - interpolation method                        I
*      a,na,ds  - tabulated Green function            I
*      wrk - work array                               D
*      ierr  - error parameter. Ierr=0 is OK          O  
***************************************************************     
      subroutine GenWK(wl,ns,nx,ny,ia,ka,rf,n,m,dx,dy,x0,y0,
     %      xs,ys,h,btru,ipol,a,na,ds,ierr,wrk)
      integer ns,nx,ny,ia,ka,n,m,ipol,na,ierr
      real wl(ns),rf,dx,dy,x0,y0,xs,ys,h,btru,a(0:na),ds,wrk(*)
*----------------------------------------------------------------
      integer i,k,ib,kb,nap,nam,ntot
      real xdis,ydis,y2,eps,bteff,rh2,sum,tmp,hg
      external hg

      ierr=0

      rh2=1.0/(h*h)
      tmp=h*sqrt(na*ds*0.5)
      if(btru.ge.tmp) then
        bteff=tmp
      else
        bteff=btru
      end if 

      ia=(xs-bteff-x0)/dx+1
      ia=max(ia,1)
      ib=(xs+bteff-x0)/dx+2
      ib=min(ib,n)

      ka=(ys-bteff-y0)/dy+1
      ka=max(ka,1)
      kb=(ys+bteff-y0)/dy+2
      kb=min(kb,m)

      eps=dx*dx*0.0001
      nap=na+1
      nam=na-1
c     nam are number of internal points for splines

      nx=ib-ia+1
      ny=kb-ka+1
      ntot=nx*ny
      if(ntot.gt.ns) then
        write(0,*)'nx*ny>ns',nx,ny,ns
        ierr=1
        return
      end if

      ydis=y0+(ka-2)*dy-ys
      do 100 k=1,ny
       ydis=ydis+dy
       y2=ydis*ydis
       xdis=x0+(ia-2)*dx-xs
       do 50 i=1,nx
        xdis=xdis+dx
        wrk(i+(k-1)*nx)=rh2*(xdis*xdis+y2)
  50  continue
 100  continue

      ntot=nx*ny

      If(ipol.le.1) then
        call  regran(a,nap,ds,0.0,wrk,wl(1),ntot,2,0.0,0.0,
     %                   eps,ierr)
cc       do 55 i=1,ntot
cc            call  hregulin(a,nap,ds,0.0,wrk(i),wl(i),1,2,0.0,0.0,
cc     %                   eps,ierr)
cc 55   continue
      else
       do 60 i=1,ntot
         wl(i)=hg(wrk(i),0.0,ds,a,nam)
 60    continue
      end if     

      sum=0.0
      do 70 i=1,ntot
         wl(i)=wl(i)*rh2
         sum=sum+wl(i)
 70   continue

      if(sum.gt.0.0) then
       sum=sum*dx*dy
       rf=1.0/sum
      else
       ierr=12
c     c       write(0,*)'Negative sum in GenWK', ntot,sum,ds,nam,a(2)
        write(0,*)'Negative sum in GenWK', rh2, nx,ny, wl(20), wrk(20)
       return
      end if

      return
      end


********************************************************************************
*
*         Sums up the volume defined by a and dx,dy. Taking half the 
*         contribution from lower and left boundaries, and a quarter from the 
*         lower-left corner point
************************************************************************
      subroutine arrint(a,n,m,dx,dy,sum)
      integer n,m
      real a(n,m),dx,dy,sum
*----------------------------------------------------
      integer i,k
      real sx,sy,sxy

      sy=0.0
      sxy=0.0
      do 100 k=2,m
        sy=sy+a(1,k)
        do 50 i=2,n
          sxy=sxy+a(i,k)
 50     continue
 100  continue
 
      sx=0.0
      do 200 i=2,n
          sx=sx+a(i,1)
 200  continue

      sum=dx*dy*(sxy+0.5*(sx+sy)+0.25*a(1,1))
      return
      end


********************************************************************************
*         Expanding a quarter array to a full one
*         Sums up the volume defined by b and dx,dy. 
************************************************************************
      subroutine arrexp(a,b,n,m,dx,dy,sum)
      integer n,m
      real a(n,m),b(1-n:n-1,1-m:m-1),dx,dy,sum
*----------------------------------------------------
      integer i,k,im,km
      real sx,sy,sxy,aval
 
      do 40 k=1,m
      km=k-1
        do 20 i=1,n
          im=i-1
          aval=a(i,k)
          b(im,km)=aval
          b(im,-km)=aval
          b(-im,km)=aval
          b(-im,-km)=aval
 20     continue
 40   continue

      sy=0.0
      sxy=0.0
      do 100 k=2,m
        sy=sy+a(1,k)
        do 50 i=2,n
          sxy=sxy+a(i,k)
 50     continue
 100  continue
 
      sx=0.0
      do 200 i=2,n
          sx=sx+a(i,1)
 200  continue

      sum=dx*dy*(sxy+0.5*(sx+sy)+0.25*a(1,1))
      return
      end


********************************************************************************
*        
************************************************************************
      subroutine asmuns(d,sm,n,m,b,dx,dy,nm,mm,rs)
      integer n,m,nm,mm
      real d(n,m),sm(n,m),b(-nm:nm,-mm:mm),dx,dy,rs
*----------------------------------------------------
      integer i,k,na,nb,ma,mb,ii,kk,im,km
      real deff
 
      do 200 k=1,m
        ma=-min(mm,k-1)
        mb=min(mm,m-k)
        do 100 i=1,n
          deff=d(i,k)*rs*dx*dy
          na=-min(nm,i-1)
          nb=min(nm,n-i)
          do 40 kk=ma,mb
           do 30 ii=na,nb
             sm(i+ii,k+kk)=sm(i+ii,k+kk)+deff*b(ii,kk)
 30        continue   
 40      continue     

 100    continue
 200  continue

      return
      end


********************************************************************************
*        dx, dy, dxg,dyg - must be Cartesian
*        if ipb=2, then h must contain spline coeff.
************************************************************************
      subroutine tranht(d,sm,hset,n,m,dx,dy,x0,y0,il,ir,kl,kr,
     % h,ng,mg,dxg,dyg,x0g,y0g,a,naa,ds,iph,ipb,b,blf,hmin,rs,wrk)
      integer n,m,il,ir,kl,kr,ng,mg,naa,iph,ipb
      real d(n,m),sm(1-il:n+ir,1-kl:m+kr),hset(n,m),b(100),h(ng,mg)
      real a(naa+2),wrk(*)
      real dx,dy,x0,y0,ds,dxg,dyg,x0g,y0g,blf,hmin,rs
*----------------------------------------------------
      integer i,k,na,nb,ma,mb,ii,kk,im,km,nx,ny,nn,mm,nsw,ipos
      integer iu,ku,nnt,mmt
      parameter(nsw=1000000)
      real deff,hv,bl,sum,eps,dxn,dyn,wl(nsw)
         

      call null(sm,1,n,1,m)
cc    compute hset
      eps=0.001*dx
      call ggstagv(h,ng,mg,dxg,dyg,x0g,y0g,hset,n,n,m,dx,dy,x0,y0
     %                  ,1.0,0.0,iph,eps,hmin)

      iu=n/2
      ku=m/2
           write(0,*)'n,m =',n,m
      do 200 k=1,m
        
        do 100 i=1,n
          hv=hset(i,k)
          if(hv.lt.hmin)hv=hmin
          dxn=dx/hv
          dyn=dy/hv
          bl=blf*hv
          call WeKaji(wl,nsw,nx,ny,dx,dy,hv,bl,ipb,a,naa,ds,wrk)
          nn=nx-1
          mm=ny-1
          call arrexp(wl,wrk,nx,ny,dx,dy,sum)
          rs=0.25*hv*hv/sum
          deff=d(i,k)*rs*dxn*dyn
          na=-min(nn,i-1+il)
          nb=min(nn,n-i+ir)
          ma=-min(mm,k-1+kl)
          mb=min(mm,m-k+kr)

          if(i.eq.iu.and.k.eq.ku) then
           nnt=1+2*nn
           mmt=1+2*mm
           write(0,*)'nnt,nnt,mmt,dx,dy', nnt,nnt,mmt,dx,dy 
cc           call gphovpri(wrk,nnt,nnt,mmt,dx,dy,0.0,0.0,0,'wf!',1,.true.)
           write(0,*)'hv,rs =',hv,rs
          end if          

          do 40 kk=ma,mb
           do 30 ii=na,nb
             ipos=ii+nn+1 +(2*nn+1)*(kk+mm)
             sm(i+ii,k+kk)=sm(i+ii,k+kk)+deff*wrk(ipos)
 30        continue   
 40      continue     

 100    continue
 200  continue

      return
      end




********************************************************************************
*        dx, dy, dxg,dyg etc - must be Cartesian
*        if ipb=2, then h must contain spline coeff.
*        uf is a factor for the source
************************************************************************
      subroutine gtraht(d,hset,n,m,dx,dy,x0,y0,
     % h,ng,mg,dxg,dyg,x0g,y0g,sm,nr,mr,dxr,dyr,x0r,y0r,
     % a,naa,ds,iph,ipb,uf,blf,hmin,irs,wl,rkv,ns,ierr)
      integer n,m,ng,mg,nr,mr,naa,irs,iph,ipb,ns,ierr
      real d(n,m),hset(n,m),sm(nr,mr),h(ng,mg)
      real a(naa+2),wl(ns),rkv(ns)
      real dx,dy,x0,y0,ds,dxg,dyg,x0g,y0g,uf,blf,hmin
      real dxr,dyr,x0r,y0r
*----------------------------------------------------
      integer i,k,ia,ka,nx,ny,ipos,ii,kk,kam,iam,ntmp
      real deff,rs,hv,eps,xs,ys,tmp,x0u,y0u
         

cc      call null(sm,1,nr,1,mr)
cc    compute hset
      eps=0.001*dx
      call ggstagv(h,ng,mg,dxg,dyg,x0g,y0g,hset,n,n,m,dx,dy,x0,y0
     %                  ,1.0,0.0,iph,eps,hmin)


       do 200 k=1,m
        ys=y0+(k-1)*dy

        do 100 i=1,n
          hv=hset(i,k)
          xs=x0+(i-1)*dx
         
          if(hv.lt.hmin)hv=hmin
cc           write(0,*)'GenWK hv=', hv
c          call WeKaji(wl,nsw,nx,ny,dx,dy,hv,bl,ipb,a,naa,ds,wrk)
          call GenWK(wl,ns,nx,ny,ia,ka,rs,nr,mr,dxr,dyr,x0r,y0r,
     %      xs,ys,hv,blf,ipb,a,naa,ds,ierr,rkv)

C#ABgsfrsdf
          iam=ia-1
          kam=ka-1
          deff=uf*d(i,k)*(1+irs*(rs-1))*dx*dy

          do 40 kk=1,ny
           do 30 ii=1,nx
             ipos=ii+(kk-1)*nx
             sm(iam+ii,kam+kk)=sm(iam+ii,kam+kk)+deff*wl(ipos)
 30        continue   
 40      continue     

 100    continue
 200  continue

      return
      end



******************************************************************
*
*     Computes a Kajiuara dipole on an independent grid
*     wl  - the response                 O
*     n,m,dx,dy,x0,y0  - specify external grid            I
*     xs, ys - location of source                         I
*     h      -   depth                                    I
*     vol - volume of moving body                         I
*      ipol  - interpolation method                        I
*      da,na,ds  - tabulated Green function differentiate with respect to r**2            I
*      wrk - work array                               D
*      ierr  - error parameter. Ierr=0 is OK          O  
***************************************************************     
      subroutine Kdip(wl,n,m,dx,dy,x0,y0,
     %      xs,ys,h,vol,ipol,da,na,ds,ierr,wrk)
      integer n,m,ipol,na,ierr
      real wl(n*m),dx,dy,x0,y0,xs,ys,h,vol,da(0:na),ds,wrk(*)
*----------------------------------------------------------------
      integer i,k,ntot,nam,nap,ipos
      real xdis,ydis,y2,eps,rh2,sum,tmp,hg
      external hg

      ierr=0

      rh2=1.0/(h*h)
      
      ntot=n*m
      eps=dx*dx*0.0001

      nam=na-1
      nap=na+1
c     nam are number of internal points for splines


      ydis=y0-dy-ys
      do 100 k=1,m
       ydis=ydis+dy
       y2=ydis*ydis
       xdis=x0-dx-xs
       do 50 i=1,n
        xdis=xdis+dx
        wrk(i+(k-1)*n)=rh2*(xdis*xdis+y2)
  50  continue
 100  continue

      ntot=n*m

      If(ipol.le.1) then
cc    will not work because wrk is not monotoneous
cc      call  hregulin(da,nap,ds,0.0,wrk,wl(1),ntot,.false.,0.0,0.0,
cc     %                   eps,ierr)
      else
       do 60 i=1,ntot
         wl(i)=hg(wrk(i),0.0,ds,da,nam)
 60    continue
      end if     

      sum=0.0
      do 75 k=1,m
      do 70 i=1,n
         xdis=x0+(i-1)*dx -xs
         ipos=i+(k-1)*n 
         wl(ipos)=2.0*vol*wl(ipos)*rh2*rh2*xdis
         sum=sum+wl(ipos)
 70   continue
 75   continue



      return
      end
******************************************************
*
*     Routine computes the parallel and normal components to
*     a line xa,ya --> xb,yb. The tangential and normal unit
*     vectors to the line is organized in counterclockwise
*     sequence.  
*        xa,ya,xb,yb  -- define the line                     I
*        xp,yp  -- define the point s,n                      I
*        slin   --  the value of parallel coordinate corresponding
*                   to xb,yb                                 O
*        s,n    --  tangential and normal components         O
*********************************************************
      subroutine linedist(xa,ya,xb,yb,xp,yp,slin,s,n)
      real xa,ya,xb,yb,xp,yp,slin,s,n
*--------------------------------------------------------------------
      real cs,ss,cn,sn,dx,dy

      slin=sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya))
      cs=(xb-xa)/slin
      ss=(yb-ya)/slin
      cn=-ss
      sn=cs

      dx=xp-xa
      dy=yp-ya

      s=cs*dx+ss*dy
      n=cn*dx+sn*dy

      return
      end

******************************************************
*
*     Routine that computes the minimum distance from a point
*     to a line xa,ya --> xb,yb. 
*        xa,ya,xb,yb  -- define the line                     I
*        xp,yp  -- define the point s,n                      I
*        dmin   --  minimum distance                         O
*        rmin    -- minimum distance relative to length of line O
*********************************************************
      subroutine mindist(xa,ya,xb,yb,xp,yp,dmin,rmin)
      real xa,ya,xb,yb,xp,yp,dmin,rmin
*--------------------------------------------------------------------
      real slin,s,n 


      call linedist(xa,ya,xb,yb,xp,yp,slin,s,n)
      dmin=sqrt(s*s+n*n)
      if(s.lt.0.0)then
        dmin=sqrt(s*s+n*n)
      else
        if(s.le.slin)then
          dmin=abs(n)
        else
         dmin=sqrt((s-slin)*(s-slin)+n*n)
        end if
      end if
 
      rmin=dmin/slin

      return
      end


******************************************************
*
*     Finds intersection point of two lines
*     Routine computes the parallel and normal components to
*     a line xa,ya --> xb,yb. The tangential and normal unit
*     vectors to the line is organized in counterclockwise
*     sequence.  
*        xa,ya,xb,yb  -- define first line                     I
*        x1,y1,x2,y2  -- define second line                      I
*        sa,sc    --  tangential components for intersection points
*                     along the line the point          O
*        dmin -- minimum distance between lines, value <=0 means that
*                lines intersect between their end-points
*        xm,ym -- if dmin=0: intersection point
*********************************************************
      subroutine intersec(xa,ya,xb,yb,x1,y1,x2,y2,sa,sc,dmin,xm,ym)
      real xa,ya,xb,yb,x1,y1,x2,y2,sa,sc,dmin,xm,ym
*--------------------------------------------------------------------
      real slin,s1,n1,s2,n2,p,tmp,sl2,rmin

      call linedist(xa,ya,xb,yb,x1,y1,slin,s1,n1)
      call linedist(xa,ya,xb,yb,x2,y2,slin,s2,n2)

      if(n2.eq.n1) then
        dmin=-100
        return
      end if

c     p is the weight of p*(x1,y1)+(1-p)*(x2,y2) that yields zero n
      p=n2/(n2-n1)
      sa= p*s1+(1.0-p)*s2
      sl2=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))
      sc=(1-p)*sl2
      if(n1*n2.le.0.0. and. sa.le.slin.and.sa.ge.0.0) then
        dmin=0.0
        xm=p*x1+(1.0-p)*x2
        ym=p*y1+(1.0-p)*y2
      else
        call mindist(xa,ya,xb,yb,x1,y1,dmin,rmin)
        call mindist(xa,ya,xb,yb,x2,y2,tmp,rmin)
        if(tmp.le.dmin)dmin=tmp
        call mindist(x1,y1,x2,y2,xa,ya,tmp,rmin)
        if(tmp.le.dmin)dmin=tmp
        call mindist(x1,y1,x2,y2,xb,yb,tmp,rmin)
        if(tmp.le.dmin)dmin=tmp
      end if

      return
      end



******************************************************
*
*     Finds intersection points between a curve (polygon defined be given
*     points) and the (extension)  of a given line.
*        x,y  -- arrays that define the curve                   I
*          n - number of points                                  I
*        x1,y1,x2,y2  -- define line                             I
*        ncmax -- maximum number of crossings allowed            I
*        ncr  -- number of crossings found for extension of line O
*        isl  --  array for segments that are crossed            O
*                 segment i corresponds to (x,y)_i .. (x,y)_i+1     
*        sl   --  the value of parallel coordinate relative
*                   to segment of curve. value =0 is x_i, 
*                    value 1 is x_i+1                          O
*        sa       --  tangential components for intersection points
*                     along the line, defined as sl; values in the interval
*                     0,1 means that an intersection is found           O
*        dm   -- minimum distance from line, value <=0 means that
*                lines intersect between their end-points               O
*        eps -- tolerance for intersection detection                    I
*********************************************************
      subroutine geosec(x,y,n,x1,y1,x2,y2,ncmax,ncr,isl,sl,sa,dm,eps)
      integer n,ncmax,ncr,isl(ncmax)
      real x(n),y(n),x1,y1,x2,y2,sl(ncmax),sa(ncmax),dm(n),eps
*--------------------------------------------------------------------
      integer iseg,nm
      real xa,ya,xb,yb,sisg,ssc,xm,ym,sseg,slin
      logical prevright

      nm=n-1
      slin=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))
      ncr=0
      xb=x(1)
      yb=y(1)
      prevright=.false.

      do 100 iseg=1,nm
        xa=xb
        ya=yb
        xb=x(iseg+1)
        yb=y(iseg+1)
        sseg=sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya))
        call  intersec(xa,ya,xb,yb,x1,y1,x2,y2,sisg,ssc,dm(iseg),xm,ym)
        sisg=sisg/sseg
        ssc=ssc/slin

cc        write(*,*)'Hit',sisg,ncr
cc
c       handle crossing close to point iseg
        if(sisg.ge.-eps.and.sisg.le.1.0+eps.and.dm(iseg).gt.-99.0) then
cc        if(sisg.ge.-eps.and.sisg.le.1.0+eps) then
cc        write(*,*)'Exept'
        if( abs(sisg).le.eps) then
          if(.not.prevright) then
            if(ncr.eq.ncmax) then
               write(0,*)'WARNING to many cross in geosec'
               return
            end if
            ncr=ncr+1
            isl(ncr)=iseg
            if(sisg.lt.0.0) sisg=0.0
            sl(ncr)=sisg
            sa(ncr)=ssc
          else
            prevright=.false.
          end if
        else
          prevright=.false.
          if( abs(sisg).le.1.0+eps) then
            if(ncr.eq.ncmax) then
               write(0,*)'WARNING to many cross in geosec'
               return
            end if
            ncr=ncr+1
            isl(ncr)=iseg
            prevright=sisg.gt.1.0-eps 
            if(sisg.gt.1.0) sisg=1.0
            sl(ncr)=sisg
            sa(ncr)=ssc
          end if
        end if
        end if
 100  continue

      return
      end
      
******************************************************
*
*     Finds intersection points between a curve (polygon defined be given
*     points) and the (extension)  of a grid.
*        x,y  -- arrays that define the curve                   I
*          n - number of points                                  I
*        ng,mg  -- size of grid                             I
*        dx,dy,x0,y0 -- increments and lower left corner of grid  I
*        rot -rotation angle in degrees (around x0,y0)            I
*        nsc  -- number of crossings found for extension of 
*                grid lines  First index 1...ng corresponds to
*                columns, ng+1...ng+mg to rows                   O
*        isg  --  array for segments that are crossed            O
*                 segment i corresponds to (x,y)_i .. (x,y)_i+1     
*        ssg   --  the value of parallel coordinate relative
*                   to segment of curve. value =0 is x_i, 
*                    value 1 is x_i+1                          O
*        aar       --  tangential components for intersection points
*                     along the grid lines. Values in the interval
*                     0,1 means that an intersection is found           O
*                lines intersect between their end-points               O
*        eps -- tolerance for intersection detection                    I
*        w  -- work space of dimension n                                D
*********************************************************
      subroutine gridsec(x,y,n,ng,mg,dx,dy,x0,y0,rot,nsc,isg,ssg,aar,
     %                 eps,w)
      integer n,ng,mg,ncal,nsc(ng+mg),isg(ng+mg,10)
      real x(n),y(n),dx,dy,x0,y0,rot,ssg(ng+mg,10)
      real aar(ng+mg,10),eps,w(n)
*----------------------------------------------------------------
      integer i,j,itot,ncr,isl(10)
      real pi,cc,ss,x1,y1,x2,y2,tv1,sv1,tv2,sv2,sl(10),sa(10)

      pi=3.14159265358979

      cc=cos(rot*pi/180)
      ss=sin(rot*pi/180)

      itot=ng+mg

      do 100 i=1,itot
        if(i.le.ng) then
          tv1=0.0
          sv1=(i-1)*dx
          tv2=(mg-1)*dy
          sv2=sv1
        else
          tv1=(i-1-ng)*dy
          sv1=0.0
          tv2=tv1
          sv2=(ng-1)*dx
        end if
  
        x1=x0+cc*sv1-ss*tv1
        y1=y0+ss*sv1+cc*tv1            
        x2=x0+cc*sv2-ss*tv2
        y2=y0+ss*sv2+cc*tv2            

        call geosec(x,y,n,x1,y1,x2,y2,10,ncr,isl,sl,sa,w,eps)
        nsc(i)=ncr
        do 50 j=1,ncr
         isg(i,j)=isl(j)
         ssg(i,j)=sl(j)
         aar(i,j)=sa(j)
        
 50     continue

 100  continue

      return
      end

********************************************************
*    NGAUSET
*
*    Returns weights and points for Gaussian integration on the unit
*    interval divided into subintervals.
*         aw, xr - weights and positions. Number: mh*mult          O
*         mh, mult - number of Guassian points per subinterval and number 
*                    of subintervals                                      I 
*
*********************************************************
      subroutine ngauset(aw,xr,mh,mult)
      integer mh,mult
      real xr(mh*mult), aw(mh*mult)
*-------------------------------------
      integer i,k,k0
      real w(40),g(40),ds,x0

      if(mh.eq.1) then
        w(1)=1.0
        g(1)=0.5
      end if
      if(mh.eq.2) then
       w(1)=0.5
       w(2)=0.5
       g(1)=0.5*(1.0-sqrt(2.0/3.0))
       g(2)=0.5*(1.0+sqrt(2.0/3.0))
      end if
      if(mh.eq.3) then
       w(1)=5.0/18.0
       w(3)=w(1)
       w(2)=8.0/18.0
       g(1)=0.5*(1.0-sqrt(3.0/5.0))
       g(3)=0.5*(1.0+sqrt(3.0/5.0))
       g(2)=0.5
      end if

      if(mh.eq.4) then
       w(1)=0.1739274225687269
       w(2)=0.3260725774312730
       w(3)=w(2)
       w(4)=w(1)

       g(1)=0.0694318442029738
       g(2)=0.3300094782075719 
       g(3)=1.0-g(2)
       g(4)=1.0-g(1)
      end if
      if(mh.eq.5) then
       w(1)=0.1184634425280945
       w(2)=0.2393143352496832
       w(3)=0.2844444444444444
       w(4)=w(2)
       w(5)=w(1)
       g(1)=0.0469100770306681
       g(2)=0.2307653449471585
       g(3)=0.5
       g(4)=1.0-g(2)
       g(5)=1.0-g(1)
      end if
      if(mh.eq.6) then
       w(1)=0.0856622461895852
       w(2)=0.1803807865240693
       w(3)=0.2339569672863455
       w(4)=w(3)
       w(5)=w(2)
       w(6)=w(1)
       g(1)=0.0337652428984240 
       g(2)=0.1693953067668678  
       g(3)=0.3806904069584015 
       g(4)=1.0-g(3)
       g(5)=1.0-g(2)
       g(6)=1.0-g(1)
      end if

      include 'gausadd.inc'
cc      do 60 i=1,mh
cc        write(*,*)'i,w,g',i,w(i),g(i)
cc 60   continue
      ds=1.0/mult

      do 100 i=1,mult
       x0=(i-1)*ds
       k0=(i-1)*mh
       do 50 k=1,mh
        xr(k0+k)=x0+ds*g(k)
        aw(k0+k)=ds*w(k)
 50    continue
 100  continue

      return
      end


************************************************************
*
*   Qaudint
*
****************************************************************
      function quadint(f,a,b,x,w,n)
      integer n
      real f,a,b,x(n),w(n),quadint
      external f
*-----------------------------------------------------------------
      integer i
      real fval,xval,sum,dx

      dx=b-a
      sum=0.0

      do 100 i=1,n
        xval=a+dx*x(i)
        sum=sum+w(i)*f(xval)
 100  continue

      quadint=sum*dx

      return
      end 
******************************************************************
      subroutine crossing(y,n0,n,ia,icross,val,up)
      integer n0,n,ia,icross
      real  y(n0:n),val
      logical up
*----------------------------------------------------------------
      integer i,sgn
      real grens

      if(up) then
        sgn=1
      else
        sgn=-1
      end if

      grens=val*sgn

      icross=ia

 100  continue
      if(sgn*y(icross).lt.grens .and. icross.lt.n) then
        icross=icross+1
        go to 100
      end if
      
      if(sgn*y(icross).lt.grens) icross=n+1
 
      return
      end 

******************************************************************
      subroutine nestemax(y,n0,n,ia,ires,ymin,maks)
      integer n0,n,ia,ires
      real  y(n0:n),ymin
      logical maks
*----------------------------------------------------------------
      integer i,sgn
      real grens,rmm,rm,rp

      if(maks) then
        sgn=1
      else
        sgn=-1
      end if

      grens=ymin*sgn
      rm=sgn*y(ia)
      rp=sgn*y(ia+1)


       do 100 i=ia+2,n
       rmm=rm
       rm=rp
       rp=sgn*y(i)
       if (rmm.lt.rm .and. rp.lt.rm .and.rm.gt.grens) then
          ires=i-1
          ymin=y(i)
          return
       end if
 100  continue

      ymin=-666
      ires=n+1

      return
      end 


*************************************************************
      subroutine crmber(r,na,n,dt,imax,rmax,tkorr)
      integer na,n,imax
      real r(na:n),dt,rmax,tkorr
*----------------------------------------------------------------
      real rmm,rm,rp,d1,d2

      if(imax.le.na.or.imax.ge.n) then
        rmax=-666
        tkorr=-666
        return
      end if   


      rmm=r(imax-1)
      rm=r(imax)
      rp=r(imax+1)
      d2=(rmm+rp-2*rm)/(dt*dt)
      d1=0.5*(rp-rmm)/dt
c         rkorr=rm-0.5*d1*d1/d2
      tkorr=-d1/d2
      return
      end 


***********************************************************************
* Finds minum and maximum of a in interval
*
***********************************************************************
      subroutine iekstrem(a,na,n,ia,ib,imin,imax,amin,amax)
      integer na,n,ia,ib,imax,imin
      real a(na:n),amax,amin
*------------------------------------------------------
      integer i

      amax=a(ia)
      amin=a(ia)
      imin=ia
      imax=ia

      do 100 i=ia+1,ib
      if(a(i).gt.amax) then
        imax=i
        amax=a(i)
      end if
      if(a(i).lt.amin) then
         imin=i
         amin=a(i)
      end if
 100  continue


      return
      end

*****************************************************************************************
      subroutine stadev(y,med,n,nber,ymid,ydev)
      integer n,nber
      real y(n),ymid,ydev
      logical med(n)
*----------------------------------------------------------------------
      integer i
      real sy2,sy

      nber=0
      sy=0.0
      sy2=0.0

      do 100 i=1,n
        if(med(i)) then
          nber=nber+1
          sy=sy+y(i)
          sy2=sy2+y(i)*y(i)         
        end if
 100  continue

      if(nber.gt.0) then
        ymid=sy/nber
        if(nber.ge.2) then
         ydev=sqrt( (sy2-ymid*ymid*nber)/(nber-1))
        else
         ydev=-666
        end if
      else
        ymid=-666
        ydev=-666
      end if      

      return
      end
***************************************************************************
*
*              D I G I T
*
*     Unders|ker om karakteren 'c' er et tall
****************************************************************************
      function digit(c)
      character c
      logical digit
*-------------------------------------------------------------------------
      integer idiff                       

      idiff=ichar(c)-ichar('0')
      digit= idiff.ge.0  .and.  idiff.le.9
      return
      end
         

********************************************************************
*
*                 i t c o n v
*
*     Konverterer heltallet 'i' til teksten 't'. 'kf' angir lengden
*     av den konverterte teksten. 
***********************************************************************
      subroutine itconv(i,kf,t)
      integer i,kf
      character*20 t
*---------------------------------------------------------------------
      integer irest,isiff,inull,ichar,isgn,j,jstop
      character c,char

      inull=ichar('0')

      if(i.ge.0) then
        isgn=1
        irest=i
      else
        isgn=-1
        irest=-i
      end if

      kf=0

 100  continue
      kf=kf+1
      isiff=mod(irest,10)
      t(kf:kf)=char(isiff+inull)
      if(kf.eq.20) go to 300
      irest=irest/10
      if(irest.eq.0) go to 200
      go to 100

 200  continue
      if(isgn.lt.0) then
        kf=kf+1
        t(kf:kf)='-'
      end if

 300  continue
 
      jstop=kf/2
      do 400 j=1,jstop
      c=t(j:j)
      t(j:j)=t(kf+1-j:kf+1-j)
      t(kf+1-j:kf+1-j)=c
 400  continue
      
      return
      end
     
*****************************************************************************
*
*                   L O K O R D
*
*    Lokaliserer foerste sekvens i strengen c(n0:n1) og angir den eventuelle
*    posisjonen c(i0:i1). Med sekvens menes en sammenhengende serie ikke-
*    blanke tegn. 
*        iflag=0   : strengen er tom.
*              1   : sekvens som begynner med bokstav funnet. 
*              2   : sekvens som begynner med tall funnet. 
*              3   : sekvens som begynner med annet tegn funnet. 
*
*
**************************************************************************
      subroutine lokord(c,n0,n1,i0,i1,iflag)
      integer n0,n1,i0,i1,iflag
      character c*80
*--------------------------------------------------------------------
      integer i
      character b

      if(n0.gt.n1) then
        iflag=0
        return
      end if

      i=n0
      b=c(i:i)                

  100 continue
      if(b.ne.' ') go to 200
      if(i.eq.n1) then
        iflag=0
        return
      end if
      i=i+1
      b=c(i:i)
      go to 100
  200 continue

      if((b.ge.'a'.and.b.le.'z').or.(b.ge.'A'.and.b.le.'Z'))then
        iflag=1
      else
        if(b.ge.'0'.and.b.le.'9')then
          iflag=2
        else
          iflag=3
        end if
      end if
      i0=i

                                 
  300 continue
      if(b.eq.' ') go to 400
      if(i.eq.n1) then
        i1=n1
        return
      end if
      i=i+1
      b=c(i:i)
      go to 300

  400 continue
      i1=i-1

      return
      end                                                

****************************************************************
*
*                 B L A N K
*
*     C(1:K) fylles med blanke
*
*****************************************************************          
      SUBROUTINE BLANK(C,K)
      INTEGER K
      CHARACTER*80 C
*--------------------------------------------------------------
      INTEGER I

      DO 100 I=1,K
  100 C(i:i)=' '

      RETURN
      END



****************************************************************             
*
*              U P C H A R
*    Gjoer om smaa bokstaver i c(k1:k2) til store.
****************************************************************
      subroutine upchar(c,k1,k2)
      character*80 c
      integer k1,k2
*------------------------------------------------------------------
      integer ila,ilz,isa,ks,iver,i
                        
      ila=ichar('a')
      ilz=ichar('z')
      isa=ichar('A')
      ks=min(80,k2)
      do 100 i=k1,k2
      iver=ichar(c(i:i))                             
      if(iver.ge.ila .and. iver.le.ilz)
     %   c(i:i)=char(isa+iver-ila)
  100 continue
      return
      end


*********************************************************************
*                 L O K I N T                                       *
*   Lokaliserer foerste heltall i strengen c(n0:n1).Ved output      *
*   angir i0 posisjon av foerste ikke-blanke tegn, et eventuelt     *
*   heltall er lokalisert til c(i0:i1). Kontroll parameter:         *
*     iflag=0    : tekst-streng er tom.                             *
*          =1    : ledende heltall etterfulgt av blank er funnet.   *
*          =2    : ledende heltall etterfulgt av  ikke-blank.       *
*          =3    : ensomt fortegn funnet i posisjon n1.             *
*          =4    : foerste ikkebl. tegn ingen mulig del av heltall. *
*                                                                   *
*********************************************************************
      subroutine lokint(c,n0,n1,i0,i1,iflag)
      integer n0,n1,i0,i1,iflag
      character c*80
*--------------------------------------------------------------------
      integer i
      character b

      if(n0.gt.n1) then
        iflag=0
        return
      end if

      i=n0
      b=c(i:i)                

  100 continue
      if(b.ne.' ') go to 200
      if(i.eq.n1) then
        iflag=0
        return
      end if
      i=i+1
      b=c(i:i)
      go to 100
  200 continue
      i0=i

      if(b.eq.'+' .or. b.eq.'-') then
        if(i.eq.n1) then
          iflag=3
          return
        end if
        i=i+1
        b=c(i:i)
      end if

      if(b.lt.'0' .or. b.gt.'9') then
        iflag=4
        return
      end if
                                 
  300 continue
      if(b.lt.'0' .or. b.gt.'9') go to 400
      if(i.eq.n1) then
        i1=n1
        iflag=1
        return
      end if
      i=i+1
      b=c(i:i)
      go to 300

  400 continue
      i1=i-1
      if(b.eq.' ') then
        iflag=1
      else
        iflag=2
      end if

      return
      end




**************************************************************************
*                                                                        *
*                   I C O N                                              *
*                                                                        *
*    Funksjon som regner om tekst t(k1:k2) til et heltall. Det m} p}     *
*  forh}nd v{re testet at alle enkeltbokstaver i sekvensen er tall       *
*  (bortsett fra fortegn) og at antallet ikke er ulovlig stort.          *
**************************************************************************

      function icon(t,k1,k2)
      integer icon,k1,k2
      character*80 t
*--------------------------------------------------------------------------
      integer i,l,n10,null,isum,kk1,isgn

      null=ichar('0')
      n10=1     
      isum=0                       
      
      if(t(k1:k1).eq.'-') then
        kk1=k1+1
        isgn=-1                          
      else                
        isgn=1
        if(t(k1:k1).eq.'+') then
          kk1=k1+1               
        else
          kk1=k1 
        end if
      end if
         
      do 100 i=kk1,k2
      l=k2 -i + kk1 
      isum=isum+(ichar(t(l:l))-null)*n10
      n10=n10*10
  100 continue
      icon=isgn*isum
      return
      end


**************************************************************************
*                                                                        *
*                   I C O NL                                              *
*                                                                        *
*    Funksjon som regner om tekst t(k1:k2) til et heltall. Det m} p}     *
*  forh}nd v{re testet at alle enkeltbokstaver i sekvensen er tall       *
*  (bortsett fra fortegn) og at antallet ikke er ulovlig stort.          *
**************************************************************************

      function iconl(t,k1,k2)
      integer*8 iconl
      integer k1,k2
      character*80 t
*--------------------------------------------------------------------------
      integer i,l,null,kk1,isgn
      integer*8 isum,n10

      null=ichar('0')
      n10=1     
      isum=0                       
      
      if(t(k1:k1).eq.'-') then
        kk1=k1+1
        isgn=-1                          
      else                
        isgn=1
        if(t(k1:k1).eq.'+') then
          kk1=k1+1               
        else
          kk1=k1 
        end if
      end if
         
      do 100 i=kk1,k2
      l=k2 -i + kk1 
      isum=isum+(ichar(t(l:l))-null)*n10
      n10=n10*10
  100 continue
      iconl=isgn*isum
      return
      end
         



*********************************************************************** 
*                                                                     *
*                 G E T I                                             *
*   Samme anvendelse som 'lokint' men returnerer med verdi av evnt.   *
*   heltall i parameter 'ires'                                        *
*                                                                     *   
*   kaller 'lokint','icon'                                                   *
***********************************************************************
      SUBROUTINE GETI(C,N0,N1,I0,I1,iRES,IFLAG)
      character*80 c
      integer n0,n1,i0,i1,ires,iflag
*------------------------------------------------
      integer icon

      call lokint(c,n0,n1,i0,i1,iflag)
      if(iflag.eq.1 .or. iflag.eq.2)then
        if(i1-i0.gt.9) then
          write(0,*)'for mange siffer >9'
          iflag=6
        else
          ires=icon(c,i0,i1)
        end if
      end if
      return
      end


*********************************************************************** 
*                                                                     *
*                 G E T L                                             *
*   Samme anvendelse som 'GETI' men gir 8 bytes ires                            *                                                                     *   
*   kaller 'lokint','icon'                                                   *
***********************************************************************
      SUBROUTINE GETL(C,N0,N1,I0,I1,iRES,IFLAG)
      character*80 c
      integer n0,n1,i0,i1,iflag
      integer*8 ires
*------------------------------------------------
      integer*8 iconl

      call lokint(c,n0,n1,i0,i1,iflag)
      if(iflag.eq.1 .or. iflag.eq.2)then
        if(i1-i0.gt.17) then
          write(0,*)'for mange siffer'
          iflag=6
        else
          ires=iconl(c,i0,i1)
        end if
      end if
      return
      end


**************************************************************************
*                                                                        *
*                   G E T R                                              *
*   Real versjon av 'geti'. 'iflag'=1,2,3,4 som for 'geti', 'iflag'=6    *
*   angir for stor eksponent.                                            *
*                                                                        *
*   kaller 'geti'                                                        *
**************************************************************************
      subroutine getr(c,n0,n1,i0,i1,res,iflag)
      character*80 c
      integer n0,n1,i0,i1,iflag
      real res
*-------------------------------------------------------------------------
      character b
      integer num,nn0,ia,i,kflag,ie,ib,isgn
      real des
                                           
      isgn=1
      call geti(c,n0,n1,i0,i1,num,iflag)
      if(iflag.eq.3 .or. iflag.eq.0) return
      b=c(i0:i0)
      if(b.eq.'-') isgn=-1
      if(iflag.eq.4) then                
        nn0=i0+1
        if(b.eq.'-'  .or.  b.eq.'+') then
          b=c(nn0:nn0)
          nn0=nn0+1   
        end if
        if(b.ne.'.') return
        res=0.0
        call geti(c,nn0,n1,ia,ib,num,kflag)
        if(kflag.gt.2. .or. kflag.eq.0 .or. ia.gt.nn0) return
        b=c(ia:ia)
      else    
        res=num*1.0
        if(iflag.ne.2) then
          return   
        end if
  
        b=c(i1+1:i1+1)
        if(b.eq.'e' .or. b.eq.'E') then
          ie=i1+1
          go to 300 
        end if

        if(b.ne.'.') return
        i1=i1+1
        nn0=i1+1
        call geti(c,nn0,n1,ia,ib,num,kflag)

        if(kflag.eq.0 .or. ia.gt.nn0) then
          iflag=1
          return
        end if  

        b=c(nn0:nn0)

        if(kflag.gt.2)  then
          if(b.eq.'e' .or. b.eq.'E') then
            ie=nn0
            go to 300
          else
            return
          end if
        end if
 
      end if

      if(b.eq.'-' .or. b.eq.'+') return

      i1=ib
      des=num*1.0*isgn
      do 100 i=ia,i1     
  100 des=des*0.1
      
      res=res+des

      if(kflag.eq.1) then
        iflag=1
        return
      end if
      
      ie=i1+1
  300 b=c(ie:ie) 
      if(b.ne.'e' .and. b.ne.'E') then
        iflag=2
        return
      end if
                               
      nn0=ie+1
      call geti(c,nn0,n1,ia,i1,num,kflag)
      if(kflag.eq.0 .or. kflag.gt.2 .or. ia.gt.nn0) then
        iflag=2
        return
      end if                      

      if(abs(num).gt.99) then
        write(0,*)'for stor eksponent'
        iflag=6
        return
      end if
   

      iflag=kflag
      res=res*10.0**(num*1.0)
      return
      end



**************************************************************************
*                                                                        *
*                   G E T R L                                              *
*   Real versjon av 'geti'. 'iflag'=1,2,3,4 som for 'geti', 'iflag'=6    *
*   angir for stor eksponent.                                            *
*                                                                        *
*   kaller 'geti'                                                        *
**************************************************************************
      subroutine getrl(c,n0,n1,i0,i1,res,iflag)
      character*80 c
      integer n0,n1,i0,i1,iflag
      real res
*-------------------------------------------------------------------------
      character b
      integer num,nn0,ia,i,kflag,ie,ib,isgn
      integer*8 numl
      real des
                                           
      isgn=1
      call geti(c,n0,n1,i0,i1,num,iflag)
      if(iflag.eq.3 .or. iflag.eq.0) return
      b=c(i0:i0)
      if(b.eq.'-') isgn=-1
      if(iflag.eq.4) then                
        nn0=i0+1
        if(b.eq.'-'  .or.  b.eq.'+') then
          b=c(nn0:nn0)
          nn0=nn0+1   
        end if
        if(b.ne.'.') return
        res=0.0
        call geti(c,nn0,n1,ia,ib,num,kflag)
        if(kflag.gt.2. .or. kflag.eq.0 .or. ia.gt.nn0) return
        b=c(ia:ia)
      else    
        res=num*1.0
        if(iflag.ne.2) then
          return   
        end if
  
        b=c(i1+1:i1+1)
        if(b.eq.'e' .or. b.eq.'E') then
          ie=i1+1
          go to 300 
        end if

        if(b.ne.'.') return
        i1=i1+1
        nn0=i1+1
        call getl(c,nn0,n1,ia,ib,numl,kflag)

        if(kflag.eq.0 .or. ia.gt.nn0) then
          iflag=1
          return
        end if  

        b=c(nn0:nn0)

        if(kflag.gt.2)  then
          if(b.eq.'e' .or. b.eq.'E') then
            ie=nn0
            go to 300
          else
            return
          end if
        end if
 
      end if

      if(b.eq.'-' .or. b.eq.'+') return

      i1=ib
      des=numl*1.0*isgn
      do 100 i=ia,i1     
  100 des=des*0.1
      
      res=res+des

      if(kflag.eq.1) then
        iflag=1
        return
      end if
      
      ie=i1+1
  300 b=c(ie:ie) 
      if(b.ne.'e' .and. b.ne.'E') then
        iflag=2
        return
      end if
                               
      nn0=ie+1
      call geti(c,nn0,n1,ia,i1,num,kflag)
      if(kflag.eq.0 .or. kflag.gt.2 .or. ia.gt.nn0) then
        iflag=2
        return
      end if                      

      if(abs(num).gt.99) then
        write(0,*)'for stor eksponent'
        iflag=6
        return
      end if
   

      iflag=kflag
      res=res*10.0**(num*1.0)
      return
      end



               

************************************************************************
*                                                                      *
*                      S T R I P                                       *
*                                                                      *
*    Finner den del, tekst(1:k), som ikke er blank. kmax angir total   *
*    lengde av tekst. Dersom tekst er tom returneres med k=0.
************************************************************************

      subroutine strip(tekst,kmax,k)
      integer k,kmax
      logical def
      character*200 tekst
*-----------------------------------------------------------------------
      integer i                                                         

      k=kmax

  100 continue
      if( tekst(k:k).ne.' ') go to 200
      if(k.eq.1) then
        k=0
        return
      end if
      k=k-1
      go to 100
  200 continue
      return
      end



************************************************************************
*                                                                      *
*                      F R O N T S T R I P                             *
*                                                                      *
*    Fjerner innledende blanke i tekst(1:kmax), kmax angir total          *
*    lengde av tekst. k er antall blanke som er funnet.
************************************************************************

      subroutine frontstrip(tekst,kmax,k)
      integer kmax
      character*200 tekst
*-----------------------------------------------------------------------
      integer k,i             

      k=1

  100 continue
      if( tekst(k:k).ne.' ') go to 200
      if(k.eq.kmax) return
      k=k+1
      go to 100
  200 continue

      do 300 i=k,kmax
      tekst(i+1-k:i+1-k)=tekst(i:i)
 300  continue

      k=k-1
      
      return
      end




************************************************************************
*                                                                      *
*                      X S T R I P                                     *
*                                                                      *
*    Finner den del, tekst(1:k), som ligger foran tegnet '#' eller '!'.*
*    def=true angir '#', def=false angir '!'. Dersom hverken '#' eller *
*    '!' finnes settes k=80 og def=true.                               *
************************************************************************

      subroutine xstrip(tekst,k,def)
      integer k
      logical def
      character*80 tekst
*-----------------------------------------------------------------------
      integer i                                                         
      character b

      k=1
      b=tekst(1:1)
  100 continue
      if( b.eq.'!' .or. b.eq.'#') go to 200
      if(k.eq.80) then
        def=.true.
        return
      end if
      k=k+1
      b=tekst(k:k)
      go to 100
  200 k=k-1
      def=.not.b.eq.'!'
      return
      end


************************************************************************
*                                                                      *
*                      X S T R I P G                                    *
*                                                                      *
*    Finner den del, tekst(1:k), som ligger foran tegnet '#' eller '!'.*
*    def=true angir '#', def=false angir '!'. Dersom hverken '#' eller *
*    '!' finnes settes k=kmax og def=true.                               *
************************************************************************

      subroutine xstripg(tekst,k,kmax,def)
      integer k,kmax
      logical def
      character*200 tekst
*-----------------------------------------------------------------------
      integer i                                                         
      character b

      k=1
      b=tekst(1:1)
  100 continue
      if( b.eq.'!' .or. b.eq.'#') go to 200
      if(k.eq.kmax) then
        def=.true.
        return
      end if
      k=k+1
      b=tekst(k:k)
      go to 100
  200 k=k-1
      def=.not.b.eq.'!'
      return
      end




************************************************************************
*
*                    D E L T E K
*
*    sekvenser i t(na:nb) lokaliseres. iant angir antall sekvenser,sekvens
*    k er lokalisert til t( l(1,k):l(2,k) ).
*
*    kaller 'lokord'
*****************************************************************************
      subroutine deltek(t,na,nb,l,iant)
      integer na,nb,l(2,20),iant
      character*80 t
*-------------------------------------------------------------------------
      integer k1,k2,iflag,n0

      n0=na
      iant=0  

  100 continue
      call lokord(t,n0,nb,k1,k2,iflag)
      if(iflag.eq.0) return
      n0=k2+1     
      iant=iant +1
      l(1,iant)=k1
      l(2,iant)=k2
      go to 100
c     return
      end        



************************************************************************     
*                      R E D
*
*    ordner teksten t(1:nb) slik at den faar form: 'sekvens1 sekvens2 ...
*     sekvensN#' ,der n er antall sekvenser.  
*
*    kaller lokord.
*************************************************************************
      subroutine red(t,nb)
      integer nb
      character*80 t
*-------------------------------------------------------------------------
      integer k1,k2,iflag,n0,kfyll,ka,iant,i

      n0=1
      kfyll=0
      ka=1                                
      iant=0

  100 continue
      call lokord(t,n0,nb,k1,k2,iflag)
      if(iflag.eq.0) then
        kfyll=kfyll+1
        if(kfyll.le.nb) t(kfyll:kfyll)='#'      
        do 75 i=kfyll+1,nb
   75   t(i:i)=' '
        return
      end if
      n0=k2+1
      iant=iant+1
      if(iant.gt.1) then
        t(kfyll+1:kfyll+1)=' '
        ka=2
      end if
      t(kfyll+ka:kfyll+ka+k2-k1)=t(k1:k2)
      kfyll=kfyll+k2-k1+ka
      go to 100
c      return
      end




********************************************************************     
*
*                 A L A C A R T
*
*   sjekker hvor mange tekster i menyen t(n) streng(na:nb) er en
*   forkortelse av. Tekstene i t(n) maa avsluttes med '#' eller '!'.
*   iant er antall treff, l(1)..l(iant) er nummerne i menyen t for
*   tekstene som er truffet. 
*
*   kaller 'fork'
***********************************************************************
      subroutine alacart(t,n,streng,na,nb,l,iant)
      integer n,na,nb,iant,l(n)
      character*80 streng,t(n)
*-------------------------------------------------------------------
      integer i   
      logical fork
      
      iant=0  
      do 100 i=1,n   
      if(fork(t(i),streng,na,nb)) then
        iant=iant+1
        l(iant)=i
      end if
  100 continue
      return
      end

********************************************************************     
*
*                 V E L G C A R T
*
*   Funksjon som alacart, men returnerer med bare ett treff dersom
*   streng(na,nb) er identisk lik en av tekstene i t.
*
*   kaller 'fork'
***********************************************************************
      subroutine velgcart(t,n,streng,na,nb,l,iant)
      integer n,na,nb,iant,l(n)
      character*80 streng,t(n)
*-------------------------------------------------------------------
      integer i,kf
      logical fork,paelm
      
      iant=0  
      do 100 i=1,n   
      if(fork(t(i),streng,na,nb)) then
        call xstrip(t(i),kf,paelm)
        if(kf.eq.(nb-na+1)) then
          iant=1
          l(1)=i
          return
        end if
        iant=iant+1
        l(iant)=i
      end if
  100 continue
      return
      end



********************************************************************     
*
*                 B A L A C
*
*   sjekker hvor mange tekster i menyen t(n) streng(na:nb) er en
*   forkortelse av. Tekstene i t(n) maa avsluttes med '#' eller '!'.
*   iant er antall treff, l(1)..l(iant) er nummerne i menyen t for
*   tekstene som er truffet. Det skilles ikke mellom sm} og store 
*   bokstaver.
*
*   kaller 'fork'.
***********************************************************************
      subroutine balac(t,n,streng,na,nb,l,iant)
      integer n,na,nb,iant,l(n)
      character*80 streng,t(n)
*-------------------------------------------------------------------
      integer i  
      logical fork
      
      iant=0  
      do 100 i=1,n   
        if(fork(t(i),streng,na,nb)) then
          iant=iant+1
          l(iant)=i
        end if
  100 continue
      return
      end




*******************************************************************          
*
*                  F O R K
*
*   tester om s(n1:n2) er en forkortelse av t(1:...). Teksten i t maa
*   avsluttes av '#' eller '!'. Testen virker uavhengig av smaa/store
*   bokstaver.
*
*   kaller 'xstrip','upchar'.
**********************************************************************
      function fork(t,s,n1,n2)
      integer n1,n2
      character*80 t,s
      logical fork
*--------------------------------------------------------------------
      integer kt,ks
      character*80 tt,ss
      logical def

      call xstrip(t,kt,def)    
      ks=n2-n1+1
      if(ks.gt.kt) then
        fork=.false.
        return
      end if

      tt(1:ks)=t(1:ks) 
      ss(1:ks)=s(n1:n2)
      call upchar(tt,1,ks)
      call upchar(ss,1,ks)
      fork= tt(1:ks).eq.ss(1:ks)

      return
      end            



****************************************************************
*                  
*                  FSOK
*
*     Leter fra gitt poisjon i en tekst til n{rmeste forekomst av
*     en gitt karakter.
*
*     parametere:
*           t(1:nmax) - tekst det s|kes i                              I
*           np - ved input: posisjon i t sIket starter fra,            I/O
*                         NB. dvs en finnner n{rmeste forekomst til
*                         utenom posisjonen selv,forekomster i endene
*                         av tekstem m} da spesialbehandles.
*                ved output: posisjon der karakter er funnet.
*                      Dersom ingen forekomst finnes returneres med
*                      np=1 eller np=nmax, avh. av ivis, og iflag=1
*                forekomst i np teller ikke
*           c - karakter som s|kes                                     I
*           ivis - =1 : forover s|k  , ellers bakover s|k              I
*           iflag - feilparameter                                      O
*                   =0 alt ok
*                    1  se ovenfor
*                    2  np er ulovlig ved input
*********************************************************************
      subroutine fsok(t,nmax,np,c,ivis,iflag)
      integer np,nmax,ivis,iflag
      character c
      character*120 t
*------------------------------------------
      integer iadd


      iflag=0

      if(np.lt.1 .or. np.gt.nmax ) then
        iflag=2
        return
      end if

      if(ivis.eq.1) then
        iadd=1
        if (np.eq.nmax) return
      else
        iadd=-1
        if (np.eq.1) return
      end if

 100  continue
      np=np+iadd
      if(np.lt.1 .or. np.gt.nmax )then
        if(np.lt.1) np=1
        if(np.gt.nmax) np=nmax
        iflag=1
        return
      end if      

      if(t(np:np).ne.c) go to 100

      return
      end

****************************************************************
*
*                 INSUBST
*
*   skifter ut foerste forekomst av en gitt karakter med et tall
*   parametere:
*            t - tekst som behandles, maa avsluttes med # el !
*            nmax - max tillat lengde av t
*            np - input: punkt der sok starter, output:
*                 forste posisjon etter insatt tal
*            c - karakter som skal byttes
*            ival - heltall som skal settes inn
*            iflag - feilparameter
*                    =0 subst. gjennomfort
*                    =1 ingen forekomst av c
*                    =2 np for stor
*                    =3 resulterende tekst for lang - intet gjort
****************************************************************
      subroutine insubst(t,nmax,np,c,ival,iflag)
      integer nmax,np,ival,iflag
      character c
      character*80 t
*----------------------------------------------------------------
      integer kf,kt,irest,i,np1
      logical hiv
      character*20 ctall
  

      call xstrip(t,kf,hiv)

      call gsok(t,kf,np,c,1,iflag)
      if(iflag.ne.0) return

      irest=kf-np
      call itconv(ival,kt,ctall)

      if( (kf+kt-1).ge.nmax) then
        iflag=3
        return
      end if

      np1=np+1

c     avsl merke maa med, derfor kf+1 som lokke-grense

      do 100 i=kf+1,np1,-1
      t(i+kt-1:i+kt-1)=t(i:i)
100   continue

      t(np:np+kt-1)=ctall(1:kt)
      np=np+kt

      return
      end
      

****************************************************************
*                  
*                  GSOK
*
*     Leter fra gitt poisjon i en tekst til n{rmeste forekomst av
*     en gitt karakter. Forskjellig fra fsok ved betyd av inp.par np
*
*     parametere:
*           t(1:nmax) - tekst det s|kes i                              I
*           np - ved input:  sIket starter fom. denne posisjon i t     I/O
*                ved output: posisjon der karakter er funnet.
*                      Dersom ingen forekomst finnes returneres med
*                      np=1 eller np=nmax, avh. av ivis, og iflag=1
*                dersom ingen forekomst finnes returneres med np uendret
*           c - karakter som s|kes                                     I
*           ivis - =1 : forover s|k  , ellers bakover s|k              I
*           iflag - feilparameter                                      O
*                   =0 alt ok
*                    1  se ovenfor
*                    2  np er ulovlig ved input
*********************************************************************
      subroutine gsok(t,nmax,np,c,ivis,iflag)
      integer np,nmax,ivis,iflag
      character c
      character*120 t
*------------------------------------------
      integer iadd,np0


      np0=np
      iflag=0

      if(np.lt.1 .or. np.gt.nmax ) then
        iflag=2
        return
      end if

      if(ivis.eq.1) then
        iadd=1
      else
        iadd=-1
      end if

      np=np-iadd

 100  continue
      np=np+iadd

      if(np.lt.1 .or. np.gt.nmax )then
        np=np0
        iflag=1
        return
      end if      

      if(t(np:np).ne.c) go to 100

      return
      end



***************************************************************************
*           K O M P R E S S
*
*     Alle blanke i t(n1:n2) som befinner seg ved siden av karakteren 
*     gitt ved c fjernes og den komprimerte teksten venstrestilles og 
*     returneres med oppdatert verdi for n2.
***************************************************************************
      subroutine kompress(t,c,n1,n2)
      integer n1,n2
      character c
      character*80 t
*-----------------------------------------------------------------------
      integer ipos,i1,iblankl,iblankr,k,i,iskip
      logical eks 


      i1=n1

 100  continue

      call lokchar(t,c,i1,n2,k,eks)
      
      if(.not.eks) return
    

      iblankl=0
 40   continue
      if(k-iblankl.eq.n1) go to 45
      if(t(k-iblankl-1:k-iblankl-1).ne.' ') go to 45
      iblankl=iblankl+1
      go to 40

 45   continue

      iblankr=0
 50   continue
      if(k+iblankr.eq.n2) go to 55
      if(t(k+iblankr+1:k+iblankr+1).ne.' ') go to 55
      iblankr=iblankr+1
      go to 50

 55   continue

      ipos= k-iblankl
      t(k:k)=' '
      t(ipos:ipos)=c

      i1=ipos+1
      iskip=iblankl+iblankr
      ipos=ipos+iskip+1
      do 70 i=ipos,n2
      t(i-iskip:i-iskip)=t(i:i)
 70   continue

      n2=n2-iskip

      go to 100

      end


******************************************************************
*     
*      Naar karakteren lik cjok forekommer, har dette samme virkning
*      som om tallet jokver sto der i stedet. Ellers virker rutinen
*      p} samme vis som geti. Det s|kes i strengen t(ia:ib)
*******************************************************************
      subroutine GETIJOK(t,cjok,jokver,ia,ib,I0,I1,ires,kFLAG)
      integer jokver,ia,ib,I0,I1,ires,kFLAG
      character cjok
      character*80 t
*----------------------------------------------------------------------

      call GETI(t,ia,ib,I0,I1,ires,kFLAG)

      if(kflag.eq.4) then
        if(t(i0:i0).eq.cjok) then
          ires=jokver
c         posisjonsparametere ma na justeres og kflag maa
c        gis verdi som om det sto jokver i stedet for cjok
          i1=i0
          if(i1.eq.ib) then
            kflag=1
          else
            if(t(i1+1:i1+1).eq.' ')then
              kflag=1
            else
              kflag=2
            end if
          end if
        end if
      end if

      return
      end


************************************************************************
*
*     Markerer tall i intervallet 1:nmax, gitt ved sekvenser i t(n1:n2). 
*     Eksempler p} sekvenser
*
*           5     -  posisjon 5 markeres
*           2:5   -  posisjoner 2,3,4 og 5 markeres. Dersom 5 er st|rre
*                    enn nmax markeres bare tom. nmax
*           2:11;2 - posisjoner 2,4,6,8 og 10 markeres, igjen bare tom. nmax
*
*     Det kan gis et vilk}rlig antall sekvenser og blanke omkring : og ; 
*     neglisjeres. Karakteren * betyr nmax og kan benyttes i alle sammen-
*     stillinger.
*     Parametere:
*            t(n1:n2) - tekst med sekvenser                              I
*            nmax    - ovre grense for tall                              I
*            med   -  Dersom med(j) er sann er j markert                 O
*            iflag -  feilparameter                                      O
*                   =0   alt i orden
*                    1   funn av sekvens som begynner ulovlig
*                    2   karakteren ':' er ikke funnet som forventet
*                    3   karakteren ':' er ikke fulgt av heltall ( eller *)
*                    6   karakteren ';' er ikke funnet som forventet
*                    7   karakteren ';' er ikke fulgt av heltall ( eller *)
******************************************************************************

      subroutine konvert(t,n1,n2,nmax,med,iflag)
      integer n1,n2,nmax,iflag
      character*80 t
      logical med(nmax)
*------------------------------------------------------------------------
      integer l(2,20),iant,ia,ib,khigh,kstep,knum,k
      integer j,i0,i1,kflag

      iflag=0

      call kompress(t,':',n1,n2)
      call kompress(t,';',n1,n2)

      call deltek(t,n1,n2,l,iant)

      if(iant.eq.0) then
        do 50 k=1,nmax
        med(k)=.true.
 50     continue
        return
      else
        do 60 k=1,nmax
        med(k)=.false.
 60     continue
      end if

      do 100 j=1,iant
      ia=l(1,j)
      ib=l(2,j)
      call GETIJOK(t,'*',nmax,ia,ib,I0,I1,knum,kFLAG)

      if(kflag.gt.2) then
        iflag=1
        return
      end if

      if(kflag.eq.1)  then
        if(knum.le.nmax) med(knum)=.true.
      else
        if(t(i1+1:i1+1).ne.':') then
          iflag=2
          return
        end if
        ia=i1+2
        call GETIJOK(t,'*',nmax,ia,ib,I0,I1,khigh,kFLAG)
        if(kflag.eq.0 .or.kflag.eq.3) then
           iflag=3
           return
        end if


        if(kflag.eq.1) then
          kstep=1
        else
          if(t(i1+1:i1+1).ne.';') then
            iflag=6
            return
          end if
  
          ia=i1+2
          call GETI(t,ia,ib,I0,I1,kstep,kFLAG)
          if(kflag.ne.1) then
            iflag=7
            return
          end if
 
        end if

        if(khigh.gt.nmax) khigh=nmax

        do 70 k=knum,khigh,kstep
        med(k)=.true.
 70     continue

      end if
 100  continue

      return
      end
          


************************************************************************
*
*     Markerer tall i intervallet n0:nmax, gitt ved sekvenser i t(n1:n2). 
*     Eksempler p} sekvenser
*
*           5     -  posisjon 5 markeres
*           2:5   -  posisjoner 2,3,4 og 5 markeres. Dersom 5 er st|rre
*                    enn nmax markeres bare tom. nmax
*           2:11;2 - posisjoner 2,4,6,8 og 10 markeres, igjen bare tom. nmax
*
*     Det kan gis et vilk}rlig antall sekvenser og blanke omkring : og ; 
*     neglisjeres. Karakteren * betyr nmax og kan benyttes i alle sammen-
*     stillinger.
*     Parametere:
*            t(n1:n2) - tekst med sekvenser                              I
*            nmax    - ovre grense for tall                              I
*            med   -  Dersom med(j) er sann er j markert                 O
*            iflag -  feilparameter                                      O
*                   =0   alt i orden
*                    1   funn av sekvens som begynner ulovlig
*                    2   karakteren ':' er ikke funnet som forventet
*                    3   karakteren ':' er ikke fulgt av heltall ( eller *)
*                    6   karakteren ';' er ikke funnet som forventet
*                    7   karakteren ';' er ikke fulgt av heltall ( eller *)
******************************************************************************

      subroutine gkonvert(t,n1,n2,n0,nmax,med,iflag)
      integer n1,n2,n0,nmax,iflag
      character*80 t
      logical med(n0:nmax)
*------------------------------------------------------------------------
      integer l(2,20),iant,ia,ib,khigh,kstep,knum,k
      integer j,i0,i1,kflag

      iflag=0

      call kompress(t,':',n1,n2)
      call kompress(t,';',n1,n2)

      call deltek(t,n1,n2,l,iant)

      if(iant.eq.0) then
        do 50 k=n0,nmax
        med(k)=.true.
 50     continue
        return
      else
        do 60 k=n0,nmax
        med(k)=.false.
 60     continue
      end if

      do 100 j=1,iant
      ia=l(1,j)
      ib=l(2,j)
      call GETIJOK(t,'*',nmax,ia,ib,I0,I1,knum,kFLAG)

      if(kflag.gt.2) then
        iflag=1
        return
      end if

      if(kflag.eq.1)  then
        if(knum.le.nmax .and. knum.ge.n0) med(knum)=.true.
      else
        if(t(i1+1:i1+1).ne.':') then
          iflag=2
          return
        end if
        ia=i1+2
        call GETIJOK(t,'*',nmax,ia,ib,I0,I1,khigh,kFLAG)
        if(kflag.eq.0 .or.kflag.eq.3) then
           iflag=3
           return
        end if


        if(kflag.eq.1) then
          kstep=1
        else
          if(t(i1+1:i1+1).ne.';') then
            iflag=6
            return
          end if
  
          ia=i1+2
          call GETI(t,ia,ib,I0,I1,kstep,kFLAG)
          if(kflag.ne.1) then
            iflag=7
            return
          end if
 
        end if

        if(khigh.gt.nmax) khigh=nmax

        do 70 k=knum,khigh,kstep
        if(k.ge.n0)med(k)=.true.
 70     continue

      end if
 100  continue

      return
      end
          


***************************************************************************
*
*                SKIPKOM
*
*    Rutinen skipper alle linjer som ikke starter med tall.
*    Den er ment for innlesning av filer i gkurv-format, og bruker
*    "getr" for aa gjenkjenne tall. 
*
*    parametere:
*             n - antall forbigaatte records                        O
*             itape - unitnummer                                    I
*             ierr -  feilparameter                                 O
*                  =0   alt klart for lesning av tallkolonner
*                  =1   EOF funnet under skipping                   
*                  =2   feil i innlesning under skipping
*                  =3   linje funnet som har feil format, men som
*                       ville blitt godkjent av gkurv
**************************************************************************
      subroutine skipkom(n,itape,ierr)
      integer n,itape,ierr
*-------------------------------------------------------------------------
      integer iflag,i0,i1
      real r
      character*80 t

      n=0

 100  continue
      read(itape,'(a80)',err=400,end=300) t(1:80)

      call getrl(t,1,80,i0,i1,r,iflag)

      if(iflag.eq.1) then
        ierr=0
        backspace itape
        return   
      end if

      if(iflag.eq.2 .or. iflag.eq.6) then
        ierr=3
        backspace itape
        return
      end if

      n=n+1
      go to 100


 300  continue

      ierr=1
      return

 400  continue
      ierr=2

      return

      end

***************************************************************************
*
*                SLOPSKIP
*
*    Rutinen skipper alle linjer som ikke starter med tall.
*    Den er ment for innlesning av filer i gkurv-format, og bruker
*    lokint for aa gjenkjenne tall. Mer sloppy en skipkom. 
*
*    parametere:
*             n - antall forbigaatte records                        O
*             itape - unitnummer                                    I
*             ierr -  feilparameter                                 O
*                  =0   alt klart for lesning av tallkolonner
*                  =1   EOF funnet under skipping                   
*                  =2   feil i innlesning under skipping
**************************************************************************
      subroutine slopskip(n,itape,ierr)
      integer n,itape,ierr
*-------------------------------------------------------------------------
      integer iflag,i0,i1
      real r
      character*80 t

      n=0

 100  continue
      read(itape,'(a80)',err=400,end=300) t(1:80)

      call lokint(t,1,80,i0,i1,iflag)

      if(iflag.eq.1 .or. iflag.eq.2) then
        ierr=0
        backspace itape
        return   
      end if


      n=n+1
      go to 100


 300  continue

      ierr=1
      return

 400  continue
      ierr=2

      return

      end

*********************************************************************
*
*             WRIBCHAR
*
*    Rutinen representerer tegnene i teksten t ved sine ascinummer som
*    saa skrives ut uformattert som 32 bits heltall. Hensikten er kunne
*    bringe tekster helskinnet gjennom byte-swap operasjoner.
*    parametere: 
*            t(1:n) - tekst som skal skrives ut                I
*            itape - unitnummer for utskrift                   I
*            iflag - feilparameter                             O
*                  = 0 alt ok.
*                  = 1 for stor n
*                  = 2 feil unders skrivning
**********************************************************************
      subroutine wribchar(t,n,itape,iflag)
      integer n,itape,iflag
      character*200 t
*----------------------------------------------------------------------
      integer i,nr(200),ichar

      if(n.gt.200) then
        iflag=1
        return
      end if


      do 100 i=1,n
      nr(i)=ichar(t(i:i))
 100  continue

      write(itape,err=200)(nr(i) , i=1,n)

      iflag=0
      return

 200  iflag=2
      return

      end

*********************************************************************
*
*             REABCHAR
*
*    Rutinen leser inn uformatterte 32 bits heltall som saa konverteres til
*    tegn i hht. til ascii-nummeret. Rutinen er beregnet for aa lese det som
*    wribchar skriver.
*    parametere: 
*            t    - tekst for plassering av kar.                  O
*            n - antall tegn                                      I
*            itape - unitnummer for lesning                       I
*            iflag - feilparameter                                O
*                  = 0 alt ok.
*                  = 1 for stor n
*                  = 2 feil under lesning
*                  = 3 slutt paa fil itape
*                  = 4 funnet ascii.nr. utenfor 0...127. Det tilhorende 
*                      tegn settes da til char(0)
**********************************************************************
      subroutine reabchar(t,n,itape,iflag)
      integer n,itape,iflag
      character*200 t
*----------------------------------------------------------------------
      integer i,nr(200),nc
      character char

      if(n.gt.200) then
        iflag=1
        return
      end if

      iflag=0

      read(itape,err=200,end=300)(nr(i) , i=1,n)


      do 100 i=1,n
      nc=nr(i)
      if(nc.lt.0 .or. nc.gt.127) then
        iflag=4
        nc=0
      end if
      t(i:i)=char(nc)
 100  continue


      return

 200  iflag=2
      return

 300  iflag=3
      return

      end












*******************************************************************************
*
*                   I T A L L 
*
*    Leser svar paa spoersmaal i teksten spor. Avsluttes dette med '#' tilbys
*    standardverdi idef som oppnaes ved returnering av tom linje. Avslutting
*    med '!' foerer til ignorering av standardverdi idef. Inngaar verken
*    '#' eller '!' i spor antas det aa ha lengde 80 og standardverdi tilbys.
*
*    kaller 'primles','xstrip' og 'geti'.
*******************************************************************************
      function itall(idef,spor)
      integer itall,idef
      character*80 spor
      include 'styr.inc'
*---------------------------------------------------------
      character*80 svar
      integer k,i1,i2,iflag,ib,ires,iverdi
      logical def        

      call xstrip(spor,k,def)
      if(k.eq.0) k=1                           

  200 continue
    
      if(def) then
        write(iustrr,*) spor(1:k) ,'/', idef ,'/'
        if(aktlog.eq.1 .and. aktkom.eq.1)
     %  write(ltape,*) charkom,spor(1:k) ,'/', idef ,'/'
      else
        write(iustrr,*) spor(1:k)
        if(aktlog.eq.1 .and. aktkom.eq.1)
     %     write(ltape,*) charkom,spor(1:k) 
      end if

      call primles(svar,iflag)

      if(iflag.eq.2) stop
      if(iflag.ne.0) go to 200

      call geti(svar,1,80,i1,i2,ires,iflag)
      if(iflag.eq.0) then
        if(def) then
          iverdi=idef
          go to 500
        else
          call primvri('No default value#')
          go to 200
        end if
      end if
      
      if(iflag.gt.1) then
        call primvri('illegal input#')
        go to 200
      else
        iverdi=ires
        ib=i2+1
        call geti(svar,ib,80,i1,i2,ires,iflag)
        if(iflag.ne.0) then
          call primvri('Please give a single integer#')
          go to 200
        end if
      end if

 500  continue

      itall=iverdi
      dispens=0

      return
      end



**************************************************************************
*
*              SETPA
*    Fyller homedir, defdir og hosttype inn i commonomraadet direk og setter
*    iustrr i com. komstyr
********************************************************************
      subroutine setpa
      include 'styr.inc'
*-------------------------------------------------------
      integer kf
      open(unit=20,file='scrp.xx',status='unknown')
      write(20,'(a)')'#!/bin/sh'
      write(20,'(a)')'echo $HOME > qrw.zzz'
      write(20,'(a)')'pwd >> qrw.zzz'
      write(20,'(a)')'hosttype >> qrw.zzz'
      close(20)
      call system('chmod +x scrp.xx ; ./scrp.xx ')
      open(unit=20,file='qrw.zzz',status='unknown')
      call blank(homedir,120)
      read(20,'(a120)') homedir(1:120)
      call strip(homedir,120,nhd)
      nhd=nhd+1
      homedir(nhd:nhd)='/'
      call blank(defdir,120)
      read(20,'(a120)') defdir(1:120)
      call strip(defdir,120,ndd)
      ndd=ndd+1
      defdir(ndd:ndd)='/'
      read(20,'(a20)') hostt(1:20)
      call strip(hostt,20,kf)
      if(hostt(1:kf).eq.'hp') then
        iustrr=7
      else
        iustrr=0
      end if

      call system('rm -f scrp.xx qrw.zzz ')

      return
      end



**************************************************************************
*
*              SETPANY
*    Fyller homedir, defdir, lognavn, prosessnummer, host og hosttype 
*    inn i commonomraadet direk og setter
*    iustrr i com. komstyr. Benytter bla .a. getenv
********************************************************************
      subroutine setpany
      include 'styr.inc'
*-------------------------------------------------------
      integer kf,ierr,getcwd,getpid

      call blank(hostt,40)
      call getenv('HOSTTYPE',hostt)
      call strip(hostt,40,kf)
      if(hostt(1:2).eq.'hp') then
        iustrr=7
      else
        iustrr=0
      end if


      call blank(homedir,120)
      call getenv('HOME',homedir)
      call strip(homedir,120,nhd)
      if(nhd.gt.0) then
        nhd=nhd+1
        homedir(nhd:nhd)='/'
      else
        write(iustrr,'(a)')'$HOME ikke satt'
      end if
      

      call blank(vert,40)
      call getenv('HOST',vert)
      
      call blank(lognav,40)
cc      call getlog(lognav)

      pronum=getpid()

      call blank(defdir,120)
      ierr=getcwd(defdir)
      if(ierr.eq.0.or.hostt(1:2).ne.'ds') then
        call strip(defdir,120,ndd)
        ndd=ndd+1
        defdir(ndd:ndd)='/'
      else
        write(iustrr,'(a)')'finner ikke def-dir'
      end if
        

      return
      end


**********************************************************************
*             FULLPATH
*
*       Gjor om lokalt navn til fullt pathnavn.
*           t - tekst med filnavn                      I/O
*           n - maks lengde av t                       I
*               dersom nodvendig trunkeres fremste part av pathnavn
**********************************************************************
      subroutine fullpath(t,n)
      integer n
      character*200 t
      include 'styr.inc'
*-------------------------------------------------------------- 
      integer kf,kadd

      call frontstrip(t,n,kf)

      if(t(1:1).eq.'/') return

      call strip(t,n,kf)

      kadd=min(ndd,n-kf)

      t(kadd+1:kadd+kf)=t(1:kf)

      t(1:kadd)=defdir(ndd-kadd+1:ndd)

      return
      end

      

***********************************************************************
*
*              FETCHENV
*
*   Henter ut opplysninger i en character*80 array.
*                cinf(1) - inneholder lognavn og pid (fra pos hhv. 1 og 41)
*                     2  -            host
*                     3  -            working dir. ( trunkeres evt. forfra)
*********************************************************************
      subroutine fetchenv(cinf)
      character*80 cinf(3)
      include 'styr.inc'
*------------------------------------------------------------------------
      integer kf,i
      character*20 cpid

      do 100 i=1,3
        call blank(cinf(i),80)
 100  continue  

      call itconv(pronum,kf,cpid)
      cinf(1)(1:40)=lognav(1:40)
      cinf(1)(41:40+kf)=cpid(1:kf)
 
      cinf(2)(1:40)=vert(1:40)
      
      if(ndd.gt.80) then
        cinf(3)(1:80)=defdir(ndd-79:ndd)
      else
        cinf(3)(1:ndd)=defdir(1:ndd)
      end if

      return
      end


**************************************************************************
*
*              SKRIVENV
*    Skriver ut homedir, defdir og hosttype fra commonomraadet direk.
********************************************************************
      subroutine skrivenv(itape)
      integer itape
      include 'styr.inc'
*-------------------------------------------------------
      integer kf,ierr

      write(itape,'(a,a)')'HOSTTYPE=',hostt(1:20)
      write(itape,'(a,a)')'HOME=',homedir(1:nhd)
      write(itape,'(a,a)')'wdir:',defdir(1:ndd)

      return
      end



**********************************************************************
      subroutine initkom
      include 'styr.inc'
*----------------------------------------------------------------------

      aktlog=0
      sattlog=0

      aktkom=0
      sattkom=0
      dispens=0

      intape=5
c     lese-enhet  settes til standard  input

      call setpany
c     def og home dir leses inn.
      return
      end

************************************************************************
*     setter lese-uniten til itape
********************************************************************
      subroutine setinp(itape)
      integer itape
      include 'styr.inc'
*------------------------------------------------------------------
      
      if(itape.le.0 .or. itape .eq.6)then
        write(iustrr,*)'ulovlig enhet for innlesning:',itape
        return
      end if

      intape=itape

      return
      end

************************************************************************
*     returnerer lese-enheten i infast
********************************************************************
      subroutine inpenh(infast)
      integer infast
      include 'styr.inc'
*------------------------------------------------------------------
      infast=intape

      return
      end


************************************************************************
*     returnerer unitnummer for standard error i iero
********************************************************************
      subroutine stdrror(iero)
      integer iero
      include 'styr.inc'
*------------------------------------------------------------------
      iero=iustrr

      return
      end


************************************************************************
*
*     Dersom fnavn inneholder navnet p} en mulig input fil, }pnes denne
*     med unit intape og innum settes som lesenhet. dersom fila ikke
*     finnes returneres uten noen endring.
**********************************************************************
      subroutine inpfil(fnavn,innum)
      integer innum
      character*80 fnavn
*-----------------------------------------------------------------------
      integer kf
      logical hiv
      integer iustrr

      call stdrror(iustrr)


      call xstrip(fnavn,kf,hiv)


      open(unit=innum,file=fnavn(1:kf),status='old',err=200)
      write(iustrr,*)'fil: ',fnavn(1:kf),' }pnet for input'
      call setinp(innum)
      return

 200  continue

      write(iustrr,*)'fil: ',fnavn(1:kf),' kan ikke }pnes'
      return
      end
      
************************************************************************
      subroutine setkomm(c)
      character c
      include 'styr.inc'
*----------------------------------------------------------------------

      sattkom=1
      aktkom=1
      charkom=c

      return
      end

************************************************************************
      subroutine stoppkomm
      include 'styr.inc'
*----------------------------------------------------------------------


      aktkom=0

      return
      end


************************************************************************
      subroutine startkomm
      include 'styr.inc'
*----------------------------------------------------------------------

      if(sattkom.eq.1) then
        aktkom=1
      else
        write(iustrr,*)'comments attempted chosen without enabling'
      end if

      return
      end


************************************************************************
      subroutine setlog(itape)
      integer itape
      include 'styr.inc'
*----------------------------------------------------------------------

      sattlog=1
      aktlog=1
      ltape=itape

      return
      end

************************************************************************
      subroutine stopplog
      include 'styr.inc'
*----------------------------------------------------------------------


      aktlog=0

      return
      end


************************************************************************
      subroutine startlog
      include 'styr.inc'
*----------------------------------------------------------------------


      if(sattlog.eq.1) then
        aktlog=1
      else
        write(iustrr,*)'logging set without enabling'
      end if

      return
      end


**************************************************************************
      subroutine primvri(t)
      character*80 t
      include 'styr.inc'

*----------------------------------------------------------------------
      integer k
      logical def        

      call xstrip(t,k,def)
      if(k.ge.78) k=78
      if(k.gt.0) then                           
        write(iustrr,*) t(1:k)
        if(aktlog.eq.1 .and. aktkom.eq.1 )
     %  write(ltape,*) charkom,t(1:k)
      else
        write(iustrr,*)' '
        if(aktlog.eq.1 .and. aktkom.eq.1)
     %  write(ltape,*) charkom
      end if

      return
      end


***************************************************************************
      subroutine primles(t,iflag)
      integer iflag
      character*80 t
      include 'styr.inc'

*----------------------------------------------------------------------
      integer i0,i1,kflag

      iflag=0

 100  continue

      if(dispens.eq.1) then
        read(5,20,end=50) t
      else
        read(intape,20,end=50) t
      end if
 20   format(a80)
 25   format(a78)

      call lokord(t,1,80,i0,i1,kflag)

      if(intape.ne.5 .and.
     %     t(i0:i0+3).eq.'????' .and.(dispens.eq.0)) then
        dispens=1
        write(iustrr,*)'input fra standard'
        read(5,20,end=50) t
      end if

      if(kflag.ne.0 .and. aktkom.eq.1) then
        if(t(i0:i0).eq.charkom) then
          if(aktlog.eq.1) write(ltape,25) t(1:78)
          go to 100
        end if
      end if

      if(aktlog.eq.1) write(ltape,25) t(1:78) 
      if(intape.ne.5 .and.(dispens.eq.0)) write(iustrr,25) t(1:78) 

      return

 50   continue
      if(dispens.eq.1 .or. intape.eq.5) then
         write(iustrr,*)'end of file funnet p} standard input'
         iflag=2
         return
      end if
      write(iustrr,*)'slutt p} input fra enhet=',intape
      intape=5
      iflag=1
      return
      
      end




*******************************************************************************
*
*                    R T A L L
*
*     Real-utgaven av itall.
*
*     kaller 'xstrip' og 'getr'

*******************************************************************************
*
*                    R T A L L
*
*     Real-utgaven av itall.
*
*     kaller 'xstrip' og 'getr'
******************************************************************************
      function rtallo(rdef,spor)
      real rtallo,rdef
      character*80 spor
      include 'styr.inc'
*---------------------------------------------------------
      character*80 svar
      integer k,i,i1,i2,iflag,ib,kflag
      real rres,rverdi
      logical def        

      call xstrip(spor,k,def)                           
      if(k.eq.0) k=1

  200 continue
   
      if(def) then
        write(iustrr,*) spor(1:k) ,'/', rdef ,'/'
        if(aktlog.eq.1 .and. aktkom.eq.1)
     %  write(ltape,*) charkom,spor(1:k) ,'/', rdef ,'/'
      else
        write(iustrr,*) spor(1:k)
        if(aktlog.eq.1 .and. aktkom.eq.1)
     %     write(ltape,*) charkom,spor(1:k) 
      end if

      call primles(svar,kflag)

      if(kflag.eq.2) stop
      if(kflag.ne.0) go to 200

      call getr(svar,1,80,i1,i2,rres,iflag)
      if(iflag.eq.0) then
        if(def) then
          rverdi=rdef
          go to 500
        else
          call primvri('No default value available#')
          go to 200
        end if
      end if
      
      if(iflag.gt.1) then
        call primvri('illegal input#')
        go to 200
      else
        rverdi=rres
        ib=i2+1
        call getr(svar,ib,80,i1,i2,rres,iflag)
        if(iflag.ne.0) then
          call primvri('Please, give a single real#')
          go to 200
        end if
      end if

 500  continue

      rtallo=rverdi
      dispens=0

      return
      end                 



******************************************************************************
      function rtall(rdef,spor)
      real rtall,rdef
      character*80 spor
      include 'styr.inc'
*---------------------------------------------------------
      character*80 svar
      integer k,i,i1,i2,iflag,ib,kflag
      real rres,rverdi
      logical def        

      call xstrip(spor,k,def)                           
      if(k.eq.0) k=1

  200 continue
   
      if(def) then
        write(iustrr,*) spor(1:k) ,'/', rdef ,'/'
        if(aktlog.eq.1 .and. aktkom.eq.1)
     %  write(ltape,*) charkom,spor(1:k) ,'/', rdef ,'/'
      else
        write(iustrr,*) spor(1:k)
        if(aktlog.eq.1 .and. aktkom.eq.1)
     %     write(ltape,*) charkom,spor(1:k) 
      end if

      call primles(svar,kflag)

      if(kflag.eq.2) stop
      if(kflag.ne.0) go to 200

      call getrl(svar,1,80,i1,i2,rres,iflag)
      if(iflag.eq.0) then
        if(def) then
          rverdi=rdef
          go to 500
        else
          call primvri('No default value available#')
          go to 200
        end if
      end if
      
      if(iflag.gt.1) then
        call primvri('illegal input#')
        go to 200
      else
        rverdi=rres
        ib=i2+1
        call getrl(svar,ib,80,i1,i2,rres,iflag)
        if(iflag.ne.0) then
          call primvri('Please, give a single real#')
          go to 200
        end if
      end if

 500  continue

      rtall=rverdi
      dispens=0

      return
      end                 


*******************************************************************************
*           
*                     L E S T E K
*
*    stiller sporsmal spor og leser svaret som teksten 'tekst'. spor angis
*    som for 'itall', standardverdien tdef avsluttes med # eller '!'(disse
*    blir ikke med i svar). Svar antas aa vaere 80 lang og etterfylles med
*    blanke. En tom tekst kan ikke vare standardverdi.                                                                 
*
*    kaller 'xstrip' og 'lokord'
*******************************************************************************
      subroutine lestek(tdef,spor,svar)
      character*80 tdef,svar, spor
      include 'styr.inc'
*---------------------------------------------------------
      integer k,i,i1,i2,k1,iflag,kflag
      logical def,daf        

      do 70 i=1,80
   70 svar(i:i)=' '

      call xstrip(spor,k,def)
      if(k.eq.0) k=1                           

      if(def) then
        call xstrip(tdef,k1,daf)
        if(k1.eq.0) def=.false.
      end if

  200 continue

      if(def) then
        write(iustrr,*) spor(1:k) ,'/', tdef(1:k1) ,'/'
        if(aktlog.eq.1 .and. aktkom.eq.1)
     %  write(ltape,*) charkom,spor(1:k) ,'/', tdef(1:k1) ,'/'
      else
        call primvri(spor)
      end if

      call primles(svar,kflag)
      if(kflag.eq.2) stop
      if(kflag.ne.0) go to 200


      call lokord(svar,1,80,i1,i2,iflag)
      if(iflag.eq.0) then
        if(def) then
          svar(1:k1)=tdef(1:k1)
          go to 500
        else
          call primvri('No default value available#')
          go to 200
        end if
      end if
      
 500  continue

      dispens=0
      
      return
      end



***********************************************************************
*
*                      L S J E K K
*
*    stiller spoersmal spor ,standardverdi tdef faaes eller ignoreres
*    paa samme vis som i 'itall'. I svaret (teksten svar) lokaliseres
*    sekvensene til l(1,k):l(2,k) k=1..iant. Som lovlig svar godtas
*    bare de hvis foerste sekvens er en entydig forkortelse til et ord
*    i menyen liste(nl); lnum angir hvilket ord. 
*
*    kaller 'alacart','deltek','red' ,'lestek' og 'upchar'.
*****************************************************************************
      subroutine lsjekk(tdef,spor,svar,l,iant,liste,nl,lnum)
      integer iant,lnum,l(2,20),nl
      character*80 tdef,spor,svar,liste(nl)
      include 'styr.inc'
*----------------------------------------------------------------------
      integer ltall,men(20),i



  100 continue
      call lestek(tdef,spor,svar)
      call red(svar,80)
      call deltek(svar,1,80,l,iant)
      i=l(2,iant)
      if(svar(i:i).eq.'#') l(2,iant)=i-1
      call alacart(liste,nl,svar,l(1,1),l(2,1),men,ltall)

      if(ltall.eq.0) then
        write(iustrr,*) 'illegal command: /',svar(l(1,1):l(2,1)),'/'
        if(aktlog.eq.1 .and. aktkom.eq.1) 
     %   write(ltape,*) charkom,'illegal command: /'
        if(aktlog.eq.1 .and. aktkom.eq.1) 
     %   write(ltape,*) charkom,svar(l(1,1):l(2,1)),'/' 
        go to 100
      end if
      if(ltall.gt.1) then
        write(iustrr,*) 'ambiguous abbreviation: /',svar(1:l(2,1)),'/'
        if(aktlog.eq.1 .and. aktkom.eq.1) 
     %   write(ltape,*) charkom,'ambiguous abbreviation: /'
        if(aktlog.eq.1 .and. aktkom.eq.1) 
     %  write(ltape,*) charkom,svar(1:l(2,1)),'/'
        go to 100
      end if

      lnum=men(1)

      return
      end


************************************************************************* 
*
*                     J A
*    leser ja/nei svar(returnerer med true dersom svar er 'ja').
*    Spoersmal spor gis som for 'itall'. Svarene 'ja','nei','yes','no'
*    og deres forkortelser godtas.
*
*    kaller: 'lestek','lokord','fork'
**************************************************************************
      function ja(def,spor)
      logical def,ja
      character*80 spor
      include 'styr.inc'
*-------------------------------------------------------------------
      character*80 cdef,svar
      integer i1,i2,i3,i4,iflag,n0      
      logical fork,lik

      if(def) then
        cdef='yes#'
      else
        cdef='no#'   
      end if
  100 continue
            
      call lestek(cdef,spor,svar)      
      call lokord(svar,1,80,i1,i2,iflag)
      n0=i2+1
      call lokord(svar,n0,80,i3,i4,iflag)
      if(iflag.ne.0) then
        call primvri(' Answer yes or no#')
        go to 100
      end if

      if( fork('ja#',svar,i1,i2) ) then
        ja=.TRUE.
        return
      end if

      if( fork('nei#',svar,i1,i2) ) then
        ja=.false.
        return
      end if
      if( fork('yes#',svar,i1,i2) ) then   
        ja=.true.
        return
      end if
      if( fork('oui#',svar,i1,i2) ) then   
        ja=.true.
        return
      end if
      if( fork('no#',svar,i1,i2) ) then 
        ja=.false.
        return
      end if
     
      call primvri('Answer yes or no#')
      go to 100
      end




*************************************************************************     *  
*
*             L E S O R D 
*
*     spoer om og leser et enkelt ord. Stopptegnene '#' og '!' brukes
*     paa vanlig maate. svaret ligger i svar(1:il).
*
*    kaller: 'lestek','lokord'.
**************************************************************************
      subroutine lesord(wdef,spor,svar,il)
      character*80 wdef,spor,svar 
      integer il
      include 'styr.inc'
*--------------------------------------------------------------------
      integer i,kb,k1,k2,k3,k4,iflag
      character*80 resp

  100 continue
      call lestek(wdef,spor,resp) 
      call lokord(resp,1,80,k1,k2,iflag) 
      kb=k2+1
      call lokord(resp,kb,80,k3,k4,iflag)
      if(iflag.ne.0) then
        call primvri('Give a single word, please#')
        go to 100
      end if
      
      il=k2-k1+1
      svar(1:il)=resp(k1:k2)   
      return
      end



**********************************************************************
*          LESFRAREG
*
*    Leser ett valg fra en liste av ord, forkortelser aksepteres
*    paa vanlig vis.
*        spor - sporsmaal                                    I
*        idef - nummer paa default valg                      I
*        reg - liste av alternativer                         I
*        n  - antall alternativer                            I
*        ires - nummer paa svar                              O
**********************************************************************
      subroutine lesfrareg(spor,idef,reg,n,ires)
      integer n,ires,idef
      character*80 reg(n),spor
*-----------------------------------------------------------------
      character*80 wdef,svar,mess
      integer il,kf,l(50),i,iant,ig
      logical kast,fork

      if(idef.gt.n .or. n.gt.50) then
        call primvri('Error in lesfrareg#')
        return
      end if

      call xstrip(reg(idef),kf,kast)
      wdef(1:kf)=reg(idef)(1:kf)
      wdef(kf+1:kf+1)='#'

 100  continue
      call lesord(wdef,spor,svar,il)
      


      call velgcart(reg,n,svar,1,il,l,iant)

      if(iant.gt.1) then
        call primvri('ambiguous abbreviation#')
        go to 100
      end if

      if(iant.eq.0) then
        if(fork('hjelp#',svar,1,il) .or. fork('help#',svar,1,il)) then
          call primvri('Available alternatives:#')
          do 50 i=1,n
          call xstrip(reg(i),kf,kast)
          mess(1:kf+1)=reg(i)(1:kf+1)
          call primvri(mess)
 50       continue
        else
          mess(1:14)='illegal answer:'
          ig=min(il,64)
          mess(15:15+ig)=svar(1:ig)
          mess(16+ig:16+ig)='#'
          call primvri(mess)
        end if
        go to 100
      end if

      ires=l(1)

      end





**********************************************************************
*          MERKIREG
*
*    Leser ett utvalg fra en liste av ord, forkortelser aksepteres
*    paa vanlig vis dersom de er entydige.
*        spor - sporsmaal                                    I
*        def - default svar, maa væere korekt                      I
*        reg - liste av alternativer                         I
*        n  - antall alternativer                            I
*        mark - gir antall angivelser av hvert ord           O
**********************************************************************
      subroutine merkireg(spor,def,reg,n,mark)
      integer n,mark(n)
      character*80 reg(n),spor,def
*-----------------------------------------------------------------
      character*80 wdef,svar,mess
      integer il,kf,l(100),i,iant,ig,lp(2,40),k,iord,ires
      logical kast,fork

      if(n.gt.100) then
        call primvri('Error in lesfrareg#')
        return
      end if


 100  continue

      do 30 i=1,n
        mark(i)=0
 30   continue
  
      call lestek(def,spor,svar)
      
      call deltek(svar,1,80,lp,iord)

      do 70 k=1,iord
      call velgcart(reg,n,svar,lp(1,k),lp(2,k),l,iant)

      if(iant.gt.1) then
        call primvri('ambiguous abbreviation#')
        go to 100
      end if

      if(iant.eq.0) then
        if(fork('hjelp#',svar,lp(1,k),lp(2,k)) .or.
     %       fork('help#',svar,lp(1,k),lp(2,k))) then
          call primvri('Available alternatives:#')
          do 50 i=1,n
          call xstrip(reg(i),kf,kast)
          mess(1:kf+1)=reg(i)(1:kf+1)
          call primvri(mess)
 50       continue
        else
          mess(1:16)='illegal answer: '
          ig=17+lp(2,k)-lp(1,k)
          mess(17:ig)=svar(lp(1,k):lp(2,k))
          mess(ig+1:ig+1)='#'
          call primvri(mess)
        end if
        go to 100
      end if

      ires=l(1)
      mark(ires)=mark(ires)+1

 70   continue


      end









************************************************************************
*
*              LESNRS
***********************************************************************
      subroutine lesnrs(def,spor,svar,med,irad)
      integer irad
      logical med(irad)
      character*80 def,spor,svar
*----------------------------------------------------------------------
      integer iflag,i2,kf
      character*80 mess

      mess(1:14)='antall rader: '
      call itconv(irad,kf,svar)
      mess(15:15+kf)=svar(1:kf)
      mess(16+kf:16+kf)='#'

 100  call primvri(mess)
      call blank(svar,80)
      call lestek(def,spor,svar)
      call strip(svar,80,i2)
      if(i2.lt.80) svar(i2+1:i2+1)='#'
      call konvert(svar,1,i2,irad,med,iflag)
      if(iflag.gt.0) then
        call primvri('ulovlig nummerangivelse#')
        go to 100
      end if

      return
      end


******************************************************************************
*
*                L E S I 2
*     virker som 'itall' men leser to heltall ir1,ir2 med standardverdier
*     id1,id2.                                                                
*
*     kaller 'lesveki'.
******************************************************************************
      subroutine lesi2(id1,id2,spor,ir1,ir2)
      integer id1,id2,ir1,ir2
      character*80 spor
*---------------------------------------------------------------------------
      integer id(2),isv(2)
      
      id(1)=id1
      id(2)=id2
      call lesveki(id,spor,isv,2)
      ir1=isv(1)
      ir2=isv(2)
      return
      end

******************************************************************************
      subroutine lesi3(id1,id2,id3,spor,ir1,ir2,ir3)
      integer id1,id2,id3,ir1,ir2,ir3
      character*80 spor
*---------------------------------------------------------------------------
      integer id(3),isv(3)
      
      id(1)=id1
      id(2)=id2
      id(3)=id3
      call lesveki(id,spor,isv,3)
      ir1=isv(1)
      ir2=isv(2)
      ir3=isv(3)
      return
      end


******************************************************************************
      subroutine lesi4(id1,id2,id3,id4,spor,ir1,ir2,ir3,ir4)
      integer id1,id2,id3,id4,ir1,ir2,ir3,ir4
      character*80 spor
*---------------------------------------------------------------------------
      integer id(4),isv(4)
      
      id(1)=id1
      id(2)=id2
      id(3)=id3   
      id(4)=id4
      call lesveki(id,spor,isv,4)
      ir1=isv(1)
      ir2=isv(2)
      ir3=isv(3)                
      ir4=isv(4)
      return
      end



******************************************************************************
*
*                  L E S V E K I
*
*     leser en array med heltall. 
*
*     kaller: 'xstrip','deltek','geti'
*****************************************************************************
      subroutine lesveki(id,spor,isv,n)
      integer n,id(n),isv(n)
      character*80 spor
      include 'styr.inc'
*----------------------------------------------------------------------------
      integer l(2,20),ires,io ,ks ,i ,iflag ,k1,k2,kflag
      logical def
      character*80 linje

      call xstrip(spor,ks,def)
      if(ks.eq.0) ks=1
   
  100 continue

      if(def) then
        write(iustrr,*) spor(1:ks),'/',(id(i) , i=1,n),'/'
        if(aktlog.eq.1 .and. aktkom.eq.1) then
          write(ltape,*) charkom,spor(1:ks),'/'
          do 105 i=1,n
          write(ltape,*)charkom,id(i) 
 105      continue
          write(ltape,*) charkom,'/'
        end if

      else
        call primvri(spor)
      end if              

      call primles(linje,kflag)
      if(kflag.eq.2) stop
      if(kflag.ne.0) go to 100

      call deltek(linje,1,80,l,io)

      if(io.eq.0) then
          if(def) then
          do 50 i=1,n
   50     isv(i)=id(i)
          go to 500
        else
          call primvri('standardverdier er ikke tilgjengelige#')
          go to 100
        end if
      end if

      if(io.ne.n)  then
        call primvri('feil antall tall#')
        go to 100
      end if

      do 200 i=1,n
      call geti(linje,l(1,i),l(2,i),k1,k2,ires,iflag)
      if(iflag.eq.1) then
        isv(i)=ires
      else
        write(iustrr,*)'ulovlig svar:',linje(l(1,i):l(2,i))
        if(aktlog.eq.1 .and. aktkom.eq.1)
     %   write(ltape,*)charkom,'ulovlig svar',linje(l(1,i):l(2,i))
        go to 100
      end if
  200 continue

 500  continue

      dispens=0

      return
      end 

******************************************************************************
      subroutine lesr2(rd1,rd2,spor,rr1,rr2)
      real rd1,rd2,rr1,rr2
      character*80 spor
*---------------------------------------------------------------------------
      real rd(2),rsv(2)
      
      rd(1)=rd1
      rd(2)=rd2
      call lesvekr(rd,spor,rsv,2)
      rr1=rsv(1)
      rr2=rsv(2)
      return
      end

******************************************************************************
      subroutine lesr3(rd1,rd2,rd3,spor,rr1,rr2,rr3)
      real rd1,rd2,rd3,rr1,rr2,rr3
      character*80 spor
*---------------------------------------------------------------------------
      real rd(3),rsv(3)
      
      rd(1)=rd1
      rd(2)=rd2
      rd(3)=rd3
      call lesvekr(rd,spor,rsv,3)
      rr1=rsv(1)
      rr2=rsv(2)
      rr3=rsv(3)
      return
      end


******************************************************************************
      subroutine lesr4(rd1,rd2,rd3,rd4,spor,rr1,rr2,rr3,rr4)
      real rd1,rd2,rd3,rd4,rr1,rr2,rr3,rr4
      character*80 spor
*---------------------------------------------------------------------------
      real rd(4),rsv(4)
      
      rd(1)=rd1
      rd(2)=rd2       
      rd(3)=rd3   
      rd(4)=rd4
      call lesvekr(rd,spor,rsv,4)
      rr1=rsv(1)
      rr2=rsv(2)
      rr3=rsv(3)                
      rr4=rsv(4)
      return
      end



******************************************************************************
      subroutine lesvekr(rd,spor,rsv,n)
      integer n
      real rd(n),rsv(n)
      character*80 spor
      include 'styr.inc'
*----------------------------------------------------------------------------
      integer l(2,20),io ,ks ,i ,iflag ,k1,k2,kflag
      real rres
      logical def
      character*80 linje

      call xstrip(spor,ks,def)   
      if(ks.eq.0) ks=1

  100 continue

      if(def) then
        write(iustrr,*) spor(1:ks),'/',(rd(i) , i=1,n),'/'
        if(aktlog.eq.1 .and. aktkom.eq.1) then
          write(ltape,*) charkom,spor(1:ks),'/'
          do 105 i=1,n
          write(ltape,*)charkom,rd(i) 
 105      continue
          write(ltape,*)charkom, '/'
        end if
      else
        call primvri(spor)
      end if                                  


      call primles(linje,kflag)
      if(kflag.eq.2) stop
      if(kflag.gt.0) go to 100

      call deltek(linje,1,80,l,io)

      if(io.eq.0) then
          if(def) then
          do 50 i=1,n
   50     rsv(i)=rd(i)
          go to 500
        else
          call primvri('standardverdier er ikke tilgjengelige#')
          go to 100
        end if
      end if

      if(io.ne.n)  then
        call primvri('feil antall tall i#')
          call primvri(linje)
         go to 100
      end if

      do 200 i=1,n
      call getr(linje,l(1,i),l(2,i),k1,k2,rres,iflag)
      if(iflag.eq.1) then
        rsv(i)=rres
      else
        write(iustrr,*)'ulovlig svar:',linje(l(1,i):l(2,i))
        if(aktlog.eq.1 .and. aktkom.eq.1)
     %  write(ltape,*)charkom,'ulovlig svar:',linje(l(1,i):l(2,i))
        go to 100
      end if
  200 continue
 500  continue

      dispens=0

      return
      end 

*************************************************************************
      subroutine beglog(blogg,ilogg)
      integer ilogg
      logical blogg
*--------------------------------------------------------
      integer kf
      character*20 lnavn
      logical ja

      call initkom
      call setkomm('!')
      blogg=ja(blogg,'skal dialogen protokolleres?#')
      if(blogg) then
        kf=20
        call utfil(ilogg,'gi log-fil#','indat#',lnavn,kf)
        if(ilogg.gt.0) then
          write(ilogg,*) 'nei'
          call setlog(ilogg)
        else
          if(ilogg.eq.0) then          
             call setlog(ilogg)
          end if
        end if
      else
        ilogg=-2
      end if

      return
      end 



*************************************************************************
      subroutine beglogv(blogg,ilogg,defn)
      integer ilogg
      logical blogg
      character*80 defn
*--------------------------------------------------------
      integer kf
      character*20 lnavn
      logical ja

      call initkom
      call setkomm('!')
      blogg=ja(blogg,'skal dialogen protokolleres?#')
      if(blogg) then
        kf=20
        call utfil(ilogg,'gi log-fil#',defn,lnavn,kf)
        if(ilogg.gt.0) then
          write(ilogg,*) 'nei'
          call setlog(ilogg)
        else
          if(ilogg.eq.0) then          
             call setlog(ilogg)
          end if
        end if
      else
        ilogg=-2
      end if

      return
      end 


*************************************************************************
*
*       Denne gj|r det samme som beglog men skriver intet p} fila og
*       kaller heller ikke initkom. Rutinen er beregnet for bruk
*       etter at dialogen er startet.
************************************************************************
      subroutine openlog(blogg,ilogg)
      integer ilogg
      logical blogg
*--------------------------------------------------------
      integer kf
      character*20 lnavn
      logical ja

      call setkomm('!')
      blogg=ja(blogg,'skal dialogen protokolleres?#')
      if(blogg) then
        kf=20
        call utfil(ilogg,'gi log-fil#','indat#',lnavn,kf)
        if(ilogg.gt.0) then
          call setlog(ilogg)
        else
          if(ilogg.eq.0) then          
             call setlog(ilogg)
          end if
        end if
      else
        ilogg=-2
      end if

      return
      end 


************************************************************************
      subroutine openlogv(blogg,ilogg,defn)
      integer ilogg
      logical blogg
      character*80 defn
*--------------------------------------------------------
      integer kf
      character*20 lnavn
      logical ja

      call setkomm('!')
      blogg=ja(blogg,'skal dialogen protokolleres?#')
      if(blogg) then
        kf=20
        call utfil(ilogg,'gi log-fil#',defn,lnavn,kf)
        if(ilogg.gt.0) then
          call setlog(ilogg)
        else
          if(ilogg.eq.0) then          
             call setlog(ilogg)
          end if
        end if
      else
        ilogg=-2
      end if

      return
      end 



****************************************************************************
*      E X P A N D
*
*   Expanderer filnavn pa c-shell vis dvs. at feks:
*
*        ~/h.dat      gir  $HOME/h.dat
*        ~rupert/foo  gir  /mn/kastor/home/u1/rupert/foo dersom home-dir
*                      er /mn/kastor/home/u1/meg
*        ../kla.bla   gir kla.bla ett direktori opp i forhold til default
*
*  'Initkom' m} v{re kalt f|r f|rste kall p} 'expand'
*  parametere:
*       t(1:nt)  - tekst som skal expanderes                          I
*       tut(1:nut) - expandert tekst                                  O
*       iflag - feilparameter, iflag >0 betyr ulovlig innhold i t     O
*
**************************************************************************
      subroutine expand(t,nt,tut,nut,iflag)
      integer nt,nut,iflag
      character*120 t,tut
      include 'styr.inc'
*--------------------------------------
      integer mg,nb,kflag

      iflag=0

      if(t(1:1).ne.'.' .and. t(1:1).ne.'~') then
        nut=nt
        tut(1:nt)=t(1:nt)
        return
      end if

      if(t(1:1).eq.'~') then
        mg=nhd
        if(t(2:2).eq.'/') then
          nb=3
          if(nt.le.2) then
            iflag=1
            return
          end if
        else
          nb=2
          call fsok(homedir,nhd,mg,'/',-1,kflag)
        end if

        nut=mg+nt+1-nb       
        tut(1:mg)=homedir(1:mg)
        tut(mg+1:nut)=t(nb:nt)
        return
      end if

c  naa er . foerste tegn


      mg=ndd
      nb=1

 100  continue

      if(t(nb:nb).ne.'.' .or. nb.eq.nt) go to 200
      if(t(nb+1:nb+1).eq.'/') then
        nb=nb+2
        if(nb.gt.nt) then
          iflag=2
          return
        end if
        go to 100
      end if
    
      if(t(nb+1:nb+1).ne.'.') then
        go to 200
      else
        if(nb+1.eq.nt) then
          iflag=3
          return
        end if
        if(t(nb+2:nb+2).eq.'/') then
          if(nb+2.eq.nt) then
            iflag=4
            return
          end if
          call  fsok(defdir,ndd,mg,'/',-1,kflag)
          nb=nb+3
        else
          go to 200
        end if
      end if
 
      go to 100

 200  continue

      nut=mg+nt+1-nb
      tut(1:mg)=defdir(1:mg)
      tut(mg+1:nut)=t(nb:nt)
      return
      end




*************************************************************************
      SUBROUTINE INNFILNY(ITAPE,SPOR,DEFN,NAVN,KF,asc)
      INTEGER ITAPE,KF
      CHARACTER*80 SPOR,DEFN,NAVN
      logical asc
      include 'styr.inc'
*-----------------------------------------------------------------------
      CHARACTER*80 HNAVN
      LOGICAL FORK
      INTEGER KF0,kfex,iflag
                 
      KF0=KF
 100  continue
      CALL BLANK(HNAVN,KF0)
      CALL LESORD(DEFN,SPOR,HNAVN,KF)

      IF(KF.GT.KF0) THEN
        WRITE(IUSTRR,11) ' ',KF,KF0
        if(aktlog.eq.1 .and. aktkom.eq.1)
     %    WRITE(ltape,11) charkom,KF,KF0
  11    FORMAT(1X,a1,'LENGDE:',I4,'  MAX. LENGDE:',I4)
        call primvri('GI NYTT OG KORTERE NAVN#')
        GO TO 100
      ELSE       
        CALL BLANK(NAVN,KF0)
        NAVN(1:KF)=HNAVN(1:KF)
      END IF

      IF(FORK('STOPP#',NAVN,1,KF)) THEN
        ITAPE=-1
        RETURN
      END IF


      IF(FORK('TTY#',NAVN,1,KF).and. asc) THEN
        ITAPE=5
        RETURN
      END IF

      IF(FORK('TERMINAL#',NAVN,1,KF) .and. asc) THEN
        ITAPE=5
        RETURN
      END IF   

      IF(FORK('STANDARD#',NAVN,1,KF) .and. asc) THEN
        ITAPE=5
        RETURN
      END IF   

      call expand(navn,kf,hnavn,kfex,iflag)
      if(iflag.gt.0) then
        call primvri('ulovlig filnavn#')
        go to 100
      end if

      if(kfex.le.kf0) then
        kf=kfex
        navn(1:kf)=hnavn(1:kf)
      end if



      if(asc) then
        open(UNIT=itape,file=navn(1:kf),ERR=117,status='old',
     %  form='formatted')
      else
        open(UNIT=itape,file=navn(1:kf),ERR=117,status='old',
     %  form='unformatted')
      end if
      return
  117 call primvri('fil ikke funnet#')
      call primvri('gi nytt navn ("stopp" dersom du gir opp)#')
      go to 100

      end


*************************************************************************
      SUBROUTINE INNFIL(ITAPE,SPOR,DEFN,NAVN,KF)
      INTEGER ITAPE,KF
      CHARACTER*80 SPOR,DEFN,NAVN
      include 'styr.inc'
*-----------------------------------------------------------------------
      CHARACTER*80 HNAVN
      LOGICAL FORK
      INTEGER KF0
                 
      KF0=KF
 100  continue
      CALL BLANK(HNAVN,KF0)
      CALL LESORD(DEFN,SPOR,HNAVN,KF)

      IF(KF.GT.KF0) THEN
        WRITE(IUSTRR,11) ' ',KF,KF0
        if(aktlog.eq.1 .and. aktkom.eq.1)
     %   WRITE(ltape,11) charkom,KF,KF0
  11    FORMAT(1X,a1,'LENGDE:',I4,'  MAX. LENGDE:',I4)
        call primvri('GI NYTT OG KORTERE NAVN#')
        GO TO 100
      ELSE       
        CALL BLANK(NAVN,KF0)
        NAVN(1:KF)=HNAVN(1:KF)
      END IF

      IF(FORK('STOPP#',NAVN,1,KF)) THEN
        ITAPE=-1
        RETURN
      END IF


      IF(FORK('TTY#',NAVN,1,KF)) THEN
        ITAPE=5
        RETURN
      END IF

      IF(FORK('TERMINAL#',NAVN,1,KF)) THEN
        ITAPE=5
        RETURN
      END IF   

      IF(FORK('STANDARD#',NAVN,1,KF)) THEN
        ITAPE=5
        RETURN
      END IF   


      open(UNIT=itape,file=navn(1:kf0),ERR=117,status='old')
      return
  117 call primvri('fil ikke funnet#')
      call primvri('gi nytt navn ("stopp" dersom du gir opp)#')
      go to 100

      end


*************************************************************************
      SUBROUTINE UTFIL(ITAPE,SPOR,DEFN,NAVN,KF)
      INTEGER ITAPE,KF
      CHARACTER*80 SPOR,DEFN,NAVN
      include 'styr.inc'
*-----------------------------------------------------------------------
      CHARACTER*80 HNAVN
      LOGICAL FORK
      INTEGER KF0
                 
      KF0=KF

  100 CONTINUE               

      CALL BLANK(HNAVN,80)                         
      CALL LESORD(DEFN,SPOR,HNAVN,KF)

      IF(KF.GT.KF0) THEN
        WRITE(IUSTRR,11) ' ',KF,KF0
        if(aktlog.eq.1 . and. aktkom.eq.1)
     %        WRITE(ltape,11)charkom,KF,KF0
  11    FORMAT(1X,a1,'LENGDE:',I4,'  MAX. LENGDE:',I4)
        call primvri('GI NYTT OG KORTERE NAVN#')
        GO TO 100
      ELSE       
        CALL BLANK(NAVN,KF0)
        NAVN(1:KF)=HNAVN(1:KF)
      END IF

      IF(FORK('INGEN#',NAVN,1,KF)) THEN
        ITAPE=-1
        RETURN
      END IF

      IF(FORK('TTY#',NAVN,1,KF)) THEN
        ITAPE=iustrr
        RETURN
      END IF

      IF(FORK('TERMINAL#',NAVN,1,KF)) THEN
        ITAPE=iustrr
        RETURN
      END IF   

      IF(FORK('ERROR#',NAVN,1,KF)) THEN
        ITAPE=iustrr
        RETURN
      END IF   

      IF(FORK('OUTPUT#',NAVN,1,KF)) THEN
        ITAPE=6
        RETURN
      END IF   

      open(UNIT=itape,file=navn(1:kf),status='unknown')
      return

      end


************************************************************
      subroutine lukk(itape)
      integer itape
      include 'styr.inc'

      if(itape.ge.0 .and. itape.ne.5 .and. itape.ne.6 
     %    .and. itape.ne.iustrr) then
        close(itape)
        itape=-1
      end if
      return
      end



***********************************************************
*   
*        NUMOPP
*  
*   Aapner fil med navn avledet av stamme pluss nummer.
*   feks. vil nummopp('gagga#',17,23) aapne fil med navn
*   gagga17 og enhet 23.
*   Parametere:
*    stem - navnestamme, avsluttet med # eller !    I
*    isyk - nummer                                    I   
*    itape - unitnummer                              I
********************************************************
      subroutine nummop(stem,isyk,itape)
      integer isyk,itape
      character*80 stem
*--------------------------------------------------------
      integer kf,kf2
      character*80 tall,utfil
      logical def

      call itconv(isyk,kf,tall)
      call xstrip(stem,kf2,def)
      utfil(1:kf2)=stem(1:kf2)
      utfil(kf2+1:kf2+kf)=tall(1:kf)
      kf=kf2+kf

      open(unit=itape,file=utfil(1:kf),status='unknown')
      return
      end

***********************************************************
*   
*        NUMMOPG
*  
*   Aapner fil med navn avledet av stamme pluss nummer.
*   feks. vil nummopp('gagga#',17,23,.true.) aapne ascii fil med navn
*   gagga17 og enhet 23.
*   Parametere:
*    stem - navnestamme, avsluttet med # eller !    I
*    isyk - nummer                                    I   
*    itape - unitnummer                              I
*    asc - true gir asci, false binaert
********************************************************
      subroutine nummopg(stem,isyk,itape,asc)
      integer isyk,itape
      character*80 stem
      logical asc
*--------------------------------------------------------
      integer kf,kf2
      character*80 tall,utfil
      logical def

      call itconv(isyk,kf,tall)
      call xstrip(stem,kf2,def)
      utfil(1:kf2)=stem(1:kf2)
      utfil(kf2+1:kf2+kf)=tall(1:kf)
      kf=kf2+kf

      if(asc) then
        open(unit=itape,file=utfil(1:kf),status='unknown')
      else
        open(unit=itape,file=utfil(1:kf),status='unknown',
     %       form='unformatted')
      end if
      return
      end




***********************************************************
*   
*        NOPPG
*  
*   Aapner fil med navn avledet av stamme pluss nummer.
*   feks. vil nummopp('gagga#',17,23,.true.) aapne ascii fil med navn
*   gagga17 og enhet 23.
*   Parametere:
*    stem - navnestamme, avsluttet med # eller !    I
*            med ! legges ikke isyk til navn
*    isyk - nummer                                    I   
*    itape - unitnummer                              I
*    asc - true gir asci, false binaert
********************************************************
      subroutine noppg(stem,isyk,itape,asc)
      integer isyk,itape
      character*80 stem
      logical asc
*--------------------------------------------------------
      integer kf,kf2
      character*80 tall,utfil
      logical def

      call xstrip(stem,kf2,def)
      utfil(1:kf2)=stem(1:kf2)

      if(def) then
        call itconv(isyk,kf,tall)
        utfil(kf2+1:kf2+kf)=tall(1:kf)
        kf=kf2+kf
      else
         kf=kf2
      end if


      if(asc) then
        open(unit=itape,file=utfil(1:kf),status='unknown')
      else
        open(unit=itape,file=utfil(1:kf),status='unknown',
     %       form='unformatted')
      end if
      return
      end














**************************************************************************
*              MENYSYSTEM.
*   De etterf|lgende rutiner styrer oppbygning og bruk av et kommando-
*   system bygd opp som et hierarki av menyer. En kommandopost ser
*   ut som   <ordre> <spes1> <spes2>...<spesN> der <ordre> er en
*   sekvens som begynner med en bokstav og <spes1> etc. er tall.
*   <ordre> legges inn med rutinen 'setordre' og spesifikasjonene
*   (spes1 etc) opprettes med 'setspes'. Alle ordrer m} h|re
*   hjemme i en meny. Menyer opprettes med 'setmen' og en meny
*   m} v{re opprettet f|r den tilegnes ordre.
*    
*   Menyene er organisert i ett eller flere hierarkier ved at hver
*   ordre i en meny tilordnes en peker til en undermeny. Kommandoene
*   kan kjedes sammen etter m|nsteret:
*        <ordre1><spes...> <ordre2><spes...> <ordre3><spes..>
*   der ordre2 tilh|rer undermenyen til ordre1 etc. Spesifikasjonene
*   kan utelates dersom standardverdier er tilgjengelige. Dialogen
*   som er styrt ved 'kom' kan settes i sekvensiell mode ved kall p}
*   'setsekvens'. Dvs. at ordrer sp|rres om og leses en etter en,
*    og at spesifikasjonene gis p} egen linje etter at ordren de tilh|rer
*   er lest.
*
*   Noen punkter:
*       $ 'initmeny' m} kalles f|r noen andre rutiner
*       $ deler av et svar eller kommandopost adskilles alltid  av
*         en eller flere blanke.
*       $ ingen ordrer kan v{re over 80 lange.
*       $ hele kommandoskjeder og ikke bare neste kommando kan velges
*         standard dersom siste linje lest fra terminalen tar slutt
*         f|r noen feil er funnet eller en ordre uten undermeny er lest.
*         Dette gjelder ikke dersom dialogen er i sekvensiell mode.
*         V{r forsiktig med standardverdier.
*       $ Systemet bygger p} dia-rutinene for lesing og behandling
*         av tekst og tall.
*       $ All behandling av systemets registere er ment } foreg}
*         gjennom kall p} styrerutiner.
***************************************************************************
*---------------------------------------------------------------------------








***************************************************************************
*
*              S E T S E K V E N S
*
*     Dersom argument 'a' er sann settes dialogen til } v{re 
*     sekvensiell, dvs. at kommandoer bare kan gis en om gangen og at
*     spesifikasjoner gis separat etter kommandoen. 'a' lik usann setter
*     dialogen til simultan.
*************************************************************************
      subroutine setsekvens(a)
      logical a
      include 'men.inc'
*---------------------------------------------------------
      sekvens=a
      return
      end


***************************************************************************
*
*                  I N I T M E N Y 
*
*     Denne initialiserer menysystemets registere og m} kalles
*     f|r noen annen rutine som ber|rer systemet.
***************************************************************************
      subroutine initmeny
      include 'men.inc'
*-------------------------------------------------------------------------
      integer i

      sekvens=.false.

      eksmen(1)='HJELP# gir denne utskrift'
      eksmen(2)='LISTE# lister alle ordrer p} dette niv}et'
      eksmen(3)='AVBRYT# kommandosekvens avbrytes og startes p} nytt'
      eksmen(4)='VIS_SPESIFIKASJON# spes. for etterf. ordre listes'

      do 100 i=1,mordre
      spes(i)='#'
      menpek(i)=0
  100 continue

      do 200 i=1,mmen
      meny(i)=0
      menl(i)=0
  200 continue

      iopptatt=0
      menant=0 
      iint=0
      ireal=0
      ifeil=0
      ikommando=0
      menstart=1
              
      return
      end    

****************************************************************************
*
*               S P E S S J E K K
*
*     Rutinen plasserer spesifikasjoner fra teksten 't' inn i registere
*     tilh|rende kommando 'ikomm'.
*     parametere:                                        
*            ikomm    -   kommandonummer (globalt)                    I 
*            t        -   tekststreng inneholdene tallspes.           I
*            l        -   array med posisjoner av ord i 't'.          I
*                         l gis verdier ved kall p} 'deltek'
*            iseq     -   antall sekvenser i 't'                      I
*            ibeg     -   f|rste sekvens som tolkes som spesifik.     I
*            iend     -   siste sekvens som er behandlet              O
*            iflag    -   feilparameter                               O
*                       =0    alt  i orden
*                       =1    ikke eksterende standardverdier fors|kt
*                             benyttet.
*                       =2    ikke nok korrekte tall
*                       =3    for mange tall
*****************************************************************************

      subroutine spessjekk(ikomm,t,l,iseq,ibeg,iend,iflag)
      include 'men.inc'
      integer ikomm,iseq,ibeg,iend,iflag
      integer l(2,iseq)
      character*80 t
*--------------------------------------------------------------------------
      integer ihelt,iflyt,ierr,ibuff(10),iant,na,nb,ia,ib
      integer ix,i,i0
      real rbuff(10),rx
      character c
      logical def,nomore,digit,bkast
                            
      iflag=0

      def=nrordre(ikomm,1).gt.0 

      call xstrip(spes(ikomm),iant,bkast)

      if(ibeg.gt.iseq) then
        nomore=.true.      
      else                
        na=l(1,ibeg)
        c=t(na:na)
        nomore=(.not.digit(c))  .and.  c.ne.'+'  .and.  c.ne.'-' 
     %          .and.  c.ne.'.'
      end if

      if(nomore) then
        if(.not.(nrordre(ikomm,1).gt.0  .or. iant.eq.0)) iflag=1
        iend=ibeg-1
        return
      end if
                 
      if(ibeg+iant-1.gt.iseq) then
        iflag=2
        iend=iseq
        return
      end if


      ihelt=0
      iflyt=0

      do 100 i=1,iant  
        c=spes(ikomm)(i:i)
        na=l(1,ibeg+i-1)
        nb=l(2,ibeg+i-1)

        if(c.eq.'I') then
          call geti(t,na,nb,ia,ib,ix,ierr)
          if(ierr.ne.1) then
            iflag=2
            iend=ibeg-1
            return
          end if
          ihelt=ihelt+1
          ibuff(ihelt)=ix
        else  
          call getr(t,na,nb,ia,ib,rx,ierr)
          if(ierr.ne.1) then
            iflag=2
            iend=ibeg-1
            return
          end if 
          iflyt=iflyt+1
          rbuff(iflyt)=rx 
        end if

  100 continue
 

      if(iseq.gt.(ibeg+iant-1)) then
        na=l(1,ibeg+iant)
        c=t(na:na)
        if(digit(c)) then
          iflag=3
          iend=ibeg-1
          return
        end if
      end if
  

      iend=ibeg+iant-1

      i0=nrordre(ikomm,2)-1
      do 200 i=1,ihelt
        ispes(i0+i)=ibuff(i)
  200 continue

      i0=nrordre(ikomm,3)-1
      do 300 i=1,iflyt
        rspes(i0+i)=rbuff(i)
  300 continue

      if(nrordre(ikomm,1).gt.-1) nrordre(ikomm,1)=nrordre(ikomm,1)+1

      
      return
      end


****************************************************************************
*
*                    S E T M E N
*
*     Definerer en ny meny.
*     parametere:
*          spor  -  sp|rsm}l som stilles ved innlesning av ordre i     I
*                   menyen. Avsluttes med: '#' dersom standardverdi
*                   tillates, med '!' ellers.                              
*          nomax -  maksimalt antall ordrer i menyen                   I
*          numm  -  menynummer                                         O
*          iflag -  dersom meny er opprettet gis denne verdi 0         O     
*****************************************************************************
      subroutine setmen(spor,nomax,numm,iflag)
      include 'men.inc'
      integer nomax,numm,iflag
      character*80 spor
*--------------------------------------------------------------------------
      integer kf
      logical def ,bkast
          
      iflag=0
        
      numm=menant+1                             
      if(numm.gt.mmen .or. (iopptatt+nomax).gt.mordre) then
        iflag=1
        return
      end if                                    

      menant=numm
      meny(numm)=iopptatt+1                      
      menl(numm)=0
      iopptatt=iopptatt+nomax
     
      call xstrip(spor,kf,def)
      if(def) then
        mendef(numm)=0
      else
        mendef(numm)=-1
      end if

      call blank(menspor(numm),80)
      menspor(numm)(1:kf+1)=spor(1:kf+1)

      return
      end        

**************************************************************************** 
*
*                 S E T O R D R E
*
*     Definerer en ordre i en meny
*     parametere:
*          spor - selve ordren. Dersom denne skal kunne brukes som       I
*                 standardverdi m} den avsluttes med '#', hvis ikke
*                 med '!' 
*          bes  - beskrivelse av ordre, avsluttes med '#'. 'spor'        I
*                 pluss 'bes' b|r ikke inneholde mer enn 80 tegn
*                 tilsammen.
*          imen - Meny ordren skal tilh|re.                              I
*          nestemen - nummer for ordrens undermeny. Dersom verdien er    I
*                     0 eller -1 har ordren ingen undermeny. 0 tillater
*                     at kommandoposten etterf|lges av en rest som taes
*                     ut som output av 'kom'. -1 krever at posten
*                     etterf|lges av blanke.
*          numm - Globalt nummer tildelt ordren                          O
*          iflag - Feilparameter.                                        O
*                =0 alt i orden
*                =1 'imen' for stor
*                =2 'imen' svarer ikke til eksisterende meny
*                =3 antall ordrer tilh|rende 'imen' blir for stort
***************************************************************************** 
*
      subroutine setordre(spor,bes,imen,nestemen,numm,iflag)
      include 'men.inc'
      integer imen,nestemen,numm,iflag
      character*80 spor,bes
*--------------------------------------------------------------------------
      integer kf,ierr,kf1,nomax
      logical def ,bkast
          
      iflag=0
         
      if(imen.gt.mmen) then
        iflag=1
        return
      end if
  
      if(meny(imen).eq.0) then
        iflag=2
        return
      end if
                               
      if(imen.eq.mmen) then
        nomax=mordre                            
      else
        if(meny(imen+1).eq.0) then
          nomax=mordre
        else
          nomax=meny(imen+1)-1
        end if
      end if
      numm= meny(imen)+menl(imen)                    
      if(numm.gt.nomax) then
        iflag=3
        return
      end if                                    

      menl(imen)=menl(imen)+1
     
      call xstrip(spor,kf,def)

      call blank(ordre(numm),80)
      ordre(numm)(1:kf+1)=spor(1:kf+1)
      call xstrip(bes,kf1,bkast)
      if(kf1+kf+1.gt.80) kf1=78-kf
      if(kf1.gt.0) ordre(numm)(kf+2:kf1+2+kf)  = bes(1:kf1)

      spes(numm)='#'                            
      menpek(numm)=nestemen

      return
      end           

****************************************************************************
*
*               S E T S P E S
*
*     Setter spesifikasjon for en kommando.
*     parametere:
*           ikomm   -  kommandonummer (globalt)                          I
*           spesi   -  angir hvordan spesifikasjoner skal leses.         I
*                      Strengen 'IRRI#' betyr feks. at tallene gis som:
*                      helt reelt reelt helt. Avsluttes med '#' vil
*                      sist innlagte tallsett v{re standardverdier,
*                      avsluttes med '!' vil standardverdier aldri bli
*                      benyttet. Maksimalt antall spesifikasjoner er 9.
*           spor    -  sp|rsm}l som stilles etter spesifikasjon.         I
*                      avsluttes med '#'
*           bes     -  beskrivelse av sp|rsm}l. lengde av spor og        I
*                      bes b|r i sum ikke overskride 80. Avsluttes med
*                      '#'
*           iflag   -  feilparameter                                     O
*                   =  0     alt i orden
*                   =  1     for lang spesifikasjon
*                   =  2     ulovlig spesifikasjon
*                   =  3     registere er fulle
*******************************************************************************
      subroutine setspes(ikomm,spesi,spor,bes,iflag)
      include 'men.inc'
      integer ikomm,iflag
      character*10 spesi
      character*80 spor,bes
*--------------------------------------------------------------------------
      integer ihelt,iflyt,kf,ierr,kf1
      logical def ,bkast
      character*10 spesil
          
      iflag=0
             
      call xstrip(spesi,kf,def)
                        
  
      if(kf.gt.9) then
        iflag=1
        return
      end if            

      spesil(1:kf+1)=spesi(1:kf+1)
      if(kf.gt.0) call upchar(spesil,1,kf)

      call suspes(spesil,ihelt,iflyt,ierr)
                                 
      if(ierr.gt.0) then
        iflag=2
        return
      end if
            
      if(ireal+iflyt.gt.mspes  .or. iint+ihelt.gt.mspes) then
        iflag=3
        return
      end if

      if(def) then
        nrordre(ikomm,1)=0
      else
        nrordre(ikomm,1)=-1
      end if
                     
                 
      spes(ikomm)(1:kf+1)=spesil(1:kf+1)
      call xstrip(spor,kf,bkast)
      call xstrip(bes,kf1,bkast)
      call blank(spesspor(ikomm),80) 
      spesspor(ikomm)(1:kf+1)=spor(1:kf+1)
      if(kf+kf1+1 .gt.80) kf1=78-kf
      if(kf1.gt.0) spesspor(ikomm)(kf+2:kf+2+kf1)=bes(1:kf1)
      nrordre(ikomm,2)=iint+1
      nrordre(ikomm,3)=ireal+1
      ireal=ireal+iflyt                                                   
      iint=iint+ihelt
 
      return
      end
                                                                       


************************************************************************       
*
*               S U S P E S                                               
*
*     Teller antall flytende og heltall svarende til en streng med 
*     spesifikasjoner.
*        spesil - streng med spesifikasjoner, avsluttes med '#'       I
*                 eller '!' Maks antall=9.
*        ihelt,iflyt - antall av hhv hele og reelle tall.             O
*        iflag  - feilparameter, verdi=1 betyr at det er feil i       O
*                 spesifikasjonene.
************************************************************************
      subroutine suspes(spesil,ihelt,iflyt,iflag)
      integer ihelt,iflyt,iflag
      character*10 spesil
*---------------------------------------------------------------------
      integer kf,i             
      character c  
      logical bkast
                
      iflag=0

      call xstrip(spesil,kf,bkast)
                   
      ihelt=0
      iflyt=0

      do 100 i=1,kf  
        c=spesil(i:i)
        if(c.eq.'I') then
          ihelt=ihelt+1
        else  
          if(c.eq.'R') then
            iflyt=iflyt+1
          else
            iflag=1
            return
          end if
        end if

  100 continue
 
      return
      end



***************************************************************************
*
*                   S P O R S P E S
*
*         Sp|r om spesifikasjoner til en kommando og plasserer svaret i
*         registere.
*         parametere:
*            ikomm   -  kommandonummer (globalt)                         I
*            iflag   -  dersom verdien er 1 har brukeren avbrutt         O
*                       innlesning og registere er ikke endret.
***************************************************************************
      subroutine sporspes(ikomm,iflag)
      include 'men.inc'
      integer ikomm,iflag
*--------------------------------------------------------------------------
      integer l(2,20),kf,kf1,iseq,ib,ie,ierr
      character*80 t,sp 
      logical bkast,fork
      integer iustrr

      call stdrror(iustrr)
                       
      iflag=0

      call xstrip(spes(ikomm),kf1,bkast)
      if(kf1.eq.0) return
      sp(1:80)=spesspor(ikomm)
      call xstrip(sp,kf,bkast)
      if(kf.lt.68) then     
        sp(kf+1:kf+1)='['
        sp(kf+2:kf+1+kf1)=spes(ikomm)(1:kf1)
        sp(kf+2+kf1:kf+2+kf1)=']'
        if(nrordre(ikomm,1).gt.0) then
          sp(kf+3+kf1:kf+3+kf1)='#'
        else
          sp(kf+3+kf1:kf+3+kf1)='!'
        end if
      end if
         
  100 continue

      call lestek('gamle verd.#',sp,t)
      call deltek(t,1,80,l,iseq) 

      if(iseq.eq.1) then
        if(fork('avbryt#',t,l(1,1),l(2,1))) then
          iflag=1
          return
        end if
      end if
                                         
      ib=1
      call spessjekk(ikomm,t,l,iseq,ib,ie,ierr)
      if(ierr.ne.0) then
        write(iustrr,*)'ukorrekt format'
        go to 100
      end if
 
      return
      end




****************************************************************************
*
*              I N U M
*
*     Funksjonen returnerer globalt kommando nummer til lokal kommando
*     'ilokal' i meny nr. 'imeny'. Returnert verdi lik null betyr at 
*     input svarer til en ikke-eksisterende kommando.
****************************************************************************
      function inum(imeny,ilokal)
      include 'men.inc'
      integer inum,imeny,ilokal
*-------------------------------------------------------------
      
      if(imeny.gt.mmen) then
        inum=0
        return
      end if

      if(ilokal.gt.menl(imeny)) then
        inum=0
        return
      end if
                    
      inum=meny(imeny)+ilokal-1
      return
      end


*************************************************************************
*
*          P U T T S P E S
*
*      Henter eller plasserer spesifikasjoner for en kommando.
*
*      parametere
*          ikomm - globalt kommandonummer                                I
*          ibuff,rbuff - arrayer som inneholder hhv. heltallsverdier     I/O
*                        og reelle verdier som hentes eller plasseres.
*          ihelt,iflyt - antall av hhv. hele og reelle tall              O
*          hent - dersom hent er sann hentes tall, hvis ikke plasseres   I
*                 tallene
****************************************************************************
      subroutine puttspes(ikomm,ibuff,rbuff,ihelt,iflyt,hent)
      include 'men.inc'
      integer ikomm,ibuff(10),ihelt,iflyt
      real rbuff(10)
      logical hent
*---------------------------------------------------------------------------
      integer i,i0,iflag
      character*10 sp
     
      sp(1:10)=spes(ikomm)(1:10)
      call suspes(sp,ihelt,iflyt,iflag)
      if(ihelt+iflyt.eq.0) return

      i0=nrordre(ikomm,2)
      do 100 i=1,ihelt
        if(hent) then
          ibuff(i)=ispes(i0+i-1)
        else
          ispes(i0+i-1)=ibuff(i)
        end if
  100 continue

      i0=nrordre(ikomm,3)
      do 200 i=1,iflyt
        if(hent) then
          rbuff(i)=rspes(i0+i-1)
        else
          rspes(i0+i-1)=rbuff(i)
        end if
  200 continue

      if((.not.hent) .and. nrordre(ikomm,1).ne.-1) 
     %     nrordre(ikomm,1)=nrordre(ikomm,1)+1

      return
      end


********************************************************************
*
*                 P R I N T S P E S
*
*    Rutinen printer spesifikasjoner for kommando med globalt nummer
*    ikomm.
*********************************************************************
      subroutine printspes(ikomm)
      include 'men.inc'
      integer ikomm
*------------------------------------------------------------------- 
      integer ihelt,iflyt,ith,rth,ibuff(10),itot,i
      integer kf
      real rbuff(10)
      character*10 sp                   
      logical bkast
      integer iustrr

      call stdrror(iustrr)

                 
      if(ikomm.gt.mordre) then
        write(iustrr,*)' henvist til for stort ordrenummer'
        return
      end if

      sp=spes(ikomm)              
      call puttspes(ikomm,ibuff,rbuff,ihelt,iflyt,.true.)

      call xstrip(ordre(ikomm),kf,bkast)

      if(kf.eq.0) then
        write(iustrr,*)' henvist til ikke-opprettet ordrenummer'
        return
      end if

      write(iustrr,*)'spes. for: /',ordre(ikomm)(1:kf),'/'

      itot=iflyt+ihelt     
      ith=0
      rth=0
      do 100 i=1,itot                         
        if(sp(i:i).eq.'I') then
          ith=ith+1
          write(iustrr,50) ibuff(ith)
        else
          rth=rth+1
          write(iustrr,60) rbuff(rth)
        end if
  100 continue
     
   50 format(1x,i10)
   60 format(14x,E15.6)
      return
      end



**************************************************************************
*
*              S P U N D E F
*
*     Rutinen fjerner defaultsetting for spesifik. til kommando 'ikomm'
**************************************************************************
      subroutine spundef(ikomm)
      include 'men.inc'
      integer ikomm
*------------------------------------------------------------------------
      integer iustrr

      call stdrror(iustrr)

      if(ikomm.gt.mordre) then
       write(iustrr,*)'fors|k p} behandling av for h|yt kommando-nr'
       return
      end if

      if(nrordre(ikomm,1).gt.0) nrordre(ikomm,1)=0

      return
      end


***************************************************************************
*
*                 defset
*
*  Rutinen setter en ordre til } v{re standardsvar i en meny.
*   parametere:
*       imen   - menynummer                                         I
*       iordre - globalt ordrenummer                                I
*       iflag  - feilparameter                                      O
*              =0   alt i orden
*              =1   menynummmer eksisterer ikke
*              =2   ordre eksisterer ikke i angitte meny
*              =3   ordre kan ikke benyttes som standard
***********************************************************************
      subroutine defset(imen,iordre,iflag)
      integer imen,iordre,iflag
      include 'men.inc'
*-----------------------------------------------------------------------
      integer ilok,kf
      logical def

      iflag=0

      if(imen.gt.menant .or. imen.le.0) then
        iflag=1
        return
      end if
 
      ilok=iordre+1-meny(imen)

      if(ilok.lt.0 .or. ilok.gt.menl(imen)) then
        iflag=2
        return
      end if

      call xstrip(ordre(iordre),kf,def)
      
      if(def) then
        mendef(imen)=ilok
      else
        iflag=3
      end if
      return
      end



**************************************************************************
*
*              M E U N D E F
*                                          
*     Rutinen fjerner defaultsetting for meny nr 'imen'
**************************************************************************
      subroutine meundef(imen)
      include 'men.inc'
      integer imen
*------------------------------------------------------------------------
      integer iustrr

      call stdrror(iustrr)

      if(imen.gt.mmen) then
       write(iustrr,*)'fors|k p} behandling av for h|yt meny-nr'
       return
      end if

      if(mendef(imen).gt.0) mendef(imen)=0

      return
      end
                                  

****************************************************************************
*
*           T O P P M E N
*
*     Setter meny nr 'imen' til } v{re |verste meny i hierarkiet
***************************************************************************
      subroutine toppmen(imen)
      include 'men.inc'
      integer imen
*------------------------------------------------------------------------
      integer iustrr

      call stdrror(iustrr)

      if(imen.gt.mmen) then
       write(iustrr,*)'fors|k p} behandling av for h|yt meny-nr'
       return
      end if

      if(menl(imen).eq.0) then     
       write(iustrr,*)'fors|k p} } sette tom meny som toppmeny'
       return
      end if                                             

      menstart=imen
      return
      end 





**************************************************************************
*
*              M E N Y L I S
*
*      Lister ut meny nr kmen
**************************************************************************
      subroutine menylis(kmen) 
      include 'men.inc'
      integer kmen
*-------------------------------------------------------------------------
      character*10 add
      integer i,istop,i0,kf,ii
      logical bkast
      integer iustrr

      call stdrror(iustrr)


      write(iustrr,*)' '
      write(iustrr,*)' '
      write(iustrr,*)' '
      write(iustrr,*)'meny nr:',kmen
      istop=menl(kmen)   
      i0=meny(kmen)

      do 100 i=1,istop
      ii=i0+i-1
      write(iustrr,*) ordre(ii)(1:78)
      call xstrip(spes(ii),kf,bkast)      
      if(kf.gt.0) then 
        kf=kf+1
        add(1:kf)=spes(ii)(1:kf)
      else
        kf=5
        add(1:kf)='ingen'
      end if
      write(iustrr,55) add(1:kf)
   55 format(1x,'          spesifikasjoner:',a10)
  100 continue
      return
      end

      
****************************************************************************
*
*                 K O M
*
*      leser en kommando sekvens og oppdaterer registere for
*      spesifikasjoner.
*      parametere:
*            iord - iord(1)...iord(ioant) inneholder ved utgang            O
*                   lokale nummere for de leste kommandoer.    
*            imen - imen(1)...imen(ioant) inneholder ved utgang            O
*                   menynummere for de leste kommandoer.    
*            ioant- antall kommandoer i sekvensen.                         O
*            rest - ubehandlet del av siste leste svar som eventuelt       O
*                   kan brukes/tolkes senere.
**************************************************************************** 
      subroutine kom(iord,imen,ioant,rest)
      integer iord(10),imen(10),ioant
      character*80 rest
      include 'men.inc'
*-------------------------------------------------------------------------
      integer kmen,i0,idef,kf,l(2,30),iseq,na,nb,ikt,i
      integer ktreff(20),iat,ibeg,iend,iflag,ikomm    
      integer inum
      character*80 sp,tdef,tk,presord
      logical bkast,fork,ja
      integer iustrr

      call stdrror(iustrr)

                  
ccccc      forklaring av tellere/variabler  ccccccccccccccccccccccccccccccc
c                                 
c     i0 -  absolutt posisjon av forste ordre i n}v{rende meny
c     kmen -   nummer p} n{v{rende meny
c     iseq - totalt antall sekvenser i tk
c     ibeg - neste sekvens som skal behandles.
c     iend - sekvens som sist er behandlet i spessjekk
c     ioant- angir fortl|pende antall lovlige kommandoer.
c     menstart - nummer for toppmeny i hierarkiet.
c     kf   - lengde av tekster
c     na,nb - grenser for sekvenser i tekster.
c     idef - nummer for standard ordre i n}v{rende meny.
c     iat  - antall treff i menyen. (er en dersom svaret er lovlig)
c     tk   - siste leste kommandolinje
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   50 continue  
c
c       Start p} ytterste l|kke. Her startes innlesning av kommando
c       fra begynnelsen. Dersom brukeren avbryter en kommando-givning 
c       g}es det til denne adressen
c
      ioant=0
      kmen=menstart


  100 continue
c
c      Start p} nest ytterste l|kke. Behandling av kommandosekvens 
c      fortsetter ved innl|esning av ny linje fra terminal.
c      Denne adressen opps|kes dersom siste linje lest fra terminal
c      er slutt eller inneholdt en feil. Delen av linjen forut for feilen
c      annuleres ikke.
c
      i0=meny(kmen)
      idef=mendef(kmen)  
      call xstrip(menspor(kmen),kf,bkast)
      if(kf.gt.0) sp(1:kf)=menspor(kmen)(1:kf)
      if(idef.gt.0) then
        tdef=ordre(i0-1+idef)
        sp(kf+1:kf+1)='#'
      else
        tdef='#'
        sp(kf+1:kf+1)='!'
      end if
  
      call lestek(tdef,sp,tk)
      call deltek(tk,1,80,l,iseq)
      if(sekvens .and. iseq.gt.1) then
        write(iustrr,*)'gi kommandoene en om gangen!'
        go to 100
      end if  
                                                         
      na=l(1,1)
      nb=l(2,1)
                        
      ibeg=1

  200 continue        
c
c       Innerste l|kke. Her behandles en kommando
c       og dennes spesifikasjoner pr gjennomgang.
c 
      idef=mendef(kmen)
      if(ibeg.gt.iseq) then
        if(idef.le.0) go to 100
        ktreff(1)=idef
      else
        na=l(1,ibeg)
        nb=l(2,ibeg)
        call balac(ordre(i0),menl(kmen),tk,na,nb,ktreff,iat)
        if(iat.eq.0 ) then
c         er den gitte kommando en hjelp-kommando? I s} fall gis hjelp.
          call balac(eksmen,4,tk,na,nb,ktreff,iat)                     
          if(iat.eq.1) then
            ibeg=ibeg+1
            if(ktreff(1).eq.4) then        
              if(ibeg.gt.iseq) then
                 call lestek('#','spesifikasjon av hva:!',tk)
                 call deltek(tk,1,80,l,iseq)
                 ibeg=1
              end if
              if(iseq.gt.ibeg) then
                call xstrip(eksmen(4),kf,bkast)
                if(kf.eq.0) kf=1
                write(iustrr,*)'/',eksmen(4)(1:kf),
     %                         '/  m} etterf|lges av en'
                write(iustrr,*)'enkelt ordre'
              else           
                na=l(1,ibeg)
                nb=l(2,ibeg)
                call balac(ordre(i0),menl(kmen),tk,na,nb,ktreff,iat)
                if(iat.ne.1) then
                  write(iustrr,*)'feil eller upresis ordre'
                else
                  ikt=inum(kmen,ktreff(1))
                  call printspes(ikt)
                end if
              end if       
              go to 100
            end if     
            if(iseq.gt.ibeg) then
              write(iustrr,*)'ulovlig spesifikasjon av hjelpe-kommando'
              go to 100
            end if
            if(ktreff(1).eq.3) go to 50
            if(ktreff(1).eq.1) then
              do 120 i=1,4
              write(iustrr,*)eksmen(i)(1:80)
 120          continue
            end if
            if(ktreff(1).eq.2) call menylis(kmen)
            go to 100
          end if
        end if                            
        
        if(iat.eq.0) then
          write(iustrr,*)'ulovlig kommando:/',tk(na:nb),'/'
          go to 100
        end if
        if(iat.gt.1) then
          write(iustrr,*)'flertydig fork.:/',tk(na:nb),'/'
          go to 100
        end if                      
              
      end if


      ioant=ioant+1
      if(ioant.gt.10) then
        write(iustrr,*)'kommando for lang'
        go to 50
      end if

      imen(ioant)=kmen
      iord(ioant)=ktreff(1)
      if(idef.gt.-1) then
        call xstrip(ordre(ktreff(1)+i0-1),kf,bkast) 
        if(bkast) mendef(kmen) = ktreff(1)
      end if
      ibeg=ibeg+1
      ikomm=i0+ktreff(1)-1
      presord=ordre(ikomm)      
      
      if(.not.sekvens)
     % call spessjekk(ikomm,tk,l,iseq,ibeg,iend,iflag)
      if(iflag.ne.0 .or. sekvens) then
        if(iflag.ne.0) then
          call xstrip(presord,kf,bkast)                  
          if(kf.eq.0) kf=1
          write(iustrr,*)'feil spesifikasjon for:/',presord(1:kf),'/'
        end if
        call sporspes(ikomm,iflag)
        if(iflag.ne.0) go to 50
      end if     

      ibeg=iend+1
      kmen=menpek(ikomm)
        
      if(kmen.gt.0) then
        i0=meny(kmen)
        if(sekvens) go to 100
        go to 200
      end if
                              
      rest='#'
      if(kmen.eq.0 .and.ibeg.le.iseq) then
        na=l(1,ibeg)
        nb=l(2,iseq)
        kf=nb-na+1
        rest(1:kf)=tk(na:nb)
        if(kf.lt.80) rest(kf+1:kf+1)='#'
      end if  

      if(kmen.lt.0 .and.ibeg.le.iseq) then
        na=l(1,ibeg)
        nb=l(2,iseq)
        write(iustrr,*) 'det er en rest av komm:/',tk(na:nb),'/'
        if(.not.ja(.true.,'ignoreres?#')) go to 50
      end if
     

      return
      end




**************************************************************************
      subroutine lokchar(t,c,i1,i2,k,eks)
      integer i1,i2,k
      character*80 t
      character c
      logical eks
*------------------------------------------------------------------------
      integer i

      if(i1.gt.i2) then
        eks=.false.
        return
      end if

      i=i1-1

 100  continue
      i=i+1

      if(t(i:i).eq.c) then
        k=i
        eks=.true.
        return
      end if

      if(i.eq.i2) then
        eks=.false.
        return
      end if

      go to 100

      end



***************************************************************************
*
*          t - i1:i2 inneholder ett ord.
*************************************************************************
      subroutine extension(t,i1,i2,ten,kf)
      integer i1,i2,kf
      character*80 t,ten
*-------------------------------------------------------------------------
      integer k,nmax,np,iflag
      logical eks
      character*80 tt

       nmax=i2-i1+1
       tt(1:nmax)=t(i1:i2)
c      call lokchar(t,'.',i1,i2,k,eks)
       np=nmax
       iflag=0
      if(tt(nmax:nmax).ne.'.') call fsok(tt,nmax,np,'.',-1,iflag)
      eks=iflag.eq.0 
      if(eks) then
        kf=nmax-np
        if(kf.gt.0) then
          ten(1:kf)=tt(np+1:nmax)
        end if
      else
        kf=-1
      end if

      return
      end



*************************************************************************
*
*                 FILGEN
*
*  Stiller sporsmaal paa standard error og leser filangivelse som spesifisert i
*  beskrivelse av rutinen 'ordfile'. Filnavn ekspanderes som beskrevet i
*   expand.
*  'Initkom' m} v{re kalt f|r f|rste kall p} filgen. I traad med vanlig 
*   praksis i dia rutiner ender inputstrenger med ! eller #.
*   parametere:
*        ITAPE - unitnummer                                        I/O
*                dersom noe gaar galt settes itape < 0
*        spor  - sporsmaal                                         I
*        defn  - standard filnavn, derom spor ender med ! neglisjeres
*                denne og svar maa gis.                            I
*        NAVN(1:kf) - filnavn                                       O
*        kf - se ovenfor                                            O
*        asc - = .true. angir asci, = .false. binaer               O
*        inp - = .true. angir inputfil, = .false. outputfil        I
*************************************************************************
      SUBROUTINE FILgen(ITAPE,SPOR,DEFN,NAVN,KF,asc,inp)
      INTEGER ITAPE,KF
      CHARACTER*80 SPOR,DEFN,NAVN
      logical asc,inp
      include 'styr.inc'
*-----------------------------------------------------------------------
      CHARACTER*80 HNAVN,FNA
      character*11 fmat
      LOGICAL FORK
      INTEGER KF0,iflag,kfex
                 
      kf0=kf
 100  continue


      CALL BLANK(HNAVN,80)
      CALL BLANK(FNA,80)
      CALL LESTEK(DEFN,SPOR,HNAVN)

      call ordfile(hnavn,1,80,fna,kf,asc,iflag)

      if(iflag.gt.0) then
        write(0,*)'iflag=',iflag
        call primvri('ulovlig eller ufullstendig fil-spes.#')
        go to 100
      end if

      IF(KF.GT.KF0) THEN
        WRITE(IUSTRR,11) ' ',KF,KF0
        if(aktlog.eq.1 .and. aktkom.eq.1)
     %     WRITE(ltape,11) charkom,KF,KF0
  11    FORMAT(1X,a1,'LENGDE:',I4,'  MAX. LENGDE:',I4)
        call primvri('GI NYTT OG KORTERE NAVN#')
        GO TO 100
      ELSE       
        CALL BLANK(NAVN,KF0)
        NAVN(1:KF)=fna(1:KF)
      END IF

      IF(FORK('STOPP#',NAVN,1,KF) .or. 
     %     fork('INGEN#',navn,1,kf) ) THEN
        ITAPE=-1
        RETURN
      END IF


      if(asc) then
        IF(FORK('TTY#',NAVN,1,KF) .OR. 
     %       FORK('TERMINAL#',NAVN,1,KF) ) THEN
          if(inp) then
            ITAPE=5
          else
            itape=iustrr
          end if
          RETURN
        END IF


        if(inp) then
          IF(FORK('STANDARD#',NAVN,1,KF) .or. 
     %        fork('input#',navn,1,kf) ) THEN
            ITAPE=5
            RETURN
          END IF  
        else
         IF(FORK('ERROR#',NAVN,1,KF)) THEN
           ITAPE=iustrr
           RETURN
         END IF   

         IF(FORK('OUTPUT#',NAVN,1,KF)) THEN
           ITAPE=6
           RETURN
          END IF   
       end if
          
      end if 

      if(asc) then
        fmat='formatted  '
      else
        fmat='unformatted'
      end if

      call expand(navn,kf,fna,kfex,iflag)
      if(iflag.gt.0) then
        call primvri('ulovlig filnavn#')
        go to 100
      end if

      if(kfex.le.kf0) then
        kf=kfex
        navn(1:kf)=fna(1:kf)
      end if

      if(inp) then
       open(UNIT=itape,file=fna(1:kfex),ERR=117,status='old',
     % form=fmat)
cc       write(0,*)'aaa, kfex=',kfex

        return
  117   call primvri('fil ikke funnet#')
      else
       open(UNIT=itape,file=fna(1:kfex),ERR=119,status='unknown'
     %      ,form=fmat)
cc       write(0,*)'bbb, kfex=',kfex
       return
 119   call primvri('fil kan ikke }pnes#')
      end if
      call primvri('gi nytt navn ("stopp" dersom du gir opp)#')
      go to 100

      end




***********************************************************************
*
*     Rutinen unders|ker en streng for } finne en filbskrivelse som
*     best}r av navn og type-angivelse. Type angivelse er en forkortelse,
*     minst 3 tegn lang, av 'ascii' eller 'binary'. Store og sm} bokstaver
*     behandles likt. 
*          Dersom to ord finnes i strengen tolkes det ene som
*     type, det andre som filnavn uavhengig av rekkef|lge. Dette betyr
*     at filnavn som er mulige typeangivelser vanligvis ikke kan benyttes.
*          Hvis strengen inneholder bare ett ord, tolkes dette som filnavn.
*     Type settes da som ascii med ett unntak: filnavn inneholder et og
*     bare ett punktum etterfulgt av type angivelse for bin{r fil.
*          eksempler
*   
*               streng              filnavn           type
*          
*             bin  a.dat             a.dat            bin{r
*             a.dat ascii            a.dat            ascii
*             bin asc                    ulovlig streng
*             a.dat                  a.dat            ascii
*             a.binary               a.binary         bin{r
*             a.bin.dat              a.bin.dat        ascii
*             .bin                    .bin            bin{r
*             asci a.bin             a.bin            ascii
*             ascii                      ulovlig streng
*
*     parametere:
*             t(i1:i2) - streng med filangivelse                     I
*               i1,i2  - angir del av tekst t som skal unders|kes    I
*             name(1:kf) - tekst som angir funnet filnavn            O
*             kf        - lengde av filnavn                          O
*             asc       - verdi .true. angir ascii, .false. bin{r    O
*             iflag     - feilparameter                              O
*                        =0  lovlig og fullstendig informasjon i streng
*                        =1  streng er tom
*                        =2  streng inneholder to ord, ingen kan tolkes 
*                            som type-betegnelse
*                        =3  streng inneholder 2 ord som begge er
*                            type-betegnelser
*                        =4,5  streng inneholder bare ett ord som er en
*                            betegnelse  for hhv ascii og bin{r.
*                        =6  Too many words
*******************************************************************************
      subroutine ordfile(t,i1,i2,name,kf,asc,iflag)
      integer i1,i2,kf,iflag
      character*80 t,name
      logical asc
*-----------------------------------------------------------------------
      integer l(2,20),iant,j,itreff(2),m,nfil,ntyp
      character*80 ten
      logical fork

      iflag=0
      call deltek(t,i1,i2,l,iant)

       if(iant.eq.0) then
        iflag=1
        return
      end if

      if(iant.gt.2) then
        iflag=6
        return
      end if

      if(iant.eq.1) then
        kf=l(2,1)-l(1,1)+1
        name(1:kf)=t(l(1,1):l(2,1))
        call extension(name,1,kf,ten,j)
        if(j.gt.0) then
          asc =.not.(fork('binary#',ten,1,j) .and. j.ge.3)
        else
          asc=.true.
          if(kf.ge.3) then
            if(fork('ascii#',name,1,kf)) iflag=4
            if(fork('binary#',name,1,kf)) iflag=5
          end if
        end if
        return
      end if


      do 100 m=1,2
      if(l(2,m)-l(1,m).lt.2) then
        itreff(m)=0
      else
        if(fork('ascii#',t,l(1,m),l(2,m)) ) then
          itreff(m)=1
        else
          if(fork('binary#',t,l(1,m),l(2,m)) ) then
            itreff(m)=2
          else
            itreff(m)=0
          end if
        end if
      end if
 100  continue

      if( itreff(1).eq.0 .and. itreff(2).eq.0) then
        write(0,'(a,a,a)')t(l(1,1):l(2,1)),' ',t(l(1,2):l(2,2))
        iflag=2
        return
      end if

      if( itreff(1).gt.0 .and. itreff(2).gt.0) then
        iflag=3
        return
      end if

      if(itreff(1).gt.0) then
        ntyp=1
        nfil=2
      else
        ntyp=2
        nfil=1
      end if

      asc=itreff(ntyp).eq.1
      kf=l(2,nfil)-l(1,nfil)+1
      name(1:kf)=t(l(1,nfil):l(2,nfil))
    
      return
      end



*************************************************************************
*
*                 FILOPN
*
*  AApner fil fra gitt filangivelse som beskrevet i expand og ordfile.
*  'Initkom' m} v{re kalt f|r f|rste kall p} filopn dersom expandering
*   skal virke.
*   parametere:
*        ITAPE - unitnummer                                        I
*        NAVN - filangivelse                                       I
*               eksempe paa lov. ang:  bin ~/hh.jj 
*        fulln(1:kf) -filnavn, dersom det er ekspandert inneholder   O
*                     fulln hele pathen.
*        kf - se ovenfor                                            O
*        asc - = .true. angir asci, = .false. binaer               O
*        inp - = .true. angir inputfil, = .false. outputfil        I
*        ierr - feilparameter                                      O
*              = 0          alt ok
*                1          gal angivelse av type ( bin/asc etc.)
*                2          feil under expandering
*                3          fil kan ikke aapnes
***********************************************************************
      SUBROUTINE FILopn(ITAPE,NAVN,fulln,kf,asc,inp,ierr)
      INTEGER ITAPE,kf,ierr
      CHARACTER*80 NAVN
      character*120 fulln
      logical asc,inp
      include 'styr.inc'
*-----------------------------------------------------------------------
      CHARACTER*120 FNA
      character*11 fmat
      character*7 fstat
      INTEGER iflag,kfa,kfn
      logical paelm
                 
      ierr=0

      call xstrip(navn,kfn,paelm)
      CALL BLANK(FNA,120)
      call ordfile(navn,1,kfn,fna,kfa,asc,iflag)

      if(iflag.gt.0) then
        ierr=1
        return
      end if



      if(asc) then
        fmat='formatted  '
      else
        fmat='unformatted'
      end if

      if(inp) then
        fstat='old    '
      else
        fstat='unknown'
      end if

      call expand(fna,kfa,fulln,kf,iflag)
      if(iflag.gt.0) then
        ierr=2
        return
      end if

      open(UNIT=itape,file=fulln(1:kf),ERR=119,status=fstat
     %      ,form=fmat)
      return
 119  continue
      ierr=3
      return
    

      end


********************************************************************
*
*                LEXYFIL
*                
*     Reads data from two-column file
*
*     fnavn  - Name on file                                            I
*     kl  -   lengtn of filename. If kl<=0 the name is assumed to be   I
*             ended by '#' or '!' and the length is calculated
*     x,y  -  arrays for data                                          O
*     nmax  - maximum number of data                                   I
*     n  -    number of data read                                      O
*     ope  -  set to false if file needs to be opened                  I
*     itape  - unitnumber                                              I
*     nskip  -  number of leading comment lines in file                O  
*     iflag  - error flag                                              O
*                0     : OK
*                1     : maximum number of data reached before end of file.
*                        The data that is read is returned in x and y
*                3     : error in comment processing ( no data ?)
*                30    : inappropriate filename 
*                50    : cannot open file
********************************************************************
      subroutine lexyfil(fnavn,kl,x,y,nmax,n,ope,itape,nskip,iflag)
      integer kl,nmax,n,itape,nskip,iflag
      real x(nmax),y(nmax)
      character*80 fnavn
      logical ope
*------------------------------------------------------------------
      integer kf,ierr
      real a,b
      logical def
 
      iflag=0

      if(kl.le.0) then
      call xstrip(fnavn,kf,def)
         if(kf.le.0) then
           iflag=30
           return
         end if         
      else
         kf=kl
      end if


      if( .not.ope) then
       open(unit=itape,file=fnavn(1:kf),err=117,status='old')
      end if


cc      call skipkom(nskip,itape,ierr)
      call slopskip(nskip,itape,ierr)
      if(ierr.gt.0) then
      write(0,*)'skipkom: nskip,ierr=',nskip,ierr
         iflag=3
         close(itape)
         return
      end if

      n=0

 100  continue

      read(itape,*,end=200)a,b
      n=n+1
      if(n.gt.nmax) then
        iflag=1
        close(itape)
        return
      end if
      x(n)=a
      y(n)=b

      go to 100

 200  continue
      close(itape)
      return

 117  continue
      iflag=50
      return

      end



********************************************************************
*
*                LEXYzFIL
*                
*     Reads data from three-column file
*
*     fnavn  - Name on file                                            I
*     kl  -   lengtn of filename. If kl<=0 the name is assumed to be   I
*             ended by '#' or '!' and the length is calculated
*     x,y,z  -  arrays for data                                          O
*     nmax  - maximum number of data                                   I
*     n  -    number of data read                                      O
*     ope  -  set to false if file needs to be opened                  I
*     itape  - unitnumber                                              I
*     nskip  -  number of leading comment lines in file                O  
*     iflag  - error flag                                              O
*                0     : OK
*                1     : maximum number of data reached before end of file.
*                        The data that is read is returned in x,y and z
*                3     : error in comment processing ( no data ?)
*                30    : inappropriate filename 
*                50    : cannot open file
********************************************************************
      subroutine lexyzfil(fnavn,kl,x,y,z,nmax,n,ope,itape,nskip,iflag)
      integer kl,nmax,n,itape,nskip,iflag
      real x(nmax),y(nmax),z(nmax)
      character*80 fnavn
      logical ope
*------------------------------------------------------------------
      integer kf,ierr
      real a,b,c
      logical def
 
      iflag=0

      if(kl.le.0) then
      call xstrip(fnavn,kf,def)
         if(kf.le.0) then
           iflag=30
           return
         end if         
      else
         kf=kl
      end if


      if( .not.ope) then
       open(unit=itape,file=fnavn(1:kf),err=117,status='old')
      end if


cc      call skipkom(nskip,itape,ierr)
      call slopskip(nskip,itape,ierr)
      if(ierr.gt.0) then
      write(0,*)'skipkom: nskip,ierr=',nskip,ierr
         iflag=3
         close(itape)
         return
      end if

      n=0

 100  continue

      read(itape,*,end=200)a,b,c
      n=n+1
      if(n.gt.nmax) then
        iflag=1
        close(itape)
        return
      end if
      x(n)=a
      y(n)=b
      z(n)=c

      go to 100

 200  continue
      close(itape)
      return

 117  continue
      iflag=50
      return

      end


********************************************************************
*
*                LEXFIL
*                
*     Reads data from one-column file
*
*     fnavn  - Name on file                                            I
*     kl  -   lengtn of filename. If kl<=0 the name is assumed to be   I
*             ended by '#' or '!' and the length is calculated
*     x  -  array for data                                          O
*     nmax  - maximum number of data                                   I
*     n  -    number of data read                                      O
*     ope  -  set to false if file needs to be opened                  I
*     itape  - unitnumber                                              I
*     nskip  -  number of leading comment lines in file                O  
*     iflag  - error flag                                              O
*                0     : OK
*                1     : maximum number of data reached before end of file.
*                        The data that is read is returned in x and y
*                3     : error in comment processing ( no data ?)
*                30    : inappropriate filename 
*                50    : cannot open file
********************************************************************
      subroutine lexfil(fnavn,kl,x,nmax,n,ope,itape,nskip,iflag)
      integer kl,nmax,n,itape,nskip,iflag
      real x(nmax)
      character*80 fnavn
      logical ope
*------------------------------------------------------------------
      integer kf,ierr
      real a
      logical def
 
      iflag=0

      if(kl.le.0) then
      call xstrip(fnavn,kf,def)
         if(kf.le.0) then
           iflag=30
           return
         end if         
      else
         kf=kl
      end if


      if( .not.ope) then
       open(unit=itape,file=fnavn(1:kf),err=117,status='old')
      end if


cc      call skipkom(nskip,itape,ierr)
      call slopskip(nskip,itape,ierr)
      if(ierr.gt.0) then
      write(0,*)'skipkom: nskip,ierr=',nskip,ierr
         iflag=3
         close(itape)
         return
      end if

      n=0

 100  continue

      read(itape,*,end=200)a
      n=n+1
      if(n.gt.nmax) then
        iflag=1
        close(itape)
        return
      end if
      x(n)=a

      go to 100

 200  continue
      close(itape)
      return

 117  continue
      iflag=50
      return

      end


********************************************************************
*
*                LENMFIL
*                
*     Reads data from two-column file
*
*     fnavn  - Name on file                                            I
*     kl  -   lengtn of filename. If kl<=0 the name is assumed to be   I
*             ended by '#' or '!' and the length is calculated
*     nx,ny  -  arrays for data                                        O
*     nmax  - maximum number of data                                   I
*     n  -    number of data read                                      O
*     ope  -  set to false if file needs to be opened                  I
*     itape  - unitnumber                                              I
*     nskip  -  number of leading comment lines in file                O  
*     iflag  - error flag                                              O
*                0     : OK
*                1     : maximum number of data reached before end of file.
*                        The data that is read is returned in x and y
*                3     : error in comment processing ( no data ?)
*                30    : inappropriate filename 
*                50    : cannot open file
********************************************************************
      subroutine lenmfil(fnavn,kl,nx,ny,nmax,n,ope,itape,nskip,iflag)
      integer kl,nmax,n,itape,nskip,iflag,nx(nmax),ny(nmax)
      character*80 fnavn
      logical ope
*------------------------------------------------------------------
      integer kf,ierr,na,nb
      logical def
 
      iflag=0

      if(kl.le.0) then
      call xstrip(fnavn,kf,def)
         if(kf.le.0) then
           iflag=30
           return
         end if         
      else
         kf=kl
      end if


      if( .not.ope) then
       open(unit=itape,file=fnavn(1:kf),err=117,status='old')
      end if


      call skipkom(nskip,itape,ierr)
      if(ierr.gt.0) then
      write(0,*)'skipkom: nskip,ierr=',nskip,ierr
         iflag=3
         close(itape)
         return
      end if

      n=0

 100  continue

      read(itape,*,end=200)na,nb
      n=n+1
      if(n.gt.nmax) then
        iflag=1
        close(itape)
        return
      end if
      nx(n)=na
      ny(n)=nb
 
      go to 100

 200  continue
      close(itape)
      return

 117  continue
      iflag=50
      return

      end


********************************************************************
*
*                LENFIL
*                
*     Reads data from a one-column file
*
*     fnavn  - Name on file                                            I
*     kl  -   lengtn of filename. If kl<=0 the name is assumed to be   I
*             ended by '#' or '!' and the length is calculated
*     nx  -  array for data                                        O
*     nmax  - maximum number of data                                   I
*     n  -    number of data read                                      O
*     ope  -  set to false if file needs to be opened                  I
*     itape  - unitnumber                                              I
*     nskip  -  number of leading comment lines in file                O  
*     iflag  - error flag                                              O
*                0     : OK
*                1     : maximum number of data reached before end of file.
*                        The data that is read is returned in x and y
*                3     : error in comment processing ( no data ?)
*                30    : inappropriate filename 
*                50    : cannot open file
********************************************************************
      subroutine lenfil(fnavn,kl,nx,nmax,n,ope,itape,nskip,iflag)
      integer kl,nmax,n,itape,nskip,iflag,nx(nmax)
      character*80 fnavn
      logical ope
*------------------------------------------------------------------
      integer kf,ierr,na,nb
      logical def
 
      iflag=0

      if(kl.le.0) then
      call xstrip(fnavn,kf,def)
         if(kf.le.0) then
           iflag=30
           return
         end if         
      else
         kf=kl
      end if


      if( .not.ope) then
       open(unit=itape,file=fnavn(1:kf),err=117,status='old')
      end if


      call skipkom(nskip,itape,ierr)
      if(ierr.gt.0) then
      write(0,*)'skipkom: nskip,ierr=',nskip,ierr
         iflag=3
         close(itape)
         return
      end if

      n=0

 100  continue

      read(itape,*,end=200)na
      n=n+1
      if(n.gt.nmax) then
        iflag=1
        close(itape)
        return
      end if
      nx(n)=na
 
      go to 100

 200  continue
      close(itape)
      return

 117  continue
      iflag=50
      return

      end





********************************************************************
*
*                PROMPTFIL
*                
*     Prompts for name and reads data from two-column file
*
*     spor   - question                                                I
*     fnavn,defn  - Name on file and def value                         O
*     kl  -   lengtn of filename. If kl<=0 the name is assumed to be   O
*             ended by '#' or '!' and the length is calculated
*     x,y  -  arrays for data                                          O
*     nmax  - maximum number of data                                   I
*     n  -    number of data read                                      O
*     nskip  -  number of leading comment lines in file                O  
*     iflag  - error flag                                              O
*                0     : OK
*                1     : maximum number of data reached before end of file.
*                        The data that is read is returned in x and y
*               11     : file not opened
*                3     : error in comment processing ( no data ?)
*                30    : inappropriate filename 
*                50    : cannot open file
********************************************************************
      subroutine promptfil(spor,fnavn,defn,kl,x,y,nmax,n,nskip,iflag)
      integer kl,nmax,n,nskip,iflag
      real x(nmax),y(nmax)
      character*80 spor,fnavn,defn
*------------------------------------------------------------------
      integer itape
      logical asc

      itape=55
      iflag=0
   
      kl=79
      call filgen(itape,spor,defn,fnavn,kl,asc,.true.)

      if(itape.le.0 .or.(.not.asc)) then
          iflag=11
          return
      end if

      call lexyfil(fnavn,kl,x,y,nmax,n,.true.,itape,nskip,iflag)
 
      return
      end

********************************************************************
*
*                gpromptf
*                
*     Prompts for name and reads data from a one- or two-column file
*
*     spor   - question                                                I
*     fnavn,defn  - Name on file and def value                         O
*     kl  -   lengtn of filename. If kl<=0 the name is assumed to be   O
*             ended by '#' or '!' and the length is calculated
*     x,y  -  arrays for data                                          O
*     nmax  - maximum number of data                                   I
*     ic  - number of columns                                          I
*     n  -    number of data read                                      O
*     nskip  -  number of leading comment lines in file                O  
*     iflag  - error flag                                              O
*                0     : OK
*                1     : maximum number of data reached before end of file.
*                        The data that is read is returned in x and y
*               11     : file not opened
*                3     : error in comment processing ( no data ?)
*                30    : inappropriate filename 
*                50    : cannot open file
********************************************************************
      subroutine gpromptf(spor,fnavn,defn,kl,x,y,nmax,ic,n,nskip,iflag)
      integer kl,nmax,ic,n,nskip,iflag
      real x(nmax),y(nmax)
      character*80 spor,fnavn,defn
*------------------------------------------------------------------
      integer itape
      logical asc

      itape=55
      iflag=0
   
      kl=79
      call filgen(itape,spor,defn,fnavn,kl,asc,.true.)
      if(itape.le.0 .or.(.not.asc)) then
          iflag=11
          return
      end if

      if(ic.eq.2) then
       call lexyfil(fnavn,kl,x,y,nmax,n,.true.,itape,nskip,iflag)
      else
       call lexfil(fnavn,kl,x,nmax,n,.true.,itape,nskip,iflag)
      end if

      return
      end

********************************************************************
*
*                npromptf
*                
*     Prompts for name and reads data from a one- or two-column file
*     with integers
*
*     spor   - question                                                I
*     fnavn,defn  - Name on file and def value                         O
*     kl  -   lengtn of filename. If kl<=0 the name is assumed to be   O
*             ended by '#' or '!' and the length is calculated
*     nx,ny  -  arrays for data                                          O
*     nmax  - maximum number of data                                   I
*     ic  - number of columns
*     n  -    number of data read                                      O
*     nskip  -  number of leading comment lines in file                O  
*     iflag  - error flag                                              O
*                0     : OK
*                1     : maximum number of data reached before end of file.
*                        The data that is read is returned in x and y
*               11     : file not opened
*                3     : error in comment processing ( no data ?)
*                30    : inappropriate filename 
*                50    : cannot open file
********************************************************************
      subroutine npromptf(spor,fnavn,defn,kl,nx,ny,nmax,ic,n,
     %  nskip,iflag)
      integer kl,nmax,ic,n,nskip,iflag
      integer nx(nmax),ny(nmax)
      character*80 spor,fnavn,defn
*------------------------------------------------------------------
      integer itape
      logical asc

      itape=55
      iflag=0
   
      kl=79
      call filgen(itape,spor,defn,fnavn,kl,asc,.true.)
      if(itape.le.0 .or.(.not.asc)) then
          iflag=11
          return
      end if

      if(ic.eq.2) then
       call lenmfil(fnavn,kl,nx,ny,nmax,n,.true.,itape,nskip,iflag)
      else
       call lenfil(fnavn,kl,nx,nmax,n,.true.,itape,nskip,iflag)
      end if

      return
      end


********************************************************************
*
*                nnpromptf
*                
*     Prompts for name and reads data from one, two or three-column file
*
*     spor   - question                                                I
*     fnavn,defn  - Name on file and def value                         O
*     kl  -   lengtn of filename. If kl<=0 the name is assumed to be   O
*             ended by '#' or '!' and the length is calculated
*     x,y,z  -  arrays for data                                          O
*     nmax  - maximum number of data                                   I
*     ic  - number of columns
*     n  -    number of data read                                      O
*     nskip  -  number of leading comment lines in file                O  
*     iflag  - error flag                                              O
*                0     : OK
*                1     : maximum number of data reached before end of file.
*                        The data that is read is returned in x and y
*               11     : file not opened
*                3     : error in comment processing ( no data ?)
*                30    : inappropriate filename 
*                50    : cannot open file
********************************************************************
      subroutine nnpromptf(spor,fnavn,defn,kl,x,y,z,nmax,ic,n,
     %nskip,iflag)
      integer kl,nmax,ic,n,nskip,iflag
      real x(nmax),y(nmax),z(nmax)
      character*80 spor,fnavn,defn
*------------------------------------------------------------------
      integer itape
      logical asc

      itape=55
      iflag=0
   
      kl=79

      call filgen(itape,spor,defn,fnavn,kl,asc,.true.)
      if(itape.le.0 .or.(.not.asc)) then
          iflag=11
          return
      end if


      if(ic.eq.3) then
       call lexyzfil(fnavn,kl,x,y,z,nmax,n,.true.,itape,nskip,iflag)
      else
         if(ic.eq.2) then
           call lexyfil(fnavn,kl,x,y,nmax,n,.true.,itape,nskip,iflag)
         else
           call lexfil(fnavn,kl,x,nmax,n,.true.,itape,nskip,iflag)
         end if
      end if

      return
      end

      SUBROUTINE LBLANK(C,K)
      INTEGER K
      CHARACTER*200 C
*--------------------------------------------------------------
      INTEGER I

      DO 100 I=1,K
  100 C(i:i)=' '

      RETURN
      END

**************************************************************
      subroutine utkopi
      integer ktape
      include 'simfil.inc'
*-----------------------------------------------------------
      integer ilen
      call strip(bilde,irstopp,ilen)
      if(iut.ge.0.and.ilen.gt.0)
     %     write(iut,'(a)') bilde(1:ilen)
      return
      end


**************************************************************
      subroutine utbilde
      integer ktape
      include 'simfil.inc'
*-----------------------------------------------------------
      if(iut.ge.0)write(iut,'(a)') utbuff(1:utpos)
      call lblank(utbuff,200)
      utpos=0
      return
      end

************************************************************
      subroutine inbilde(slt,err)
      logical slt,err
      include 'simfil.inc'
*-----------------------------------------------------------
      integer k,ichar
      logical eks

      err=.false.
      slt=.false.

      read(itape,'(a)',err=200,end=300) bilde(1:irec)
      if(cm.eq.' ') then
        irstopp=irec
      else
        call lokchar(bilde,cm,1,irec,k,eks)
        if(eks) then
          irstopp=k-1
        else
          irstopp=irec
        end if
      end if

c    bytter kontroll-tegn og tabulatorer med blanke

      do 50 k=1,irstopp
      if(ichar(bilde(k:k)).lt.32) bilde(k:k)=' '
 50   continue


      ipos=1

      return
 200  continue
      err=.true.
      write(0,*)'feil ved linje av unit=',itape
      return

 300  continue
      slt=.true.
      return
      end



*************************************************************
      subroutine opsfil(itapez,irecz,iutz,iautz,cms,err)
      integer itapez,irecz,iutz,iautz
      character cms
      logical err
      include 'simfil.inc'
*-----------------------------------------------------------
      logical eks,slutt

      itape=itapez
      irec=irecz
      iut=iutz
      iaut=iautz
      cm=cms
      err=.false.

      call  lblank(utbuff,200)

      utpos=0

      call inbilde(slutt,err)

      return
      end


*********************************************************************
      subroutine utord(w,kw)
      integer kw
      Character*200 w
      include 'simfil.inc'
*--------------------------------------------------------------------
      integer nypos,nstart

      if(utpos.eq.0)  then
        nstart=1
      else
        nstart=utpos+2
      end if

      nypos=nstart+kw-1
      if(nypos.gt.200) then
        call utbilde
        nstart=1
        nypos=kw
      end if

      utbuff(nstart:nypos)=w(1:kw)
      utpos=nypos

      return
      end
        

*********************************************************************
      subroutine simord(w,kw,ctc)
      integer kw
      character*200 w
      character ctc
      include 'simfil.inc'
*--------------------------------------------------------------------
      integer i0,i1,iflag,k0,k1,kflag,k
      real tmp
      logical eks,slt,err

 100  continue

      call lokord(bilde,ipos,irstopp,i0,i1,iflag)

      if(iflag.eq.0) then
        if(iaut.eq.1 .and.utpos.gt.0) call utbilde
        call inbilde(slt,err)
        if(slt) then
          ctc='e'
          return
        end if
        go to 100
      end if

      kw=i1-i0+1
      w(1:kw)=bilde(i0:i1)
      call getr(bilde,i0,i1,k0,k1,tmp,kflag)
      if(kflag.eq.1) then
          ctc='n'
      else
       if(iflag.eq.1) then
          ctc='o'
        else
          ctc='c'
        end if
      end if

      ipos=i1+1
    
      return
      end



















***************************************************************************
*
*     GPRI
*
* printer array i plotxy-format p} fil
* parametere:
*         y(n)  - array                   I
*         n  - array-grense                  I
*         dx - gitterinkrement (betydning bare nar noun=.false.)          I
*         x - posisjoner  av y -verdier, naar noun=.false.                  I
*             benyttes bare x(1) som er pos av y(1)
*         noun - .false. angir uniformt gitter             I
*         isyk - nummer paa fil              I
*         stem - navnestamme, feks stem='eta#', isyk=11 gir datafil eta11 I
*                hvis stem ender paa ! ignoreres nummer
*         fak - skaleringsfaktor for y, slik at y/skal skrives ut     I
*         fakx - skaleringsfaktor for x, slik at x/fakx skrives ut    I
****************************************************************************
      subroutine gpri(y,n,dx,x,noun,isyk,stem,fak,fakx)
      integer isyk,n
      real y(n),dx,x(n),fak,fakx
      logical noun
      character*80 stem
*--------------------------------------------------------
      integer i,kf,kf2
      real xver,yver,faki,faxi
      character*80 tall,utfil
      logical def

      faki=1.0/fak
      faxi=1.0/fakx

      call xstrip(stem,kf2,def)
      utfil(1:kf2)=stem(1:kf2)

      if(def) then
        call itconv(isyk,kf,tall)
        utfil(kf2+1:kf2+kf)=tall(1:kf)
        kf=kf2+kf
      else
         kf=kf2
      end if

      open(unit=30,file=utfil(1:kf))


      do 100 i=1,n
      if(noun) then
        xver=faxi*x(i)
      else
        xver=faxi*(x(1)+(i-1)*dx)
      end if
      yver=faki*y(i)
      write(30,*)xver,yver
 100  continue

      close(30)

      return
      end

****************************************************************************
*
*     As gpri but y is an integer array. No stretch of y then
********************************************************************
      subroutine igpri(y,n,dx,x,noun,isyk,stem,fakx)
      integer isyk,n
      integer y(n)
      real dx,x(n),fak,fakx
      logical noun
      character*80 stem
*--------------------------------------------------------
      integer i,kf,kf2,yver
      real xver,faxi
      character*80 tall,utfil
      logical def


      faxi=1.0/fakx

      call xstrip(stem,kf2,def)
      utfil(1:kf2)=stem(1:kf2)

      if(def) then
        call itconv(isyk,kf,tall)
        utfil(kf2+1:kf2+kf)=tall(1:kf)
        kf=kf2+kf
      else
         kf=kf2
      end if

      open(unit=30,file=utfil(1:kf))


      do 100 i=1,n
      if(noun) then
        xver=faxi*x(i)
      else
        xver=faxi*(x(1)+(i-1)*dx)
      end if
      yver=y(i)
      write(30,*)xver,yver
 100  continue

      close(30)

      return
      end



***************************************************************************
*
*     CPRI
*
* printer complex array i plotxy-format p} fil.
*       Kolonner:  x Re Im abs arg
* parametere:
*         y(n)  - array                   I
*         n  - array-grense                  I
*         dx - gitterinkrement (betydning bare nar noun=.false.)          I
*         x - posisjoner  av y -verdier, naar noun=.false.                  I
*             benyttes bare x(1) som er pos av y(1)
*         noun - .false. angir uniformt gitter             I
*         isyk - nummer paa fil              I
*         stem - navnestamme, feks stem='eta#', isyk=11 gir datafil eta11 I
*                hvis stem ender paa ! ignoreres nummer
*         fak - skaleringsfaktor for y, slik at y/skal skrives ut     I
*         fakx - skaleringsfaktor for x, slik at x/fakx skrives ut    I
****************************************************************************
      subroutine cpri(y,n,dx,x,noun,isyk,stem,fak,fakx)
      integer isyk,n
      real dx,x(n),fak,fakx
      complex y(n)
      logical noun
      character*80 stem
*--------------------------------------------------------
      integer i,kf,kf2
      real xver,faki,faxi,rdel,imdel,abver,arg
      complex yver
      character*80 tall,utfil
      logical def

      faki=1.0/fak
      faxi=1.0/fakx

      call xstrip(stem,kf2,def)
      utfil(1:kf2)=stem(1:kf2)

      if(def) then
        call itconv(isyk,kf,tall)
        utfil(kf2+1:kf2+kf)=tall(1:kf)
        kf=kf2+kf
      else
         kf=kf2
      end if

      open(unit=30,file=utfil(1:kf))


      do 100 i=1,n
      if(noun) then
        xver=faxi*x(i)
      else
        xver=faxi*(x(1)+(i-1)*dx)
      end if
      yver=faki*y(i)
      rdel=real(yver)
cc      imdel=aimag(yver)
cc      abver=cabs(yver)
cc      arg=aimag(log(yver))
      imdel=imag(yver)
      abver=abs(yver)
      arg=imag(log(yver))
      write(30,'(5e16.6)')xver,rdel,imdel,abver,arg
 100  continue

      close(30)

      return
      end



 
***************************************************************************
*
*     IPRI
*
* printer integer array i plotxy-format p} fil
* parametere:
*         ny(n0:n)  - array                   I
*         n0,n  - array-grenser                  I
*         noun - .false. angir uniformt gitter             I
*         isyk - nummer paa fil              I
*         stem - navnestamme, feks stem='eta#', isyk=11 gir datafil eta11 I
****************************************************************************
      subroutine ipri(ny,n0,n,isyk,stem)
      integer isyk,n0,n,ny(n0:n)
      character*80 stem
*--------------------------------------------------------
      integer i,kf,kf2
      character*80 tall,utfil
      logical def


      call xstrip(stem,kf2,def)
      utfil(1:kf2)=stem(1:kf2)

      if(def) then
        call itconv(isyk,kf,tall)
        utfil(kf2+1:kf2+kf)=tall(1:kf)
        kf=kf2+kf
      else
         kf=kf2
      end if

      open(unit=30,file=utfil(1:kf))


      do 100 i=n0,n

      write(30,'(2i10)')i,ny(i)
 100  continue

      close(30)

      return
      end
 

***************************************************************************
*
*     YPRI
*
* prints a single array in a one-column ascii file
* parametere:
*         y(n)  - array                   I
*         n  - array-grense                  I
*         isyk - nummer paa fil              I
*         stem - navnestamme, feks stem='eta#', isyk=11 gir datafil eta11 I
*                hvis stem ender paa ! ignoreres nummer
*         fak - skaleringsfaktor for y, slik at y/skal skrives ut     I
****************************************************************************
      subroutine ypri(y,n,isyk,stem,fak)
      integer isyk,n
      real y(n),dx,x(n),fak,fakx
      logical noun
      character*80 stem
*--------------------------------------------------------
      integer i,kf,kf2
      real yver,faki
      character*80 tall,utfil
      logical def

      faki=1.0/fak

      call xstrip(stem,kf2,def)
      utfil(1:kf2)=stem(1:kf2)

      if(def) then
        call itconv(isyk,kf,tall)
        utfil(kf2+1:kf2+kf)=tall(1:kf)
        kf=kf2+kf
      else
         kf=kf2
      end if

      open(unit=30,file=utfil(1:kf))


      do 100 i=1,n
      yver=faki*y(i)
      write(30,*)yver
 100  continue

      close(30)

      return
      end

***************************************************************************
*
*     TPRI
*
* printer array i plotxy-format p} fil
* parametere:
*         y(nc,n)  - array                   I
*         nc - number of ordinate columns             I
*         n  - array-grense                  I
*         dx - gitterinkrement (betydning bare nar noun=.false.)          I
*         x - posisjoner  av y -verdier, naar noun=.false.                  I
*             benyttes bare x(1) som er pos av y(1)
*         noun - .false. angir uniformt gitter             I
*         isyk - nummer paa fil              I
*         stem - navnestamme, feks stem='eta#', isyk=11 gir datafil eta11 I
*                hvis stem ender paa ! ignoreres nummer
*         fak - skaleringsfaktor for y, slik at y/skal skrives ut     I
*         fakx - skaleringsfaktor for x, slik at x/fakx skrives ut    I
****************************************************************************
      subroutine tpri(y,nc,n,dx,x,noun,isyk,stem,fak,fakx)
      integer isyk,n,nc
      real y(nc,n),dx,x(n),fak,fakx
      logical noun
      character*80 stem
*--------------------------------------------------------
      integer i,kf,kf2,k
      real xver,faki,faxi,row(100)
      character*80 tall,utfil
      logical def

      faki=1.0/fak
      faxi=1.0/fakx

      call xstrip(stem,kf2,def)
      utfil(1:kf2)=stem(1:kf2)

      if(def) then
        call itconv(isyk,kf,tall)
        utfil(kf2+1:kf2+kf)=tall(1:kf)
        kf=kf2+kf
      else
         kf=kf2
      end if

      open(unit=30,file=utfil(1:kf))


      do 100 i=1,n
      if(noun) then
        xver=faxi*x(i)
      else
        xver=faxi*(x(1)+(i-1)*dx)
      end if
      do 80 k=1,nc
        row(k)=faki*y(k,i)
 80   continue
      write(30,*)xver,(row(k),k=1,nc)
 100  continue

      close(30)

      return
      end

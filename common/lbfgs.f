      SUBROUTINE VA15AD(N,M,X,F,G,DIAGCO,DIAG,IPRINT,EPS,S,Y,           
     *          POINT,W,IFLAG,ITER)                                     
C                                                                       
C                                                                       
      DOUBLE PRECISION X(N),G(N),S(M*N),Y(M*N),DIAG(N),W(N+2*M)         
      DOUBLE PRECISION FTOL,GTOL,XTOL,STPMIN,STPMAX,STP,F,YS,SQ,        
     *              YR,BETA,ONE,ZERO,EPS,XNORM,GNORM,YY,DDOT,STP1       
C                                                                       
      INTEGER BOUND,LP,ITER,NFUN,NFEV,IPRINT(2),POINT,CP,IFLAG          
      LOGICAL FINISH,DIAGCO                                             
C      COMMON /VA15DD/MP,LP, GTOL                                        
C MP, LP, GTOL ARE GIVEN HERE TO PREVENT COMMON IN MAIN 7/31/02 Y.HONDA
C      DATA   MP, LP /6, 6/
      COMMON /VA15DD/ MP, LP
      SAVE                                                              
      DATA ONE,ZERO/1.0D+0,0.0D+0/                                      
C                                                                       
C     ------------------------------------------------------------      
C     INITIALIZE                                                        
C     ------------------------------------------------------------      
C                                                                       
      IF(IFLAG.EQ.0) GO TO 1                                            
      GO TO (72,10) IFLAG                                               
   1  ITER= 0                                                           
C MP, LP, GTOL ARE GIVEN HERE TO PREVENT COMMON IN MAIN 7/31/02 Y.HONDA
      CALL VA15CD
C     
      IF(N.LE.0.OR.M.LE.0) GO TO 96                                     
      IF(GTOL.LE.1.D-04) THEN                                           
        IF(LP.GT.0) WRITE(LP,145)                                       
        GTOL=1.D-02                                                     
      ENDIF                                                             
      NFUN= 1                                                           
      POINT= 0                                                          
      FINISH= .FALSE.                                                   
      IF(DIAGCO) THEN                                                   
         DO 3 I=1,N                                                     
 3       IF (DIAG(I).LE.ZERO) GO TO 95                                  
      ELSE                                                              
         DO 4 I=1,N                                                     
 4       DIAG(I)= 1.0D0                                                 
      ENDIF                                                             
      DO 5 I=1,N                                                        
 5    S(I)= -G(I)*DIAG(I)                                               
      GNORM= DSQRT(DDOT(N,G,1,G,1))                                     

      STP1= ONE/GNORM                                                   
C                                                                       
C     PARAMETERS FOR LINE SEARCH ROUTINE                                
C     ----------------------------------                                
      FTOL= 1.0D-4                                                      
      GTOL= 9.0D-1
      XTOL= 1.0D-17                                                     
      STPMIN= 1.0D-20                                                   
      STPMAX= 1.0D+20                                                   
      MAXFEV= 20                                                        
C                                                                       
       CALL VA15BD(IPRINT,ITER,NFUN,                                    
     *                     N,M,X,F,G,STP,FINISH)                        
C                                                                       
C    ------------------------------------------------------------       
C     MAIN ITERATION LOOP                                               
C    --------------------------------------------------------           
C                                                                       
 8    ITER= ITER+1                                                      
      INFO=0                                                            
      BOUND=ITER-1                                                      
      IF (ITER .GT. M)BOUND=M                                           
      IF(ITER.EQ.1) GO TO 65                                            
C                                                                       
C     ------------------------------------------------------------      
C     COMPUTE -HG AND THE DIAGONAL SCALING MATRIX IN DIAG               
C     ------------------------------------------------------------      
C                                                                       
      IF(.NOT.DIAGCO) THEN                                              
CZOU         PRINT*,'form scaling ITER=', ITER                          
         DO 9 I=1,N                                                     
   9     DIAG(I)= YS/YY                                                 
      ELSE                                                              
         IFLAG=2                                                        
         RETURN                                                         
      ENDIF                                                             
  10  CONTINUE                                                          
      DO 11 I=1,N                                                       
  11  IF (DIAG(I).LE.ZERO) GO TO 95                                     
C                                                                       
      CP= POINT                                                         
      IF (POINT.EQ.0) CP=M                                              
      W(N+CP)= ONE/YS                                                   
      DO 12 I=1,N                                                       
  12  W(I)= -G(I)                                                       
      CP= POINT                                                         
      DO 25 II= 1,BOUND                                                 
         CP=CP-1                                                        
         IF (CP.EQ. -1)CP=M-1                                           
         SQ= DDOT(N,S(CP*N+1),1,W,1)                                    
         W(N+M+CP+1)= W(N+CP+1)*SQ                                      
         DO 20 K=1,N                                                    
  20     W(K)= W(K)-W(N+M+CP+1)*Y(CP*N+K)                               
  25  CONTINUE                                                          
C                                                                       
      DO 30 I=1,N                                                       
  30  W(I)=DIAG(I)*W(I)                                                 
      DO 45 II=1,BOUND                                                  
         YR= DDOT(N,Y(CP*N+1),1,W,1)                                    
         BETA= W(N+CP+1)*YR                                             
         DO 40 K=1,N                                                    
  40     W(K)= W(K)+S(CP*N+K)*(W(N+M+CP+1)-BETA)                        
         CP=CP+1                                                        
         IF (CP.EQ.M)CP=0                                               
  45  CONTINUE                                                          
C                                                                       
C     ------------------------------------------------------------      
C     STORE THE NEW DIRECTION IN S                                      
C     ------------------------------------------------------------      
C                                                                       
       DO 60 J=1,N                                                      
  60   S(POINT*N+J)= W(J)                                               
C                                                                       
C     ------------------------------------------------------------      
C     OBTAIN THE MINIMIZER OF THE FUNCTION ALONG THE                    
C     DIRECTION S BY USING THE LINE SEARCH ROUTINE OF VD05AD            
C     ------------------------------------------------------------      
  65  NFEV=0                                                            
      STP=ONE                                                           
      IF (ITER.EQ.1) STP=STP1                                           
      DO 70 I=1,N                                                       
  70  W(I)=G(I)                                                         
  72  CONTINUE                                                          
C                                                                       
      CALL VD05AD(N,X,F,G,S(POINT*N+1),STP,FTOL,GTOL,                   
     *            XTOL,STPMIN,STPMAX,MAXFEV,INFO,NFEV,DIAG)             
C                                                                       
      IF (INFO .EQ. -1) THEN                                            
        IFLAG=1                                                         
        RETURN                                                          
      ENDIF                                                             
      IF (INFO .NE. 1) GO TO 90                                         
      NFUN= NFUN + NFEV                                                 
C                                                                       
C     ------------------------------------------------------------      
C     COMPUTE THE NEW S AND Y                                           
C     ------------------------------------------------------------      
C                                                                       
      NPT=POINT*N                                                       
      DO 75 I=1,N                                                       
      S(NPT+I)= STP*S(NPT+I)                                            
  75  Y(NPT+I)= G(I)-W(I)                                               
      YS= DDOT(N,Y(NPT+1),1,S(NPT+1),1)                                 
      YY= DDOT(N,Y(NPT+1),1,Y(NPT+1),1)                                 
      POINT=POINT+1                                                     
      IF (POINT.EQ.M)POINT=0                                            
C                                                                       
C     ------------------------------------------------------------      
C     CONVERGENCE CHECK                                                 
C     ------------------------------------------------------------      
C                                                                       
      GNORM= DDOT(N,G,1,G,1)                                            
      GNORM=DSQRT(GNORM)                                                
      XNORM= DDOT(N,X,1,X,1)                                            
      XNORM=DSQRT(XNORM)                                                
      XNORM= DMAX1(1.0D0,XNORM)                                         
C                                                                       
CZOU start                                                              
C      XIAO=GNORM/XNORM                                                 
C      PRINT*,'GNORM/(1, XNORM)=', XIAO                                 
CZOU end                                                                
      IF (GNORM/XNORM .LE. EPS) FINISH=.TRUE.                           
C                                                                       
      CALL VA15BD(IPRINT,ITER,NFUN,                                     
     *               N,M,X,F,G,STP,FINISH)                              
      IF (FINISH) THEN                                                  
         IFLAG=0                                                        
         RETURN                                                         
      ENDIF                                                             
      GO TO 8                                                           
C                                                                       
C     ------------------------------------------------------------      
C     END OF MAIN ITERATION LOOP. ERROR EXITS.                          
C     ------------------------------------------------------------      
C                                                                       
  90  IF(LP.LE.0) RETURN                                                
      IF (INFO.EQ.0) THEN                                               
           IFLAG= -1                                                    
           WRITE(LP,100)IFLAG                                           
      ELSE IF (INFO.EQ.2) THEN                                          
           IFLAG= -2                                                    
           WRITE(LP,105)IFLAG                                           
      ELSE IF (INFO.EQ.3) THEN                                          
           IFLAG= -3                                                    
           WRITE(LP,110)IFLAG                                           
      ELSE IF (INFO.EQ.4) THEN                                          
           IFLAG= -4                                                    
           WRITE(LP,115)IFLAG                                           
      ELSE IF (INFO.EQ.5) THEN                                          
           IFLAG= -5                                                    
           WRITE(LP,120)IFLAG                                           
      ELSE IF (INFO.EQ.6) THEN                                          
           IFLAG= -6                                                    
           WRITE(LP,125)IFLAG                                           
      ENDIF                                                             
      RETURN                                                            
C                                                                       
  95  IFLAG= -7                                                         
      IF(LP.GT.0) WRITE(LP,135)IFLAG,I                                  
      RETURN                                                            
  96  IFLAG= -8                                                         
      IF(LP.GT.0) WRITE(LP,140)IFLAG                                    
C                                                                       
C     ------------------------------------------------------------      
C     FORMATS                                                           
C     ------------------------------------------------------------      
C                                                                       
 100  FORMAT(/' IFLAG= ',I2,/' IMPROPER INPUT PARAMETERS DURING',       
     .       ' THE LINE SEARCH.')                                       
 105  FORMAT(/' IFLAG= ',I2,/' RELATIVE WIDTH OF THE INTERVAL OF',      
     .       ' UNCERTAINTY IN THE LINE SEARCH'/'IS OF THE ORDER OF      
     .          MACHINE ROUNDOFF.')                                     
 110  FORMAT(/' IFLAG= ',I2,/' NUMBER OF CALLS TO FUNCTION IN THE',     
     .       ' LINE SEARCH HAS REACHED 20.')                            
 115  FORMAT(/' IFLAG= ',I2,/' THE STEP IN THE LINE SEARCH IS',         
     .       ' TOO SMALL.')                                             
 120  FORMAT(/' IFLAG= ',I2,/' THE STEP IN THE LINE SEARCH IS',         
     .       ' TOO LARGE.')                                             
 125  FORMAT(/' IFLAG= ',I2,/' ROUNDING ERRORS PREVENT FURTHER',        
     .       ' PROGRESS IN THE LINE SEARCH.')                           
 135  FORMAT(/' IFLAG= ',I2,/' THE',I5,'-TH DIAGONAL ELEMENT OF THE',   
     .       ' INVERSE HESSIAN APPROXIMATION IS NOT POSITIVE')          
 140  FORMAT(/' IFLAG= ',I2,/' IMPROPER INPUT PARAMETERS (N OR M',      
     .       ' ARE NOT POSITIVE)')                                      
 145  FORMAT(/'  GTOL IS LESS THAN OR EQUAL TO 1.D-04',                 
     .       / 'IT HAS BEEN RESET TO 1.D-02')                           
      RETURN                                                            
      END                                                               
C                                                                       
C                                                                       
C                                                                       
C                                                                       
      SUBROUTINE VA15BD(IPRINT,ITER,NFUN,                               
     *                     N,M,X,F,G,STP,FINISH)                        
C                                                                       
C     ------------------------------------------------------------------
C     THIS ROUTINE PRINTS MONITORING INFORMATION. THE FREQUENCY AND AMOU
C     OF OUTPUT ARE SPECIFIED AS FOLLOWS:                               
C                                                                       
C     IPRINT(1) < 0 : NO OUTPUT IS GENERATED                            
C     IPRINT(1) = 0 : OUTPUT ONLY AT FIRST AND LAST ITERATION           
C     IPRINT(1) > 0 : OUTPUT EVERY IPRINT(1) ITERATION                  
C     IPRINT(2) = 0 : ITERATION COUNT, FUNCTION VALUE, NORM OF THE GRADI
C                     ,NUMBER OF FUNCTION CALLS AND STEP LENGTH         
C     IPRINT(2) = 1 : + VECTOR OF VARIABLES AND GRADIENT VECTOR AT THE  
C                       INITIAL POINT                                   
C     IPRINT(2) = 2 : + VECTOR OF VARIABLES                             
C     IPRINT(2) = 3 : + GRADIENT VECTOR                                 
C     ------------------------------------------------------------------
C                                                                       
      DOUBLE PRECISION X(N),G(N),F,GNORM,STP,FACTOR,DDOT,GTOL           
      INTEGER IPRINT(2),ITER,NFUN,PROB,LP                               
      LOGICAL FINISH                                                    
      COMMON /SET/ FACTOR,PROB                                          
      COMMON /VA15DD/MP,LP, GTOL                                        
C                                                                       
      IF (IPRINT(1).LT.0)RETURN                                         
      GNORM= DDOT(N,G,1,G,1)                                            
      GNORM= DSQRT(GNORM)                                               

      IF (ITER.EQ.0)THEN                                                
           WRITE(MP,10)                                                 
           WRITE(MP,20) PROB,N,M                                        
           WRITE(MP,30)F,GNORM                                          
                 IF (IPRINT(2).GE.1)THEN                                
                     WRITE(MP,40)                                       
                     WRITE(MP,50) (X(I),I=1,N)                          
                     WRITE(MP,60)                                       
                     WRITE(MP,50) (G(I),I=1,N)                          
                  ENDIF                                                 
           WRITE(MP,10)                                                 
           WRITE(MP,70)                                                 
      ELSE                                                              
          IF ((IPRINT(1).EQ.0).AND.(ITER.NE.1.AND..NOT.FINISH))RETURN   
              IF (IPRINT(1).NE.0)THEN                                   
                   IF(MOD(ITER-1,IPRINT(1)).EQ.0.OR.FINISH)THEN         
                         WRITE(MP,80)ITER,NFUN,F,GNORM,STP              
                   ELSE                                                 
                         RETURN                                         
                   ENDIF                                                
              ELSE                                                      
                   WRITE(MP,80)ITER,NFUN,F,GNORM,STP                    
              ENDIF                                                     
              IF (IPRINT(2).EQ.2.OR.IPRINT(2).EQ.3)THEN                 
                    IF (FINISH)THEN                                     
                        WRITE(MP,90)                                    
               WRITE(MP,50) (X(I),I=1,N)                                
        WRITE(MP,60)                                                    
        WRITE(MP,50) (G(I),I=1,N)                                       
                    ENDIF                                               
              ENDIF                                                     
            IF (FINISH) WRITE(MP,100)                                   
      ENDIF                                                             
C                                                                       
C10   FORMAT('*************************************************')       
 10   FORMAT(/'*************************************************')      
C20   FORMAT(' PROB=',I3,'   N=',I5,'   NUMBER OF CORRECTIONS=',I2)     
 20   FORMAT(' PROB=',I3,'   N=',I7,'   NUMBER OF CORRECTIONS=',I2)     
 30   FORMAT('    ',E15.8,'          ',E15.8)                           
 40   FORMAT(' VECTOR X= ')                                             
 50   FORMAT(6(2X,1E15.8))                                              
 60   FORMAT(' GRADIENT VECTOR G= ')                                    
C70   FORMAT(/'I    NFN',4X,'FUNC',8X,'GNORM',7X,'STEPLENGTH'/)         
 70   FORMAT(/3X,'I',2X,'NFN',11X,'FUNC',12X,'GNORM',9X,'STEPLENGTH'/)  
 80   FORMAT(2(I4,1X),3X,3(E15.8,2X))                                   
 90   FORMAT(' FINAL POINT X= ')                                        
 100  FORMAT(/' THE MINIMIZATION TERMINATED WITHOUT DETECTING ERRORS.', 
     .       /' IFLAG = 0')                                             
C                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C   ----------------------------------------------------------          
C   DATA BLOCK                                                          
C   ----------------------------------------------------------          
C                                                                       
CTT   BLOCK DATA VA15CD                                                 
      SUBROUTINE VA15CD                                                 
      COMMON /VA15DD/MP,LP, GTOL                                        
      INTEGER LP                                                        
      DOUBLE PRECISION GTOL                                             
CTT   DATA MP,LP,GTOL/6,6,9.0D-01/                                      
c     MP   = 6                                                          
c     LP   = 6                                                          
c     GTOL = 9.0D-01                                                    
      MP   = 66                                                         
      LP   = 66                                                         
      GTOL = 9.0D-01                                                    
CTS                                                                     
      RETURN                                                            
CTE                                                                     
      END                                                               
C                                                                       
C   -------------------------------------------------------------       
C                                                                       
      SUBROUTINE VD05AD(N,X,F,G,S,STP,FTOL,GTOL,XTOL,                   
     *                  STPMIN,STPMAX,MAXFEV,INFO,NFEV,WA)              
      INTEGER N,MAXFEV,INFO,NFEV                                        
      DOUBLE PRECISION F,STP,FTOL,GTOL,XTOL,STPMIN,STPMAX               
      DOUBLE PRECISION X(N),G(N),S(N),WA(N)                             
      SAVE                                                              
C     **********                                                        
C                                                                       
C     SUBROUTINE VD05AD                                                 
C                                                                       
C     THE PURPOSE OF VD05AD IS TO FIND A STEP WHICH SATISFIES           
C     A SUFFICIENT DECREASE CONDITION AND A CURVATURE CONDITION.        
C     THE USER MUST PROVIDE A SUBROUTINE WHICH CALCULATES THE           
C     FUNCTION AND THE GRADIENT.                                        
C                                                                       
C     AT EACH STAGE THE SUBROUTINE UPDATES AN INTERVAL OF               
C     UNCERTAINTY WITH ENDPOINTS STX AND STY. THE INTERVAL OF           
C     UNCERTAINTY IS INITIALLY CHOSEN SO THAT IT CONTAINS A             
C     MINIMIZER OF THE MODIFIED FUNCTION                                
C                                                                       
C          F(X+STP*S) - F(X) - FTOL*STP*(GRADF(X)'S).                   
C                                                                       
C     IF A STEP IS OBTAINED FOR WHICH THE MODIFIED FUNCTION             
C     HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE DERIVATIVE,      
C     THEN THE INTERVAL OF UNCERTAINTY IS CHOSEN SO THAT IT             
C     CONTAINS A MINIMIZER OF F(X+STP*S).                               
C                                                                       
C     THE ALGORITHM IS DESIGNED TO FIND A STEP WHICH SATISFIES          
C     THE SUFFICIENT DECREASE CONDITION                                 
C                                                                       
C           F(X+STP*S) .LE. F(X) + FTOL*STP*(GRADF(X)'S),               
C                                                                       
C     AND THE CURVATURE CONDITION                                       
C                                                                       
C           ABS(GRADF(X+STP*S)'S)) .LE. GTOL*ABS(GRADF(X)'S).           
C                                                                       
C     IF FTOL IS LESS THAN GTOL AND IF, FOR EXAMPLE, THE FUNCTION       
C     IS BOUNDED BELOW, THEN THERE IS ALWAYS A STEP WHICH SATISFIES     
C     BOTH CONDITIONS. IF NO STEP CAN BE FOUND WHICH SATISFIES BOTH     
C     CONDITIONS, THEN THE ALGORITHM USUALLY STOPS WHEN ROUNDING        
C     ERRORS PREVENT FURTHER PROGRESS. IN THIS CASE STP ONLY            
C     SATISFIES THE SUFFICIENT DECREASE CONDITION.                      
C                                                                       
C     THE SUBROUTINE STATEMENT IS                                       
C                                                                       
C        SUBROUTINE VD05AD(N,X,F,G,S,STP,FTOL,GTOL,XTOL,                
C                          STPMIN,STPMAX,MAXFEV,INFO,NFEV,WA)           
C     WHERE                                                             
C                                                                       
C       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER        
C         OF VARIABLES.                                                 
C                                                                       
C       X IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE         
C         BASE POINT FOR THE LINE SEARCH. ON OUTPUT IT CONTAINS         
C         X + STP*S.                                                    
C                                                                       
C       F IS A VARIABLE. ON INPUT IT MUST CONTAIN THE VALUE OF F        
C         AT X. ON OUTPUT IT CONTAINS THE VALUE OF F AT X + STP*S.      
C                                                                       
C       G IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE         
C         GRADIENT OF F AT X. ON OUTPUT IT CONTAINS THE GRADIENT        
C         OF F AT X + STP*S.                                            
C                                                                       
C       S IS AN INPUT ARRAY OF LENGTH N WHICH SPECIFIES THE             
C         SEARCH DIRECTION.                                             
C                                                                       
C       STP IS A NONNEGATIVE VARIABLE. ON INPUT STP CONTAINS AN         
C         INITIAL ESTIMATE OF A SATISFACTORY STEP. ON OUTPUT            
C         STP CONTAINS THE FINAL ESTIMATE.                              
C                                                                       
C       FTOL AND GTOL ARE NONNEGATIVE INPUT VARIABLES. TERMINATION      
C         OCCURS WHEN THE SUFFICIENT DECREASE CONDITION AND THE         
C         DIRECTIONAL DERIVATIVE CONDITION ARE SATISFIED.               
C                                                                       
C       XTOL IS A NONNEGATIVE INPUT VARIABLE. TERMINATION OCCURS        
C         WHEN THE RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY        
C         IS AT MOST XTOL.                                              
C                                                                       
C       STPMIN AND STPMAX ARE NONNEGATIVE INPUT VARIABLES WHICH         
C         SPECIFY LOWER AND UPPER BOUNDS FOR THE STEP.                  
C                                                                       
C       MAXFEV IS A POSITIVE INTEGER INPUT VARIABLE. TERMINATION        
C         OCCURS WHEN THE NUMBER OF CALLS TO FCN IS AT LEAST            
C         MAXFEV BY THE END OF AN ITERATION.                            
C                                                                       
C       INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS:              
C                                                                       
C         INFO = 0  IMPROPER INPUT PARAMETERS.                          
C                                                                       
C         INFO =-1  A RETURN IS MADE TO COMPUTE THE FUNCTION AND GRADIEN
C                                                                       
C         INFO = 1  THE SUFFICIENT DECREASE CONDITION AND THE           
C                   DIRECTIONAL DERIVATIVE CONDITION HOLD.              
C                                                                       
C         INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY       
C                   IS AT MOST XTOL.                                    
C                                                                       
C         INFO = 3  NUMBER OF CALLS TO FCN HAS REACHED MAXFEV.          
C                                                                       
C         INFO = 4  THE STEP IS AT THE LOWER BOUND STPMIN.              
C                                                                       
C         INFO = 5  THE STEP IS AT THE UPPER BOUND STPMAX.              
C                                                                       
C         INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS.           
C                   THERE MAY NOT BE A STEP WHICH SATISFIES THE         
C                   SUFFICIENT DECREASE AND CURVATURE CONDITIONS.       
C                   TOLERANCES MAY BE TOO SMALL.                        
C                                                                       
C       NFEV IS AN INTEGER OUTPUT VARIABLE SET TO THE NUMBER OF         
C         CALLS TO FCN.                                                 
C                                                                       
C       WA IS A WORK ARRAY OF LENGTH N.                                 
C                                                                       
C     SUBPROGRAMS CALLED                                                
C                                                                       
C       HARWELL-SUPPLIED...VD05BD                                       
C                                                                       
C       FORTRAN-SUPPLIED...ABS,MAX,MIN                                  
C                                                                       
C     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983           
C     JORGE J. MORE', DAVID J. THUENTE                                  
C                                                                       
C     **********                                                        
      INTEGER INFOC,J                                                   
      LOGICAL BRACKT,STAGE1                                             
      DOUBLE PRECISION DG,DGM,DGINIT,DGTEST,DGX,DGXM,DGY,DGYM,          
     *       FINIT,FTEST1,FM,FX,FXM,FY,FYM,P5,P66,STX,STY,              
     *       STMIN,STMAX,WIDTH,WIDTH1,XTRAPF,ZERO                       
      DATA P5,P66,XTRAPF,ZERO /0.5D0,0.66D0,4.0D0,0.0D0/                
      IF(INFO.EQ.-1) GO TO 45                                           
      INFOC = 1                                                         
C                                                                       
C     CHECK THE INPUT PARAMETERS FOR ERRORS.                            
C                                                                       
      IF (N .LE. 0 .OR. STP .LE. ZERO .OR. FTOL .LT. ZERO .OR.          
     *    GTOL .LT. ZERO .OR. XTOL .LT. ZERO .OR. STPMIN .LT. ZERO      
     *    .OR. STPMAX .LT. STPMIN .OR. MAXFEV .LE. 0) RETURN            
C                                                                       
C     COMPUTE THE INITIAL GRADIENT IN THE SEARCH DIRECTION              
C     AND CHECK THAT S IS A DESCENT DIRECTION.                          
C                                                                       
      DGINIT = ZERO                                                     
      DO 10 J = 1, N                                                    
         DGINIT = DGINIT + G(J)*S(J)                                    
   10    CONTINUE                                                       
      IF (DGINIT .GE. ZERO) RETURN                                      
C                                                                       
C     INITIALIZE LOCAL VARIABLES.                                       
C                                                                       
      BRACKT = .FALSE.                                                  
      STAGE1 = .TRUE.                                                   
      NFEV = 0                                                          
      FINIT = F                                                         
      DGTEST = FTOL*DGINIT                                              
      WIDTH = STPMAX - STPMIN                                           
      WIDTH1 = WIDTH/P5                                                 
      DO 20 J = 1, N                                                    
         WA(J) = X(J)                                                   
   20    CONTINUE                                                       
C                                                                       
C     THE VARIABLES STX, FX, DGX CONTAIN THE VALUES OF THE STEP,        
C     FUNCTION, AND DIRECTIONAL DERIVATIVE AT THE BEST STEP.            
C     THE VARIABLES STY, FY, DGY CONTAIN THE VALUE OF THE STEP,         
C     FUNCTION, AND DERIVATIVE AT THE OTHER ENDPOINT OF                 
C     THE INTERVAL OF UNCERTAINTY.                                      
C     THE VARIABLES STP, F, DG CONTAIN THE VALUES OF THE STEP,          
C     FUNCTION, AND DERIVATIVE AT THE CURRENT STEP.                     
C                                                                       
      STX = ZERO                                                        
      FX = FINIT                                                        
      DGX = DGINIT                                                      
      STY = ZERO                                                        
      FY = FINIT                                                        
      DGY = DGINIT                                                      
C                                                                       
C     START OF ITERATION.                                               
C                                                                       
   30 CONTINUE                                                          
C                                                                       
C        SET THE MINIMUM AND MAXIMUM STEPS TO CORRESPOND                
C        TO THE PRESENT INTERVAL OF UNCERTAINTY.                        
C                                                                       
         IF (BRACKT) THEN                                               
            STMIN = MIN(STX,STY)                                        
            STMAX = MAX(STX,STY)                                        
         ELSE                                                           
            STMIN = STX                                                 
            STMAX = STP + XTRAPF*(STP - STX)                            
            END IF                                                      
C                                                                       
C        FORCE THE STEP TO BE WITHIN THE BOUNDS STPMAX AND STPMIN.      
C                                                                       
         STP = MAX(STP,STPMIN)                                          
         STP = MIN(STP,STPMAX)                                          
C                                                                       
C        IF AN UNUSUAL TERMINATION IS TO OCCUR THEN LET                 
C        STP BE THE LOWEST POINT OBTAINED SO FAR.                       
C                                                                       
         IF ((BRACKT .AND. (STP .LE. STMIN .OR. STP .GE. STMAX))        
     *      .OR. NFEV .GE. MAXFEV-1 .OR. INFOC .EQ. 0                   
     *      .OR. (BRACKT .AND. STMAX-STMIN .LE. XTOL*STMAX)) STP = STX  
C                                                                       
C        EVALUATE THE FUNCTION AND GRADIENT AT STP                      
C        AND COMPUTE THE DIRECTIONAL DERIVATIVE.                        
C                                                                       
         DO 40 J = 1, N                                                 
            X(J) = WA(J) + STP*S(J)                                     
   40       CONTINUE                                                    
         INFO=-1                                                        
         RETURN                                                         
C                                                                       
   45    INFO=0                                                         
         NFEV = NFEV + 1                                                
         DG = ZERO                                                      
         DO 50 J = 1, N                                                 
            DG = DG + G(J)*S(J)                                         
   50       CONTINUE                                                    
         FTEST1 = FINIT + STP*DGTEST                                    
C                                                                       
C        TEST FOR CONVERGENCE.                                          
C                                                                       
         IF ((BRACKT .AND. (STP .LE. STMIN .OR. STP .GE. STMAX))        
     *      .OR. INFOC .EQ. 0) INFO = 6                                 
         IF (STP .EQ. STPMAX .AND.                                      
     *       F .LE. FTEST1 .AND. DG .LE. DGTEST) INFO = 5               
         IF (STP .EQ. STPMIN .AND.                                      
     *       (F .GT. FTEST1 .OR. DG .GE. DGTEST)) INFO = 4              
         IF (NFEV .GE. MAXFEV) INFO = 3                                 
         IF (BRACKT .AND. STMAX-STMIN .LE. XTOL*STMAX) INFO = 2         
         IF (F .LE. FTEST1 .AND. ABS(DG) .LE. GTOL*(-DGINIT)) INFO = 1  
C                                                                       
C        CHECK FOR TERMINATION.                                         
C                                                                       
         IF (INFO .NE. 0) RETURN                                        
C                                                                       
C        IN THE FIRST STAGE WE SEEK A STEP FOR WHICH THE MODIFIED       
C        FUNCTION HAS A NONPOSITIVE VALUE AND NONNEGATIVE DERIVATIVE.   
C                                                                       
         IF (STAGE1 .AND. F .LE. FTEST1 .AND.                           
     *       DG .GE. MIN(FTOL,GTOL)*DGINIT) STAGE1 = .FALSE.            
C                                                                       
C        A MODIFIED FUNCTION IS USED TO PREDICT THE STEP ONLY IF        
C        WE HAVE NOT OBTAINED A STEP FOR WHICH THE MODIFIED             
C        FUNCTION HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE      
C        DERIVATIVE, AND IF A LOWER FUNCTION VALUE HAS BEEN             
C        OBTAINED BUT THE DECREASE IS NOT SUFFICIENT.                   
C                                                                       
         IF (STAGE1 .AND. F .LE. FX .AND. F .GT. FTEST1) THEN           
C                                                                       
C           DEFINE THE MODIFIED FUNCTION AND DERIVATIVE VALUES.         
C                                                                       
            FM = F - STP*DGTEST                                         
            FXM = FX - STX*DGTEST                                       
            FYM = FY - STY*DGTEST                                       
            DGM = DG - DGTEST                                           
            DGXM = DGX - DGTEST                                         
            DGYM = DGY - DGTEST                                         
C                                                                       
C           CALL CSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY            
C           AND TO COMPUTE THE NEW STEP.                                
C                                                                       
            CALL VD05BD(STX,FXM,DGXM,STY,FYM,DGYM,STP,FM,DGM,           
     *                 BRACKT,STMIN,STMAX,INFOC)                        
C                                                                       
C           RESET THE FUNCTION AND GRADIENT VALUES FOR F.               
C                                                                       
            FX = FXM + STX*DGTEST                                       
            FY = FYM + STY*DGTEST                                       
            DGX = DGXM + DGTEST                                         
            DGY = DGYM + DGTEST                                         
         ELSE                                                           
C                                                                       
C           CALL VD05BD TO UPDATE THE INTERVAL OF UNCERTAINTY           
C           AND TO COMPUTE THE NEW STEP.                                
C                                                                       
            CALL VD05BD(STX,FX,DGX,STY,FY,DGY,STP,F,DG,                 
     *                 BRACKT,STMIN,STMAX,INFOC)                        
            END IF                                                      
C                                                                       
C        FORCE A SUFFICIENT DECREASE IN THE SIZE OF THE                 
C        INTERVAL OF UNCERTAINTY.                                       
C                                                                       
         IF (BRACKT) THEN                                               
            IF (ABS(STY-STX) .GE. P66*WIDTH1)                           
     *         STP = STX + P5*(STY - STX)                               
            WIDTH1 = WIDTH                                              
            WIDTH = ABS(STY-STX)                                        
            END IF                                                      
C                                                                       
C        END OF ITERATION.                                              
C                                                                       
         GO TO 30                                                       
C                                                                       
C     LAST CARD OF SUBROUTINE VD05AD.                                   
C                                                                       
      END                                                               
      SUBROUTINE VD05BD(STX,FX,DX,STY,FY,DY,STP,FP,DP,BRACKT,           
     *                 STPMIN,STPMAX,INFO)                              
      INTEGER INFO                                                      
      DOUBLE PRECISION STX,FX,DX,STY,FY,DY,STP,FP,DP,STPMIN,STPMAX      
      LOGICAL BRACKT,BOUND                                              
C     **********                                                        
C                                                                       
C     SUBROUTINE VD05BD                                                 
C                                                                       
C     THE PURPOSE OF VD05BD IS TO COMPUTE A SAFEGUARDED STEP FOR        
C     A LINESEARCH AND TO UPDATE AN INTERVAL OF UNCERTAINTY FOR         
C     A MINIMIZER OF THE FUNCTION.                                      
C                                                                       
C     THE PARAMETER STX CONTAINS THE STEP WITH THE LEAST FUNCTION       
C     VALUE. THE PARAMETER STP CONTAINS THE CURRENT STEP. IT IS         
C     ASSUMED THAT THE DERIVATIVE AT STX IS NEGATIVE IN THE             
C     DIRECTION OF THE STEP. IF BRACKT IS SET TRUE THEN A               
C     MINIMIZER HAS BEEN BRACKETED IN AN INTERVAL OF UNCERTAINTY        
C     WITH ENDPOINTS STX AND STY.                                       
C                                                                       
C     THE SUBROUTINE STATEMENT IS                                       
C                                                                       
C       SUBROUTINE VD05BD(STX,FX,DX,STY,FY,DY,STP,FP,DP,BRACKT,         
C                        STPMIN,STPMAX,INFO)                            
C                                                                       
C     WHERE                                                             
C                                                                       
C       STX, FX, AND DX ARE VARIABLES WHICH SPECIFY THE STEP,           
C         THE FUNCTION, AND THE DERIVATIVE AT THE BEST STEP OBTAINED    
C         SO FAR. THE DERIVATIVE MUST BE NEGATIVE IN THE DIRECTION      
C         OF THE STEP, THAT IS, DX AND STP-STX MUST HAVE OPPOSITE       
C         SIGNS. ON OUTPUT THESE PARAMETERS ARE UPDATED APPROPRIATELY.  
C                                                                       
C       STY, FY, AND DY ARE VARIABLES WHICH SPECIFY THE STEP,           
C         THE FUNCTION, AND THE DERIVATIVE AT THE OTHER ENDPOINT OF     
C         THE INTERVAL OF UNCERTAINTY. ON OUTPUT THESE PARAMETERS ARE   
C         UPDATED APPROPRIATELY.                                        
C                                                                       
C       STP, FP, AND DP ARE VARIABLES WHICH SPECIFY THE STEP,           
C         THE FUNCTION, AND THE DERIVATIVE AT THE CURRENT STEP.         
C         IF BRACKT IS SET TRUE THEN ON INPUT STP MUST BE               
C         BETWEEN STX AND STY. ON OUTPUT STP IS SET TO THE NEW STEP.    
C                                                                       
C       BRACKT IS A LOGICAL VARIABLE WHICH SPECIFIES IF A MINIMIZER     
C         HAS BEEN BRACKETED. IF THE MINIMIZER HAS NOT BEEN BRACKETED   
C         THEN ON INPUT BRACKT MUST BE SET FALSE. IF THE MINIMIZER      
C         IS BRACKETED THEN ON OUTPUT BRACKT IS SET TRUE.               
C                                                                       
C       STPMIN AND STPMAX ARE INPUT VARIABLES WHICH SPECIFY LOWER       
C         AND UPPER BOUNDS FOR THE STEP.                                
C                                                                       
C       INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS:              
C         IF INFO = 1,2,3,4,5, THEN THE STEP HAS BEEN COMPUTED          
C         ACCORDING TO ONE OF THE FIVE CASES BELOW. OTHERWISE           
C         INFO = 0, AND THIS INDICATES IMPROPER INPUT PARAMETERS.       
C                                                                       
C     SUBPROGRAMS CALLED                                                
C                                                                       
C       FORTRAN-SUPPLIED ... ABS,MAX,MIN,SQRT                           
C                                                                       
C     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983           
C     JORGE J. MORE', DAVID J. THUENTE                                  
C                                                                       
C     **********                                                        
      DOUBLE PRECISION GAMMA,P,Q,R,S,SGND,STPC,STPF,STPQ,THETA          
      INFO = 0                                                          
C                                                                       
C     CHECK THE INPUT PARAMETERS FOR ERRORS.                            
C                                                                       
      IF ((BRACKT .AND. (STP .LE. MIN(STX,STY) .OR.                     
     *    STP .GE. MAX(STX,STY))) .OR.                                  
     *    DX*(STP-STX) .GE. 0.0 .OR. STPMAX .LT. STPMIN) RETURN         
C                                                                       
C     DETERMINE IF THE DERIVATIVES HAVE OPPOSITE SIGN.                  
C                                                                       
      SGND = DP*(DX/ABS(DX))                                            
C                                                                       
C     FIRST CASE. A HIGHER FUNCTION VALUE.                              
C     THE MINIMUM IS BRACKETED. IF THE CUBIC STEP IS CLOSER             
C     TO STX THAN THE QUADRATIC STEP, THE CUBIC STEP IS TAKEN,          
C     ELSE THE AVERAGE OF THE CUBIC AND QUADRATIC STEPS IS TAKEN.       
C                                                                       
      IF (FP .GT. FX) THEN                                              
         INFO = 1                                                       
         BOUND = .TRUE.                                                 
         THETA = 3*(FX - FP)/(STP - STX) + DX + DP                      
         S = MAX(ABS(THETA),ABS(DX),ABS(DP))                            
         GAMMA = S*SQRT((THETA/S)**2 - (DX/S)*(DP/S))                   
         IF (STP .LT. STX) GAMMA = -GAMMA                               
         P = (GAMMA - DX) + THETA                                       
         Q = ((GAMMA - DX) + GAMMA) + DP                                
         R = P/Q                                                        
         STPC = STX + R*(STP - STX)                                     
         STPQ = STX + ((DX/((FX-FP)/(STP-STX)+DX))/2)*(STP - STX)       
         IF (ABS(STPC-STX) .LT. ABS(STPQ-STX)) THEN                     
            STPF = STPC                                                 
         ELSE                                                           
           STPF = STPC + (STPQ - STPC)/2                                
           END IF                                                       
         BRACKT = .TRUE.                                                
C                                                                       
C     SECOND CASE. A LOWER FUNCTION VALUE AND DERIVATIVES OF            
C     OPPOSITE SIGN. THE MINIMUM IS BRACKETED. IF THE CUBIC             
C     STEP IS CLOSER TO STX THAN THE QUADRATIC (SECANT) STEP,           
C     THE CUBIC STEP IS TAKEN, ELSE THE QUADRATIC STEP IS TAKEN.        
C                                                                       
      ELSE IF (SGND .LT. 0.0) THEN                                      
         INFO = 2                                                       
         BOUND = .FALSE.                                                
         THETA = 3*(FX - FP)/(STP - STX) + DX + DP                      
         S = MAX(ABS(THETA),ABS(DX),ABS(DP))                            
         GAMMA = S*SQRT((THETA/S)**2 - (DX/S)*(DP/S))                   
         IF (STP .GT. STX) GAMMA = -GAMMA                               
         P = (GAMMA - DP) + THETA                                       
         Q = ((GAMMA - DP) + GAMMA) + DX                                
         R = P/Q                                                        
         STPC = STP + R*(STX - STP)                                     
         STPQ = STP + (DP/(DP-DX))*(STX - STP)                          
         IF (ABS(STPC-STP) .GT. ABS(STPQ-STP)) THEN                     
            STPF = STPC                                                 
         ELSE                                                           
            STPF = STPQ                                                 
            END IF                                                      
         BRACKT = .TRUE.                                                
C                                                                       
C     THIRD CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE            
C     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DECREASES.         
C     THE CUBIC STEP IS ONLY USED IF THE CUBIC TENDS TO INFINITY        
C     IN THE DIRECTION OF THE STEP OR IF THE MINIMUM OF THE CUBIC       
C     IS BEYOND STP. OTHERWISE THE CUBIC STEP IS DEFINED TO BE          
C     EITHER STPMIN OR STPMAX. THE QUADRATIC (SECANT) STEP IS ALSO      
C     COMPUTED AND IF THE MINIMUM IS BRACKETED THEN THE THE STEP        
C     CLOSEST TO STX IS TAKEN, ELSE THE STEP FARTHEST AWAY IS TAKEN.    
C                                                                       
      ELSE IF (ABS(DP) .LT. ABS(DX)) THEN                               
         INFO = 3                                                       
         BOUND = .TRUE.                                                 
         THETA = 3*(FX - FP)/(STP - STX) + DX + DP                      
         S = MAX(ABS(THETA),ABS(DX),ABS(DP))                            
C                                                                       
C        THE CASE GAMMA = 0 ONLY ARISES IF THE CUBIC DOES NOT TEND      
C        TO INFINITY IN THE DIRECTION OF THE STEP.                      
C                                                                       
         GAMMA = S*SQRT(MAX(0.0D0,(THETA/S)**2 - (DX/S)*(DP/S)))        
         IF (STP .GT. STX) GAMMA = -GAMMA                               
         P = (GAMMA - DP) + THETA                                       
         Q = (GAMMA + (DX - DP)) + GAMMA                                
         R = P/Q                                                        
         IF (R .LT. 0.0 .AND. GAMMA .NE. 0.0) THEN                      
            STPC = STP + R*(STX - STP)                                  
         ELSE IF (STP .GT. STX) THEN                                    
            STPC = STPMAX                                               
         ELSE                                                           
            STPC = STPMIN                                               
            END IF                                                      
         STPQ = STP + (DP/(DP-DX))*(STX - STP)                          
         IF (BRACKT) THEN                                               
            IF (ABS(STP-STPC) .LT. ABS(STP-STPQ)) THEN                  
               STPF = STPC                                              
            ELSE                                                        
               STPF = STPQ                                              
               END IF                                                   
         ELSE                                                           
            IF (ABS(STP-STPC) .GT. ABS(STP-STPQ)) THEN                  
               STPF = STPC                                              
            ELSE                                                        
               STPF = STPQ                                              
               END IF                                                   
            END IF                                                      
C                                                                       
C     FOURTH CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE           
C     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DOES               
C     NOT DECREASE. IF THE MINIMUM IS NOT BRACKETED, THE STEP           
C     IS EITHER STPMIN OR STPMAX, ELSE THE CUBIC STEP IS TAKEN.         
C                                                                       
      ELSE                                                              
         INFO = 4                                                       
         BOUND = .FALSE.                                                
         IF (BRACKT) THEN                                               
            THETA = 3*(FP - FY)/(STY - STP) + DY + DP                   
            S = MAX(ABS(THETA),ABS(DY),ABS(DP))                         
            GAMMA = S*SQRT((THETA/S)**2 - (DY/S)*(DP/S))                
            IF (STP .GT. STY) GAMMA = -GAMMA                            
            P = (GAMMA - DP) + THETA                                    
            Q = ((GAMMA - DP) + GAMMA) + DY                             
            R = P/Q                                                     
            STPC = STP + R*(STY - STP)                                  
            STPF = STPC                                                 
         ELSE IF (STP .GT. STX) THEN                                    
            STPF = STPMAX                                               
         ELSE                                                           
            STPF = STPMIN                                               
            END IF                                                      
         END IF                                                         
C                                                                       
C     UPDATE THE INTERVAL OF UNCERTAINTY. THIS UPDATE DOES NOT          
C     DEPEND ON THE NEW STEP OR THE CASE ANALYSIS ABOVE.                
C                                                                       
      IF (FP .GT. FX) THEN                                              
         STY = STP                                                      
         FY = FP                                                        
         DY = DP                                                        
      ELSE                                                              
         IF (SGND .LT. 0.0) THEN                                        
            STY = STX                                                   
            FY = FX                                                     
            DY = DX                                                     
            END IF                                                      
         STX = STP                                                      
         FX = FP                                                        
         DX = DP                                                        
         END IF                                                         
C                                                                       
C     COMPUTE THE NEW STEP AND SAFEGUARD IT.                            
C                                                                       
      STPF = MIN(STPMAX,STPF)                                           
      STPF = MAX(STPMIN,STPF)                                           
      STP = STPF                                                        
      IF (BRACKT .AND. BOUND) THEN                                      
         IF (STY .GT. STX) THEN                                         
            STP = MIN(STX+0.66*(STY-STX),STP)                           
         ELSE                                                           
            STP = MAX(STX+0.66*(STY-STX),STP)                           
            END IF                                                      
         END IF                                                         
      RETURN                                                            
C                                                                       
C     LAST CARD OF SUBROUTINE VD05BD.                                   
C                                                                       
      END                                                               
                                                                        
                                                                        
                                                                        
                                                                        
       DOUBLE PRECISION FUNCTION DDOT(N,D,I1,S,I2)                      
C                                                                       
C      -------------------------------------------------------          
C      THIS FUNCTION COMPUTES THE INNER PRODUCT OF TWO VECTORS          
C      -------------------------------------------------------          
C                                                                       
       DOUBLE PRECISION D(N),S(N),PROD                                  
       INTEGER I1,I2                                                    
C                                                                       
        PROD=0.0D0                                                      
        DO 10 I=1,N                                                     
 10     PROD= PROD+D(I)*S(I)                                            
C                                                                       
        DDOT= PROD                                                      
C                                                                       
       RETURN                                                           
       END                                                              

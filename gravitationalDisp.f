      program gravitationalDisp
      implicit none

          !grav accel
          real, parameter :: g = 9.81
    
          !variables
          real :: s !displacement
          real :: t !time 
          real :: u !initial speed
            
          !assigning values
          t = 5.0 
          u = 50 
            
          !displacement calculation 
          s = u*t - g*(t**2) / 2
            
          print *, "Time = ", t 
          print *, "Displacement = ", s

      end program gravitationalDisp